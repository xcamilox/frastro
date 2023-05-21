from os import listdir
import os.path as path
import os
from frastro import FITSFile, Config, LOGUtil, ImageUtils
import shutil
from frastro import SwarpUtil,Utils, PSFexUtil,SextractorUtil,CatalogUtil
import datetime
import ntpath
from frastro.core.database.mongodb.mongodb_manager import MongodbManager
from frastro.core.pipeline.pyDIA_pipeline import PYDIAPipeline
from frastro.core.pipeline.reports.generate_reports import Generatereports
from frastro.external.galfit_util import GalfitUtils
from frastro.external.psfex_util import PSFexUtil
from frastro.external.la_cosmics import CosmicsRayRemove
import numpy as np

from astropy import units as u
from astropy import coordinates as coords
from astropy import wcs
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import figure, cm
from matplotlib.gridspec import GridSpec
import pdfkit

class LiverpollSnPipeline():



    def __init__(self):
        self.project_path = Config.getPath("liverpool_sn")

        self.mgdb = MongodbManager()
        self.mgdb.setCollection("liverpool_sn_v3")


    def run(self,data_path):

        files=self.getData(data_path)
        print("all files:",files)
        for path_object in files:
            LOGUtil.log("-----------Liverpool pipeline---------"+path_object,path_object)
            stack_file = self.stackFiles(path_object)
            if stack_file != None:
                catalog_file = self.createCatalogs(path_object)
                self.moveToScienceDirectory(path_object,[stack_file,catalog_file])
            print("complite single: ",path_object)

        """
        for path_allfile in allfiles:
            LOGUtil.log("-----------Liverpool pipeline all---------"+path_allfile, path_allfile)
            stack_file = self.stackFiles(path_allfile)
            if stack_file != None:
                catalog_file = self.createCatalogs(path_allfile)
                self.moveToScienceDirectory(path_allfile, [stack_file, catalog_file])
            print("complite all: ", path_allfile)
        """


    def moveToScienceDirectory(self,path_object,files):
        science_path=Utils.createPath(path_object+"/science")
        swarp_image = files[0]
        sex_catalog = files[1]
        current_date= datetime.datetime.now().strftime("%d_%m_%Y")
        if os.path.exists(swarp_image):
            shutil.copy2(swarp_image,science_path+"/stacking_"+current_date+".fits")

        if type(sex_catalog) is list:
            for idx,cat in enumerate(sex_catalog):
                if os.path.exists(cat):
                    shutil.copy2(cat,science_path+"/stacking_catalog_"+str(idx)+"_"+current_date+".fits")
        else:
            if os.path.exists(sex_catalog):
                shutil.copy2(sex_catalog,science_path+"/stacking_catalog_"+current_date+".fits")




    def readHeaderInfo(self,img_path):
        header = FITSFile.header(img_path)[0]
        head, file = ntpath.split(img_path)
        # save image info
        date = header["DATE"]  # epoach
        epoc_data = {}

        epoc_data["date"] = date
        epoc_data["id"] = header["OBJECT"]
        epoc_data["file"] = os.path.splitext(file)[0]
        epoc_data["file_path"] = img_path
        epoc_data["date_obs"] = header["DATE-OBS"]
        epoc_data["group_id"] = header["GROUPID"]
        epoc_data["mjd"] = header["MJD"]
        epoc_data["exptime"] = header["EXPTIME"]
        epoc_data["airmass"] = header["AIRMASS"]
        epoc_data["seeing_img"] = header["SEEING"]*0.30

        return epoc_data


    def generateModels(self):
        galfit = GalfitUtils()
        mgdb=self.mgdb
        psfex = PSFexUtil()
        all_data = mgdb.getData({'epoach': {'$exists': 'true'}})
        for target in all_data:
            try:
                ra = target["ra"]
                dec = target["dec"]
                for epoach in target["epoach"]:
                    stack = target["epoach"][epoach]["stacking"]
                    if "error" not in stack.keys():
                        psf_best = stack["best"]["psf_path"]
                        science_path = os.path.dirname(psf_best)
                        image = science_path + "/stack_best.fits"
                        image_all = science_path + "/stack_all.fits"

                        # detection data
                        detections_best = stack["best"]["detection"]["best"]
                        magnitud = detections_best["mag"]
                        axis_ratio = detections_best["ELONGATION"]
                        angle = detections_best["THETA_IMAGE"]
                        zp = stack["best"]["zp"]

                        # create model from sextrctor detection data
                        model_path = galfit.createModel(image, ra, dec, psfex_path=psf_best, magnitud=magnitud,
                                                        axis_ratio=axis_ratio, angle=angle, zp=zp)

                        # update database with the new model path
                        if model_path != "":
                            mgdb.update({'id': target["id"]},
                                        {"$set": {"epoach." + epoach + ".stacking.best.model": model_path}})



                        #model_path = galfit.createModel(image, ra, dec, psfex_path=psf_best)

                            #mgdb.update({'id': target["id"]}, {
                            #    "$set": {"epoach." + epoach + ".stacking.best.model": model_path,
                            #             "epoach." + epoach + ".stacking.best.file_path": image,
                            #             "epoach." + epoach + ".stacking.all.file_path": image_all}})
            except KeyError:
                print("error with ", target["id"])

    def coordinateFile(self,target):
        id=target["id"]
        ra=target["ra"]
        dec = target["dec"]
        coordinates_cat = self.project_path + id + "/"+id+"_coordinates.dat"

        if not path.exists(coordinates_cat):
            coordinates_file = open(coordinates_cat, "w+")
            coordinates_file.write("#id ra dec\n" + id + " " + str(ra) + " " + str(dec) + "\n")
            coordinates_file.close()
            self.mgdb.update({"id":id}, {"$set": {"coordinate_file": coordinates_cat}})


    def getData(self,path_folder,force=False,all=False,targetid=""):
        paths = {}
        onlyfiles = [f for f in os.listdir(path_folder) if path.isfile(path_folder+"/"+f)]


        for file in onlyfiles:
            if path.splitext(file)[1]== ".fits" or path.splitext(file)[1]== ".fit":
                header=FITSFile.header(path_folder+"/"+file)[0]
                target = header["OBJECT"]
                if "_offset" in target:
                    target = str(target).replace("_offset","")

                if "_of" in target:
                    target = str(target).replace("_of","")
                band = str(header["FILTER1"]).replace(" ","")

                if targetid != "" and targetid != target:
                    continue

                obs_date = header["DATE"]

                science_folder = self.project_path + target + "/science/"
                ext_catalogs = self.project_path + target + "/ext_catalogs/"
                observations_folder = self.project_path + target + "/observations/"
                target_path = observations_folder + obs_date + "/"


                # create science folder
                # create ext_catalogs folder
                # create observations folder

                if all:
                    if target not in paths.keys():
                        paths[target] = {obs_date: target_path}
                    elif obs_date not in paths[target].keys():
                        paths[target][obs_date] = target_path


                #check data base
                file_name = path.splitext(file)[0]
                curret_data=self.mgdb.getData({"id": target,'epoach.'+obs_date+'.obs.file':file_name})
                if len(curret_data)>0 and force==False:
                    print("ready ",target,obs_date,file_name)
                    continue
                else:

                    if all==False:
                        if target not in paths.keys():
                            paths[target] = {obs_date: target_path}
                        elif obs_date not in paths[target].keys():
                            paths[target][obs_date] = target_path


                    if not os.path.exists(science_folder):
                        os.makedirs(science_folder)

                    if not os.path.exists(ext_catalogs):
                        os.makedirs(ext_catalogs)

                    if not os.path.exists(observations_folder):
                        os.makedirs(observations_folder)


                    newlocation=target_path+file

                    if not os.path.exists(target_path):
                        os.makedirs(target_path)


                    #create files list to each apoach
                    files_list_eapoach=target_path+"files.lst"
                    if not os.path.exists(files_list_eapoach):
                        file_list_epoach_file =open(files_list_eapoach,"w+")
                    else:
                        file_list_epoach_file = open(files_list_eapoach, "a+")

                    # create files list for all files
                    files_list_all = observations_folder + "files.lst"
                    if not os.path.exists(files_list_all):
                        file_list_all_file = open(files_list_all, "w+")
                    else:
                        file_list_all_file = open(files_list_all, "a+")


                    #shutil.move(path_folder + "/" + file, newlocation)
                    shutil.copy2(path_folder+"/"+file,newlocation)


                    file_list_epoach_file.write(newlocation + "\n")
                    file_list_all_file.write(newlocation + "\n")

                    file_list_epoach_file.close()
                    file_list_all_file.close()


                    cat_info=self.createCatalog(newlocation)

                    header_info = self.readHeaderInfo(newlocation)
                    object_info={**cat_info,**header_info}

                    filer={"id":target}



                    query={ "$addToSet": {"epoach."+obs_date+".obs": object_info }}

                    self.mgdb.update(filer,query)
                    items=self.mgdb.getData({"id":target},{'id':1,"epoach":1})
                    if len(items)>0:
                        items=items[0]
                    if "epoach" in items.keys():
                        self.mgdb.update(items["id"], {"$set": {"n_epoach": len(items["epoach"])}})
        return paths

    def stackAllObsByTarget(self,target):
        print("stacking: ",target)
        swarp = SwarpUtil()
        good_seeing = []
        base_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/" + target + "/observations/"
        curret_data = self.mgdb.getData({"id": target})[0]
        ra_source = curret_data["ra"]
        dec_source = curret_data["dec"]
        best_list = []
        log_files_list = base_path + "log_files.lst"
        f_log = open(log_files_list, "w+")
        for epoachs in curret_data["epoach"]:
            for epoach_data in curret_data["epoach"][epoachs]["obs"]:
                file_path = epoach_data["file_path"]
                f_log.write(str(epoach_data["seeing_cal"])+"\t"+file_path + "\n" )
                if float(epoach_data["seeing_cal"])< 2.4:
                    good_seeing.append([float(epoach_data["seeing_cal"]), file_path])
            good_seeing.sort(key=lambda x: x[0])
        f_log.close()
        min_seeing = min(good_seeing)[0]

        for image in good_seeing:
            print(min_seeing + (min_seeing * 0.5), image[0])
            if image[0] <= min_seeing + (min_seeing * 0.5):
                best_list.append(image[1])

        if len(good_seeing) > 1:
            all_files_list = base_path + "all_list.lst"
            f = open(all_files_list, "w+")
            for fil in good_seeing:
                f.write(fil[1] + "\n")
            f.close()

        outpaht = base_path + "stacking/swarp/"
        if not os.path.exists(outpaht):
            os.makedirs(outpaht)

        if len(best_list) >= 1:

            if len(best_list) == 1:
                # remove cosmic ray
                best_list[0] = CosmicsRayRemove().cleanImage(best_list[0])

            best_files_list = base_path + "best_list.lst"

            if os.path.exists(best_files_list):
                os.remove(best_files_list)

            f = open(best_files_list, "w+")
            for fil in best_list:
                f.write(fil + "\n")
            f.close()

            kargs = {}
            kargs["-SUBTRACT_BACK"] = "Y"
            kargs["-RESAMPLING_TYPE"] = "NEAREST"  # No iterpolate
            kargs["-IMAGEOUT_NAME"] = "coadd_nearest_sapling.fits"  # No iterpolate
            kargs["-WEIGHTOUT_NAME"] = "weightout_coadd_nearest_sapling.fits"  # No iterpolate
            kargs["ra"] = ra_source
            kargs["dec"] = dec_source
            kargs["size"] = (1900, 1900)

            stack_file = swarp.run(best_files_list, outpaht, **kargs)
            science_dir = base_path + "science/"
            if not os.path.exists(science_dir):
                os.makedirs(science_dir)
            newlocation_best = science_dir + "stack_best.fits"
            shutil.copy2(stack_file, newlocation_best)

            best_cat = self.createCatalog(stack_file)

            newlocation_catdeep_best = science_dir + "catalog_best_deep.fits"
            shutil.copy2(best_cat["cat_path_deep"], newlocation_catdeep_best)
            best_cat["cat_path_deep"] = newlocation_catdeep_best

            newlocation_cat_best = science_dir + "catalog_best.fits"
            shutil.copy2(best_cat["cat_path"], newlocation_cat_best)

            base_psf_path = path.dirname(best_cat["psf_path"])
            if path.exists(base_psf_path + "/psfex_out.cat"):
                newlocation_psf_best_cat = science_dir + "psf_stars.fits"
                shutil.copy2(base_psf_path + "/psfex_out.cat", newlocation_psf_best_cat)
                best_cat["psf_stars"] = newlocation_psf_best_cat


            best_cat["cat_path"] = newlocation_cat_best

            newlocation_psf_best = science_dir + "psf_best.fits"
            shutil.copy2(best_cat["psf_path"], newlocation_psf_best)
            best_cat["psf_path"] = newlocation_psf_best
            best_cat["file_path"] = newlocation_best
            header = FITSFile.header(newlocation_best)[0]
            best_cat["exp_time"] = header["EXPTIME"]

            if len(good_seeing) > len(best_list):

                all_files_list = outpaht + "all_list.lst"
                if os.path.exists(all_files_list):
                    os.remove(all_files_list)

                f = open(all_files_list, "w+")
                tostack = 0
                for fil in good_seeing:
                    if fil[0] < 900:
                        tostack += 1
                        f.write(fil[1] + "\n")
                f.close()

                if tostack >= 1:

                    stack_file = swarp.run(all_files_list, outpaht, **kargs)
                    newlocation_all = science_dir + "stack_all.fits"
                    shutil.copy2(stack_file, newlocation_all)

                    stack_file = swarp.run(all_files_list, outpaht, **kargs)

                    all_cat = self.createCatalog(stack_file)

                    newlocation_catdeep_all = science_dir + "catalog_all_deep.fits"
                    shutil.copy2(all_cat["cat_path_deep"], newlocation_catdeep_all)
                    all_cat["cat_path_deep"] = newlocation_catdeep_all

                    newlocation_cat_all = science_dir + "catalog_all.fits"
                    shutil.copy2(all_cat["cat_path"], newlocation_cat_all)
                    all_cat["cat_path"] = newlocation_cat_all

                    newlocation_psf_all = science_dir + "psf_all.fits"
                    base_all_psf_path = path.dirname(all_cat["psf_path"])
                    if path.exists(base_all_psf_path + "/psfex_out.cat"):
                        newlocation_psf_all_cat = science_dir + "psf_stars_all.fits"
                        shutil.copy2(base_all_psf_path + "/psfex_out.cat", newlocation_psf_all_cat)
                        all_cat["psf_stars"] = newlocation_psf_all_cat

                    shutil.copy2(all_cat["psf_path"], newlocation_psf_all)
                    all_cat["psf_path"] = newlocation_psf_all

                    newlocation_all = science_dir + "stack_all.fits"
                    shutil.copy2(stack_file, newlocation_all)

                    all_cat["file_path"] = newlocation_all
                    header_all = FITSFile.header(newlocation_all)[0]
                    all_cat["exp_time"] = header_all["EXPTIME"]

                else:
                    all_cat = best_cat
            else:
                all_cat = best_cat

            self.mgdb.update({"id": target},
                             {'$set': {'stacking': {'best': best_cat, 'all': all_cat}}})
        else:
            self.mgdb.update({"id": target},
                             {'$set': {'stacking': {"error": "Single observation"}}})


    def stackEpoach(self,files,force=False):
        swarp = SwarpUtil()
        for targets in files:
            for epoachs in files[targets]:
                print(targets,epoachs)
                curret_data = self.mgdb.getData({"id": targets,"epoach."+epoachs:{"$exists":"true"}})[0]
                if len(curret_data)>0 and force==False:
                    print("ready stacking: ",targets, epoachs)
                    continue
                else:

                    good_seeing=[]

                    file = files[targets][epoachs]
                    base_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/" + targets + "/observations/" + epoachs + "/"
                    outpaht = base_path + "reduction/stacking/swarp/"
                    if not os.path.exists(outpaht):
                        os.makedirs(outpaht)
                    curret_data = self.mgdb.getData({"id": targets, 'epoach.'+epoachs : {'$exists': 'true'}},{ 'epoach.'+epoachs:1,'ra':1,'dec':1})[0]
                    ra_source=curret_data["ra"]
                    dec_source = curret_data["dec"]
                    for epoach_data in curret_data["epoach"][epoachs]["obs"]:
                        file_path=epoach_data["file_path"]
                        good_seeing.append([float(epoach_data["seeing_cal"]),file_path])
                    good_seeing.sort(key=lambda x:x[0])

                    min_seeing=min(good_seeing)[0]

                    best_list=[]
                    for image in good_seeing:
                        print(min_seeing+(min_seeing*0.3),image[0])
                        if image[0] <= min_seeing+(min_seeing*0.3):
                            best_list.append(image[1])

                    if len(best_list) >= 1 and min_seeing < 99:

                        if len(best_list) == 1:
                            #remove cosmic ray
                            best_list[0] = CosmicsRayRemove().cleanImage(best_list[0])


                        best_files_list=outpaht+"best_list.lst"

                        if os.path.exists(best_files_list):
                            os.remove(best_files_list)

                        f=open(best_files_list,"w+")

                        for fil in best_list:
                            f.write(fil+"\n")
                        f.close()


                        kargs={}
                        kargs["-SUBTRACT_BACK"] = "Y"
                        kargs["-RESAMPLING_TYPE"] = "NEAREST"  # No iterpolate
                        kargs["-IMAGEOUT_NAME"] = "coadd_nearest_sapling.fits"  # No iterpolate
                        kargs["-WEIGHTOUT_NAME"] = "weightout_coadd_nearest_sapling.fits"  # No iterpolate
                        kargs["ra"] = ra_source
                        kargs["dec"] = dec_source
                        kargs["size"] = (1900,1900)



                        stack_file=swarp.run(best_files_list,outpaht,**kargs)
                        science_dir=base_path+"science/"
                        if not os.path.exists(science_dir):
                            os.makedirs(science_dir)
                        newlocation_best=science_dir+"stack_best.fits"

                        shutil.copy2(stack_file, newlocation_best)

                        best_cat = self.createCatalog(stack_file)

                        newlocation_catdeep_best = science_dir + "catalog_best_deep.fits"
                        shutil.copy2(best_cat["cat_path_deep"], newlocation_catdeep_best)
                        best_cat["cat_path_deep"]=newlocation_catdeep_best

                        newlocation_cat_best = science_dir + "catalog_best.fits"
                        shutil.copy2(best_cat["cat_path"], newlocation_cat_best)
                        best_cat["cat_path"] = newlocation_cat_best

                        newlocation_psf_best = science_dir + "psf_best.fits"
                        shutil.copy2(best_cat["psf_path"], newlocation_psf_best)

                        base_psf_path=path.dirname(best_cat["psf_path"])
                        if path.exists(base_psf_path+"/psfex_out.cat"):
                            newlocation_psf_best_cat = science_dir + "psf_stars.fits"
                            shutil.copy2(base_psf_path+"/psfex_out.cat", newlocation_psf_best_cat)
                            best_cat["psf_stars"] = newlocation_psf_best_cat



                        best_cat["psf_path"] = newlocation_psf_best
                        best_cat["file_path"] = newlocation_best
                        header=FITSFile.header(newlocation_best)[0]
                        best_cat["exp_time"] = header["EXPTIME"]

                        if len(good_seeing) > len(best_list):

                            all_files_list = outpaht + "all_list.lst"
                            if os.path.exists(all_files_list):
                                os.remove(all_files_list)
                            f = open(all_files_list, "w+")
                            tostack=0
                            for fil in good_seeing:
                                if fil[0]<900:
                                    tostack+=1
                                    f.write(fil[1] + "\n")
                            f.close()

                            if tostack>=1:


                                stack_file = swarp.run(all_files_list, outpaht, **kargs)
                                newlocation_all = science_dir + "stack_all.fits"
                                shutil.copy2(stack_file, newlocation_all)


                                stack_file = swarp.run(all_files_list, outpaht, **kargs)

                                all_cat = self.createCatalog(stack_file)

                                newlocation_catdeep_all = science_dir + "catalog_all_deep.fits"
                                shutil.copy2(all_cat["cat_path_deep"], newlocation_catdeep_all)
                                all_cat["cat_path_deep"] = newlocation_catdeep_all

                                newlocation_cat_all = science_dir + "catalog_all.fits"
                                shutil.copy2(all_cat["cat_path"], newlocation_cat_all)
                                all_cat["cat_path"] = newlocation_cat_all

                                newlocation_psf_all = science_dir + "psf_all.fits"
                                shutil.copy2(all_cat["psf_path"], newlocation_psf_all)

                                base_all_psf_path = path.dirname(all_cat["psf_path"])
                                if path.exists(base_all_psf_path + "/psfex_out.cat"):
                                    newlocation_psf_all_cat = science_dir + "psf_stars_all.fits"
                                    shutil.copy2(base_all_psf_path + "/psfex_out.cat", newlocation_psf_all_cat)
                                    all_cat["psf_stars"] = newlocation_psf_all_cat

                                all_cat["psf_path"] = newlocation_psf_all



                                newlocation_all = science_dir + "stack_all.fits"
                                shutil.copy2(stack_file, newlocation_all)

                                all_cat["file_path"] = newlocation_all
                                header_all = FITSFile.header(newlocation_all)[0]
                                all_cat["exp_time"] = header_all["EXPTIME"]
                            else:
                                all_cat = best_cat
                        else:
                            all_cat=best_cat


                        self.mgdb.update({"id": targets},{'$set':{'epoach.'+epoachs+'.stacking':{'best':best_cat,'all':all_cat}}})
                    else:
                        self.mgdb.update({"id": targets},
                                         {'$set': {'epoach.' + epoachs + '.stacking': {"error":"Bad seeing"}}})




    def createCatalog(self,file_path):

        cat_util = CatalogUtil()
        zeropoint, zp_from_cat = cat_util.getZeroPoint(file_path)
        sex_cat = cat_util.generateCatalogNormalized(file_path,zeropoint=zeropoint,good_cat=zp_from_cat)

        respond={"zp":zeropoint,"ref_cat":zp_from_cat[0],"ref_cat_path":zp_from_cat[1],"cat_path":sex_cat[0],"cat_path_deep":sex_cat[1],"seeing_cal":sex_cat[2]["seeing"],"psf_path":sex_cat[2]["psf_file"],"cat_deep_maglimit":sex_cat[2]["cat_deep_maglimit"],"cat_maglimit":sex_cat[2]["cat_maglimit"]}
        return respond

    def stackFiles(self,folder_path):
        LOGUtil.log("liverpool-pipeline Creating staking image for "+ folder_path, folder_path)
        swarp=SwarpUtil()
        return swarp.runSwarp(folder_path,quality="best")


    def createCatalogs(self,folder_path):
        LOGUtil.log("liverpool-pipeline Creating catalog for " + folder_path, folder_path)
        sex = SextractorUtil()
        folder_path = folder_path if folder_path[-1] == "/" else folder_path + "/"
        path_files = folder_path+"reduction/swarp/coadd.fits"
        output_path = folder_path+"reduction/sextractor"
        cata_util=CatalogUtil()
        sex_cat=cata_util.generateCatalogNormalized(path_files,output_path)
        return sex_cat
        #sex.runSextractor(path_files,output_path)



    def generateThumbnail(self,file,ra,dec,size=(200,200)):
        file_fits = fits.open(file)
        file_fits = file_fits[1]
        print(file)
        skycoor = coords.SkyCoord(ra=ra, dec=dec, unit="deg")

        orig_file_wcs = wcs.WCS(header=file_fits.header,fobj=file_fits)





        pixel_file = orig_file_wcs.wcs_world2pix(np.array([[skycoor.ra.deg, skycoor.dec.deg]], np.float_), 1)
        origin_file = (pixel_file[0][0], pixel_file[0][1])
        cutout_file = Cutout2D(file_fits.data, origin_file, size, wcs=orig_file_wcs)
        data_file = cutout_file.data
        min_file = data_file.min()
        max_file = data_file.max()
        mean_file = data_file.mean()
        std_file = data_file.std()
        plt.clf()
        plt.imshow(data_file, cmap=cm.Oranges, vmin=mean_file - std_file, vmax=mean_file + std_file, origin="lower")
        plt.axis('off')
        image_output = file[:-4] + "png"
        plt.savefig(image_output)

        return image_output

    def genereateThumbnails(self,file,ra,dec,size=(200,200)):


        ref_file=file["ref"]
        image_file=file["img"]
        diferents_file=file["dif"]

        ref_file_fits =fits.open(ref_file)[0]
        img_file_fits = fits.open(image_file)[0]
        dif_file_fits = fits.open(diferents_file)[0]

        skycoor = coords.SkyCoord(ra=ra, dec=dec, unit="deg")

        #ref image data
        orig_ref_wcs = wcs.WCS(ref_file_fits.header)

        pixel_ref = orig_ref_wcs.wcs_world2pix(np.array([[skycoor.ra.deg, skycoor.dec.deg]], np.float_), 1)
        origin_ref = (pixel_ref[0][0], pixel_ref[0][1])
        cutout_ref = Cutout2D(ref_file_fits.data, origin_ref, size, wcs=orig_ref_wcs)
        data_ref =cutout_ref.data
        min_ref = data_ref.min()
        max_ref=data_ref.max()
        mean_ref = data_ref.mean()
        std_ref = data_ref.std()

        """
        plt.hist(data_ref.flatten(),50)
        plt.axvline(mean_ref, color='k', linestyle='dashed', linewidth=1)
        plt.axvline(mean_ref + std_ref, color='r', linestyle='dashed', linewidth=1)
        plt.axvline(mean_ref - std_ref, color='r', linestyle='dashed', linewidth=1)
        plt.show()
        """
        #image epoch data
        orig_img_wcs = wcs.WCS(img_file_fits.header)
        pixel_img = orig_img_wcs.wcs_world2pix(np.array([[skycoor.ra.deg, skycoor.dec.deg]], np.float_), 1)
        origin_img = (pixel_img[0][0], pixel_img[0][1])
        cutout_img = Cutout2D(img_file_fits.data, origin_img, size, wcs=orig_img_wcs)
        data_img = cutout_img.data



        min_img = data_img.min()
        max_img = data_img.max()
        mean_img = data_img.mean()
        std_img = data_img.std()

        """
        plt.hist(data_img.flatten(), 50)
        plt.axvline(mean_img, color='k', linestyle='dashed', linewidth=1)
        plt.axvline(mean_img+std_img, color='r', linestyle='dashed', linewidth=1)
        plt.axvline(mean_img - std_img, color='r', linestyle='dashed', linewidth=1)
        plt.show()
        """
        #difference data

        orig_dif_wcs = wcs.WCS(dif_file_fits.header)

        pixel_dif = orig_dif_wcs.wcs_world2pix(np.array([[skycoor.ra.deg, skycoor.dec.deg]], np.float_), 1)
        origin_dif = (pixel_dif[0][0], pixel_dif[0][1])
        cutout_dif = Cutout2D(dif_file_fits.data, origin_dif, size, wcs=orig_dif_wcs)
        data_dif = cutout_dif.data



        min_dif = data_dif.min()
        max_dif = data_dif.max()
        mean_dif = data_dif.mean()
        std_dif = data_dif.std()
        """
        plt.hist(data_dif.flatten(), 50)
        plt.axvline(mean_dif, color='k', linestyle='dashed', linewidth=1)
        plt.axvline(mean_dif + std_dif, color='r', linestyle='dashed', linewidth=1)
        plt.axvline(mean_dif - std_dif, color='r', linestyle='dashed', linewidth=1)

        plt.show()
        """
        fig = plt.figure()

        #show reference

        gs1 = GridSpec(1, 3)
        gs1.update(wspace=0.05)

        ax1=plt.subplot(1, 3, 1)


        plt.imshow(data_ref, cmap=cm.Oranges, vmin=mean_ref-std_ref, vmax=mean_ref+std_ref, origin="lower")
        plt.axis('off')

        ax2=plt.subplot(1,3,2)
        plt.subplots_adjust(wspace=0.5)
        plt.imshow(data_img, cmap=cm.Oranges,vmin=mean_img-std_img, vmax=mean_img+std_img,origin="lower")
        plt.axis('off')


        ax3=plt.subplot(1, 3, 3)
        plt.subplots_adjust(wspace=0.01)
        plt.imshow(data_dif, cmap=cm.Oranges, vmin=mean_dif-std_dif, vmax=mean_dif+std_dif, origin="lower")

        #ax1.coords['ra'].set_axislabel('Right Ascension')
        #ax1.coords['dec'].set_axislabel('Declination')

        #ax2.coords['ra'].set_axislabel('Right Ascension')
        #ax2.coords['dec'].set_axislabel('')
        #ax2.set_yticks([])

        #ax3.coords['ra'].set_axislabel('Right Ascension')
        #ax3.coords['dec'].set_axislabel('')
        #ax3.set_yticks([])
        plt.axis('off')

        #plt.colorbar()
        image_output=image_file[:-4]+"png"
        plt.savefig(image_output)

        print(image_output)
        return image_output


    def createReport(self):
        html = """
                <html>
                  <head>
                    <meta name="pdfkit-page-size" content="Legal"/>
                    <meta name="pdfkit-orientation" content="Landscape"/>
                  </head>
                  <body>
                  """

        self.mgdb.setCollection("liverpool_sn_v3")
        data = self.mgdb.getData({"epoach": {"$exists": "true"},"pyDIA": {"$exists": "true"}})

        hst_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/lensed_sn_hst"
        hst_images=Utils.getFilesList(hst_path)



        for target in data:
            pydia = target["pyDIA"]
            html += "<h2>" + target["id"] + "</h2>"
            print(target["id"])

            hst_image = [f for f in hst_images if target["id"][-9:] in f]
            if len(hst_image)>0:
                hst_image=hst_path+"/"+hst_image[0]
                """
                try:
                    hst_thumbnail = self.generateThumbnail(hst_image, target["ra"], target["dec"])
                except:
                    hst_thumbnail= ""
                """
            epochs = list(pydia)
            epochs.sort()

            #for epoch_key in epochs:

            epoch = pydia
            ref_file = epoch["ref_file"]

            for files in epoch["files"]:

                if "results" in files and "ref" in epoch:
                    ref_file_path = epoch["ref"]
                    image_to_thumbnail = {"ref": ref_file_path, "img": files["file_path"]}



                    for result in files["results"]:

                        if os.path.basename(result)[0:2] == "d_":
                            image_to_thumbnail["dif"] = result
                            break

                    # os.path.basename(files[])


                    if "dif" in image_to_thumbnail:

                        if image_to_thumbnail["ref"] != image_to_thumbnail["img"]:

                            img_url = lvpipeline.genereateThumbnails(image_to_thumbnail, target["ra"], target["dec"],size=(100,100))
                            html += "<div ><h3>Epoch: " + os.path.basename(files["file_path"])[:10] + "</h3>"

                            html += "<img src='" + img_url + "' width='1600'>"

                            html += "<p> reference: seeing "+str(ref_file["seeing"])[:4]+ " exptime "+str(ref_file["exptime"]/60)+"min</br>image: seeing"+ str(files["seeing"])[:4] +" exptime "+ str(files["exptime"]/60)+"min</br>Difference:</br>HST image:</br>load " + image_to_thumbnail["ref"] + "</br>load " + image_to_thumbnail[
                                "img"] + "<br>load " + image_to_thumbnail["dif"] + "<br>load "+hst_image+"</p></div>"

                            if "refim" in epoch.keys():
                                refimageslst = epoch["refim"]
                                file_refimg = open(refimageslst, "r")
                                cont_ref = file_refimg.read()
                                print(cont_ref)

                                html +="Ref list Pydia: name fwhm sky signal<br>"
                                html += str(cont_ref).replace("\n","<br>")

                    else:
                        print("any results")
                else:
                    html += "<h3>Pending for pydia results " + epoch["epoach"] + "</h3></div>"
                # skycoor = coords.SkyCoord(ra=target["ra"], dec=target["dec"], unit="deg")
                # ext = fits.open(ref_file_path)[0]
                # orig_wcs = wcs.WCS(ext.header)
                #
                # pixel = orig_wcs.wcs_world2pix(np.array([[skycoor.ra.deg, skycoor.dec.deg]], np.float_), 1)
                #
                # ImageUtils.imcopy(ref_file_path, ref_file_path + "_crop.fits",
                #                   (pixel[0][0], pixel[0][1]), (100, 100))
                # print("ok")

        html += """
                  </body>
              </html>
            """
        f = open("/Users/cjimenez/Documents/PHD/data/liverpool_lens/report.html", "w+")
        f.write(html)
        f.close()
        pdfkit.from_string(html, "/Users/cjimenez/Documents/PHD/data/liverpool_lens/report.pdf")

    def getDetection(self,target):
        cat_util = CatalogUtil()

        id = target["id"]
        ra = target["ra"]
        dec = target["dec"]
        for epoach in target["epoach"]:
            try:
                print(id, epoach)
                stack = target["epoach"][epoach]["stacking"]
                if "error" not in stack.keys():
                    query = cat_util.searchDetectionbyTarget(id,ra,dec,stack,epoach)
                    self.mgdb.update({"id":id},query)
                else:
                    print("Observation error")
            except KeyError:
                    print("Stacking not found ",id,epoach)


    ####WARNING####
    #if you run this function, you will be deleate all pydia results from your local files
    def clean(self):
        mgdb = MongodbManager()
        mgdb.setCollection("liverpool_sn_v3")
        data = mgdb.getData({"epoach": {"$exists": "true"}})
        #remove all pydia results
        for target in data:
            try:
                path_to_delete=self.project_path+target["id"]+"/pydia"
                shutil.rmtree(path_to_delete)
                print("deleted "+self.project_path+target["id"]+"/pydia")
            except FileNotFoundError:
                print("not pydia files" +self.project_path+target["id"]+"/pydia")




if __name__ == "__main__":


    img_path="/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181110_1"
    #img_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181003_1"
    img_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181115_1"
    img_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181119_1"
    img_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181127_1"
    img_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/alldatareduction"
    #img_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181201_1"



    lvpipeline = LiverpollSnPipeline()
    mgdb = MongodbManager()
    mgdb.setCollection("liverpool_sn_v3")

    allFiles=True

    files = lvpipeline.getData(img_path,all=allFiles,force=False)


    #files={}
    #files["BELLSJ1349+3612"] = file["BELLSJ1349+3612"]
    #files["BELLSJ1601+2138"] = file["BELLSJ1601+2138"]
    #files["BELLSJ1337+3620"] = file["BELLSJ1337+3620"]
    """
    files={}

    files["SLACSJ1432+6317"]= file["SLACSJ1432+6317"]
    files["BELLSJ1215+0047"] = file["BELLSJ1215+0047"]
    files["BELLSJ0830+5116"] = file["BELLSJ0830+5116"]
    files["BELLSJ1611+1705"] = file["BELLSJ1611+1705"]
    files["SLACSJ1621+3931"] = file["SLACSJ1621+3931"]
    files["SLACSJ1525+3327"] = file["SLACSJ1525+3327"]
    files["SLACSJ1538+5817"] = file["SLACSJ1538+5817"]


    files["SLACSJ1430+4105"]= file["SLACSJ1430+4105"]
    files["SLACSJ1416+5136"] = file["SLACSJ1416+5136"]
    files["BELLSJ1159-0007"] = file["BELLSJ1159-0007"]
    files["S4TMJ1543+2202"] = file["S4TMJ1543+2202"]
    files["S4TMJ1550+2020"] = file["S4TMJ1550+2020"]
    files["BELLSJ1352+3216"] = file["BELLSJ1352+3216"]
    files["BELLSJ1601+2138"] = file["BELLSJ1601+2138"]
    files["BELLSJ1337+3620"] = file["BELLSJ1337+3620"]
    files["BELLSJ1349+3612"] = file["BELLSJ1349+3612"]
    files["BELLSJ1318-0104"] = file["BELLSJ1318-0104"]
    files["BELLSJ1221+3806"] = file["BELLSJ1221+3806"]
    """

    # files["SLACSJ1112+0826"] = file["SLACSJ1112+0826"]
    #files["SLACSJ1103+5322"] = file["SLACSJ1103+5322"]
    # files["SLACSJ1416+5136"] = file["SLACSJ1416+5136"]
    # files["S4TMJ1550+2020"] = file["S4TMJ1550+2020"]
    # files["BELLSJ1611+1705"] = file["BELLSJ1611+1705"]
    # files["SLACSJ1621+3931"] = file["SLACSJ1621+3931"]
    # files["S4TMJ1051+4439"] = file["S4TMJ1051+4439"]
    #files["BELLSJ0801+4727"] = file["BELLSJ0801+4727"]

    # files["BELLSJ1215+0047"] = file["BELLSJ1215+0047"]
    # files["BELLSJ1352+3216"] = file["BELLSJ1352+3216"]
    # files["SLACSJ1430+4105"] = file["SLACSJ1430+4105"]

    # files["S4TMJ1048+1313"] = file["S4TMJ1048+1313"]
    # files["SLACSJ1143-0144"] = file["SLACSJ1143-0144"]
    # files["S4TMJ1051+4439"] = file["S4TMJ1051+4439"]
    # files["BELLSJ1542+1629"] = file["BELLSJ1542+1629"]
    # files["BELLSJ1631+1854"] = file["BELLSJ1631+1854"]
    # files["BELLSJ1601+2138"] = file["BELLSJ1601+2138"]
    #files["SLACSJ1032+5322"] = file["SLACSJ1032+5322"]
    #files["S4TMJ1051+4439"] = file["S4TMJ1051+4439"]
    #files["BELLSJ1349+3612"] = file["BELLSJ1349+3612"]

    #files["BELLSJ1337+3620"] = file["BELLSJ1337+3620"]
    allFiles=False

    runProcess=["stacking","pydia"] #"copyresults","stacking","models","detection","pydia","report","list_stack","clean"]

    if allFiles:
        filters = {"epoach": {"$exists": "true"}}
    else:
        ids=files.keys()
        filters = {"epoach": {"$exists": "true"},"$or":[]}
        for id in ids:
            filters["$or"].append({"id":id})


    data = mgdb.getData(filters)

    if "list_stack" in runProcess:
        for target in data:
            print(target["id"])
            for epoch in target["epoach"]:
                if "error" in target["epoach"][epoch]["stacking"]:
                    print(epoch+" bad observation")
                else:
                    file_path=target["epoach"][epoch]["stacking"]["best"]["file_path"]
                    seeing_path = target["epoach"][epoch]["stacking"]["best"]["seeing_cal"]

                    #print("load "+file_path + " "+str(seeing_path))
                    print("load " + file_path)


    if "copyresults" in runProcess:
        pyDIA = PYDIAPipeline()
        for target in data:
            pyDIA.checkResults(target)
            pyDIA.checkResidues(target)

    if "coordinatefile" in runProcess:
        for target in data:
            lvpipeline.coordinateFile(target)


    if "stacking" in runProcess:

        lvpipeline.stackEpoach(files,force=True)

        for target in data:
            lvpipeline.stackAllObsByTarget(target["id"])
            lvpipeline.getDetection(target)


    if "detection" in runProcess:
        for target in data:
            lvpipeline.getDetection(target)

    if "pydia" in runProcess:
        pyDIA = PYDIAPipeline()
        for target in data:
            pyDIA.run(target)

    if "models" in runProcess:
        lvpipeline.generateModels()

    if "report" in runProcess:
        lvpipeline.createReport()
        Generatereports.galfitReport()

    if "clean" in runProcess:
        lvpipeline.clean()










