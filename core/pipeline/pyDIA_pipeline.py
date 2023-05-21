
from frastro.core.database.mongodb.mongodb_manager import MongodbManager
from frastro.core.fits.FITSFIle import FITSFile
from frastro.core.utils.image_util import ImageUtils
import numpy as np
from matplotlib import pyplot as plt
import os
from frastro import Config, StiltsUtil
import shutil
import math
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
from frastro.external.sextractor_util import SextractorUtil
from frastro.external.psfex_util import PSFexUtil
from frastro.external.stilts_util import StiltsUtil
from frastro.core.data.catalog.catalog_util import CatalogUtil

class PYDIAPipeline():

    command="python run_pyDIA.py --n_parallel=20 --gain=1.9 --readnoise=5.0 --name_pattern=*.fits --pdeg=2 --sdeg=2 --bdeg=2 --datekey=None  --loc_data=liverpool_lens/{0}/images_stack --loc_output=liverpool_lens/{0}/output_{1} --use_GPU=True --min_ref_images={2} --psf_fit_radius=3.0 --iterations=3" #--star_file=liverpool_lens/{0}/stars.txt --ref_image_list=ref_file.txt

    def __init__(self):
        self.db = MongodbManager()
        collection = "liverpool_sn_v3"
        self.db.setCollection(collection)
        self.project_path = Config.getPath("liverpool_sn")
        self.base_server="/Volumes/scratch/cjimenez/data/"
        self.base_server_path = "/Volumes/scratch/cjimenez/data/liverpool_lens/{0}/"
        self.server_path="/Volumes/scratch/cjimenez/data/liverpool_lens/{0}/images_stack/"
        if os.path.exists(self.project_path + "run_pydia.sh"):
            os.remove(self.project_path + "run_pydia.sh")



    def removeSky(self,image):

        output_file = image["epoach"] + "_" + str(image["seeing"])[:4] + "_single_nosky" + ".fits"

        file=image["file"]
        base_path = os.path.dirname(file)+"/"

        image_row=FITSFile.open(file)
        data=image_row[0].data
        header=image_row[0].header
        mean = data.mean()
        new_data= data #-mean

        saveFile=base_path+output_file

        FITSFile.saveFile(saveFile,new_data,header)

        return saveFile

    def removeBorders(self,image):
        file=image["file"]

        output_file = image["epoach"] + "_"+str(image["seeing"])[:4]+"_stack"+".fits"
        image=FITSFile.open(file)
        data=image[0].data.shape
        print("data",data)
        width = data[0]-300
        height = data[1] - 300
        base_path = os.path.dirname(file)
        saveFile = base_path +"/"+ output_file
        ImageUtils.imcopy(file, saveFile, origin=(data[0]/2,data[1]/2),crop=(width,height),addcounts=100)
        return saveFile



    def cropAllFiles(self,file,ra,dec,size=(100,100),addcounts=0):
        pixel = ImageUtils.getPixelFromCoordinate(file, ra, dec)
        ImageUtils.imcopy(file, file, (pixel[0][0], pixel[0][1]), size,addcounts=addcounts)


    def copyFilesToServer(self,files_list):


        if os.path.exists("/Volumes/scratch/"):

            for files in files_list:
                path = self.server_path.format(files["target"])
                print(path)
                if not os.path.exists(path):
                    os.makedirs(path,mode=0o777)
                name_file=os.path.basename(files["file_path"])
                print(files["file_path"])
                shutil.copy2(files["file_path"],path+name_file)
                ref_file = files["file_path"][:files["file_path"].find("obs")]

                shutil.copy2(ref_file+"ref_file.lst",self.base_server_path.format(files["target"])+"ref_file.txt")
#                shutil.copy2(ref_file + "stars.txt", self.base_server_path.format(files["target"]) + "stars.txt")


        else:
            print("please connect to server!")

    def checkResults(self,target):
        target_id = target["id"]
        #results={}
        if "pyDIA" in target:
            pydia=target["pyDIA"]
            print(target_id)
            epoach=pydia["epoach"]
            #results[target_id]={}

            #for epoach in pydia:

            #print(epoach)
            #status = "pending"

            epoach_data = pydia
            epoach_server_path=self.base_server_path.format(target_id)+"output_"+epoach

            #if "status" not in epoach_data or epoach_data["status"] == "pending":

            if os.path.exists(epoach_server_path):
                os.chmod(epoach_server_path + "/", mode=0o777)


                onlyfiles = [f for f in os.listdir(epoach_server_path) if os.path.isfile(epoach_server_path + "/" + f)]
                filters = ["d_", "m_", "n_","r_"]
                files = ["psf.fits", "ref.fits", "seeing","ref.images"]
                save_path = self.project_path + target_id + "/pydia/" + epoach + "/results/"
                if not os.path.exists(save_path):

                    os.makedirs(save_path,mode=0o777)
                for file in onlyfiles:
                    if file in files or file[:2] in filters:
                        print("try to copy "+file)
                        try:

                            shutil.copy2(epoach_server_path+"/"+file,save_path+file)
                            print("copied",save_path+file)
                            #results[target_id][epoach].append(save_path+file)
                            if file in files:
                                if len(file)==8:
                                    epoach_data[file[:3]]=save_path+file
                                    if file[:3] == "ref":
                                        ref_file = FITSFile.open(epoach_data["ref_file"]["file_path"])[0]
                                        result_file = FITSFile.open(save_path + file)[0]
                                        FITSFile.saveFile(save_path + file, npa_data=result_file.data,
                                                          header=ref_file.header,overwrite=True)

                                else:
                                    idx=str(file[:6]).replace(".","")
                                    epoach_data[idx] = save_path + file



                            else:
                                for files_item in epoach_data["files"]:
                                    name=os.path.basename(files_item["file_path"])
                                    print(name,file[2:])
                                    if name==file[2:]:
                                        ref_file = FITSFile.open(files_item["file_path"])[0]
                                        result_file = FITSFile.open(save_path + file)[0]
                                        FITSFile.saveFile(save_path + file,npa_data=result_file.data, header=ref_file.header)

                                        if "results" in files_item:
                                            if save_path + file not in files_item["results"]:
                                                files_item["results"].append(save_path + file)
                                        else:
                                            files_item["results"]=[save_path+file]
                                        status = "done"

                        except PermissionError:
                            print("error copy files "+epoach_server_path+"/"+file,save_path+file)
                            f=open(self.project_path+"errors.sh","a+")
                            f.write("chmod 777 "+target_id+"/output_"+epoach+"/* \n")
                            f.close()
            else:
                f = open(self.project_path + "log_server.log", "a+")
                f.write("no found "+target_id+" "+epoach+" "+epoach_server_path+"\n")
                f.close()
            #else:
            #    print("cant to copy "+epoach_data)
            self.db.update({"id":target_id},{"$set":{"pyDIA":target["pyDIA"]}})

    def generateStarsFile(self,img_file_path,stars_path,deep_cat_path,pydia_path):

        st = StiltsUtil()

        file_stars = fits.open(stars_path, mode="update")

        if len(file_stars) >= 3:
            del file_stars[1]
            file_stars.flush()
        file_stars.close()

        out_cat_stars = pydia_path + "stars.fits"

        pos_cat = "values1=X_IMAGE Y_IMAGE"
        pos_reference = "values2=ra dec"
        separation = '1'  # in arcsec
        ouput_type = "ofmt=fits-basic"

        pos_reference = "values2=X_IMAGE Y_IMAGE"

        params = [pos_cat, pos_reference, ouput_type]

        r = st.crossMatchTwoCatalogsXY(stars_path, deep_cat_path, out_cat_stars, separation, list_str_params=params)
        print(r)

        stars_list = fits.open(out_cat_stars)

        data = stars_list[1].data
        corx = data["ALPHA_J2000"]
        cory = data["DELTA_J2000"]
        mag = data["MAG_AUTO"]

        file_check_stars_coordinate = img_file_path

        ref_file_stars = fits.open(file_check_stars_coordinate)
        header = ref_file_stars[0].header
        img_w = header["NAXIS1"]
        img_h = header["NAXIS2"]
        wcs_img = wcs.WCS(header)

        w = open(pydia_path + "stars.txt", "w+")
        w.write("#id,x_pix,y_pix,mag\n")

        index = 0
        for idx in range(len(corx)):
            pixcrd2 = wcs_img.wcs_world2pix([[corx[idx], cory[idx]]], 1)
            if pixcrd2[0][0] > 0 and pixcrd2[0][0] < img_w and pixcrd2[0][1] > 0 and pixcrd2[0][1] < img_h:
                index += 1
                w.write(str(index) + " " + str(pixcrd2[0][0]) + " " + str(pixcrd2[0][1]) + " " + str(mag[idx]) + "\n")
        w.close()

    def run(self,target):

        files_list=[]
        epoach_list=[]
        print(target["id"])
        files_list_to_copy=[]
        for epoach_idx in target["epoach"]:
            epoach_list.append(epoach_idx)
            print(epoach_idx)
            epoach=target["epoach"][epoach_idx]
            if "stacking" in epoach:
                if "error" in epoach["stacking"]:
                    if len(epoach["obs"]) > 1:
                        print("bad observations")
                        continue
                    else:
                        seeing=epoach["obs"][0]["seeing_cal"]
                        if seeing < 999:
                            file = epoach["obs"][0]["file_path"]
                            new_file=self.removeSky({"seeing":seeing,"file":file,"epoach":epoach_idx})
                            files_list.append([seeing,new_file])
                else:
                    best=epoach["stacking"]["best"]
                    seeing=best["seeing_cal"]

                    if seeing<999:
                        crop_image=self.removeBorders({"seeing":seeing,"file":best["file_path"],"epoach":epoach_idx})
                        n_point_source = best["cat_maglimit"]["n_point_sources"]
                        deep = best["cat_maglimit"]["point_deep"]
                        files_list.append([seeing,crop_image,n_point_source,deep])


        if len(files_list) > 1:
            files_list.sort(key=lambda x: x[0])
            epoach_list.sort(reverse=True)
            reffile=0
            if len(files_list) > 1:
                if files_list[1][2] > files_list[0][2]:
                    reffile = 1



            pydia_path = self.project_path+target["id"]+"/pydia/"+epoach_list[0]+"/"
            if not os.path.exists(pydia_path+"obs/"):
                os.makedirs(pydia_path+"obs/")

            seeing_images=[]
            min_seeing=9999
            pydia_files=[]
            for file in files_list:
                files_pydia_path = pydia_path +"obs/"+ file[1][file[1].rfind("/")+1:]
                shutil.move(file[1],files_pydia_path)
                header=FITSFile.header(files_pydia_path)[0]
                seeing_images.append(file[0])
                if file[0]<min_seeing:
                    min_seeing = file[0]
                pydia_files.append({"seeing": file[0], "file_path": files_pydia_path,"exptime":header["EXPTIME"]})
                files_list_to_copy.append({"target":target["id"],"file_path":files_pydia_path})
            pydia_ref =  pydia_files[reffile]["file_path"]

            f=open(pydia_path+"ref_file.lst","w+")
            f.write(os.path.basename(pydia_ref))
            print("ref ",pydia_ref, pydia_files[reffile]["seeing"])
            f.close()



            good_ref = np.where(np.array(seeing_images) <=min_seeing+(min_seeing*0.3),1,0)

            num_image_ref=sum(good_ref)

            print(num_image_ref)
            ref_images=1
            if (num_image_ref/3) > 1:
                ref_images = 3

            print(ref_images)
            command=self.command.format(target["id"],epoach_list[0],ref_images)

            # genereate stars cat to pydia from psexf stars selection
            # out_stars = best["psf_stars"]
            # deep_cat = best["cat_path_deep"]
            # img_file_path = files_list_to_copy[reffile]["file_path"]
            # self.generaeteStarsFile(img_file_path,out_stars,deep_cat,pydia_path)
            #command+=" --star_file_has_magnitudes=True --star_file_is_one_based=True"

            #add commad to run single target

            f = open(pydia_path + "run_pydia.sh", "w+")
            f.write(command)
            f.close()



            items={"epoach":epoach_list[0],"files":pydia_files,"ref_file":pydia_files[0],"command":command}
            self.db.update({"id":target["id"]},{"$set":{"pyDIA":items}})




            for file_to_crop in pydia_files:
                self.cropAllFiles(file_to_crop["file_path"],target["ra"],target["dec"],(1400,1400),addcounts=100)

            if len(files_list_to_copy) == 2:
                copy_file = files_list_to_copy[1]["file_path"]
                target_copy=files_list_to_copy[1]["target"]

                if pydia_ref in copy_file:
                    copy_file = files_list_to_copy[0]["file_path"]
                    target_copy = files_list_to_copy[0]["target"]

                copy_file_to=copy_file[:-5]+"_copy.fits"
                shutil.copy2(copy_file,copy_file_to)
                files_list_to_copy.append({"target": target_copy, "file_path": copy_file_to})



            self.copyFilesToServer(files_list_to_copy)

            log= open(self.project_path+"pydia_log.log","a+")
            log.write("================================\n")
            log.write(target["id"]+" "+epoach_list[0]+"\n")
            log.write(command+"\n")
            log.close()

            # add commad to run all targets

            log = open(self.project_path + "run_pydia.sh", "a+")
            log.write(command + "\n")
            log.close()

            shutil.copy2(self.project_path + "run_pydia.sh", self.base_server + "bash_pydia.sh")

    def checkResidues(self,target):
        sextractor =SextractorUtil()
        cat_utils = CatalogUtil()
        stilts = StiltsUtil()
        target_id = target["id"]
        print(target_id)
        dbmongo = MongodbManager()
        collection = "liverpool_sn_v3"
        dbmongo.setCollection(collection)
        # results={}
        if "pyDIA" in target:
            pydia = target["pyDIA"]

            if "ref_cat" not in pydia:
                if "ref" not in pydia:
                    return
                ref = pydia["ref"]
                ref_catlogs = cat_utils.generateCatalogNormalized(ref)
                ref_catlog=ref_catlogs[1]
                header_psf=fits.open(ref)[0].header
                exptime=header_psf["EXPTIME"]
                ref_catalog_tosave = {"cat_path": ref_catlogs[0], "cat_path_deep": ref_catlogs[1],"exptime":exptime }
                input_db = {**ref_catalog_tosave, **ref_catlogs[2]}
                input_db["cat_ref"] = {"archive": input_db["cat_ref"][0], "file_path": input_db["cat_ref"][1]}
                dbmongo.update({"id": target_id}, {"$set": {"pyDIA.ref_cat": input_db}})
            else:
                ref_catlog = pydia["ref_cat"]["cat_path_deep"]


            for indx,files in enumerate(pydia["files"]):
                if "results" in files:

                    for dif in files["results"]:

                        difname=os.path.basename(dif)
                        if difname[0:1]=="d":
                            difference_img=dif
                            file_path=files["file_path"]
                            epoch=os.path.basename(file_path)[0:10]
                            zp=target["epoach"][epoch]["stacking"]["best"]["zp"]
                            psf=target["epoach"][epoch]["stacking"]["best"]["psf_path"]
                            seeing = target["epoach"][epoch]["stacking"]["best"]["seeing_cal"]
                            cat_ref=target["epoach"][epoch]["stacking"]["best"]["cat_path_deep"]

                            kargs={}
                            kargs["-MAG_ZEROPOINT"] = str(zp)
                            kargs["-SEEING_FWHM"] = str(seeing)
                            kargs["-PSF_NAME"] = str(psf)
                            kargs["-DETECT_THRESH"] = str(1.5)
                            kargs["-ANALYSIS_THRESH"] = str(13)
                            kargs["-CATALOG_NAME"] = os.path.basename(difference_img)[2:12]+"_dualmode_cat.fits"
                            images=[file_path,difference_img]
                            outpath=os.path.dirname(difference_img)+"/"
                            sex_result=sextractor.runSextractor(images,output_path=outpath,dualmode=True,calc_psf=False,**kargs)

                            kargs["-DETECT_THRESH"] = str(2.5)
                            kargs["-ANALYSIS_THRESH"] = str(13)
                            kargs["-CATALOG_NAME"] = os.path.basename(difference_img)[2:12] + "_residues_cat.fits"
                            sex_residues = sextractor.runSextractor(difference_img, output_path=outpath, dualmode=False,
                                                                  calc_psf=False, **kargs)

                            print("pyDIA.files."+str(indx)+"residues_cat",sex_result)


                            pos_cat = "values1=ALPHA_J2000 DELTA_J2000"
                            pos_reference = "values2=ALPHA_J2000 DELTA_J2000"
                            separation = '1'  # in arcsec
                            ouput_type = "ofmt=fits-basic"
                            join="join=1not2"



                            #search detection around the target
                            residues_file=fits.open(sex_residues)
                            data = residues_file[1].data
                            ra_detections = data["ALPHA_J2000"]
                            dec_detections = data["DELTA_J2000"]
                            catalog = SkyCoord(ra=ra_detections * u.degree, dec=dec_detections * u.degree)
                            target_position=SkyCoord(ra=np.array([target["ra"]]) * u.degree, dec=np.array(target["dec"]) * u.degree)

                            #search around 10 arcsec
                            idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(target_position, 0.00277778 * u.deg)
                            print(d2d,d2d.to("arcsec"))
                            detections = {"detections":0}
                            if idxcatalog.size>0:
                                detections={"detections":idxcatalog.size,"ra":data[idxcatalog]["ALPHA_J2000"][0],"dec":data[idxcatalog]["DELTA_J2000"][0],"mag":float(data[idxcatalog]["MAG_AUTO"][0]),"distance":str(d2d[0].to("arcsec"))}

                            params = [pos_cat, pos_reference, ouput_type,join]



                            stilts = StiltsUtil()
                            stilts.crossMatchTwoCatalogs(sex_residues, ref_catlog, outpath+os.path.basename(difference_img)[2:12]+"_residual_detections.fits", separation, params)

                            dbmongo.update({"id":target_id},{ "$set": {"pyDIA.files."+str(indx)+".residuesforce_cat": sex_result,"pyDIA.files."+str(indx)+".residues_cat": sex_residues,"pyDIA.files."+str(indx)+".residues_detections":outpath+os.path.basename(difference_img)[2:12]+"_residual_detections.fits","pyDIA.files."+str(indx)+".detections":detections}})
                            break



if __name__ == "__main__":

    db= MongodbManager()
    pyDIA=PYDIAPipeline()
    collection = "liverpool_sn_v3"
    db.setCollection(collection)
    targets=db.getData({"id":"S4TMJ1051+4439","epoach":{"$exists":"true"}})
    path="/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/astrometry/good/"
    files=[{"seeing": 1, "file":path+"swarp_2018-11-16_1.fits", "epoach":"2018-11-16_1"},
        {"seeing": 1, "file":path+"swarp_2019-01-10_1.fits", "epoach":"2019-01-10_1"},
        {"seeing": 1, "file":path+"swarp_2019-04-04_1.fits", "epoach":"2019-04-04_1"},
        {"seeing": 1, "file":path+"swarp_2018-11-20_1.fits", "epoach":"2018-11-20_1"},
        {"seeing": 1, "file":path+"swarp_2019-02-03_1.fits", "epoach":"2019-02-03_1"},
        {"seeing": 1, "file":path+"swarp_2019-04-23_1.fits", "epoach":"2019-04-23_1"},
        {"seeing": 1, "file":path+"swarp_2018-12-17_1.fits", "epoach":"2018-12-17_1"},
        {"seeing": 1, "file":path+"swarp_2019-02-26_1.fits", "epoach":"2019-02-26_1"},
        {"seeing": 1, "file":path+"swarp_2019-04-24_1.fits", "epoach":"2019-04-24_1"}]

    for file in files:
        pyDIA.removeBorders(file)

    #for target in targets:
        #pyDIA.run(target)
        #pyDIA.checkResults(target)
        #pyDIA.checkResidues(target)
        """
        if "pyDIA" in target:
            pydia = target["pyDIA"]
            if "ref_cat" in pydia:
                ref = pydia["ref"]
                header_psf = fits.open(ref)[0].header
                exptime = header_psf["EXPTIME"]
                db.update({"id":target["id"]},{"$set":{"pyDIA.ref_cat.exptime":int(exptime)}})

        #pyDIA.checkResidues(target)
        """






