#!/bin/bash
import os
import shutil

from astropy.io import fits

from frastro import ImageUtils, Config, Utils, CatalogUtil, SextractorUtil, MongodbManager
from frastro.core.pipeline.astrometry.gaia_astrometry import GaiaAstrometry
from astropy.wcs._wcs import InvalidTransformError
import pwd
import getpass
import numpy as np
#Do astrometry from liris files using astrometry.net api and gaia catalogs
#if astrometry.net can't find a solution do manual usin gaia program and astrometry from fix star base on gaia catalog
from frastro.external.swarp_util import SwarpUtil
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
sns.set(style="whitegrid")

class LirisAstrometry():

    #from file path, do astrometry

    def __init__(self,file_path=""):
        if file_path!="":
            self.init(file_path)
        self.db = MongodbManager()
        collection = "unlens_liris"
        self.db.setCollection(collection)

        # plt.ion()
        # print(matplotlib.get_backend())


    def init(self,file_path):
        self.__folder_path = os.path.dirname(file_path)+"/"
        self.__image_path = os.path.basename(file_path)
        self.__file_path = file_path
        self.gc = GaiaAstrometry()

    def fixHeaderError(self,file_path=""):
        if file_path != "":
            file= file_path
        else:
            file = self.__file_path

        with fits.open(file, mode='update') as hdul:
            hdr = hdul[0].header

            if "CUNIT1" in hdr:
                hdr["CUNIT1"] = 'deg'

            if "CUNIT2" in hdr:
                hdr["CUNIT2"] = 'deg'

            if "LIRLAMP1" in hdr:
                del hdr["LIRLAMP1"]

            if "LIRLAMP2" in hdr:
                del hdr["LIRLAMP2"]

            if "OBJECT" in hdr:
                if str(hdr["OBJECT"]).find(":") > -1:
                    hdr["OBJECT"] = str(hdr["OBJECT"])[str(hdr["OBJECT"]).find(":")+1:]

            hdul.flush()
            hdul.close()
            print("update header",hdr["OBJECT"])


    def doAstrometry(self):
        self.fixHeaderError()
        catalog_path = self.__createGaiaCatalog()
        self.gc.generateIndex(catalog_path)
        file_sh=self.__createSHCommands()
        return file_sh
        #self.runShCommand(file_sh)


    # check or donwload gaia catalog from center ra dec image in a radio of 2 degree
    def __createGaiaCatalog(self):
        gaia_file_path = self.__folder_path + "gaia/" + self.__image_path[0:-5] + "_gaia2deg.fits"
        if not os.path.exists(self.__folder_path+"gaia/"):


            #the ra dec dont be perfect due the cover range will be larger thant the fild of view of image
            skycoor = self.__getRaDec()

            #find and download gaia catalog from tapquery services and save in local fits in a especific path, by default this only return a astropy table
            #the file only will be have ra,dec, mag columns for astrometry.net if you need a full catalog use gaia_archive_cp.py getCatalog instead gaia_astrometry

            catalog = self.gc.getGaiaCatalog(skycoor.ra.deg, skycoor.dec.deg,save=True,output=gaia_file_path)
        return gaia_file_path

    #find ra dec from center of image, !warning if the wcs is wrong the ra dec will be wrong
    def __getRaDec(self):
        try:
            skycoor = ImageUtils.getCenterCoordinate(self.__file_path)
        except InvalidTransformError:
            skycoor = self.__getRaDec()
        return skycoor


    def __createSHCommands(self):

        sh_file_path = self.__file_path[0:-5]+".sh"
        sh_file = open(sh_file_path, "w+")

        solved_command = Config.getCommand("astrometry_solved_file")
        sip_tav_path_command = Config.getCommand("run_sip_tav")
        conf_backend = self.__folder_path + "gaia/backend.cfg"
        low_scale = 0.20
        high_scale = 0.30
        command = solved_command.format(low_scale, high_scale, self.__file_path, conf_backend)+"\n"
        command += Config.getCommand("sip_tav_path")+"\n"
        command += sip_tav_path_command.format(self.__file_path[0:-5]+".new",self.__file_path[0:-5]+"_as.fits")+"\n"

        sh_file.write(command)
        sh_file.close()

        return sh_file_path

    def runShCommand(self,sh_file):
        r = Utils.runCommand("source ~/.bash_profile")
        os.chdir(self.__folder_path)
        r=Utils.runCommand(". "+sh_file)
        print(r)


    def listUserCM(self):

        user_id=getpass.getuser()
        pw_record = pwd.getpwnam(user_id)
        user_name = pw_record.pw_name
        user_home_dir = pw_record.pw_dir
        user_uid = pw_record.pw_uid
        user_gid = pw_record.pw_gid

        print(user_name)


    def evaluateAsytrometry(self,file_path):

        file_path = file_path[0:-5]+"_as.fits"
        cat_outpath = file_path[0:-5]+"_cat.fits"
        match_output = file_path[0:-5]+"_match.fits"
        path_folder=os.path.dirname(file_path)
        gaiapath = path_folder+"/gaia/"

        gaiapath_file = ""
        onlyfiles = [f for f in os.listdir(gaiapath) if os.path.isfile(gaiapath + "/" + f)]
        for file in onlyfiles:
            if os.path.splitext(file)[1]== ".fits" or os.path.splitext(file)[1]== ".fit":
                if file.find("_gaia2deg")!= -1:
                    gaiapath_file = gaiapath+file


        if gaiapath_file != "":
            ca=CatalogUtil()
            cat=ca.generateCatalogFromFitsPath(file_path)
            r=ca.crossMathTables(cat,gaiapath_file,match_output,ref_table2="values2=RA DEC",separation="10")
            cat_math_file=fits.open(r)
            data = cat_math_file[1].data
            stats_file = open(file_path[0:-5] + "_stats.txt", "w+")
            if len(data) > 0:
                mean=data["separation"].mean()
                std = data["separation"].std()
                average = np.average(data["separation"])
                median =  np.median(data["separation"])
                min=data["separation"].min()
                max = data["separation"].max()
                stats_file.write("mean: "+str(mean)+"\n")
                stats_file.write("mean pixel: "+str(mean/0.25)+"\n")
                stats_file.write("std: " + str(std) + "\n")
                stats_file.write("std pixel: " + str(std/0.25) + "\n")
                stats_file.write("average: " + str(average) + "\n")
                stats_file.write("median: " + str(median) + "\n")
                stats_file.write("min: " + str(min) + "\n")
                stats_file.write("max: " + str(max) + "\n")
                print(mean,std,str(mean/0.25),str(std/0.25),average,median,min,max,file_path)
            else:
                stats_file.write("any match found \n")
                print("any match found",file_path)
            stats_file.close()


            return r
        else:
            print("gaia catalog no found",file_path)

    def scamp(self,file_path):
        SextractorUtil().scampCatalog(file_path)


    def validateAsFile(self,file_path):
        errors="/Users/cjimenez/Documents/PHD/data/unlens_liris/erros.lst"
        errors=open(errors,"w+")
        if not os.path.exists(file_path[0:-5]+"_as.fits"):
            errors.write(file_path+ "\n")
            print("error",file_path)
        errors.close()

    def saveOnDataBase(self,query):



        items = self.db.getData({"id": query["id"]})
        if len(items) > 0:
            query_str = {"$set": {"band."+query["band"]: query["obs"]}}
            self.db.update(query["id"], query_str)
        else:
            self.db.saveData({"id":query["id"],"band":{query["band"]: query["obs"]}})

    def checkZeroPoints(self):
        base_path="/Users/cjimenez/Documents/PHD/data/unlens_liris/report"
        dates=["2019-06-14","2019-06-15","2019-06-16","2019-06-17","2019-06-18","2019-06-19"]
        filters=["mag_h","mag_j","mag_ks"]
        for date in dates:

            for filter in filters:
                labels = []
                chart = []
                x=[]
                y=[]
                query={"band."+filter+"."+date:{"$exists":1}}
                collection = "unlens_liris"
                self.db.setCollection(collection)
                data=self.db.getData(filter=query)
                items_in_cat = {}
                for items in data:
                    if items["id"] in labels:
                        idx=labels.index(items["id"])
                    else:
                        labels.append(items["id"])
                        idx = labels.index(items["id"])

                    obj=items["band"][filter][date]
                    photometry=obj["photometry"]
                    for photo in photometry:


                        zp = np.float32(photo["zp_mean"])
                        std = float(photo["zp_std"])
                        n_stars=int(photo["n_stars"])
                        cat=photo["catalog"]
                        x.append(idx)
                        y.append(zp)
                        chart.append([idx,obj["seeing"],obj["airmass"],zp,std,n_stars,cat])
                        if cat in items_in_cat.keys():
                            items_in_cat[cat].append([zp,n_stars])
                        else:
                            items_in_cat[cat]=[]
                            items_in_cat[cat].append([zp,n_stars])

                chart=np.array(chart)
                print("==========")
                print(date, filter)
                collection = "unlens_liris_zp"
                self.db.setCollection(collection)
                for cat_idx in items_in_cat.keys():
                    mean_cat=np.array(items_in_cat[cat_idx])
                    print(cat_idx,mean_cat[:,0].mean(),mean_cat[:,0].std(),mean_cat.size)
                    print("items in db",self.db.query({"id":date}).count())
                    if self.db.query({"id":date}).count() > 0:
                        query = {"$set":{str(date)+"."+str(filter)+"."+str(cat_idx):{"zp_mean": float(mean_cat[:,0].mean()), "std": float(mean_cat[:,0].std()),
                                           "nmeasures": int(mean_cat.size),"stars_mean":float(mean_cat[:,1].mean()), "stars_std": float(mean_cat[:,1].std())}}}
                        self.db.update({"id":str(date)},query)
                    else:
                        query={"id":str(date),str(date):{str(filter):{str(cat_idx):{"zp_men":float(mean_cat[:,0].mean()),"std":float(mean_cat[:,0].std()),"nmeasures":int(mean_cat.size),"stars_mean":float(mean_cat[:,1].mean()), "stars_std": float(mean_cat[:,1].std())}}}}
                        self.db.saveData(query)

                f, ax = plt.subplots(figsize=(9, 9))
                sns.despine(f, left=True, bottom=True)
                columns=['obs','seeing','airmass','zp','std','nstars','cat']
                dataframe=pd.DataFrame(np.array(chart),columns=columns)
                dataframe['obs'] = pd.to_numeric(dataframe['obs'])
                dataframe['zp'] = pd.to_numeric(dataframe['zp'])
                dataframe['seeing'] = pd.to_numeric(dataframe['seeing'])
                dataframe['airmass'] = pd.to_numeric(dataframe['airmass'])
                dataframe['std'] = pd.to_numeric(dataframe['std'])
                dataframe['nstars'] = pd.to_numeric(dataframe['nstars'])
                markers = ["o","*",">","D"]
                sns.scatterplot(x="obs", y="zp",
                                hue="airmass", size="nstars",
                                palette="rainbow", #"ch:r=-.2,d=.3_r",
                                sizes=(30, 50), linewidth=0,
                                data=dataframe, ax=ax,markers=markers,style="cat")


                for index, row in dataframe.iterrows():
                    ax.annotate(row["cat"]+"("+str(row["nstars"])+")", (row["obs"], row["zp"]))




                plt.savefig(base_path + "/" + date + filter + ".png")




    def genereatePhotometry(self,file_path):
        cat = CatalogUtil()
        path_as = file_path[0:-5] + "_as.fits"
        outfile = os.path.basename(path_as)
        dir_name= os.path.dirname(path_as)+"/"


        sw = SwarpUtil()
        stack_file = outfile[0:-5]+"_sw.fits"
        kargs={"-IMAGEOUT_NAME":stack_file}

        file_out = sw.runSwarp(path_as, **kargs)

        outfile_stack = os.path.basename(path_as)
        print(file_out)
        shutil.copy2(file_out,dir_name+stack_file)

        if os.path.exists(dir_name+stack_file):
           self.fixHeaderError(dir_name+stack_file)
           output_path=stack_file[:-5]+"_crop.fits"
           output_crop=ImageUtils.removeBorders(dir_name+stack_file,output_name=output_path,border_pix=300)
           cat_file=cat.createCatalogNormReport(output_crop)
           self.saveOnDataBase(cat_file)
           print(cat_file)
           return cat_file

    def createCatalogWithZP(self):
        collection = "unlens_liris"
        self.db.setCollection(collection)
        query = {"id": {"$exists":1}}
        alldata = self.db.getData(filter=query)

        collection = "unlens_liris_zp"
        self.db.setCollection(collection)
        query = {"id": {"$exists": 1}}
        allzp = self.db.getData(filter=query)

        base_path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/"
        ca = CatalogUtil()
        for objects in alldata:
            id = objects["id"]
            bands = objects["band"]
            for band in bands:
                date=list(bands[band].keys())[0]
                print(id,date,band)
                zp_idx = [zp for zp in allzp if date in zp][0]
                archives=zp_idx[date][band]
                seeing=bands[band][date]["catalog"][2]["seeing"]
                psf_file=bands[band][date]["catalog"][2]["psf_file"]
                #generate catalog from current zp basend on mean aff all observations for same archive (twomaass, ukidds)
                for archive in archives:

                    try:
                        zp = archives[archive]["zp_men"]
                        print(str(band)[band.find("_")+1:],band)
                        fits_file_image = base_path+id+"/"+id+"_"+str(band)[band.find("_")+1:]+"_as_sw_crop.fits"
                        output_path = base_path+id+"/"
                        filename=id+"_"+str(band)[band.find("_")+1:]+"_cat_"+archive
                        kargs = {}
                        kargs["-MAG_ZEROPOINT"] = str(zp)

                        kargs["-SEEING_FWHM"] = str(seeing)
                        kargs["-PSF_NAME"] = str(psf_file)
                        kargs["-CATALOG_NAME"] = filename+"_norm.fits"

                        sex_catalog_path = ca.generateCatalogFromFitsPath(fits_file_image, output_path=output_path,**kargs)

                        kargs["-DETECT_THRESH"] = str(1.5)
                        kargs["-ANALYSIS_THRESH"] = str(13)
                        kargs["-CATALOG_NAME"] = filename+"_norm_deep.fits"

                        sex_max_catalog_path = ca.generateCatalogFromFitsPath(fits_file_image, output_path=output_path,**kargs)
                        self.db.setCollection("unlens_liris")
                        query={"$set":{"band."+band+"."+date+".cat_phot."+archive:{"norm":sex_catalog_path,"deep":sex_max_catalog_path}}}
                        self.db.update({"id":id},query)
                    except Exception:
                        print(str(band)[band.find("_")+1:],band)
                #[band]
                print(date)
                print(zp_idx)
            print(id,band,date)

if __name__ == "__main__":
    path="/Users/cjimenez/Documents/PHD/data/unlens_liris/files.lst"
    path_sh = "/Users/cjimenez/Documents/PHD/data/unlens_liris/global.sh"



    path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/S1419.lst"
    files_list = np.loadtxt(path, dtype="str")
    lirisAs = LirisAstrometry()
    # lirisAs.createCatalogWithZP()
    # lirisAs.checkZeroPoints()
    # file="/Users/cjimenez/Documents/PHD/data/unlens_liris/L1715+2909/L1715+2909_J.fits"
    # catalog = lirisAs.genereatePhotometry(file)
    # sh_global_file = open(path_sh, "w+")
    for file in files_list:
        print("start astrometry for ",file)
        lirisAs.init(file)
        sh_file=lirisAs.doAstrometry()
        #catalog=lirisAs.genereatePhotometry(file)
        #print(catalog)

        #lirisAs.evaluateAsytrometry(file)
        # lirisAs.createCatalogWithZP(file)
        # sh_file=lirisAs.doAstrometry()
        # sh_global_file.write(sh_file+"\n")
    # sh_global_file.close()

