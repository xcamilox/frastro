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

class LasCumbresPipeline():
    def __init__(self):
        self.project_path = "/Volumes/scratch/tmp_LCO/"
        files=self.getData(self.project_path)
        print(files)
        for path_file in files:
            #LOGUtil.log("-----------Liverpool pipeline---------"+path_object,path_object)
            path_object = files[path_file]
            stack_file = self.stackFiles(path_object)
            if stack_file != None:
                catalog_file = self.createCatalogs(path_object)
                self.moveToScienceDirectory(path_object,[stack_file,catalog_file])

    def moveToScienceDirectory(self,path_object,files):
        science_path=Utils.createPath(path_object+"/science")
        swarp_image = files[0]
        sex_catalog = files[1]

        if os.path.exists(swarp_image):
            shutil.copy2(swarp_image,science_path+"/stacking_img.fits")

        shutil.copy2(sex_catalog[0],science_path+"/stacking_catalog.fits")
        shutil.copy2(sex_catalog[1],science_path+"/stacking_catalog_deep.fits")


        if path.exists(science_path+"/report.json"):
            os.remove(science_path+"/report.json")

        f=open(science_path+"/report.json","w+")
        f.write("zeropoint:"+str(sex_catalog[2]["zp"])+" cat: "+str(sex_catalog[2]["cat_ref"])+"\n")
        f.close()


    def getData(self,path_folder):
        paths = {}
        onlyfiles = [f for f in os.listdir(path_folder) if path.isfile(path_folder + "/" + f)]

        for file in onlyfiles:
            if path.splitext(file)[1] == ".fits" or path.splitext(file)[1] == ".fit":
                header = FITSFile.header(path_folder + "/" + file)[0]
                target = str(header["OBJECT"]).replace(" ","").replace(".","").replace("-","")
                obs_date = header["DATE"]

                science_folder = self.project_path + target + "/science/"
                ext_catalogs = self.project_path + target + "/ext_catalogs/"
                observations_folder = self.project_path + target + "/observations/"
                target_path = observations_folder

                # create science folder
                # create ext_catalogs folder
                # create observations folder


                # check data base
                file_name = path.splitext(file)[0]
                if target not in paths.keys():
                    paths[target] = target_path


                if not os.path.exists(science_folder):
                    os.makedirs(science_folder)

                if not os.path.exists(ext_catalogs):
                    os.makedirs(ext_catalogs)

                if not os.path.exists(observations_folder):
                    os.makedirs(observations_folder)

                newlocation = target_path + file

                if not os.path.exists(target_path):
                    os.makedirs(target_path)

                # create files list to each apoach
                files_list_eapoach = target_path + "files.txt"
                if not os.path.exists(files_list_eapoach):
                    file_list_epoach_file = open(files_list_eapoach, "w+")
                else:
                    file_list_epoach_file = open(files_list_eapoach, "a+")

                # create files list for all files
                files_list_all = observations_folder + "files.txt"
                if not os.path.exists(files_list_all):
                    file_list_all_file = open(files_list_all, "w+")
                else:
                    file_list_all_file = open(files_list_all, "a+")

                # shutil.move(path_folder + "/" + file, newlocation)
                file=fits.open(path_folder + "/" + file)[0]
                data=file.data
                header=file.header
                FITSFile.saveFile(newlocation,npa_data=data,header=header)
                #shutil.copy2(path_folder + "/" + file, newlocation)

                file_list_epoach_file.write(newlocation + "\n")
                file_list_all_file.write(newlocation + "\n")

                file_list_epoach_file.close()
                file_list_all_file.close()

                #cat_info = self.createCatalog(newlocation)

                #header_info = self.readHeaderInfo(newlocation)
                #object_info = {**cat_info, **header_info}



        return paths

    def createCatalog(self,file_path):

        cat_util = CatalogUtil()
        zeropoint, zp_from_cat = cat_util.getZeroPoint(file_path)
        sex_cat = cat_util.generateCatalogNormalized(file_path,zeropoint=zeropoint,good_cat=zp_from_cat)

        respond={"zp":zeropoint,"ref_cat":zp_from_cat[0],"ref_cat_path":zp_from_cat[1],"cat_path":sex_cat[0],"cat_path_deep":sex_cat[1],"seeing_cal":sex_cat[2]["seeing"],"psf_path":sex_cat[2]["psf_file"],"cat_deep_maglimit":sex_cat[2]["cat_deep_maglimit"],"cat_maglimit":sex_cat[2]["cat_maglimit"]}
        return respond

    def stackFiles(self,folder_path):
        LOGUtil.log("liverpool-pipeline Creating staking image for "+ folder_path, folder_path)
        swarp=SwarpUtil()
        return swarp.runSwarp(folder_path,quality="all")


    def createCatalogs(self,folder_path):
        LOGUtil.log("liverpool-pipeline Creating catalog for " + folder_path, folder_path)
        sex = SextractorUtil()
        folder_path = folder_path if folder_path[-1] == "/" else folder_path + "/"
        path_files = folder_path+"reduction/swarp/coadd.fits"
        output_path = folder_path+"reduction/sextractor/"
        cata_util=CatalogUtil()
        sex_cat=cata_util.generateCatalogNormalized(path_files,output_path)
        return sex_cat
        #sex.runSextractor(path_files,output_path)


if __name__ == "__main__":
    LasCumbresPipeline()