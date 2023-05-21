import ntpath

from frastro import FITSFile, ImageUtils, Telescopes, WaveLenghtCover, Config, VOTableUtil, LOGUtil, MongodbManager

from frastro import SDSSArchiveCP
from frastro import UkidssArchiveCP
from frastro import PanSTARRSArchiveCP
from frastro import DecalsSurveyArchiveCP
from frastro import CFHTAchiveCP
from frastro import DesAchiveCP
from frastro import WiseAllWiseArchiveCP
from frastro import SpitzerArcvhiveCP
from frastro import VHSAchiveCP
from frastro import TwoMassArchiveCP
from frastro.core.data.archive.apass_archive_cp import APASSArchiveCP
import json

from external.sextractor_util import SextractorUtil
from external.stilts_util import StiltsUtil
from external.psfex_util import PSFexUtil
from astropy.coordinates import SkyCoord
from astropy import units as u

import datetime
import numpy as np
import os
import matplotlib.pyplot as plt

from telescopes.hst import HSTelescope


class CatalogUtil():

    object_item={"ra":"","dec":"","telescope":"","band":"","photometry":[]}


    def __init__(self):
        self.object_item = {"ra": "", "dec": "", "telescope": "","instrument":"", "band": "", "photometry": [],"airmass":""}


    catlist={
        "sdss":SDSSArchiveCP(),
        "decals":DecalsSurveyArchiveCP(),
        "cfht":CFHTAchiveCP(),
        "ukidss":UkidssArchiveCP(),
        "panstarrs":PanSTARRSArchiveCP(),
        "des":DesAchiveCP(),
        "wise_unwise":WiseAllWiseArchiveCP(),
        "spitzer":SpitzerArcvhiveCP(),
        "two_mass":TwoMassArchiveCP(),
        "vhs":VHSAchiveCP(),
        "apass":APASSArchiveCP()
    }
    __band=""




    """
    return a list of catalogs where is center in center of img ref and band of image
    """
    def findCatalogFromImage(self,fits_image_path):

        output_path = os.path.dirname(fits_image_path)
        LOGUtil.log("========== Check Image Information : ==========\n", output_path)

        filters = ImageUtils.getFilterImage(fits_image_path)
        tele_info = Telescopes.getTelescopeData(fits_image_path)
        #if "FILTER" in filters[0]:
        #    filter=

        #filter = "SDSS-R"
        filter = filters[0][1]

        band = tele_info["instance"].getBand(filter)

        self.object_item["band"] = band


        header = FITSFile.header(fits_image_path)[0]
        object = str(header["OBJECT"]).replace(" ", "")
        name = object
        if object.find("_") > -1:
            name = object[:object.find("_")]
        self.__object = name
        self.__telescope = tele_info["telescope"]
        self.object_item["telescope"] = tele_info["telescope"]
        self.object_item["object"] = name
        self.object_item["date"] = str(header["DATE-OBS"])[:10]
        self.object_item["airmass"] = str(header["AIRMASS"])
        self.object_item["instrument"] = tele_info["instrument"]
        self.__band = band

        LOGUtil.log("Telescope Info: \n", output_path)

        range=WaveLenghtCover.getWavelenghtByTelescope(tele_info["telescope"],filter)

        LOGUtil.log("Telescope: "+tele_info["telescope"]+ " filter: "+filter, output_path)

        archives = WaveLenghtCover.waveLenghtCoverByFilter(range)



        archives=[str(archive).lower() for archive in archives["archive"].tolist()]

        if "gaia" in archives:
            archives.remove("gaia")

        #if "panstarrs" in archives:
        #    archives.remove("panstarrs")

        #if "sdss" in archives:
        #    archives.remove("sdss")


        print("archives",archives)
        if tele_info["telescope"] in archives:
            archives.remove(tele_info["telescope"])

        if "liverpool" in archives:
            archives.remove("liverpool")


        LOGUtil.log("Band: " + band, output_path)

        gap=2 #get 2 arcmin more for references catalogs
        radii_arcmin=ImageUtils.getImageRadiiSizeArcmin(fits_image_path)+gap
        centerImage_coordinate=ImageUtils.getCenterCoordinate(fits_image_path)

        self.object_item["ra"] = centerImage_coordinate.ra.deg
        self.object_item["dec"] = centerImage_coordinate.dec.deg

        LOGUtil.log("Find catalogs from : " + centerImage_coordinate.to_string('hmsdms')+ " in " +str(radii_arcmin) + "arcmin" ,output_path )



        dir_path=os.path.dirname(fits_image_path)
        tmp_path = dir_path +"/"+ object + "/ext_catalogs/"
        if not os.path.exists(tmp_path):
            os.makedirs(tmp_path)


        current_catalog=[]
        available_cat=[]
        for archive in archives:
            cat_ext_path=tmp_path+"catalog_to_" + archive + ".xml"
            cat_ext_path_empty = tmp_path + "catalog_to_" + archive + "empty.xml"
            if os.path.exists(cat_ext_path):
                current_catalog.append([archive,cat_ext_path])
            elif os.path.exists(cat_ext_path_empty):
                print("any detections",archive)
            else:
                available_cat.append(archive)

        if len(available_cat)>0:
            catalogs=self.searchCatalogInArchiveList(available_cat,centerImage_coordinate,radii_arcmin)

            for cat in catalogs:
                if len(cat[1])>0:
                    ref_cat_path = tmp_path + "catalog_to_" + cat[0] + ".xml"
                    VOTableUtil.saveFromTable(cat[1], ref_cat_path)
                    current_catalog.append([cat[0],ref_cat_path])
                else:
                    ref_cat_path = tmp_path + "catalog_to_" + cat[0] + "empty.xml"
                    f = open(ref_cat_path, "w+")
                    f.write("no founds!!")
                    f.close()


        str_cat = "Found Catalogs: \n"
        for cat in current_catalog:
            str_cat+=cat[0]+ " : "+str(len(cat[1]))+" sources \n"
        LOGUtil.log(str_cat, output_path)

        return current_catalog,band





    def generateCatalogNormalized(self,fits_file_image,output_path="",zeropoint="",good_cat=""):

        if output_path == "":
            output_path = os.path.dirname(fits_file_image)
            file_name=os.path.basename(fits_file_image)
            file_name=os.path.splitext(file_name)[0]
            output_path = output_path+"/reduction/"+file_name+"/"

        if zeropoint == "":
            zero_point,good_cat=self.getZeroPoint(fits_file_image,output_path)

            print("best catalog zp:", good_cat,zero_point)
        else:
            zero_point=zeropoint

        kargs={}
        kargs["-MAG_ZEROPOINT"]=str(zero_point)

        psf = PSFexUtil()
        seeing, psf_file = psf.calcSeeing(fits_file_image)
        kargs["-SEEING_FWHM"] = str(seeing)
        kargs["-PSF_NAME"] = str(psf_file)
        kargs["-CATALOG_NAME"] = "catalog_norm.fits"

        sex_catalog_path = self.generateCatalogFromFitsPath(fits_file_image,output_path=output_path,**kargs)

        kargs["-DETECT_THRESH"] = str(1.5)
        kargs["-ANALYSIS_THRESH"] = str(13)
        kargs["-CATALOG_NAME"] = "catalog_norm_deep.fits"

        sex_max_catalog_path = self.generateCatalogFromFitsPath(fits_file_image,output_path=output_path, **kargs)


        deep_cat=self.getDeepCatalog(sex_catalog_path)
        deep_max_cat = self.getDeepCatalog(sex_max_catalog_path)
        print(type(good_cat))
        if type(good_cat) == np.ndarray:
            self.createCatalogNormPlot(good_cat,sex_catalog_path,self.__object +" "+self.__band+" cat vs "+good_cat[0])
            self.createCatalogNormPlot(good_cat,sex_max_catalog_path,self.__object +" "+self.__band+" deep cat vs "+good_cat[0])
            ref_cat_to_save = good_cat.tolist()
        else:
            ref_cat_to_save = good_cat

        data_result = [sex_catalog_path,sex_max_catalog_path,{"seeing":seeing,"psf_file":psf_file,"cat_maglimit":deep_cat,"cat_deep_maglimit":deep_max_cat,"zp":zero_point,"cat_ref":ref_cat_to_save}]
        data_result = np.array(data_result)

        result_json=json.dumps(data_result.tolist())
        print(result_json)
        log_file = sex_catalog_path[0:-5]+"_cat_norm.log"
        log = open(log_file,"w+")
        log.write(result_json)
        log.close()

        item={}
        item["id"] = self.object_item["object"]
        item["band"] = self.object_item["band"]
        item["obs"]={self.object_item["date"]:{
                "photometry":self.object_item["photometry"],
                "catalog":data_result.tolist(),
                "telescope":self.object_item["telescope"],
                "seeing":seeing,
                "cat_ref":good_cat.tolist()[0],
                "airmass":self.object_item["airmass"],
                "instrument":self.object_item["instrument"]
            }}
        self.__item_report =item
        return data_result

    def createCatalogNormReport(self,fits_file_image,output_path="",zeropoint="",good_cat=""):
        self.generateCatalogNormalized(fits_file_image,output_path,zeropoint,good_cat)
        return self.__item_report


    def createCatalogNormPlot(self,cat_ref, cat_target,title):


        outpath_cat = self.crossMathTables(cat_target,cat_ref[1],cat_target[:-5]+"_refmacth.fits",cat_ref[0])


        hdul = FITSFile.open(outpath_cat)
        catalog = hdul[1].data


        cat_ref_archive = self.catlist[cat_ref[0]]
        band = "MAG_AUTO"
        ref1 = "MAG_AUTO"
        ref2 = "MAG_AUTO"
        if self.__band != "":
            band = cat_ref_archive.getBand(self.__band)
            ref1 = band
            index = np.where(catalog[ref1]>0,True,False)
            catalog = catalog[index]


        if ref1 == ref2:
            ref1="MAG_AUTO1"
            ref2 = "MAG_AUTO2"

        index = np.where(catalog[ref1] > 40, False, True)
        catalog = catalog[index]

        index = np.where(catalog[ref2] > 40, False, True)
        catalog = catalog[index]

        mag_auto_ref = catalog[ref1]
        mag_auto_target = catalog[ref2]
        plt.clf()
        plt.plot(mag_auto_ref,mag_auto_target,"o")
        plt.xlabel(ref2)
        plt.ylabel(ref1)
        plt.title(title)
        plt.savefig(cat_target[:-4]+"_macth.png")





    def getZeroPoint(self,fits_file_image,output_path =""):


        if output_path == "":
            output_path = os.path.dirname(fits_file_image)
            file_name = os.path.basename(fits_file_image)
            file_name = os.path.splitext(file_name)[0]
            match_path = output_path + "/reduction/" + file_name + "/stilts/"


        else:

            output_path = output_path if output_path[-1] == "/" else output_path + "/"
            tmp_path = output_path + "tmp"
            now = datetime.datetime.now()
            tmp_path = tmp_path +"/"+ str(now.day) + str(now.month) + str(now.year)
            match_path = tmp_path

        if not os.path.exists(match_path):
                os.makedirs(match_path)

        tele_info = Telescopes.getTelescopeData(fits_file_image)
        if str(tele_info["telescope"]).lower() == "hst":
            self.setTargetParams(fits_file_image)
            zp = HSTelescope.getZeroPoint(fits_file_image)["ABmag"]
            return zp, "hst"
        else:

            catalogs,band = self.findCatalogFromImage(fits_file_image)



            print(catalogs,band)


            sex_catalog_path = self.generateCatalogFromFitsPath(fits_file_image)
            zero_points = []
            good_cat=[]
            for cat in catalogs:
                out_math = match_path + "cross_match_sex_" + cat[0] + ".fits"
                path_math_catalogs = self.crossMathTables(sex_catalog_path, cat[1], out_math, cat[0])
                band_key = self.catlist[cat[0]].getBand(band)
                print("band_key",band_key,cat[0])
                filter = {}
                if cat[0] == "sdss":
                    filter = {"type": "==6"} # get zero point only from stars

                if cat[0] == "vhs":
                    filter = {band_key: ">0"}  # get zero point only from stars

                mean_zp, std_zp,items,cat_fit = self.calcZeroPoint(path_math_catalogs, band_key,cat[0],filter)
                print(mean_zp, std_zp,path_math_catalogs)
                if mean_zp > 0 and std_zp >= 0:
                    zero_points.append([mean_zp, std_zp,items])
                    good_cat.append([str(cat[0]),str(cat[1])])

            print("zero_points:",zero_points)


            if len(zero_points) > 0:
                zero_points = np.array(zero_points)

                min=zero_points[:, 1].min()
                max = zero_points[:, 2].max()


                idx = np.where(zero_points[:, 1]==min,True,False)

                zero_point = zero_points[idx][0][0]
                amount_data = zero_points[idx][0][2]

                if max/amount_data > 0.7:
                    idx = np.where(zero_points[:, 2] == max, True, False)
                    zero_point = zero_points[idx][0][0]

                good_cat = np.array(good_cat)[idx][0]

                return zero_point,np.array(good_cat).astype(str)
            else:
                return 0,"any"


    def setTargetParams(self,fits_image_path):
        header = FITSFile.header(fits_image_path)[0]
        if "TARGNAME" in header:
            object = str(header["TARGNAME"]).replace(" ", "")
        else:
            object = str(header["OBJECT"]).replace(" ", "")
        self.__object = object

    def generateCatalogFromFitsPath(self,fits_path,output_path="",**kargs):
        sex = SextractorUtil()
        return sex.runSextractor(fits_path,output_path=output_path, **kargs)

    def crossMathTables(self,table1,table2,output_path,catalog="",ref_table2="values2=ra dec",separation="1"):

        print("generating crossmath catalog",table1,table2)
        pos_cat="values1=ALPHA_J2000 DELTA_J2000"
        pos_reference=ref_table2
        #separation='1' # in arcsec
        ouput_type="ofmt=fits-basic"

        if catalog=="sdss":
            pos_reference = "values2=ra dec"
        elif catalog=="panstarrs":
            pos_reference = 'values2=hmsToDegrees(raMean) dmsToDegrees(decMean)'
        elif catalog == "cfht" or catalog == "wise_unwise" or catalog == "ukidss" or catalog == "apass":
            pos_reference = "values2=RAJ2000 DEJ2000"

        if catalog == "two_mass":
            separation= '2'




        params = [pos_cat, pos_reference, ouput_type]

        stilts=StiltsUtil()
        stilts.crossMatchTwoCatalogs(table1,table2,output_path,separation,params)



        return output_path

    #get zero porin from a differences between a reference catalog(ej:sdss,panstars) and sextractor catalog from image
    # the key band is a column in reference catalog, and compare with mag_auto sextractor column
    # the catalog is a fits file with data in 1 index of hdulist, this could be generated by stilts
    """
        return the zero poin base on difference mangnitud from sextractor catalog to reference catalog
       :returns mean, std  
    """
    def calcZeroPoint(self,master_cat_fits,band,cat,filter={}):
        print("band",band)
        print("calculating zero point",master_cat_fits,band,filter)
        hdul=FITSFile.open(master_cat_fits)
        catalog=hdul[1].data

        if len(filter.keys()) > 0:
            try:
                for key in filter.keys():
                    filter_column = catalog[key]
                    index=np.where(eval("filter_column "+filter[key]),True,False)
                    catalog=catalog[index]
            except KeyError:
                print("No found key in catalog")

        print("filters", filter.keys(),len(filter.keys()))


        #If no found data in each band, skip this catalog
        try:
            #select stars from class_star > 0.7

            band_in_cat_ref = catalog[band]

            mag_auto = catalog["MAG_AUTO"]
            index = np.where(band_in_cat_ref > 0, True, False)
            band_in_cat_ref = band_in_cat_ref[index]
            catalog = catalog[index]

            stars = catalog["CLASS_STAR"]
            index_stars = np.where(stars > 0.7, True, False)
            if index_stars.min() == index_stars.max():
                print("no found stars, try to use all data")
                index_stars = [True]*len(index_stars)
            mag_auto = mag_auto[index]
            good_stars = stars[index_stars]
            mag_auto = mag_auto[index_stars]
            band_in_cat_ref = band_in_cat_ref[index_stars]
            catalog = catalog[index_stars]
            index_stars = np.where(mag_auto < 90., True, False)

            mag_auto = mag_auto[index_stars]
            band_in_cat_ref = band_in_cat_ref[index_stars]
            catalog = catalog[index_stars]


            difference=band_in_cat_ref-mag_auto
            mean=difference.mean()
            std=difference.std()
            if band_in_cat_ref.size > 0:
                self.object_item["photometry"].append({"catalog":cat,"file":master_cat_fits,"n_stars":band_in_cat_ref.size,"zp_mean":float(np.abs(mean)),"zp_std":float(std),"stars_ra":catalog["ALPHA_J2000"].tolist(),"stars_dec":catalog["ALPHA_J2000"].tolist(),"stars_mag":catalog["MAG_AUTO"].tolist()})
            #zeropoint = 31.000208598926257-1.74
            #return zeropoint, std, band_in_cat_ref.size, catalog
            return np.abs(mean),std,band_in_cat_ref.size,catalog
        except KeyError:
            return -99,-99,0,"any"




    def searchCatalogInArchiveList(self,cat_list,coordinates,radio):

        catalogs=[]
        print("catlist",cat_list)
        for cat in cat_list:
            catalog = self.catlist[cat]
            try:
                print(coordinates.ra.deg, coordinates.dec.deg,radio)
                r = catalog.getCatalog(coordinates.ra.deg, coordinates.dec.deg, radio)
                good_cat = 0
                if type(r) is tuple:
                    r = r[0]
                if len(r) > 2 and len(catalogs) < 5:
                    catalogs.append([cat, r])
                    good_cat += 1
                    if good_cat >= 5:
                        break
                else:
                    catalogs.append([cat, []])
            except:
                print("service unabled "+cat)
        return catalogs


    def getDeepCatalog(self, cat_path):
        hdul = FITSFile.open(cat_path)
        data=hdul[1].data

        remove_bad_detections = data["MAG_AUTO"] < 50.0
        data = data[remove_bad_detections]

        filter=data["CLASS_STAR"] > 0.6
        good = data[filter]
        best = good["MAG_AUTO"] > 16.0

        good = good[best]

        mag = good["MAG_AUTO"]



        return {"point_deep":float(mag.mean()) , "point_std":float(mag.std()), "n_point_sources":int(mag.shape[0]),"sources_deep":float(data["MAG_AUTO"].mean()),"sources_std":float(data["MAG_AUTO"].std()),"n_sources":int(data.shape[0])}



    def searchAroundCoordinates(self,catalog_path,ra,dec,distance_arc=5,ra_val="ALPHA_J2000",dec_val="DELTA_J2000"):
        hdul = FITSFile.open(catalog_path)
        data = hdul[1].data
        ra_detections= data[ra_val]
        dec_detections = data[dec_val]
        # search detection around the target

        catalog = SkyCoord(ra=ra_detections * u.degree, dec=dec_detections * u.degree)
        target_position = SkyCoord(ra=np.array([ra]) * u.degree, dec=np.array([dec]) * u.degree)

        # search around 10 arcsec
        idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(target_position, distance_arc * u.arcsec)
        print(d2d, d2d.to("arcsec"))
        detections = {"detections": 0}
        if idxcatalog.size > 0:
            detections = {"detections": idxcatalog.size, "ra": data[idxcatalog]["ALPHA_J2000"][0],
                          "dec": data[idxcatalog]["DELTA_J2000"][0], "mag": float(data[idxcatalog]["MAG_AUTO"][0]),
                          "mag_err":float(data[idxcatalog]["MAGERR_AUTO"][0]),
                          "flux":float(data[idxcatalog]["FLUX_AUTO"][0]),"flux_err":float(data[idxcatalog]["FLUXERR_AUTO"][0]),
                          "distance": str(d2d[0].to("arcsec")),"ELLIPTICITY":float(data[idxcatalog]["ELLIPTICITY"][0]),"CLASS_STAR":float(data[idxcatalog]["CLASS_STAR"][0]),"ELONGATION":float(data[idxcatalog]["ELONGATION"][0]),"THETA_IMAGE":float(data[idxcatalog]["THETA_IMAGE"][0])}



        return detections

    @staticmethod
    def saveTable(self,table,output_path,format):
        table.write(output_path,format=format)


    def searchDetectionbyTarget(self,id,ra,dec,stack,epoach):
        cat = CatalogUtil()
        try:

            #for epoach in target["epoach"]:
            #stack = epoach
            if "error" not in stack.keys():
                cat_best = stack["best"]["cat_path"]
                cat_best_deep = stack["best"]["cat_path_deep"]

                cat_best_all = stack["all"]["cat_path"]
                cat_best_deep_all = stack["all"]["cat_path_deep"]

                best_cat_detections = cat.searchAroundCoordinates(cat_best,ra,dec)
                best_cat_deep_detections = cat.searchAroundCoordinates(cat_best_deep, ra, dec)

                all_cat_detections = cat.searchAroundCoordinates(cat_best_all, ra, dec)
                all_cat_deep_detections = cat.searchAroundCoordinates(cat_best_deep_all, ra, dec)
                query = {"$set":{"epoach."+epoach+".stacking.best.detection.best":best_cat_detections,"epoach."+epoach+".stacking.best.detection.best_deep":best_cat_deep_detections,"epoach."+epoach+".stacking.all.detection.best":all_cat_detections,"epoach."+epoach+".stacking.all.detection.best_deep":all_cat_deep_detections}}
                return query


        except KeyError:
            print("error with ", id)
            return None


if __name__ == "__main__":
    cat=CatalogUtil()
    """
    #search current detection
    mgdb = MongodbManager()
    mgdb.setCollection("liverpool_sn_v3")
    all_data = mgdb.getData({'epoach': {'$exists': 'true'}})
    for target in all_data:
        try:
            id = target["id"]
            ra = target["ra"]
            dec = target["dec"]
            for epoach in target["epoach"]:
                stack = target["epoach"][epoach]["stacking"]
                if "error" not in stack.keys():
                    cat_best = stack["best"]["cat_path"]
                    cat_best_deep = stack["best"]["cat_path_deep"]

                    cat_best_all = stack["all"]["cat_path"]
                    cat_best_deep_all = stack["all"]["cat_path_deep"]

                    best_cat_detections = cat.searchAroundCoordinates(cat_best,ra,dec)
                    best_cat_deep_detections = cat.searchAroundCoordinates(cat_best_deep, ra, dec)

                    all_cat_detections = cat.searchAroundCoordinates(cat_best_all, ra, dec)
                    all_cat_deep_detections = cat.searchAroundCoordinates(cat_best_deep_all, ra, dec)

                    mgdb.update({"id":id},{"$set":{"epoach."+epoach+".stacking.best.detection.best":best_cat_detections,"epoach."+epoach+".stacking.best.detection.best_deep":best_cat_deep_detections,"epoach."+epoach+".stacking.all.detection.best":all_cat_detections,"epoach."+epoach+".stacking.all.detection.best_deep":all_cat_deep_detections}})
        except KeyError:
            print("error with ", target["id"])
            
            
            



    #generate catalog normalize
    path_folder="/Users/cjimenez/Documents/PHD/data/unlens_liris/astrometria_ismael/"

    onlyfiles = [f for f in os.listdir(path_folder) if os.path.isfile(path_folder + "/" + f)]

    # path=path_folder+"SDSSJ0850+1549_wht_liris_h_20190123_ipf_astro_gaia_dr2_0p03rms.fits"
    # path="/Users/cjimenez/Documents/PHD/data/unlens_liris/astrometria_ismael/J020947_wht_liris_ks_20190119_ipf_astro_gaia_dr2_0p06rms.fits"
    # cat.generateCatalogNormalized(path)
    
    

    for file in onlyfiles:
      if os.path.splitext(file)[1] == ".fits" or os.path.splitext(file)[1] == ".fit":
          print(file)
          cat_norm = cat.generateCatalogNormalized(path_folder+file)
          logfile=path_folder+file[:-5]+".log"

          file=open(logfile,"w+")
          file.write(str(cat_norm))
          file.close()
          print(logfile)


    path="/Users/cjimenez/Documents/PHD/data/adfs27/hst/idt106010_drz.fits"
    #image=HSTelescope.getSingleImage(path,path[:-5]+"_single.fits")
    #image ='/Users/cjimenez/Documents/PHD/data/adfs27/hst/idt106010_drz_single.fits'
    file1 = "/Users/cjimenez/Documents/PHD/data/adfs27/science/adfs_27_fors2_r_band_crop.fits"
    file2 = "/Users/cjimenez/Documents/PHD/data/adfs27/science/adfs_27_fors2_z_band_crop.fits"

    images=[file1,file2]
    images = [file2]
    for image in images:
        cat_norm=cat.generateCatalogNormalized(image)
        print(cat_norm)
        logfile=image[:-5]+".log"
        print(logfile)
        file=open(logfile,"w+")
        file.write(str(cat_norm))
        file.close()
    #HSTelescope.getZeroPoint(image)
    #cat.generateCatalogNormalized(path)

    """

    path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/S1419+4341_ks_stack.fits"
    path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/S1419+4341_has.fits"
    #path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/S1419+4341_jas.fits"
    #path="/Users/cjimenez/Documents/PHD/data/liverpool_sn_2020/2020ank/images/h_e_20200130_55_1_1_1.fits" #band u

    #path="/Users/cjimenez/Documents/PHD/data/liverpool_sn_2020/2020ank/images/h_e_20200130_56_1_1_1.fits" #band g
    #path="/Users/cjimenez/Documents/PHD/data/liverpool_sn_2020/2020ank/images/h_e_20200130_57_1_1_1.fits" #band R
    #path="/Users/cjimenez/Documents/PHD/data/liverpool_sn_2020/2020ank/images/h_e_20200130_58_1_1_1.fits" #band I
    path="/Users/cjimenez/Documents/PHD/data/liverpool_sn_2020/2020ank/images/h_e_20200130_59_1_1_1.fits" #band Z
    path="/Users/camilojimenez/Desktop/fred/h_e_20200130_56_1_1_1.fits"


    cat = cat.generateCatalogNormalized(path)
    print(cat)