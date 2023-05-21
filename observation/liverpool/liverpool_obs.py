import os

from frastro import AstroSource, FITSFile
from frastro.external.psfex_util import PSFexUtil
from frastro.core.data.catalog.catalog_util import CatalogUtil
from astropy.time import Time
from frastro import MongodbManager
import numpy as np
from astropy.io import fits

class LiverPoolObsservation():
    __target = AstroSource
    __group = ""
    __texp = 0.0
    __seeing= 0.0
    __date=0
    __mjd=0
    __collection="liverpool_sn_v3"
    __mgdb = MongodbManager()


    def __init__(self,**kargs):
        if len(kargs.keys())>0:
            self.setObject(**kargs)

        self.__mgdb.setCollection(self.__collection)

    def setObject(self,**kargs):
        if "source" in kargs:
            self.__target=kargs["source"]

        if "group" in kargs:
            self.__group=kargs["group"]

        if "texp" in kargs:
            self.__texp==kargs["texp"]

        if "seeing" in kargs:
            self.__seeing==kargs["seeing"]

        if "date" in kargs:
            self.__date = kargs["date"]
            date = Time(self.__date, format='isot', scale='utc')
            self.__mjd = date.mjd

    def getObservations(self,id=None):
        if id==None:
            list = self.__mgdb.getData()#projection={"id":1,"ra":1,"dec":1,"magnification":1,"n_epoach":1,"z_lens":1,"z_source":1})
        else:
            list = self.__mgdb.getData(filter={"id":id})
        return list




    def saveObservation(self):

        if self.__collection != "":

            self.__mgdb.saveData(self.__dict__())


    def __dict__(self):
        return {"target":"","group":self.__group,"texp":self.__texp,"seeing":self.__seeing,"date":self.__date,"mjd":self.__mjd}






if __name__ == "__main__":
    liverpoolObs= LiverPoolObsservation()
    mgdb = MongodbManager()
    collection="liverpool_sn_v3"
    mgdb.setCollection(collection=collection)
    """
    data=mgdb.getData(filter={},projection={"id":1,"apoach":1})

    for items in data:
        if "apoach" in items.keys():
            print(len(items["apoach"]),items["id"])
            #mgdb.update(items["id"],{"$set":{"n_epoach":len(items["apoach"])}})


    #psf = PSFexUtil()
    #cat = CatalogUtil()
    
    """
    # redfile
    base_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/"
    location_file = base_path + "liverpool_lens_radec_name.fits"

    file = fits.open(location_file)
    lens_list = file[1]
    data= lens_list.data
    objects=[]

    

    for idx,source in enumerate(data.field('name')):
        path_source = base_path+source
        ra=data["RADeg"][idx]
        dec = data["DECDeg"][idx]
        coordinates=str(data["ra"][idx]) + " "+str(data["dec"][idx])
        z_lent=str(data["z_l"][idx])
        z_source = str(data["z_s"][idx])
        magnification = str(data["magnification"][idx])
        item={"id":source,"ra":ra,"dec":dec,"z_lens":z_lent,"z_source":z_source,"magnification":magnification,"coordinates":coordinates,"n_epoach":0}
        mgdb.saveData(item, collection=collection)

        """
        if not os.path.exists(path_source):
            print("new source",source)
            item["n_epoach"] = 0
            os.makedirs(path_source)
        else:
            all_observations = path_source+"/liverpool_obs"
            epoachs = {}
            if os.path.exists(all_observations):
                dirnames = [dI for dI in os.listdir(all_observations) if os.path.isdir(os.path.join(all_observations, dI))]

                for epoach in dirnames:

                    path_folder = all_observations+"/"+epoach
                    onlyfiles = [f for f in os.listdir(path_folder) if os.path.isfile(path_folder + "/" + f)]
                    all_files_path = []
                    for file in onlyfiles:
                        if os.path.splitext(file)[1] == ".fits" or os.path.splitext(file)[1] == ".fit":
                            img_path = path_folder + "/" + file

                            header = FITSFile.header(img_path)[0]
                            #save image info
                            date = header["DATE"] #epoach
                            epoc_data={}
                            if date not in epoachs.keys():
                                epoachs[date]=[]

                            epoc_data["date"]=date
                            epoc_data["id"] = os.path.splitext(file)[0]
                            epoc_data["date_obs"] = header["DATE-OBS"]
                            epoc_data["object"] = header["OBJECT"]

                            epoc_data["group_id"] = header["GROUPID"]
                            epoc_data["mjd"] = header["MJD"]
                            epoc_data["exptime"] = header["EXPTIME"]
                            epoc_data["airmass"] = header["AIRMASS"]

                            seeing,psfpath = psf.calcSeeing(img_path)
                            output_psf=psfpath[:-11]+"psf.fits"
                            zeropoint,zp_from_cat = cat.getZeroPoint(img_path)
                            psf.saveSinglePSFFile(psfpath, output_psf)
                            epoc_data["seeing_cal"] = seeing
                            epoc_data["seeing_img"] = header["SEEING"]
                            epoc_data["img_path"] = img_path
                            epoc_data["psf_path"] = output_psf
                            epoc_data["zeropoint"] = {"zp":zeropoint,"cat":zp_from_cat[0],"cat_file":zp_from_cat[1]}

                            epoachs[date].append(epoc_data)


                print(source,"epoach:",len(dirnames))
                item["n_epoach"]=len(dirnames)
            print("already take",source)
            item["epoachs"]=epoachs

        objects.append(item)
        mgdb.saveData(item,collection=collection)
        print(objects)
    """










