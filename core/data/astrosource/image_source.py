import json
from frastro import Utils
from frastro import ImageUtils
import datetime
class ImageSource(object):
    __name=""
    __cutouts=[]
    __files=[]
    __provider=""
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/{1}"

    def __init__(self,name,provider,save_path=""):
        self.__name = name
        self.__provider=provider
        self.__cutouts=[]
        self.__files = []
        if save_path !="":
            self.__save_path = save_path




    def addCutout(self,url,size,band,magnitude=-99.,link=""):

        self.__cutouts.append({"url":url,"size":size,"band":str(band),"magnitude":str(magnitude),"link":link})

    def addFile(self,name,url,type="fits",download=True,local_path="",uncompress=False,thumbnail=True,external=False):

        if download:
            if local_path == "":
                local_path = self.__save_path = self.__save_path.format(self.__provider,datetime.datetime.now().timestamp())

            dowload_path = Utils.saveLocalFile(url, local_path, uncompress,external=external)

            if dowload_path is not None:
                local_path = dowload_path
                if thumbnail:
                    pos = self.__name.split("_")
                    mark = Utils.getPixelPositionFromRaDec(float(pos[0]), float(pos[1]), local_path)
                    thumbnail=ImageUtils.getBase64FromFitsURL(local_path,mark=mark)
                else:
                    thumbnail = ""

        self.__files.append({"url": url, "name": name, "type": type,"local_path":local_path,"thumbnail":thumbnail})

    def getFiles(self):
        return self.__files
    def __dict__(self):
        image={"name":self.__name,"provider":self.__provider,"cutouts":self.__cutouts,"files":self.__files}
        return image