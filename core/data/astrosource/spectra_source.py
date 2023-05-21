import json
class SpectraSource(object):
    __name=""
    __cutouts=[]
    __files=[]
    __provider=""
    __cache_date = False

    def __init__(self,name,provider):
        self.__name = name
        self.__provider=provider
        self.__cutouts = []
        self.__files = []

    def addCutout(self,url,size,band):
        self.__cutouts.append({"url":url,"size":size,"band":band})

    def addFile(self,name,url,type):
        self.__files.append({"url": url, "name": name, "type": type})

    def __dict__(self):
        spectra = {"name": self.__name, "provider": self.__provider, "cutouts": self.__cutouts, "files": self.__files}
        return spectra
