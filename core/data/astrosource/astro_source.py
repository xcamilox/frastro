import time
import json
class AstroSource(object):
    __id=""
    __catalogs=[] # list of catalogs source
    __images=[] #list of image source
    __spectra=[] #list of spectral
    __sumary={} #usfully data like: magnitud, redshift, name, id..etc
    __date_request=000000 #timestamp(utc) with the last request to this ra-dec


    def __init__(self,coordinates):
        ra=coordinates.ra.degree
        dec=coordinates.dec.degree
        position=coordinates.to_string('hmsdms').replace("h",":").replace("m",":").replace("d",":").replace("s","")

        self.__sumary={
            'id':str(str(ra)+"_"+str(dec)).replace(" ",""),
            'ra':ra,
            'dec':dec,
            'pos':position
        }
        self.__id = str(str(ra)+"_"+str(dec)).replace(" ","")
        self.__date_request = time.time()
        self.__catalogs = []
        self.__images = []
        self.__spectra = []

    def addImage(self,sourceImage):
        self.__images.append(sourceImage)

    def addCatalog(self,sourceCatalog):
        self.__catalogs.append(sourceCatalog)

    def addSpectra(self,sourceSpectra):
        self.__spectra.append(sourceSpectra)

    def addSummaryParams(self,key,value):
        self.__sumary[key]=value

    def getCatalogs(self):
        return self.__catalogs

    def getImages(self):
        return self.__images

    def getSpectra(self):
        return self.__spectra

    def getSummary(self):
        return self.__sumary

    def getId(self):
        return self.__id

    def __dict__(self):
        astrosource = {"id":self.__id,"catalogs": [], "images": [], "spectra": [], "summary": self.__sumary,
                       "request_date": self.__date_request}
        for cat in self.__catalogs:
            astrosource["catalogs"].append(cat.__dict__())

        for img in self.__images:
            astrosource["images"].append(img.__dict__())

        for spc in self.__spectra:
            astrosource["spectra"].append(spc.__dict__())

        return astrosource

    def toJson(self):
        return json.dumps(self.__dict__())