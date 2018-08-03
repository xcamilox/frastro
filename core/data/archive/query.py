from astroquery.utils import TableList
from astropy.table import Table
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy import coordinates as coords
from astropy.coordinates import Angle
from django.conf.urls.static import static
from frastro import CoordinateParser
from frastro import SDSSArchiveCP
from frastro import UkidssArchiveCP
from frastro import PanSTARRSArchiveCP
from frastro import DecalsSurveyArchiveCP
from frastro import CFHTAchiveCP
from frastro import DesAchiveCP
from frastro import WiseAllWiseArchiveCP
from frastro import SpitzerArcvhiveCP
from frastro import VHSAchiveCP


import json

from astropy import units as u
import time

class Query():

    __coordinates="0,0"
    __results={}
    __procedureId=00
    __radius=1.0
    __id=0
    __tags=[]


    def __init__(self):
        self.__procedureId = time.time()
        self.__results={"procedureId":self.__procedureId,"archives":[],"ra":0,"dec":0,"id":0,"pos":"","tag":[],"summary":{}}

    def __addResult(self,archive,status,result,tags):



        self.__results["summary"][archive] = {"det": 1 if len(result["summary"]) > 4 else 0,
                                              "img": 1 if len(result["images"]) >= 1 else 0,
                                              "spec": 1 if len(result["spectra"]) >= 1 else 0}
        if len(tags)>0:
            self.__results["tag"]=tags
        if "error" in result["summary"]:
            status = 300
            self.__results["summary"][archive]["det"] = 0


        self.__results["archives"].append({"archive":archive,"status":status,"data":result})



    def search(self,coordinates,radius=1,tags=[]):


        self.__pos = CoordinateParser.validateCoordinates(coordinates)
        self.__results["ra"]=self.__pos.ra.degree
        self.__results["dec"] = self.__pos.dec.degree
        self.__results["pos"] = self.__pos.to_string('hmsdms').replace("h",":").replace("m",":").replace("d",":").replace("s","")
        self.__id = str(self.__pos.ra.degree)+"_"+str(self.__pos.dec.degree)
        self.__results["id"] = self.__id
        self.__radius = radius
        self.__addResult("SDSS",200, self.searchSDSS(),tags)
        self.__addResult("PanSTARRS", 200, self.searchPanStarrs(),tags)
        self.__addResult("Ukidss", 200, self.serchUkidss(),tags)
        self.__addResult("CFHT",200,self.searchCanadianCFHT(),tags)
        self.__addResult("decals",200, self.searchLegacySurvey(),tags)
        self.__addResult("wise_unwise", 200, self.searchWiseAllWise(),tags)
        self.__addResult("spitzer", 200, self.searchSpitzer(), tags)
        self.__addResult("DES", 200, self.searchDesSurvey(), tags)
        self.__addResult("vhs", 200, self.searchVHS(), tags)


        #self.__addResult("vizier", self.searchVizier())
        #self.__addResult("simbad",self.searchSimbad())

        return self.__results



    def serchUkidss(self):
        search_catalog_cp = UkidssArchiveCP()
        respond = search_catalog_cp.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree,radius=self.__radius)
        return respond.__dict__()

    def searchPanStarrs(self):
        pns = PanSTARRSArchiveCP()
        respond = pns.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree,radius=self.__radius)
        return respond.__dict__()


    def searchSimbad(self):
        return Simbad.query_region(self.__pos, radius=self.__radius)


    def searchVizier(self):

        return Vizier.query_region(self.__pos, radius=self.__radius)

    def getId(self):
        return self.__id

    def searchSDSS(self):
        sdss = SDSSArchiveCP()
        respond = sdss.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree, release=14,radius=self.__radius)
        return respond.__dict__()

    def searchLegacySurvey(self):
        lgsurveryAC = DecalsSurveyArchiveCP()
        r = lgsurveryAC.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree, radius=self.__radius)
        return r.__dict__()

    def searchCanadianCFHT(self):
        canadianSurvey = CFHTAchiveCP()
        r = canadianSurvey.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree, radius=self.__radius)
        return r.__dict__()

    def searchDesSurvey(self):
        deSurveryAC = DesAchiveCP()
        r = deSurveryAC.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree, radius=self.__radius)
        return r.__dict__()

    def searchWiseAllWise(self):
        wiseallwise=WiseAllWiseArchiveCP()
        r = wiseallwise.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree, radius=self.__radius)
        return r.__dict__()

    def searchSpitzer(self):
        spitzer=SpitzerArcvhiveCP()
        r = spitzer.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree, radius=self.__radius)
        return r.__dict__()

    def searchVHS(self):
        search_catalog_cp = VHSAchiveCP()
        r = search_catalog_cp.query(ra=self.__pos.ra.degree, dec=self.__pos.dec.degree, radius=self.__radius)
        return r.__dict__()