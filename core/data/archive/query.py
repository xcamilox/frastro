from astroquery.utils import TableList
from astropy.table import Table
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy import coordinates as coords
from astropy.coordinates import Angle
from django.conf.urls.static import static
from frastro import CoordinateParser, VOTableUtil
from frastro import SDSSArchiveCP
from frastro import UkidssArchiveCP
from frastro import PanSTARRSArchiveCP
from frastro import DecalsSurveyArchiveCP
from frastro import CFHTAchiveCP
from frastro import DesAchiveCP
from frastro import WiseAllWiseArchiveCP
from frastro import SpitzerArcvhiveCP
from frastro import VHSAchiveCP
from frastro import Config

import numpy as np
import pandas as pd


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
    __config=Config()


    def __init__(self):
        self.__procedureId = time.time()
        self.__results={"procedureId":self.__procedureId,"archives":[],"ra":0,"dec":0,"id":0,"pos":"","tag":[],"summary":{}}
        self.__tmpPath = self.__config.getPath("tmp")

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
        self.createSed()
        return self.__results


    def createSed(self):
        x_wavelenght=[]
        y_magnitudes=[]
        labels=[]
        item = self.__results
        id = self.__results["id"]
        for archive in self.__results["archives"]:
            for prop in archive["data"]["summary"]:
                sumary=archive["data"]["summary"][prop]
                if isinstance(sumary, dict):
                    if "lambda" in sumary:

                        x_wavelenght.append(sumary["lambda"])
                        y_magnitudes.append(sumary["ab"])
                        labels.append(archive["archive"]+"_"+prop)

        dtype = [('label', 'U20'), ('wavelength', float), ('abmagnitude', float)]

        for inx in range(len(y_magnitudes)):
            if y_magnitudes[inx] == "--" or float(y_magnitudes[inx]) < 0:
                y_magnitudes[inx] = 0

        tosort = []
        for ind in range(len(labels)):
            tosort.append((labels[ind], x_wavelenght[ind], y_magnitudes[ind]))

        sortval = np.array(tosort, dtype=dtype)
        sed = np.sort(sortval, order='wavelength')
        typesval = []
        typesval = [label.replace(" ", "_") for label in sed["label"].tolist()]

        vals = sed["abmagnitude"].tolist()

        # typesval.append(("ra",float))
        typesval.append("ra")
        vals.append(item["ra"])
        # typesval.append(("dec", float))
        typesval.append("dec")
        vals.append(item["dec"])

        setvals = np.array(vals).reshape(1, len(vals))

        df = pd.DataFrame(setvals, columns=typesval, dtype=float)
        catTable = Table.from_pandas(df)

        path_cat = self.__tmpPath + id + "/catalog.xml"

        df = pd.DataFrame(sed)
        sedTable = Table.from_pandas(df)

        SEDpath = self.__tmpPath + id + "/sed_" + id + ".xml"


        VOTableUtil.saveFromTable(catTable, path_cat)
        VOTableUtil.saveFromTable(sedTable, SEDpath)

        self.__results["sed"] = {"labels": sed["label"].tolist(), 'wavelength': sed["wavelength"].tolist(),
                        'abmagnitude': sed["abmagnitude"].tolist(), "local_path": SEDpath}

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