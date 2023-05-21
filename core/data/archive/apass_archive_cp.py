from frastro import ContentProvider
from frastro import TAPManager
from frastro import AstroSource
from frastro import CoordinateParser
import astropy.units as u
from frastro import ImageSource
from frastro import CatalogSource
import numpy as np
from frastro.core.utils.html_parse import HtmlParser
import requests
import re
from bs4.element import Tag
from frastro import ImageUtils
from frastro import VOTableUtil
from frastro import Utils
from frastro import WaveLenghtCover

class APASSArchiveCP(ContentProvider):
    __catalog_provider = "CFHT"

    __service_provider = {
        "tap_url": "http://tapvizier.u-strasbg.fr/TAPVizieR/tap",
        "width_cutout":"http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/megapipe/imc.pl" #"?lang=en&object=&ra={0}&dec={1}&size=256"
    }

    __simple_rec_query = 'SELECT * FROM "II/336/apass9" as w WHERE '+"1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {0},{1}, {2} ))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/apass/"
    #all unit in micras
    __wavelenght={"g":475*u.nm,"r":640*u.nm,"i":776*u.nm}
    __key_band = { "mag_g": "g'mag", "mag_r": "r'mag", "mag_i": "i'mag"}


    def __init__(self):
        pass

    def getBand(self, band):

        key_band = ""
        if band in self.__key_band:
            key_band = self.__key_band[band]
        return key_band


    def query(self, **kwargs):
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"]) + "," + str(kwargs["dec"])
        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))


        radius_degree = CoordinateParser.getMinToDegree(radius)

        query = self.__simple_rec_query.format(self.__coordinates.ra.degree, self.__coordinates.dec.degree, radius_degree)
        respond = self.getTapRequest(query=query)


        result = AstroSource(self.__coordinates)
        if len(respond)>0:
            Utils.createPath(self.__save_path)

            catalog = CatalogSource("apass", "Search in radius " + str(radius))
            VOTableUtil.saveFromTable(respond, self.__save_path + "catalog.xml")
            catalog.addFile("apass", self.__save_path + "catalog.xml", "vo")

            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          respond["RAJ2000"], respond["DEJ2000"])



        return result





    def getCatalog(self,ra,dec,radius):

        self.__coordinates = str(ra) + "," + str(dec)
        radius = radius  # arcmin
        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))


        radius_degree = CoordinateParser.getMinToDegree(radius)

        query = self.__simple_rec_query.format(self.__coordinates.ra.degree, self.__coordinates.dec.degree, radius_degree)
        respond = self.getTapRequest(query=query)
        return respond


    def getTapRequest(self, query=""):

        tap = TAPManager()
        tap_url = self.__service_provider["tap_url"]
        tap.connect(url=tap_url)

        result=tap.sync_query(query=query)
        r=result.get_results()
        return r


    def contentUrl(self):
        pass
    def delete(self):
        pass
    def getType(self):
        pass
    def insert(self):
        pass
    def onCreate(self):
        pass
    def update(self):
        pass

if __name__ == "__main__":
    ap=APASSArchiveCP()

    cat=ap.getCatalog(207.7481064334807,-51.208771664529074,23.51855172386045)
    print(cat)