from frastro import ContentProvider
from frastro import TAPManager
from frastro import AstroSource
from frastro import CoordinateParser
import astropy.units as u
from frastro import ImageSource
import numpy as np
from frastro import NOAOArchiveCP
from frastro import CatalogSource
from frastro import VOTableUtil
from frastro import Utils, WaveLenghtCover
from astropy.table import Table

class GAIAAchiveCP(ContentProvider):
    __catalog_provider = "des"

    __service_provider = {
        "tap_url": "http://tapvizier.u-strasbg.fr/TAPVizieR/tap"
    }

    __simple_rec_query = 'SELECT * FROM "I/345/gaia2" WHERE ' + "1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {0},{1}, {2} ))"
    __save_tmp_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/gaia/"

    __wavelenght = {"g": 673.54 * u.nm}
    __key_band = {"mag_g": "phot_g_mean_mag"}

    def __init__(self):
        pass

    def getBand(self, band):

        key_band = ""
        if band in self.__key_band:
            key_band = self.__key_band[band]
        return key_band

    def query(self, **kwargs):
        # check for data release to use
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"]) + "," + str(kwargs["dec"])
        else:
            raise ValueError("Not valid coordinates found. Used coordinates key or ra dec keys")

        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)

        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        self.__radius_degree = CoordinateParser.getMinToDegree(radius)

        query = self.__simple_rec_query.format(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                               self.__radius_degree)
        respond = self.getTapRequest(query=query)

        return respond



    #return astropy table
    def getCatalog(self,ra,dec,radius):

        self.__coordinates = str(ra) + "," + str(dec)
        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.__radius_degree = CoordinateParser.getMinToDegree(radius)
        query = self.__simple_rec_query.format(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                               self.__radius_degree)
        respond = self.getTapRequest(query=query)


        return respond


    def getTapRequest(self, query=""):

        tap = TAPManager()
        tap_url = self.__service_provider["tap_url"]
        tap.connect(url=tap_url)

        result=tap.async_query(query=query)
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
