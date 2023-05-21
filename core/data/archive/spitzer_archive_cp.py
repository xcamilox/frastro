from frastro import ContentProvider
import astropy.units as u
from frastro import CoordinateParser
from frastro import AstroSource
from frastro import ImageSource
from frastro import SpectraSource
from frastro import TAPManager
from astroquery import sha
import numpy as np

from frastro import CatalogSource
from frastro import VOTableUtil
from frastro import Utils
from astropy.table import Table

class SpitzerArcvhiveCP(ContentProvider):
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/spitzer/"

    __wavelenght = {"irac1": 3557.26 * u.nm, "irac2": 4504.93 * u.nm, "irac3": 5738.57 * u.nm, "irac4": 7927.37 * u.nm, "mips24": 23843.31 * u.nm,"mips70": 72555.53 * u.nm,"mips160": 156962.71 * u.nm}

    def __init__(self):
        pass
    def query(self, **kwargs):
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"]) + "," + str(kwargs["dec"])
        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]
        radius_degree = CoordinateParser.getMinToDegree(radius)
        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)

        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        try:
            result=sha.query(self.__coordinates, size=radius_degree)
        except IndexError:
            result = []
        except ValueError:
            result=[]
        respond = AstroSource(self.__coordinates)
        if len(result)>0:
            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          result["ra"], result["dec"])
            row = result[index]


            Utils.createPath(self.__save_path)
            catalog = CatalogSource("spitzer", "Search in radius " + str(radius))
            VOTableUtil.saveFromTable(result, self.__save_path + "catalog.xml")
            catalog.addFile("spitzer", self.__save_path + "catalog.xml", "vo")

            respond.addCatalog(catalog)

            respond.addSummaryParams("ra", row["ra"])
            respond.addSummaryParams("dec", row["dec"])
            respond.addSummaryParams("distance", loss)
            respond.addSummaryParams("filetype",row["filetype"])
            respond.addSummaryParams("modedisplayname", row["modedisplayname"])
            respond.addSummaryParams("wavelength", row["wavelength"])
            respond.addSummaryParams("minwavelength", row["minwavelength"])
            respond.addSummaryParams("maxwavelength", row["maxwavelength"])

        return respond


    def onCreate(self):
        pass
    def getTapRequest(self, query=""):
        tap = TAPManager()
        tap_url = self.__service_provider["tap_url"]
        tap.connect(url=tap_url)

        result = tap.sync_query(query=query)
        r = result.get_results()
        return r

    def contentUrl(self):
        pass

    def delete(self):
        pass

    def getType(self):
        pass

    def insert(self):
        pass

    def update(self):
        pass