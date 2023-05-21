from frastro import ContentProvider
from astroquery.mast import Mast,Observations,Catalogs
from astroquery.esa_hubble import ESAHubble
from frastro import CoordinateParser
from frastro import AstroSource
from frastro import ImageSource
from frastro import CatalogSource
from frastro import SpectraSource
import astropy.units as u
from frastro import VOTableUtil
from frastro import Utils, WaveLenghtCover
import requests

class MASTArchiveCP(ContentProvider):

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


        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)

        obsTable = Observations.query_region(self.__coordinates,radius=radius*u.arcmin)

        print(obsTable[0:10])


    def getCatalog(self,ra,dec,radius,catalog):

        self.__coordinates = str(ra) + "" + str(dec)
        self.__radius_degree = CoordinateParser.getMinToDegree(radius)

        #respond=Catalogs.query_object(self.__coordinates,radius=self.__radius_degree, catalog=catalog)

        respond = ESAHubble.cone_search(self.__coordinates, 7, "cone_search_m31_5.vot")

        return str(respond)

    def missionList(self):
        return Observations.list_missions()




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
    mast=MASTArchiveCP()
    list_obs=mast.missionList()
    print(list_obs)
    ra="14:16:38.54"
    dec="+51:55:25.5"
    mission="HST"
    catalog=mast.getCatalog(ra,dec,5,mission)
    print(catalog)