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
from astroquery.irsa import Irsa

class TwoMassArchiveCP(ContentProvider):
    __wavelenght = {"j": 1235.00 * u.nm, "h":1662.00 * u.nm, "ks":2159.00 * u.nm}

    __key_band = {"mag_j": "j_m", "mag_h": "h_m", "mag_ks": "k_m"}

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

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]

        r=self.getCatalog(self.__coordinates.ra.degree,self.__coordinates.dec.degree,radius)
        return r

    def getCatalog(self,ra,dec,radio):

        self.__coordinates = str(ra) + "," + str(dec)
        radius = radio  # arcmin
        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)


        """
        fp_psc: 2MASS All-Sky Point Source Catalog (PSC)
        fp_xsc: 2MASS All-Sky Extended Source Catalog (XSC)
        """

        #search on point source catalog
        r = Irsa.query_region(self.__coordinates, catalog='fp_psc',spatial="Cone",
                          radius=radius * u.arcmin)
        if len(r)<=0:
            #search in extended source catalog
            r = Irsa.query_region(self.__coordinates, catalog='fp_xsc', spatial="Cone",
                                  radius=radius * u.arcmin)

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



if __name__=="__main__":
    twomass=TwoMassArchiveCP()
    ra=185.16966
    dec=8.710594
    radii=10
    r=twomass.query(ra=ra,dec=dec,radius=radii)
    print(r)