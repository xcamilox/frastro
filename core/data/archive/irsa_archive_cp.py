import requests
import re
from frastro import ContentProvider
from astropy.io.votable import parse_single_table
from frastro.core.utils.html_parse import HtmlParser
from bs4.element import Tag
from frastro import ImageSource
from frastro import AstroSource
from frastro import CoordinateParser
from astroquery.irsa import Irsa
import astropy.units as u

class IRSAArchiveCP(ContentProvider):
    __service_provider = {
        "finding_chart": "http://irsa.ipac.caltech.edu/cgi-bin/bgTools/nph-bgExec"
    }

    __coordinates = 000


    def __init__(self,catalog_provider='IRSA'):
        self.__catalog_provider=catalog_provider

        self.onCreate()

    def onCreate(self):
        pass

    def query(self, **kwargs):


        # check for data release to use
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"])+","+str(kwargs["dec"])
        else:
            raise ValueError("Not valid coordinates found. Used coordinates key or ra dec keys")

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]

        params = {
            "bgApp": "/FinderChart/nph-finder",
            "romeserver": "ROMEDEV",
            "srchsize": 12.0,
            "outsize": 200,
            "colortbl": 1,
            "nthread": 10,
            "markercolor": "red",
            "markersize": 10,
            "mode": "cgi",
            "outtype": "single",
            "locstr": "20.48371,0.4223",
            "subsetsize": 5.0,
            "survey": "sdss",
            "survey": "dss",
            "survey": "2mass",
            "markervis_shrunk": "true"
        }
        print(Irsa.list_catalogs())
        r=Irsa.query_region(self.__coordinates, catalog='fp_psc',
                          radius=radius*u.arcmin)


        return r





    def irsa_request(self, url, params):
        r = requests.post(url, data=params)
        return r

    def contentUrl(self):
        pass
    def delete(self, **kwargs):
        pass
    def getType(self):
        pass
    def insert(self, **kwargs):
        pass
    def update(self, **kwargs):
        pass