from frastro import ContentProvider
from astroquery.ukidss import Ukidss
import astropy.units as u
from frastro import CoordinateParser
from frastro import AstroSource
from frastro import ImageSource
from frastro import SpectraSource
from frastro import TAPManager
import numpy as np

from frastro import CatalogSource
from frastro import VOTableUtil
from frastro import Utils
from astropy.table import Table
from frastro import WaveLenghtCover

class UkidssArchiveCP(ContentProvider):
    __service_provider = {
        "tap_url": "http://tapvizier.u-strasbg.fr/TAPVizieR/tap"
    }

    __simple_rec_query = 'SELECT * FROM "{0}" as w WHERE ' + "1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {1},{2}, {3} ))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/ukidss/"

    __wavelenght = {"z": 883.07 * u.nm, "y": 1031.87 * u.nm, "j": 1251.02 * u.nm, "h": 1637.71 * u.nm, "k": 2208.28 * u.nm}

    def __init__(self,catalog_provider='Ukidss'):
        self.__catalog_provider=catalog_provider
        self.onCreate()

    def onCreate(self):

        all_databases = ("UKIDSSDR9PLUS", "UKIDSSDR8PLUS", "UKIDSSDR7PLUS",
                         "UKIDSSDR6PLUS", "UKIDSSDR5PLUS", "UKIDSSDR4PLUS",
                         "UKIDSSDR3PLUS", "UKIDSSDR2PLUS", "UKIDSSDR1PLUS",
                         "UKIDSSDR1", "UKIDSSEDRPLUS", "UKIDSSEDR", "UKIDSSSV",
                         "WFCAMCAL08B", "U09B8v20120403", "U09B8v20100414","UKIDSSDR10PLUS","UKIDSSDR11PLUSUDSONLY")

        ukidss_programmes_short = {'ALL':'all',
                                   'LAS': 101,
                                   'GPS': 102,
                                   'GCS': 103,
                                   'DXS': 104,
                                   'UDS': 105, }

        ukidss_programmes_long = {
                                'ALL Survey':'all',
                                'Large Area Survey': 101,
                                'Galactic Plane Survey': 102,
                                'Galactic Clusters Survey': 103,
                                'Deep Extragalactic Survey': 104,
                                'Ultra Deep Survey': 105}


        Ukidss.all_databases = all_databases
        Ukidss.ukidss_programmes_long=ukidss_programmes_long
        Ukidss.ukidss_programmes_short=ukidss_programmes_short


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

        large_suvey="II/319/las9"
        galactic_cluster_suvey = "II/319/gcs9"
        deep_suvey = "II/319/dxs9"
        galatic_plane_survey="II/316/gps6"

        radius_degree = CoordinateParser.getMinToDegree(radius)
        program_id="ALL"

        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        query = self.__simple_rec_query.format(large_suvey, self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                               radius_degree)
        respond = self.getTapRequest(query=query)

        result = AstroSource(self.__coordinates)
        if len(respond) <= 0:
            query = self.__simple_rec_query.format(galactic_cluster_suvey, self.__coordinates.ra.degree,
                                                   self.__coordinates.dec.degree, radius_degree)
            respond = self.getTapRequest(query=query)
            program_id = "GCS"

        if len(respond) <= 0:
            query = self.__simple_rec_query.format(deep_suvey, self.__coordinates.ra.degree,
                                                   self.__coordinates.dec.degree, radius_degree)
            respond = self.getTapRequest(query=query)
            program_id = "DXS"

        if len(respond) <= 0:
            query = self.__simple_rec_query.format(galatic_plane_survey, self.__coordinates.ra.degree,
                                                   self.__coordinates.dec.degree, radius_degree)
            respond = self.getTapRequest(query=query)
            program_id = "GPS"

        if len(respond)>0:


            Utils.createPath(self.__save_path)
            catalog = CatalogSource("ukidss", "Search in radius " + str(radius))
            VOTableUtil.saveFromTable(respond, self.__save_path + "catalog.xml")
            catalog.addFile("ukidss", self.__save_path + "catalog.xml", "vo")

            result.addCatalog(catalog)

            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          respond["RAJ2000"], respond["DEJ2000"])

            ukidss = WaveLenghtCover.ukidss()
            result.addSummaryParams("ra", respond["RAJ2000"][index])
            result.addSummaryParams("dec", respond["DEJ2000"][index])
            result.addSummaryParams("distance", loss)

            magnitudes={"Y":-99,"J":-99,"H":-99,"K":-99,"Z":-99}

            try:
                data={"lambda": ukidss["h"], "ab": str(respond["Hmag"][index]), "err": str(respond["e_Hmag"][index])}
                result.addSummaryParams("mag_h", data)
                magnitudes["H"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_h found")

            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["Jmag"][index]), "err": str(respond["e_Jmag"][index])}
                result.addSummaryParams("mag_j", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_j found")

            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["Jmag1"][index]), "err": str(respond["e_Jmag1"][index])}
                result.addSummaryParams("mag_j1", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_j1 found")

            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["Jmag2"][index]),"err": str(respond["e_Jmag2"][index])}
                result.addSummaryParams("mag_j2", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_j2 found")

            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["pJmag"][index]),"err": str(respond["e_pJmag"][index])}
                result.addSummaryParams("pmag_j", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No pmag_j found")

            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["pKmag"][index]),"err": str(respond["e_pKmag"][index])}
                result.addSummaryParams("pmag_k", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No pmag_k found")

            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["Kmag"][index]),"err": str(respond["e_Kmag"][index])}
                result.addSummaryParams("mag_k", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_k found")

            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["Kmag1"][index]), "err": str(respond["e_Kmag1"][index])}
                result.addSummaryParams("mag_k1", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_k1 found")
            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["Kmag2"][index]),"err": str(respond["e_Kmag2"][index])}
                result.addSummaryParams("mag_k2", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_k2 found")

            try:
                data = {"lambda": ukidss["y"], "ab": str(respond["Ymag"][index]),"err": str(respond["e_Ymag"][index])}
                result.addSummaryParams("mag_y", data)
                magnitudes["Y"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_y found")
            try:
                data = {"lambda": ukidss["z"], "ab": str(respond["Zmag"][index]), "err": str(respond["e_Zmag"][index])}
                result.addSummaryParams("mag_z", data)
                magnitudes["Z"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_z found")

            images = Ukidss.get_image_list(self.__coordinates,waveband='all', frame_type='stack',
                                            image_width=radius * u.arcmin,
                                            database="UKIDSSDR10PLUS",programme_id=program_id)

            display_urls = [link.replace("getFImage", "getImage")for link in images]
            image_urls = [link.replace("getFImage", "getJImage") for link in images]

            image = ImageSource(result.getSummary()["id"], self.__catalog_provider)


            for index in range(len(images)):
                img= image_urls[index]
                band = str(img[img.find("band")+5:img.find("band")+6]).replace(" ","")
                mag=magnitudes[band]
                image.addCutout(image_urls[index],300, band, magnitude=mag,link=display_urls[index])

                local_path = self.__save_path + "ukidss_band_" + band + ".fits"
                image.addFile(name="band-" + band, url=images[index], download=True,
                              local_path=local_path, uncompress=False, thumbnail=True,external=False)

                #image.addFile("band-" + band,images[index] , "fits")

            result.addImage(image)

        return result


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
    def update(self):
        pass