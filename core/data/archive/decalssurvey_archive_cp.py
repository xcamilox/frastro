from frastro import ContentProvider
from frastro import TAPManager
from frastro import AstroSource
from frastro import CoordinateParser
import astropy.units as u
from frastro import ImageSource
from frastro import CatalogSource
from frastro import VOTableUtil
from frastro import Utils,WaveLenghtCover

import numpy as np

class DecalsSurveyArchiveCP(ContentProvider):
    __service_provider = {
        "tap_url": "http://datalab.noao.edu/tap"
    }
    __catalog_provider="decals"
    __simple_rec_query= "SELECT * FROM ls_dr6.tractor WHERE((ra between {0} and {1}) and (dec between {2} and {3}))"
    __simple_rec_query_dr7 = "SELECT * FROM ls_dr7.tractor WHERE((ra between {0} and {1}) and (dec between {2} and {3}))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/decals/"

    __wavelenght = {"g": 475 * u.nm, "r": 635 * u.nm, "i": 775 * u.nm, "z": 925 * u.nm}
    __key_band = {"mag_g": "mag_g", "mag_r": "mag_r", "mag_i": "mag_i", "mag_z": "mag_z"}

    __data_release = 6
    #flux in nanomagies

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
        self.__save_path=self.__save_path.format(str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))



        degrees_range=CoordinateParser.getMinToDegree(radius)

        ra_min = self.__coordinates.ra.degree - (degrees_range/2)
        ra_max = self.__coordinates.ra.degree + (degrees_range/2)
        dec_min = self.__coordinates.dec.degree - (degrees_range/2)
        dec_max = self.__coordinates.dec.degree + (degrees_range/2)
        query = self.__simple_rec_query_dr7.format(ra_min,ra_max,dec_min,dec_max)
        respond=self.getTapRequest(query=query)
        result = AstroSource(self.__coordinates)

        if len(respond)<=0:
            query = self.__simple_rec_query.format(ra_min, ra_max, dec_min, dec_max)
            respond = self.getTapRequest(query=query)

        if len(respond)>0:

            Utils.createPath(self.__save_path)
            catalog = CatalogSource("decals", "Search in radius " + str(radius))
            VOTableUtil.saveFromTable(respond, self.__save_path + "catalog.xml")
            catalog.addFile("decals", self.__save_path + "catalog.xml", "vo")

            result.addCatalog(catalog)

            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          respond["ra"], respond["dec"])
            result.addSummaryParams("ra", respond["ra"][index])
            result.addSummaryParams("dec", respond["dec"][index])
            result.addSummaryParams("distance", loss)

            decals = WaveLenghtCover.decals()

            data={"lambda": decals["g"], "ab": str(respond["mag_g"][index]),"err":-99}
            result.addSummaryParams("mag_g", data)

            data = {"lambda": decals["r"], "ab": str(respond["mag_r"][index]), "err": -99}
            result.addSummaryParams("mag_r", data)

            data = {"lambda": decals["z"], "ab": str(respond["mag_z"][index]), "err": -99}
            result.addSummaryParams("mag_z", data)


            result.addSummaryParams("flux_g", str(respond["flux_g"][index]))
            result.addSummaryParams("flux_r", str(respond["flux_r"][index]))
            result.addSummaryParams("flux_z", str(respond["flux_z"][index]))
            result.addSummaryParams("type", str(respond["type"][index]))

            image = ImageSource(result.getSummary()["id"], self.__catalog_provider)

            bands = ["grz", "gz", "g"]
            mag_band = [result.getSummary()["mag_r"], result.getSummary()["mag_z"], result.getSummary()["mag_g"]]

            for index in range(len(bands)):
                band=bands[index]
                mag=str(mag_band[index]["ab"])+"+/-"+str(mag_band[index]["err"])

                img=self.getImage(self.__coordinates,band)
                link = img["link"]
                image.addCutout(img["jpeg"],256,band,mag,link)
                #image.addFile(band,img["fits"])

                local_path = self.__save_path + "decals_band_" + band + ".fits"
                image.addFile(band, img["fits"], "fits", download=True, local_path=local_path,
                              uncompress=False,thumbnail=False,external=False)

            result.addImage(image)


        return result



    def getCatalog(self,ra,dec,radius):
        self.__coordinates = str(ra) + "," + str(dec)
        radius = radius  # arcmin

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        degrees_range = CoordinateParser.getMinToDegree(radius)

        ra_min = self.__coordinates.ra.degree - (degrees_range / 2)
        ra_max = self.__coordinates.ra.degree + (degrees_range / 2)
        dec_min = self.__coordinates.dec.degree - (degrees_range / 2)
        dec_max = self.__coordinates.dec.degree + (degrees_range / 2)
        query = self.__simple_rec_query_dr7.format(ra_min, ra_max, dec_min, dec_max)
        respond = self.getTapRequest(query=query)

        return respond


    def getTapRequest(self, query=""):

        tap = TAPManager()
        tap_url = self.__service_provider["tap_url"]
        tap.connect(url=tap_url)

        result=tap.sync_query(query=query)
        r=result.get_results()
        return r

    def getImage(self,coordinates,band):

        grz="http://legacysurvey.org/viewer/{0}-cutout?ra={1}&dec={2}&layer=mzls+bass-dr"+self.__data_release+"&pixscale=0.27&bands={3}"

        display="http://legacysurvey.org/viewer?ra={0}&dec={1}&zoom=15&layer=mzls+bass-dr"+self.__data_release

        image={'fits':grz.format("fits",coordinates.ra.degree,coordinates.dec.degree,band),'jpeg':grz.format("jpeg",coordinates.ra.degree,coordinates.dec.degree,band),'link':display.format(coordinates.ra.degree,coordinates.dec.degree)}

        return image



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
