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

class DesAchiveCP(ContentProvider):
    __catalog_provider = "des"

    __service_provider = {
        "tap_url": "http://datalab.noao.edu/tap"
    }

    __simple_rec_query = "SELECT * FROM des_dr1.main WHERE((ra between {0} and {1}) and (dec between {2} and {3}))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/des/"

    __wavelenght = {"g": 475 * u.nm, "r": 635 * u.nm, "i": 775 * u.nm, "z": 925 * u.nm, "y": 1000 * u.nm}
    __key_band = {"mag_g": "mag_auto_g", "mag_r": "mag_auto_r", "mag_i": "mag_auto_i", "mag_z": "mag_auto_z",
                  "mag_y": "mag_auto_y", }

    def __init__(self):
        self.__noao = NOAOArchiveCP()

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

        self.__arcmin_range = radius
        self.__degrees_range = CoordinateParser.getMinToDegree(radius)
        result = AstroSource(self.__coordinates)
        table,df =self.getCatalogFromSIAS()

        # GET catalogs from SIAS DB, the object is selected from the near source using Pythagorean theorem
        if len(df.index)>0:
            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree, df["ra"].values, df["dec"].values)
            self.loss=loss
            self.setSummary(result,df.iloc[index])

            #table = Table.from_pandas(df)
            Utils.createPath(self.__save_path)
            catalog = CatalogSource("des", "Search in radius " + str(self.__arcmin_range))
            VOTableUtil.saveFromTable(table, self.__save_path + "catalog.xml")
            catalog.addFile("des", self.__save_path + "catalog.xml", "vo")

            result.addCatalog(catalog)

            #GET IMAGES From SIAS by ra dec radio in arcmin(the service recieve radii in degree)
            img=self.getImageFromSIAS()
            if img!=None:
                result.addImage(img)

        return result

    def setSummary(self,result,respond):

        result.addSummaryParams("ra", respond["ra"])
        result.addSummaryParams("dec", respond["dec"])
        result.addSummaryParams("distance", self.loss)

        des = WaveLenghtCover.des()

        data={"lambda": des["g"], "ab": str(respond["mag_auto_g"]),"err":str(respond["magerr_auto_g"])}

        result.addSummaryParams("mag_g", data)

        data = {"lambda": des["r"], "ab": str(respond["mag_auto_r"]), "err": str(respond["magerr_auto_r"])}
        result.addSummaryParams("mag_r", data)

        data = {"lambda": des["i"], "ab": str(respond["mag_auto_i"]), "err": str(respond["magerr_auto_i"])}
        result.addSummaryParams("mag_i", data)

        data = {"lambda": des["z"], "ab": str(respond["mag_auto_z"]), "err": str(respond["magerr_auto_z"])}
        result.addSummaryParams("mag_z", data)

        data = {"lambda": des["y"], "ab": str(respond["mag_auto_y"]), "err": str(respond["magerr_auto_y"])}
        result.addSummaryParams("mag_y", data)

        """
        result.addSummaryParams("magerr_g", str(respond["magerr_auto_g"]))
        result.addSummaryParams("magerr_r", str(respond["magerr_auto_r"]))
        result.addSummaryParams("magerr_i", str(respond["magerr_auto_i"]))
        result.addSummaryParams("magerr_z", str(respond["magerr_auto_z"]))
        result.addSummaryParams("magerr_y", str(respond["magerr_auto_y"]))
        """

        result.addSummaryParams("flux_g", str(respond["flux_auto_g"]))
        result.addSummaryParams("flux_r", str(respond["flux_auto_r"]))
        result.addSummaryParams("flux_i", str(respond["flux_auto_i"]))
        result.addSummaryParams("flux_z", str(respond["flux_auto_z"]))
        result.addSummaryParams("flux_y", str(respond["flux_auto_y"]))

        result.addSummaryParams("fluxerr_g", str(respond["fluxerr_auto_g"]))
        result.addSummaryParams("fluxerr_r", str(respond["fluxerr_auto_r"]))
        result.addSummaryParams("fluxerr_i", str(respond["fluxerr_auto_i"]))
        result.addSummaryParams("fluxerr_z", str(respond["fluxerr_auto_z"]))
        result.addSummaryParams("fluxerr_y", str(respond["fluxerr_auto_y"]))


    def getCatalog(self,ra,dec,radius):

        self.__coordinates = str(ra) + "," + str(dec)
        radius = radius  # arcmin
        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        self.__arcmin_range = radius
        self.__degrees_range = CoordinateParser.getMinToDegree(radius)
        result = AstroSource(self.__coordinates)
        df =self.getCatalogFromSIAS()
        try:
            table = Table.from_pandas(df)
        except AttributeError:
            table = []
        return table,df

    def getCatalogFromSIAS(self):

        df = self.__noao.desQuery(self.__coordinates.ra.degree, self.__coordinates.dec.degree, radius_arcmin=self.__arcmin_range)
        return df

    def getCatalagoFromTAP(self):

        degrees_range=self.__degrees_range
        ra_min = self.__coordinates.ra.degree - (degrees_range / 2)
        ra_max = self.__coordinates.ra.degree + (degrees_range / 2)
        dec_min = self.__coordinates.dec.degree - (degrees_range / 2)
        dec_max = self.__coordinates.dec.degree + (degrees_range / 2)
        query = self.__simple_rec_query.format(ra_min, ra_max, dec_min, dec_max)

        respond = self.getTapRequest(query=query)

        result = AstroSource(self.__coordinates)
        if len(respond)>0:
            self.setSummary(result,respond)


            """
            image = ImageSource(result.getSummary()["id"], self.__catalog_provider)

            bands = ["grz", "gz", "g"]
            mag_band = [result.getSummary()["mag_r"], result.getSummary()["mag_z"], result.getSummary()["mag_g"]]

            for index in range(len(bands)):
                band = bands[index]
                mag = mag_band[index]
                img = self.getImage(self.__coordinates, band)
                image.addCutout(img["jpeg"], 256, band, mag)
                image.addFile(band, img["fits"])

            result.addImage(image)
            """

        return result

    def getImageFromSIAS(self):
        img=self.__noao.download_deepest_image(self.__coordinates.ra.degree, self.__coordinates.dec.degree, radius_arcmin=self.__arcmin_range,path=self.__save_path)
        return img


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

    ra=75.757204
    dec=-39.762352
    des=DesAchiveCP()
    cat=des.getCatalog(ra,dec,2)
    print(cat)