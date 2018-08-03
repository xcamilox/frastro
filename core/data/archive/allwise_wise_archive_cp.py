from frastro import ContentProvider
from frastro import TAPManager
from frastro import AstroSource
from frastro import CoordinateParser
import astropy.units as u
from frastro import ImageSource
from frastro import CatalogSource
import numpy as np
from frastro import Utils
from frastro import VOTableUtil,WaveLenghtCover

class WiseAllWiseArchiveCP(ContentProvider):
    ##WIse and allwise from IRSA

    # __service_provider = {
    #     "tap_url": "http://irsa.ipac.caltech.edu/TAP/",
    #     "wise_band2": "allsky_2band_p1bs_psd",
    #     "wise_band3": "allsky_3band_p1bs_psd",
    #     "wise_band4": "allsky_4band_p1bs_psd",
    #     "all_wise":"allwise_p3as_psd"
    # }
    #Catalog from vizier

    __service_provider = {
        "tap_url": "http://tapvizier.u-strasbg.fr/TAPVizieR/tap",
        "wise": "II/311/wise",
        "all_wise":"II/328/allwise"
    }

    __simple_rec_query = 'SELECT * FROM "{0}" WHERE ' + "1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {1},{2}, {3} ))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/wise/"

    __wavelenght = {"w1": 3352.60 * u.nm, "w2": 4602.80 * u.nm, "w3": 11560.80 * u.nm, "w4": 22088.30 * u.nm}

    def __init__(self, catalog_provider='wiseAllwise'):
        self.__catalog_provider = catalog_provider

    def onCreate(self):
        pass

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

        all_wise=self.getData(self.__service_provider["all_wise"]);
        wise = self.getData(self.__service_provider["wise"]);


        result = AstroSource(self.__coordinates)
        catalog=None
        allWiseWavelenght = WaveLenghtCover.wise()

        if len(all_wise) > 0:
            Utils.createPath(self.__save_path)
            catalog = CatalogSource("allwise","Search in radius "+str(radius))
            VOTableUtil.saveFromTable(all_wise, self.__save_path + "allwise_catalog.xml")
            catalog.addFile("allwise",self.__save_path + "allwise_catalog.xml","vo")

            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          all_wise["RAJ2000"], all_wise["DEJ2000"])



            result.addSummaryParams("allwise_ra", all_wise["RAJ2000"][index])
            result.addSummaryParams("allwise_dec", all_wise["DEJ2000"][index])
            result.addSummaryParams("distance", loss)

            data = {"lambda": allWiseWavelenght["w1"], "ab": str(all_wise["W1mag"][index]), "err": str(all_wise["e_W1mag"][index])}

            result.addSummaryParams("allwise_w1mag", data)

            data = {"lambda": allWiseWavelenght["w2"], "ab": str(all_wise["W2mag"][index]),
                    "err": str(all_wise["e_W2mag"][index])}
            result.addSummaryParams("allwise_w2mag", data)

            data = {"lambda": allWiseWavelenght["w3"], "ab": str(all_wise["W3mag"][index]),
                    "err": str(all_wise["e_W3mag"][index])}

            result.addSummaryParams("allwise_w3mag", data)

            data = {"lambda": allWiseWavelenght["w4"], "ab": str(all_wise["W4mag"][index]),
                    "err": str(all_wise["e_W4mag"][index])}

            result.addSummaryParams("allwise_w4mag", data)

            data = {"lambda": allWiseWavelenght["j"], "ab": str(all_wise["Jmag"][index]),
                    "err": str(all_wise["e_Jmag"][index])}


            result.addSummaryParams("allwise_Jmag", data)

            data = {"lambda": allWiseWavelenght["h"], "ab": str(all_wise["Hmag"][index]),
                    "err": str(all_wise["e_Hmag"][index])}
            result.addSummaryParams("allwise_Hmag", data)

            data = {"lambda": allWiseWavelenght["k"], "ab": str(all_wise["Kmag"][index]),
                    "err": str(all_wise["e_Kmag"][index])}
            result.addSummaryParams("allwise_Kmag", data)

            result.addCatalog(catalog)




        if len(wise) > 0:
            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          wise["RAJ2000"], wise["DEJ2000"])

            if catalog is None:
                Utils.createPath(self.__save_path)
                catalog = CatalogSource("wise","Search in radius(arcmin) "+str(radius))


            VOTableUtil.saveFromTable(wise, self.__save_path + "wise_catalog.xml")
            catalog.addFile("wise", self.__save_path + "wise_catalog.xml", "vo")



            result.addSummaryParams("wise_ra", wise["RAJ2000"][index])
            result.addSummaryParams("wise_dec", wise["DEJ2000"][index])
            result.addSummaryParams("wise_distance", loss)

            data = {"lambda": allWiseWavelenght["w1"], "ab": str(wise["W1mag"][index]), "err": str(wise["e_W1mag"][index])}

            result.addSummaryParams("wise_w1mag", data)

            data = {"lambda": allWiseWavelenght["w2"], "ab": str(wise["W2mag"][index]),
                    "err": str(wise["e_W2mag"][index])}

            result.addSummaryParams("wise_w2mag", data)

            data = {"lambda": allWiseWavelenght["w3"], "ab": str(wise["W3mag"][index]),
                    "err": str(wise["e_W3mag"][index])}

            result.addSummaryParams("wise_w3mag", data)

            data = {"lambda": allWiseWavelenght["w4"], "ab": str(wise["W4mag"][index]),
                    "err": str(wise["e_W4mag"][index])}

            result.addSummaryParams("wise_w4mag", data)

            data = {"lambda": allWiseWavelenght["j"], "ab": str(wise["Jmag"][index]),
                    "err": str(wise["e_Jmag"][index])}

            result.addSummaryParams("wise_Jmag", data)

            data = {"lambda": allWiseWavelenght["h"], "ab": str(wise["Hmag"][index]),
                    "err": str(wise["e_Hmag"][index])}

            result.addSummaryParams("wise_Hmag", data)

            data = {"lambda": allWiseWavelenght["k"], "ab": str(wise["Kmag"][index]),
                    "err": str(wise["e_Kmag"][index])}

            result.addSummaryParams("wise_Kmag", data)

            result.addCatalog(catalog)

        return result



    def getData(self,table):

        query = self.__simple_rec_query.format(table, self.__coordinates.ra.degree, self.__coordinates.dec.degree,self.__radius_degree)
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
    def update(self):
        pass
