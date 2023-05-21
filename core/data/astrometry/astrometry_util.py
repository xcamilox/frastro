from frastro import FITSFile, Config, LOGUtil, ImageUtils, TAPManager, CoordinateParser


class AstrometryUtil():

    def __init__(self):
        self.query = Config.getQuery("gaia-astrometry")
        self.service_url = Config.getServices("vizier")

    def getCatalog(self, ra, dec, radius=5):
        self.query

        self.__coordinates = str(ra) + "," + str(dec)

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.query=self.query.format(self.__coordinates.ra.degree, self.__coordinates.dec.degree, radius)
        return self.getTable()

    def getTable(self):
        tap = TAPManager()
        tap_url = self.service_url
        tap.connect(url=tap_url)
        print(self.query)
        result = tap.sync_query(query=self.query)
        r = result.get_results()
        print("results", len(r))
        return r

if __name__ == "__main__":
    astrometry=AstrometryUtil()
    table=astrometry.getCatalog(162.78920833333333, 44.65236388888888, 0.5)
    print(table)
