
class HipsSkyMaps():

    __sky_map_url={
        "CFHT":"CDS/P/CFHTLS/W/Color/ugi"
    }




    @staticmethod
    def getMap(archive):
        self = HipsSkyMaps()
        if archive in self.__sky_map_url:
            return self.__sky_map_url[archive]
        else:
            return None