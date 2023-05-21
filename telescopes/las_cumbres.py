class LasCumbres():
    __params = {
        "KB98":{
            "gain": 1.6,#e-/ADU
            "readnoise": 14.5,#e-
            "readtime":13.0, #s
            "darkcurrent":0.03, #e-/s
            "pixelscale":0.571, #0.2637 / 0.2626 arcsec/pixel,
            "telescope": "lascumbres"
        },
        "KB27": {
            "gain": 1.6,  # e-/ADU
            "readnoise": 14.5,  # e-
            "readtime": 13.0,  # s
            "darkcurrent": 0.03,  # e-/s
            "pixelscale": 0.571,  # 0.2637 / 0.2626 arcsec/pixel,
            "telescope": "lascumbres"
        },
        "FA12":{
            "gain": 1.0,  # e-/ADU
            "readnoise": 7.6,  # e-
            "readtime": 60.0,  # s
            "darkcurrent": 0.00,  # e-/s
            "pixelscale": 0.38,  # 0.2637 / 0.2626 arcsec/pixel,
            "telescope": "lascumbres"
        }

    }

    __key_band = {"up":"mag_u","gp":"mag_g","rp":"mag_r","ip":"mag_i","sdss-r":"mag_r"}

    def getBand(self,band):
        band = str(band).lower()
        if band in self.__key_band.keys():
            return self.__key_band[band]
        elif band in self.__key_band.values():
            for mag, filter in self.__key_band.items():
                if filter == band:
                    return mag
        else:
            return None


    @staticmethod
    def getParam(key):
        lv = LasCumbres()
        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError()

    @staticmethod
    def getParams(instrument="all"):
        lv = LasCumbres()
        if instrument == "all":
            return lv.__params
        else:
            if instrument.upper() in lv.__params:
                return lv.__params[instrument.upper()]
            else:
                raise KeyError("No found " + instrument)