class LivepoolTelescope():

    __instrument_default = "IO:O"

    __params={
        "IO:O":{
            "gain":1.5, #e-/ADU
            "readnoise": 15.0, #15-,
            "readtime": 20.0,  # 1.45 S(Non-destructive)
            "pixelscale": 0.15*2, #0.184 arcsec/pixel (unbinned)
            "telescope":"liverpool"
        },
        "IOO": {
            "gain": 1.5,  # e-/ADU
            "readnoise": 15.0,  # 15-,
            "readtime": 20.0,  # 1.45 S(Non-destructive)
            "pixelscale": 0.15 * 2,  # 0.184 arcsec/pixel (unbinned)
            "telescope": "liverpool"
        }
    }

    __key_band = {"mag_u": "sdss-u", "mag_g": "sdss-g", "mag_r": "sdss-r", "mag_i": "sdss-i", "mag_z": "sdss-z"}

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

        lv=LivepoolTelescope()

        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError("No found "+key)

    @staticmethod
    def getParams(instrument="all"):
        lv = LivepoolTelescope()
        if instrument == "all":
            return lv.__params
        else:
            instrument= instrument.replace(":","")
            if instrument in lv.__params:
                return lv.__params[instrument]
            else:
                raise KeyError("No found "+instrument)