class GTCTelescope():

    __instrument_default = "OSIRIS"

    __params={
        "OSIRIS":{
            "gain":0.95, #e-/ADU
            "readnoise": 4.5, #15-,
            "readtime": 21.0,  # 1.45 S(Non-destructive)
            "pixelscale": 0.127*2, #0.184 arcsec/pixel (unbinned)
            "telescope":"gtc"
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

        lv=GTCTelescope()

        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError("No found "+key)

    @staticmethod
    def getParams(instrument="all"):
        lv = GTCTelescope()
        if instrument == "all":
            return lv.__params
        else:
            instrument= instrument.replace(":","")
            if instrument in lv.__params:
                return lv.__params[instrument]
            else:
                raise KeyError("No found "+instrument)