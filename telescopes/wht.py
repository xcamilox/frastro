from astropy import units as u
class WHTTelescope():

    __instrument_default = "LIRIS"

    __params={
        "LIRIS":{
            "gain":4.0, #e-/ADU
            "readnoise": 17.0, #15-,
            "readtime": 0.9,  # 1.45 S(Non-destructive)
            "pixelscale": 0.25, #0.184 arcsec/pixel (unbinned)
            "telescope":"wht",
            "fov":[4.27*u.arcmin,4.27*u.arcmin]
        }
    }

    __key_band = {"mag_y": "y", "mag_j": "j", "mag_h": "h", "mag_ks": "ks", "mag_z": "z"}

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

        lv=WHTTelescope()

        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError("No found "+key)

    @staticmethod
    def getParams(instrument="all"):
        lv = WHTTelescope()
        if instrument == "all":
            return lv.__params
        else:
            instrument= instrument.replace(":","")
            if instrument in lv.__params:
                return lv.__params[instrument]
            else:
                raise KeyError("No found "+instrument)