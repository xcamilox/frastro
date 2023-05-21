class DecalsTelescope():
    __params = {
        "DECAM":{
            "gain": 4.0,#e-/ADU
            "readnoise": 6.0,#e-
            "readtime":20.0, #s
            "darkcurrent":0.17, #e-/s
            "pixelscale":0.26 #0.2637 / 0.2626 arcsec/pixel
        }
    }

    @staticmethod
    def getParam(key):
        lv = DecalsTelescope()
        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError()

    @staticmethod
    def getParams(instrument="all"):
        lv = DecalsTelescope()
        if instrument == "all":
            return lv.__params
        else:
            if instrument.upper() in lv.__params:
                return lv.__params[instrument.upper()]
            else:
                raise KeyError("No found " + instrument)