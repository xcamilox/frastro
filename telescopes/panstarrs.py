class PanstarrsTelescope():

    __instrument_default = "GPC1"

    __params={
        "GPC1":{
            "gain":1.011, #e-/ADU
            "readnoise": 14.167, #15-,
            "readtime": 8.0,  # 1.45 S(Non-destructive)
            "pixelscale": 0.258, #0.184 arcsec/pixel (unbinned)
            "telescope":"panstarrs"
        }
    }

    @staticmethod
    def getParam(key):
        lv = PanstarrsTelescope()
        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError()

    @staticmethod
    def getParams(instrument="all"):
        lv = PanstarrsTelescope()
        if instrument == "all":
            return lv.__params
        else:
            if instrument.upper() in lv.__params:
                return lv.__params[instrument.upper()]
            else:
                raise KeyError("No found " + instrument)
