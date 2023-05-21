from astropy.io import fits
import numpy as np
from core.fits.FITSFIle import FITSFile


class VLTelescope():
    __instrument_default = "FORS2"

    __params = {
        "FORS2": {
            "gain": 2.5,  # e-/ADU
            "readnoise": 5,  # 15-,
            "readtime": 28.0,  # 1.45 S(Non-destructive)
            "pixelscale": 0.25,  # 0.125 arcsec/pixel (unbinned)
            "telescope": "vlt"
        }
    }

    __key_band = {"mag_r": "r_special", "mag_z": "z_gunn"}

    def getBand(self, band):
        band = str(band).lower()
        if band in self.__key_band.keys():
            return self.__key_band[band]
        elif band in self.__key_band.values():
            for mag, filter in self.__key_band.items():
                if filter == band:
                    return mag
        else:
            return None

    # telescope="ESO-VLT-U1"
    # filter="R_SPECIAL"
    # instrument = "FORS2"
    @staticmethod
    def getParam(key):
        lv = VLTelescope()
        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError()

    @staticmethod
    def getParams(instrument="all"):
        lv = VLTelescope()
        if instrument == "all":
            return lv.__params
        else:
            if instrument.upper() in lv.__params:
                return lv.__params[instrument.upper()]
            else:
                raise KeyError("No found " + instrument)


