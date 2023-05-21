from astropy.io import fits
import numpy as np
from core.fits.FITSFIle import FITSFile


class HSTelescope():

    __instrument_default = "WFC3"

    __params={
        "WFC3":{
            "gain":2.5, #e-/ADU
            "readnoise": 21, #15-,
            "readtime": 8.0,  # 1.45 S(Non-destructive)
            "pixelscale": 0.13, #0.184 arcsec/pixel (unbinned)
            "telescope":"hst"
        }
    }

    __key_band = {"mag_y": "y", "mag_j": "j", "mag_h": "h", "mag_ks": "ks", "mag_z": "z"}

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

    @staticmethod
    def getParam(key):
        lv = HSTelescope()
        if key in lv.__params:
            return lv.__params[key]
        else:
            raise KeyError()

    @staticmethod
    def getParams(instrument="all"):
        lv = HSTelescope()
        if instrument == "all":
            return lv.__params
        else:
            if instrument.upper() in lv.__params:
                return lv.__params[instrument.upper()]
            else:
                raise KeyError("No found " + instrument)

    """
    
    http://www.stsci.edu/hst/wfc3/phot_zp_lbn
    Photometric Systems
        The STmag and ABmag systems define an equivalent flux density for a source,
        corresponding to the flux density of a source of predefined spectral shape that
        would produce the observed count rate, and convert this equivalent flux to a magnitude.
        The conversion is chosen so that the magnitude in V corresponds roughly to that in the Johnson system.
        
        In the STmag system, the flux density is expressed per unit wavelength, and the reference spectrum is flat in Fλ.
        An object with Fλ = 3.63 x 10-9 erg cm-2 s-1 Å-1 will have STmag=0 in every filter, and its zero point is 21.10.
        
        STmag = -2.5 log Fλ -21.10
        
        In the ABmag system, the flux density is expressed per unit frequency, and the reference spectrum is flat in Fν.
        Its zero point is 48.6.
        
        ABmag = -2.5 log Fν - 48.6
        
        ABmag = STmag - 5 log (PHOTPLAM) + 18.6921
        
        where Fν is expressed in erg cm-2 s-1 Hz-1, and Fλ in erg cm-2 s-1 Å-1. An object with Fν = 3.63 x 10-20 erg cm-2 s-1 Hz-1
        will have magnitude AB =0 in every filter.
        
        Formally, the VEGAmag system is defined such that  Vega (Alpha Lyra) by definition has magnitude 0 at all wavelengths.
        The magnitude of a star with flux F relative to Vega is
        
        mvega= -2.5 log10 (F/Fvega)
        
        where Fvega is the absolute CALSPEC flux of Vega; for photometry the fluxes must be averaged over the band pass.
        See Bohlin 2014 (AJ, 147,127, "Hubble Space Telescope CALSPEC Flux Standards: Sirius and Vega") for the equations
        that define the average flux.
        
         
    """

    @staticmethod
    def getZeroPoint(file_path):
        zeropoint={}
        file = fits.open(file_path)
        header = file[0].header

        "STmag = -2.5 log Fλ -21.10"
        STmag=-2.5*np.log10(header["PHOTFLAM"])+header["PHOTZPT"]

        "ABmag = STmag - 5 log (PHOTPLAM) + 18.6921"
        ABmag = STmag - 5 * np.log10(header["PHOTPLAM"]) + 18.6921

        zeropoint["STmag"] = STmag
        zeropoint["ABmag"] = ABmag

        return zeropoint

    @staticmethod
    def getSingleImage(file_path,output):
        file=fits.open(file_path)
        main_header=file[0].header
        header =file[1].header
        science_data=file[1].data


        header.append(('EXPTIME', main_header['EXPTIME'], main_header.comments['EXPTIME']))
        header.append(('TARGNAME', main_header['TARGNAME'], main_header.comments['TARGNAME']))
        header.append(('TELESCOP', main_header['TELESCOP'], main_header.comments['TELESCOP']))
        header.append(('FILTER', main_header['FILTER'], main_header.comments['FILTER']))
        header.append(('INSTRUME', main_header['INSTRUME'], main_header.comments['INSTRUME']))
        header.append(('PHOTMODE', main_header['PHOTMODE'], main_header.comments['PHOTMODE']))
        header.append(('PHOTFLAM', main_header['PHOTFLAM'], main_header.comments['PHOTFLAM']))
        header.append(('PHOTFNU', main_header['PHOTFNU'], main_header.comments['PHOTFNU']))
        header.append(('PHOTZPT', main_header['PHOTZPT'], main_header.comments['PHOTZPT']))
        header.append(('PHOTPLAM', main_header['PHOTPLAM'], main_header.comments['PHOTPLAM']))
        header.append(('PHOTBW', main_header['PHOTBW'], main_header.comments['PHOTBW']))

        FITSFile.saveFile(output, science_data, header)

        return output


if __name__ == "__main__":
    path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/hst/SLACSJ1051+4439_12210_55695_F814W_1.fits"
    HSTelescope.getSingleImage(path,path[:-5]+"single.fits")