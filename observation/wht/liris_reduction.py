from pyraf import iraf
import os
from astropy.io import fits


class LirisReduction():
    __lirisdr = "/Users/cjimenez/Documents/PHD/software/liris_iraf/lirisdr/"

    def __init__(self):
        # iraf.fitsio()
        pass

    def imageReduction(self, file_path):
        pass

    def spectroscopyReduction(self, file_path):
        pass

    def listFiles(self, path_folder):
        onlyfiles = [f for f in os.listdir(path_folder) if os.path.isfile(path_folder + "/" + f) and f[0:2] == "cr"]
        files = {"A": [], "B": []}
        for file in onlyfiles:
            if os.path.splitext(file)[1] == ".fit":
                f = fits.open(path_folder + file)
                h = f[0].header

                files[h["OBJECT"][0:1]].append(file)

        filea = open(path_folder + "fileson.lst", "w+")
        fileb = open(path_folder + "filesoff.lst", "w+")
        filea.write("\n".join(files["A"]))
        fileb.write("\n".join(files["B"]))
        filea.close()
        fileb.close()

    def skyFlatFields(self):
        pass

    def headerIspec(self, file_path):
        return iraf.hselect(file_path + "[1]", "$I,OBJECT,RA,DEC,LIRSLNAM,LIRGRNAM,\
          LIRLAMP1,LIRLAMP2,EXPTIME,AIRMASS,DATE-OBS,UTOBS", "yes", Stdout=1)[0].split("\t")

    def pixelScramblingCorrection(self, files_path):
        iraf.lcpixmap(files_path + "r*.fit", output="c", outlist="")

    def wavelengthCalibration(self, ref_a, ref_b):
        path = os.path.dirname(ref_a) + "/"
        iraf.imdelete(path + "oh_skylines.fits", verify="no")

        iraf.imcopy(ref_a + "[1]", path + "oh_skylines")
        iraf.imcopy(ref_b + "[1][*,600:643]", path + "oh_skylines[*,600:643]")
        iraf.display(path + "oh_skylines.fits[1]", 1)
        return path + "oh_skylines.fits"

    def reidentify(self, file_path):
        path = os.path.dirname(file_path)
        os.chdir(path)
        iraf.reidentify("oh_skylines", "oh_skylines", section="line 616", newaps="yes", overrid="yes", refit="yes",
                        trace="yes", step=20, nsum=20, nlost=3, cradius=5., match=-6.0, maxfeat=50, minsep=2.,
                        verbose="yes")

    def emisionLinesMap(self, file_path):

        path = os.path.dirname(file_path)
        os.chdir(path)

        iraf.identify("oh_skylines", section="line 616", coordlist=self.__lirisdr + "std/ohlines.dat",
                      match=-2., maxfeat=50, zwidth=100, fwidth=4.0, cradius="5.",
                      minsep=2.0, function="legendre", order=4, low_rej=3.0, high_rej=3.0)

    def checkCalibration(self, file_path):
        path = os.path.dirname(file_path)
        os.chdir(path)
        iraf.fitcoords("oh_skylines", fitname="oh_skylines", interactive="yes", combine="yes", function="legendre",
                       xorder=5, yorder=4)
        # iraf.imdelete("oh_skylines_cal", verify="no")
        # iraf.transform("oh_skylines", "oh_skylines_cal", interpt="spline3", xlog="no",
        #               ylog="no", flux="yes", database="database", fitname="oh_skylines")

    def starLiris(self):
        iraf.lirisdr()
        iraf.lspect()


if __name__ == "__main__":
    lr = LirisReduction()

    flat_path = "/Users/cjimenez/Documents/PHD/data/liris_june2019/20190619/FRED/FLAT/hkspec/lr_hk/Nflat.fits"
    science_path = "/Users/cjimenez/Documents/PHD/data/liris_june2019/20190619/L1608+4259/TARGET/hkspec/"
    a_ref = "cr2818596.fit"
    b_ref = "cr2818597.fit"
    # calibration_sky_lines = lr.wavelengthCalibration(science_path+a_ref,science_path+b_ref)

    # lr.emisionLinesMap(calibration_sky_lines)
    lr.reidentify(science_path+"oh_skylines.fits")
    # lr.checkCalibration(science_path+"oh_skylines.fits")
    # lr.starLiris()
    #lr.listFiles(science_path)
