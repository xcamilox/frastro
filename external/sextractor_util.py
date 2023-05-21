import subprocess

from core.fits.FITSFIle import FITSFile
from telescopes.telescopes import Telescopes
from core.utils.config import Config
from core.utils.log import LOGUtil
from core.utils.image_util import ImageUtils

import os

class DefaultParmasSextractor():
    basic_params=["NUMBER",
            "X_IMAGE",
            "Y_IMAGE",
            "ALPHA_J2000",
            "DELTA_J2000",
            "FLAGS",
            "ISOAREAF_IMAGE",
            "ISOAREA_IMAGE",
            "A_IMAGE",
            "B_IMAGE",
            "THETA_IMAGE",
            "ELLIPTICITY",
            "ELONGATION",
            "KRON_RADIUS",
            "FWHM_IMAGE",
            "CLASS_STAR",
            "FLUX_RADIUS",
            "BACKGROUND",
            "FLUX_MAX",
            "MAG_ISO",
            "MAGERR_ISO",
            "FLUX_ISO",
            "FLUXERR_ISO",
            "FLUX_RADIUS              #Fraction-of-light radii                                   [pixel]",
            "MAG_AUTO",
            "MAGERR_AUTO",
            "FLUX_AUTO",
            "FLUXERR_AUTO",
            "MAG_APER(5)",
            "MAGERR_APER(5)",
            "FLAGS                    #Extraction flags    ",
            "FLUX_APER(5)",
            "FLUXERR_APER(5)"]

    basic_config=["#-------------------------------- Catalog ------------------------------------",
            "CATALOG_NAME     catalog.cat       # name of the output catalog",
            "CATALOG_TYPE     FITS_1.0     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,",
            "                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC",
            "PARAMETERS_NAME  default.param  # name of the file containing catalog contents",
            "#------------------------------- Extraction ----------------------------------",
            "DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)",
            "DETECT_MINAREA   5.0              # min. # of pixels above threshold",
            "DETECT_MAXAREA   0              # max. # of pixels above threshold (0=unlimited)",
            "THRESH_TYPE      RELATIVE       # threshold type: RELATIVE (in sigmas) or ABSOLUTE (in ADUs)",
            "DETECT_THRESH    2.0           # <sigmas> or <threshold>,<ZP> in mag.arcsec-2",
            "ANALYSIS_THRESH  13.            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2",
            "FILTER           Y              # apply filter for detection (Y or N)?",
            "FILTER_NAME      /Applications/astro/sextractor/config/default.conv   # name of the file containing the filter",
            "FILTER_THRESH                   # Threshold[s] for retina filtering",
            "DEBLEND_NTHRESH  13          # Number of deblending sub-thresholds",
            "DEBLEND_MINCONT  0.001          # Minimum contrast parameter for deblending",
            "CLEAN            Y              # Clean spurious detections? (Y or N)?",
            "CLEAN_PARAM      1.0            # Cleaning efficiency",
            "MASK_TYPE        CORRECT        # type of detection MASKing: can be one of NONE, BLANK or CORRECT",
            "#-------------------------------- WEIGHTing ----------------------------------",
            "WEIGHT_TYPE      NONE           # type of WEIGHTing: NONE, BACKGROUND, MAP_RMS, MAP_VAR or MAP_WEIGHT",
            "RESCALE_WEIGHTS  Y              # Rescale input weights/variances (Y/N)?",
            "WEIGHT_IMAGE     coadd.weight.fits    # weight-map filename",
            "WEIGHT_GAIN      N,N              # modulate gain (E/ADU) with weights? (Y/N)",
            "WEIGHT_THRESH                   # weight threshold[s] for bad pixels",
            "#------------------------------ Photometry -----------------------------------",
            "PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>",
            "PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,",
            "PHOT_APERTURES  5.0,8.0,12.0,20.0      # MAG_APER aperture diameter(s) in pixels",
            "PHOT_AUTOPARAMS 2.5, 3.5        # MAG_AUTO parameters: <Kron_fact>,<min_radius>",
            "PHOT_FLUXFRAC   0.2, 0.5, 0.8   # Fraction of FLUXAUTO defining each element",
            "SATUR_LEVEL      50000.0        # level (in ADUs) at which arises saturation",
            "SATUR_KEY        SATURANT       # keyword for saturation level (in ADUs)",
            "MAG_ZEROPOINT    29          # magnitude zero-point AB=2.5*(23-math.log10(3631000000))-48.6 ",
            "MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)",
            "GAIN             0.95            # detector gain in e-/ADU ",
            "GAIN_KEY         GAIN           # keyword for detector gain in e-/ADUs",
            "PIXEL_SCALE      0.25            # size of pixel in arcsec (0=use FITS WCS info) spitzer irac band 2",
            "#------------------------- Star/Galaxy Separation ----------------------------",
            "SEEING_FWHM      1.0            # stellar FWHM in arcsec",
            "STARNNW_NAME     /Applications/astro/sextractor/config/default.nnw    # Neural-Network_Weight table filename",
            "#------------------------------ Background -----------------------------------",
            "BACK_TYPE        AUTO           # AUTO or MANUAL",
            "BACK_VALUE       0.0            # Default background value in MANUAL mode",
            "BACK_SIZE        256             # Background mesh: <size> or <width>,<height>",
            "BACK_FILTERSIZE  9              # Background filter: <size> or <width>,<height>",
            "BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL",
            "BACKPHOTO_THICK  100             # thickness of the background LOCAL annulus",
            "BACK_FILTTHRESH  0.0            # Threshold above which the background- map filter operates",
            "#------------------------------ Check Image ----------------------------------",
            "CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,",
            "                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,",
            "                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,",
            "                                # or APERTURES",
            "CHECKIMAGE_NAME  check.fits     # Filename for the check-image",
            "#--------------------- Memory (change with caution!) -------------------------",
            "MEMORY_OBJSTACK  3000           # number of objects in stack",
            "MEMORY_PIXSTACK  300000         # number of pixels in stack",
            "MEMORY_BUFSIZE   1024           # number of lines in buffer",
            "#------------------------------- ASSOCiation ---------------------------------",
            "ASSOC_NAME       sky.list       # name of the ASCII file to ASSOCiate",
            "ASSOC_DATA       2,3,4          # columns of the data to replicate (0=all)",
            "ASSOC_PARAMS     2,3,4          # columns of xpos,ypos[,mag]",
            "ASSOCCOORD_TYPE  PIXEL          # ASSOC coordinates: PIXEL or WORLD",
            "ASSOC_RADIUS     2.0            # cross-matching radius (pixels)",
            "ASSOC_TYPE       NEAREST        # ASSOCiation method: FIRST, NEAREST, MEAN,",
            "                                # MAG_MEAN, SUM, MAG_SUM, MIN or MAX",
            "ASSOCSELEC_TYPE  MATCHED        # ASSOC selection type: ALL, MATCHED or -MATCHED",
            "#----------------------------- Miscellaneous ---------------------------------",
            "VERBOSE_TYPE     NORMAL         # can be QUIET, NORMAL or FULL",
            "HEADER_SUFFIX    .head          # Filename extension for additional headers",
            "WRITE_XML        N              # Write XML file (Y/N)?",
            "XML_NAME         sex.xml        # Filename for XML output",
            "XSL_URL          file:///usr/local/share/sextractor/sextractor.xsl",
            "                                # Filename for XSL style-sheet",
            "NTHREADS         1              # 1 single thread",
            "FITS_UNSIGNED    N              # Treat FITS integer values as unsigned (Y/N)?",
            "INTERP_MAXXLAG   16             # Max. lag along X for 0-weight interpolation",
            "INTERP_MAXYLAG   16             # Max. lag along Y for 0-weight interpolation",
            "INTERP_TYPE      ALL            # Interpolation type: NONE, VAR_ONLY or ALL",
            "#--------------------------- Experimental Stuff -----------------------------",

            "PSF_NAME         default.psf    # File containing the PSF model",
            "PSF_NMAX         1              # Max.number of PSFs fitted simultaneously",
            "PATTERN_TYPE     RINGS-HARMONIC # can RINGS-QUADPOLE, RINGS-OCTOPOLE,",
            "                                # RINGS-HARMONICS or GAUSS-LAGUERRE",
            "SOM_NAME         default.som    # File containing Self-Organizing Map weights"]

    sextractor_psf_params = ["NUMBER                   Running object number                                    ",
                             "X_IMAGE                  #Object position along x                                   [pixel]",
                             "Y_IMAGE                  #Object position along y                                   [pixel]",
                             "VIGNET(45,45)            #Pixel data around detection                               [count]",
                             "FLUX_RADIUS              #Fraction-of-light radii                                   [pixel]",
                             "SNR_WIN                  #Gaussian-weighted SNR                                    ",
                             "FLUX_APER(1)             #Flux vector within fixed circular aperture(s)             [count]",
                             "FLUXERR_APER(1)          #RMS error vector for aperture flux(es)                    [count]",
                             "ELONGATION               #A_IMAGE/B_IMAGE                                          ",
                             "MAG_AUTO",
                             "MAGERR_AUTO",
                             "CLASS_STAR",
                             "FLAGS                    #Extraction flags    "]
    sextractor_psf_config = [
        "PHOT_APERTURES     #5arcsec"
    ]

    sextractor_scamp_params = [
        "NUMBER                   #Running object number                                    ",
        "X_IMAGE                  #Object position along x                                   [pixel]",
        "Y_IMAGE                  #Object position along y                                   [pixel]",
        "ALPHA_J2000",
        "DELTA_J2000",
        "MAG_AUTO",
        "MAGERR_AUTO",
        "CLASS_STAR",
        "XWIN_IMAGE",
        "YWIN_IMAGE",
        "ERRAWIN_IMAGE",
        "ERRBWIN_IMAGE",
        "ERRTHETAWIN_IMAGE",
        "FLUX_AUTO",
        "FLUXERR_AUTO",
        "FLAGS",
        "FLUX_RADIUS"


    ]
    dulamode= {}#{"-WEIGHT_IMAGE":"detect_rms.fits,measure_rms.fits",
            #"-WEIGHT_TYPE":"MAP_RMS,MAP_RMS",
            #"-WEIGHT_GAIN":"N,N"}

class SextractorUtil():

    sex_app=Config.getExtApp("sextractor")
    DEFAULT_CONFIG_CM = [sex_app, "-dd"]
    DEFAULT_PARAMS_CM = [sex_app, "-dp"]


    default_settings={
        "-FILTER_NAME":"/Applications/astro/sextractor/config/default.conv",
        "-STARNNW_NAME": "/Applications/astro/sextractor/config/default.nnw"
    }

    APP_PATH = sex_app

    def __init__(self):
        pass


    def createParams(self,params,path,name_file=""):

        if name_file=="":
            name_file = "default.param"

        f = open(path + name_file, "w+")
        f.write(" \n ".join(params))
        f.close()

        return path + name_file

    def createPSFexBasicParmas(self,path):
        path = path if path[-1] == "/" else path + "/"
        default = DefaultParmasSextractor().sextractor_psf_params

        return self.createParams(default, path,name_file="default_psf.param")



    def createScampBasicParmas(self,path):
        path = path if path[-1] == "/" else path + "/"
        default = DefaultParmasSextractor().sextractor_scamp_params

        return self.createParams(default, path,name_file="default_scamp.param")



    def createBasicParams(self,path):
        path = path if path[-1] == "/" else path + "/"
        default = DefaultParmasSextractor().basic_params
        return self.createParams(default,path)

    def createBasicConfigfile(self, path):
        path = path if path[-1] == "/" else path + "/"
        default = DefaultParmasSextractor().basic_config
        f = open(path + "default.sex", "w+")
        f.write(" \n ".join(default))
        f.close()
        return path + "default.sex"

    def createDefaultParamsfile(self, path):
        default = self.runCommand(self.DEFAULT_PARAMS_CM)
        path = path if path[-1] == "/" else path + "/"
        f = open(path + "default.param", "w+")

        f.write(default.decode("utf-8"))
        f.close()
        return path + "default.param"

    def createDefaultConfigfile(self, path):
        default = self.runCommand(self.DEFAULT_CONFIG_CM)
        path = path if path[-1] == "/" else path + "/"
        f = open(path + "default.sex", "w+")
        f.write(default.decode("utf-8"))
        f.close()
        return path + "default.sex"

    def getPixelScale(self,file):
        header=FITSFile.header(file)

        try:
            pixscale = header["PIXSCALE"]
        except Exception:
             tele= Telescopes.getTelescopeData(file)
             pixscale=tele["pixelscale"]
        return pixscale

    def runCommand(self, command):
        r = subprocess.run(command, stdout=subprocess.PIPE)
        return r.stdout

    def scampCatalog(self,file_path,output_path="",**kargs):

        if output_path == "":
            output_path = os.path.dirname(file_path)
            file_name=os.path.basename(file_path)
            file_name=os.path.splitext(file_name)[0]
            output_path = output_path+"/reduction/"+file_name+"/sextractor/"

        if not os.path.isdir(output_path):
            output_path = os.path.dirname(output_path)


        if not os.path.exists(output_path):
            os.makedirs(output_path)

        params_path = self.createScampBasicParmas(output_path)

        kargs["-CATALOG_TYPE"] = "FITS_LDAC"
        kargs["-PARAMETERS_NAME"] = "default_scamp.param"
        return self.runSextractor(file_path, output_path, dualmode=False, config_path=None, calc_psf=False,
                                  params_path=params_path, **kargs)

    def psfexCatalalog(self,file_path,output_path="",**kargs):


        if output_path == "":
            output_path = os.path.dirname(file_path)
            file_name=os.path.basename(file_path)
            file_name=os.path.splitext(file_name)[0]
            output_path = output_path+"/reduction/"+file_name+"/sextractor/"

        if not os.path.exists(output_path):
            os.makedirs(output_path)



        params_path = self.createPSFexBasicParmas(output_path)
        pixelScale = self.getPixelScale(file_path)

        kargs["-PHOT_APERTURES"] = str(5)
        kargs["-PHOT_FLUXFRAC"] = "0.5"
        kargs["-CATALOG_TYPE"] = "FITS_LDAC"
        kargs["-PARAMETERS_NAME"]="default_psf.param"
        kargs["-CATALOG_NAME"] = "cat_to_psf.fits"
        return self.runSextractor(file_path,output_path,dualmode=False,config_path=None,calc_psf=False,params_path=params_path,**kargs)


    def runSextractor(self,files_path,output_path="", dualmode=False,config_path=None,params_path=None,calc_psf=True, **kargs):


        commands = []




        #if ImageUtils.getFirstImageHDUIndex(files_path) > 0:
        #    image_tmp=ImageUtils.getFirstImageFromFitsFile(files_path)



        current_path = self.runCommand(["pwd"])

        if output_path == "":

            if dualmode:
                path = files_path[0]
            else:
                path = files_path
            output_path = os.path.dirname(path)
            file_name=os.path.basename(path)
            file_name=os.path.splitext(file_name)[0]
            output_path = output_path+"/reduction/"+file_name+"/sextractor/"

        if not os.path.exists(output_path):
            os.makedirs(output_path)

        os.chdir(output_path)

        LOGUtil.log(" ======= Run Sextractor ==========", output_path)
        LOGUtil.log("Config files in:", output_path)
        commands.append(self.APP_PATH)

        if dualmode:
            commands.append(",".join(files_path))
            kargs = {**kargs,**DefaultParmasSextractor().dulamode}
        else:
            commands.append(files_path)

        commands.append("-c")

        if params_path == None:
            #default_path = self.createDefaultParamsfile(output_path)
            self.createBasicParams(output_path)

        if config_path == None:
        #config_path = self.createDefaultConfigfile(output_path)
            config_path = self.createBasicConfigfile(output_path)

        commands.append(config_path)

        kargs.update(self.default_settings)


        for arg in kargs:
            commands.append(arg)
            commands.append(kargs[arg])

        LOGUtil.log(" ".join(commands), output_path)


        r = self.runCommand(commands)
        LOGUtil.log("Process: \n" + str(r), output_path)

        if "-CATALOG_NAME" in kargs:
            file_name = kargs["-CATALOG_NAME"]
        else:
            file_name = "/catalog.cat"

        LOGUtil.log("Output: " + output_path +file_name, output_path)
        return output_path+file_name
        # os.chdir(current_path.decode("utf-8").replace("\n",""))


if __name__ == "__main__":
    path_files = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/liverpoolb/reduction/swarp/swarp_S4TMJ0143-1006.fits"
    output_path= "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/liverpoolb/reduction/sextractor"
    sex = SextractorUtil()
    kargs={}
    sex.runSextractor(path_files,output_path, **kargs)