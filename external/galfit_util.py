import os
from os import path
import subprocess
import shutil
from astropy.io import fits
import json
#from frastro import ImageUtils, Config, MongodbManager,PSFexUtil

from external.psfex_util import PSFexUtil
from core.database.mongodb.mongodb_manager import MongodbManager
from core.utils.image_util import ImageUtils
from core.utils.config import Config

class GalfitUtils():

    default=["===============================================================================",
            "# IMAGE and GALFIT CONTROL PARAMETERS",
            "A) h_e_20181105_27_4_1_1_crop.fits            # Input data image (FITS file)",
            "B) model.fits       # Output data image block",
            "C) none                # Sigma image name (made from data if blank or none) ",
            "D) none   #        # Input PSF image and (optional) diffusion kernel",
            "E) 1                   # PSF fine sampling factor relative to data ",
            "F) none                # Bad pixel mask (FITS image or ASCII coord list)",
            "G) none                # File with parameter constraints (ASCII file) ",
            "H) 800    1100   800    1100   # Image region to fit (xmin xmax ymin ymax)",
            "I) 50    50          # Size of the convolution box (x y)",
            "J) 0              # Magnitude photometric zeropoint ",
            "K) 0.3054  0.3054        # Plate scale (dx dy)    [arcsec per pixel]",
            "O) regular             # Display type (regular, curses, both)",
            "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps",
            "",
            "# INITIAL FITTING PARAMETERS",
            "#",
            "#   For object type, the allowed functions are: ",
            "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, ",
            "#       ferrer, powsersic, sky, and isophote. ",
            "#  ",
            "#   Hidden parameters will only appear when they're specified:",
            "#       C0 (diskyness/boxyness), ",
            "#       Fn (n=integer, Azimuthal Fourier Modes),",
            "#       R0-R10 (PA rotation, for creating spiral structures).",
            "# ",
            "# -----------------------------------------------------------------------------",
            "#   par)    par value(s)    fit toggle(s)    # parameter description ",
            "# -----------------------------------------------------------------------------",
            "",
            "# Object number: 1",
            " 0) sersic                 #  object type",
            " 1) 51.55  51.55  1 1  #  position x, y",
            " 3) 20.0     1          #  Integrated magnitude	",
            " 4) 4.      1          #  R_e (half-light radius)   [pix]",
            " 5) 4.      1          #  Sersic index n (de Vaucouleurs n=4) ",
            " 9) 1.0      1          #  axis ratio (b/a)  ",
            "10) 0.0    1          #  position angle (PA) [deg: Up=0, Left=90]",
            " Z) 0                      #  output option (0 = resid., 1 = Don't subtract) ",
            "================================================================================"]

    """
    # Component number: 2",
    "0) sky                    #  Component type",
    "1) 0      1       #  Sky background at center of fitting region [ADUs]",
    "2) 0.0      0       #  dsky/dx (sky gradient in x)     [ADUs/pix]",
    "3) 0.0      0       #  dsky/dy (sky gradient in y)     [ADUs/pix]",
    "Z) 0                      #  Skip this model in output image?  (yes=1, no=0)",
    """



    def __init__(self):
        self.galfit_app = Config.getExtApp("galfit")
        self.psfex = PSFexUtil()


    def createModel(self,file_path_image, ra,dec,psfex_path="",magnitud=20,axis_ratio=1,angle=0,zp=0,crop=None):
        output_path = os.path.dirname(file_path_image)
        file = os.path.basename(file_path_image)
        reduction_path=output_path +"/model/galfit/"
        output_crop_file=reduction_path+ path.splitext(file)[0] + "_crop.fits"

        print("reductionpath",reduction_path)

        if not os.path.exists(reduction_path):
            os.makedirs(reduction_path)


        if psfex_path!="":
            self.psfex.saveSinglePSFFile(psfex_path,reduction_path+"psf.fits")

        shutil.copy2(file_path_image, reduction_path+os.path.basename(file_path_image))


        pixel = ImageUtils.getPixelFromCoordinate(file_path_image, ra, dec)


        #crop image
        if crop!=None:
            ImageUtils.imcopy(file_path_image, output_crop_file,(pixel[0][0], pixel[0][1]), crop)
            pixel = ImageUtils.getPixelFromCoordinate(output_crop_file, ra, dec)
            self.default[9] = "H) 1    "+str(crop[0])+"   1    "+str(crop[1])+"   # Image region to fit (xmin xmax ymin ymax)"
            file_path_image=output_crop_file


        self.default[2] = "A) "+os.path.basename(file_path_image)
        self.default[3] = "B) "+path.splitext(os.path.basename(file_path_image))[0]+"_model.fits"

        self.default[5] = "D) psf.fits"
        self.default[11] = "J) "+str(zp)
        self.default[33] = "1) "+str(pixel[0][0])+"  "+str(pixel[0][1])+"  1 1"
        self.default[34] = "3) " + str(magnitud) + "  1"
        self.default[37] = "9) " + str(axis_ratio) + " 1"
        self.default[38] = "10) " + str(angle) + " 1"
        feedme_path=reduction_path+path.splitext(file)[0]+".feedme"
        os.chdir(reduction_path)

        f = open(feedme_path, "w")



        f.write(" \n ".join(self.default))
        f.close()
        print("galfit ",os.path.basename(feedme_path))
        feedme_file=os.path.basename(feedme_path)
        commands=[self.galfit_app,feedme_file]

        model_path=reduction_path+path.splitext(os.path.basename(file_path_image))[0] + "_model.fits"

        if os.path.exists(model_path):
            os.remove(model_path)

        r = subprocess.run(commands, stdout=subprocess.PIPE)
        print(r)


        try:
            shutil.copy2(model_path, output_path+"/model_lens_2.fits")
            file_fits = fits.open(output_path+"/model_lens_2.fits", mode="update")
            for item in file_fits[1].header:
                print(file_fits[1].header[item])
                try:
                    file_fits[3].header[item] = file_fits[1].header[item]
                except:
                    print("Error key ", item)
            file_fits.flush()
            return output_path+"/model_lens_2.fits"
        except FileNotFoundError:
            return ""


if __name__ == "__main__":
    epoachs=["2018-11-16",
     "2019-01-10",
     "2019-04-04",
     "2019-04-25",
     "2019-04-28",
     "2019-05-01",
     "2018-11-20",
     "2019-02-03",
     "2019-04-23",
     "2019-04-26",
     "2019-04-29",
     "2019-05-02",
     "2018-12-17",
     "2019-02-26",
     "2019-04-24",
     "2019-04-27",
     "2019-04-30"]

    galfit=GalfitUtils()
    ra= "162.78920833333333"
    dec="44.652363888888888"
    magnitud = "17.52"
    zp = 30.25
    angle = -67.63
    axis_ratio = 1.58

    #GTC
    ra="162.7795830833333"
    dec="44.64684252777777"
    #psf_best = f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/science/psf_best.fits'
    image = f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/gtc/S4TMJ1051+4439_gtc_r_29apr2019.fits'
    zp=29.033
    model_path = galfit.createModel(image, ra, dec, magnitud=magnitud, axis_ratio=axis_ratio,angle=angle, zp=zp,crop=(300,300))
    print("GTC",model_path)

    jsonfile="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/liverpool_sn_v3.json"
    """"
    f = open(jsonfile, "r")

    # Reading from file
    data = json.loads(f.read())
    epocas=data[0]["epoach"]
    model_residuals="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/model_residuals/"
    for epoch in epoachs:
        psf_best=f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/science/psf_best.fits'
        image=f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/science/stack_best.fits'
        zp=epocas[epoch]["stacking"]["best"]["zp"]
        seeing = epocas[epoch]["stacking"]["best"]["seeing_cal"]
        model_path = galfit.createModel(image, ra, dec, psfex_path=psf_best, magnitud=magnitud, axis_ratio=axis_ratio,
                                    angle=angle, zp=zp)
#        print(model_path)


#Generate sex tractro catalog of the model residueal

        path=f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/science/'
        commands=["cd",path]
        r = subprocess.call(commands, shell=True)
        os.chdir(path)
        file_path_image=path+"model_lens_2.fits"
        output_crop_file=path+f'img_model_residual_{epoch}.fits'
        psffile=path+"psf.fits"
        #pixel = ImageUtils.getPixelFromCoordinate(file_path_image, ra, dec)


        shutil.copy2(file_path_image, output_crop_file)
        file_fits = fits.open(output_crop_file, mode="update")
        del file_fits[0]
        del file_fits[0]
        del file_fits[0]
        file_fits.flush()
        file_fits.close()



        #ImageUtils.imcopy(file_path_image, output_crop_file, (pixel[0][0], pixel[0][1]),n_ext=3)
        commands=["sex", "model_lens_residual.fits", "-c", "default.sex","-CATALOG_NAME",f'cat_residual_model_{epoch}.fits',"-PSF_NAME",psffile,"-CHECKIMAGE_TYPE","BACKGROUND","-PHOT_APERTURES","5.0,8.0,12.0,20.0","-DETECT_MINAREA","5","-DETECT_THRESH","5.0","-ANALYSIS_THRESH","7","-MAG_ZEROPOINT",str(zp),"-SEEING_FWHM",str(seeing)]
        r = subprocess.run(commands, stdout=subprocess.PIPE)

        model_img = model_residuals + f'img_model_residual_{str(seeing)[0:5]}_{epoch}.fits'
        shutil.copy2(output_crop_file, model_img)
        output_cat=path+f'cat_residual_model_{epoch}.fits'
        model_cat=model_residuals+f'cat_model_residual_{str(seeing)[0:5]}_{epoch}.fits'
        shutil.copy2(output_cat, model_cat)


"""
"""
    mgdb = MongodbManager()
    mgdb.setCollection("liverpool_sn_v3")
    all_data = mgdb.getData({"id":"BELLSJ1337+3620",'epoach':{'$exists': 'true'}})
    for target in all_data:
        #try:
        ra=target["ra"]
        dec = target["dec"]
        for epoach in target["epoach"]:
            stack= target["epoach"][epoach]["stacking"]
            if "error" not in stack.keys():

                #psf and image path from database for best stacking
                psf_best=stack["best"]["psf_path"]
                science_path = os.path.dirname(psf_best)
                image=science_path+"/stack_best.fits"

                #detection data
                detections_best = stack["best"]["detection"]["best"]
                magnitud=detections_best["mag"]
                axis_ratio = detections_best["ELONGATION"]
                angle = detections_best["THETA_IMAGE"]
                zp =stack["best"]["zp"]

                #create model from sextrctor detection data
                model_path = galfit.createModel(image,ra,dec,psfex_path=psf_best,magnitud=magnitud,axis_ratio=axis_ratio,angle=angle,zp=zp)

                #update database with the new model path
                if model_path != "":
                    mgdb.update({'id':target["id"]},{"$set":{"epoach."+epoach+".stacking.best.model":model_path}})
                else:
                    print("error,",target["id"],epoach)

        #except KeyError:
        #    print("error with ",target["id"])
"""


