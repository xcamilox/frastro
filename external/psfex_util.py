import subprocess


from core.utils.config import Config
from core.utils.log import LOGUtil
from telescopes.telescopes import Telescopes
from core.fits.FITSFIle import FITSFile


from external.sextractor_util import SextractorUtil
import os
from astropy.io.votable import parse
from datetime import datetime
import os.path as path
from astropy.io import fits
import numpy as np

class PSFexUtil():
    app_path = Config.getExtApp("psfex")
    default_parmas = [app_path,"-dd"]


    def createDefaultParms(self,path):
        params=self.runCommands(self.default_parmas)
        f = open(path + "/default.psfex", "w+")
        f.write(params.decode("utf-8"))
        f.close()
        return path + "/default.psfex"

    def runCommands(self,command):
        r = subprocess.run(command, stdout=subprocess.PIPE)
        return r.stdout


    def runPSFEx(self,fits_image_path,output_path, **kargs):

        commands=[self.app_path]
        main_path=output_path

        if output_path == "":
            output_path = os.path.dirname(fits_image_path)
            file_name=os.path.basename(fits_image_path)
            file_name=os.path.splitext(file_name)[0]
            output_path = output_path+"/reduction/"+file_name+"/psf/"

        if not os.path.exists(output_path):
            os.makedirs(output_path)



        sex=SextractorUtil()
        kargs={}
        kargs["-MAG_ZEROPOINT"]="29.03"
        sex_cat_path = sex.psfexCatalalog(fits_image_path,main_path,**kargs)

        # select only stars

        cat=fits.open(sex_cat_path,mode="update")
        filter=np.where(cat[2].data["CLASS_STAR"] >= 0.7,True,False)

        cat[2].data=cat[2].data[filter]
        cat.flush()
        cat.close()


        file_name=os.path.basename(sex_cat_path)
        file_name = os.path.splitext(file_name)[0]
        commands.append(sex_cat_path)
        os.chdir(output_path)

        commands.append("-c")

        commands.append(self.createDefaultParms(output_path))
        pixelscale=Telescopes.getPixelScale(fits_image_path)

        kargs["-PSF_DIR"]=output_path
        kargs["-PSF_SAMPLING"] = "1.0"
        kargs["-PSF_PIXELSIZE"] = str(pixelscale)
        kargs["-OUTCAT_TYPE"] = "FITS_LDAC"
        kargs["-SAMPLE_MAXELLIP"] = "0.5"
        kargs["-PSF_SIZE"] = "15"

        for arg in kargs:
            commands.append(arg)
            commands.append(kargs[arg])

        LOGUtil.log(" ======= Run PSFEX ==========", output_path)
        LOGUtil.log(" ".join(commands), output_path)
        print(" ".join(commands))
        LOGUtil.log(" Output: "+output_path+"/"+file_name+".psf", output_path)

        r=self.runCommands(commands)

        return output_path,file_name+".psf"


    def getMeanFWHM(self,path):
        return self.readXMLOutput(path,field="FWHM_Mean")

    def readXMLOutput(self,path,field="FWHM_Mean"):

        votable = parse(path)
        table=votable.get_first_table()
        return table.array[field][0]


    def getPSF(self,file):
        hdu=fits.open(file)
        hdu[1].data
        data=np.array(hdu[1].data[0].array)
        return data[0][0][0],hdu[1].header


    def saveSinglePSFFile(self,file,output):
        data,header=self.getPSF(file)
        print(file)
        print("psf:", data.min(),data.max(),data.mean())
        FITSFile.saveFile(output,data,header)


    #return seeing from a fits file, is necesary pixel scale in header file
    #this method will be read the fwhm from psfex.xml file
    def calcSeeing(self, fits_file_image,output_path=""):

        if output_path == "":
            output_path = os.path.dirname(fits_file_image)
            file_name=os.path.basename(fits_file_image)
            file_name=os.path.splitext(file_name)[0]
            output_path = output_path+"/reduction/"+file_name+"/psf/"

        if not os.path.exists(output_path):
            os.makedirs(output_path)


        outpath_psf,file_psf=self.runPSFEx(fits_file_image, output_path)
        sex = SextractorUtil()
        outcat=sex.runSextractor(fits_file_image)



        cat=fits.open(outcat)
        data=cat[1].data
        sources=len(data)
        flux_radius=data["FLUX_RADIUS"].mean()


        fwhm_mean = self.getMeanFWHM(outpath_psf + "psfex.xml")

        pixel = Telescopes.getPixelScale(fits_file_image)

        if fwhm_mean <=0.0 or flux_radius<0.7:
            seeing=9999
        else:
            seeing =fwhm_mean * pixel

        print(outpath_psf+"psfex_out.cat",path.exists(outpath_psf+"psfex_out.cat"))

        return seeing,outpath_psf+"/"+file_psf




if __name__=="__main__":
    psf=PSFexUtil()
    fits_file="/Users/cjimenez/Documents/PHD/data/tmp/349.4015416666667_12.582483333333334/panstarss/panstars_band_g.fits"
    fits_file = "/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/liris/SDSSJ122040.72+084238.1/SDSSJ122040.72+084238.1_bandH.fits"
    fits_file = "/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/liris/SDSSJ122040.72+084238.1/SDSSJ122040.72+084238.1_bandH_scale.fits"


    file="/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/reduction/galfit/catalog.psf"
    output="/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/reduction/galfit/psf.fits"

    output_scale = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/reduction/galfit/psf_scale.fits"

    file = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/reduction/sextractor/tmp/12112018/reduction/psfex/psf.psf"
    output = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/reduction/sextractor/tmp/12112018/reduction/psfex/psf_single.fits"

    file="/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ0143-1006/observations/2018-08-10/reduction/stacking/swarp/coadd_nearest_sapling.fits"
    file="/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ0143-1006/observations/2018-10-04/reduction/stacking/swarp/coadd_nearest_sapling.fits"
    #file="/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ0143-1006/observations/2018-10-04/science/stack_best.fits"

    file="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/gtc/model/galfit/S4TMJ1051+4439_gtc_r_29apr2019_crop.fits"
    file="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/gtc/model/galfit/reduction/S4TMJ1051+4439_gtc_r_29apr2019_crop/psf/catalog.psf"
    output = "/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/gtc/model/galfit/reduction/S4TMJ1051+4439_gtc_r_29apr2019_crop/psf/psf.psf"
    psf.saveSinglePSFFile(file, output)
    #seeing=psf.calcSeeing(file)
    #print(seeing)
    file = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/science/reduction/psfex/stacking_12_11_2018/2018-11-26/catalog.psf"
    output = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/science/reduction/psfex/stacking_12_11_2018/2018-11-26/psf_single.fits"
    file ="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/2019-04-04/science/psf_best.fits"
    output = "/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/2019-04-04/science/psf_best_img.fits"
    output = "/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/2019-04-04/science/psf_best_img.fits"

    file = "/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/2019-04-25/science/psf_best.fits"
    output = "/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/2019-04-25/science/psf_best_img.fits"


    path="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/2019-04-28/"
    path="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/2019-04-26/"
    epochs = ["2018-11-16",
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
    for epoch in epochs:
        path=f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/'
        file=path+"science/psf_best.fits"
        output=path+"science/psf.fits"

        psf.saveSinglePSFFile(file,output)


    #FITSFile.scaleImage(output,output_scale,0,30000,deg=2)

    """

    path_folder = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/"
    onlyfiles = [f for f in os.listdir(path_folder) if path.isfile(path_folder + "/" + f)]
    result=[]
    f = open(path_folder + "seeing.log", "w+")
    for file in onlyfiles:
        if path.splitext(file)[1] == ".fits" or path.splitext(file)[1] == ".fit":
            fwhm,psf_file=psf.calcSeeing(path_folder+file)
            result.append([file,fwhm])
            f.write(file+" , "+str(fwhm)+"\n")

    print(result)
    f.close()

    """
    """
    output="/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181005_1"
    psf.runPSFEx(fits_file,output)
    fwhm_mean = psf.getMeanFWHM(output + "/reduction/psfex/psfex.xml")
    pixel = Telescopes.getPixelScale(fits_file)
    print(fwhm_mean * pixel)
    """