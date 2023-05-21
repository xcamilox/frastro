import shutil

from astropy.io import fits

from core.data.archive.gaia_archive_cp import GAIAAchiveCP
from core.utils.config import Config
from core.utils.utils import Utils
from astropy.table import Table
import os

class GaiaAstrometry():

    build_index="build-astrometry-index -i {0} -e1 -o {1} -P {2} -S MAG -E"
    solved_file="/usr/local/astrometry/bin/solve-field --scale-units arcsecperpix --scale-low {0} --scale-high {1} {2} --overwrite --backend-config {3}"


    def __init__(self):
        pass
    #@return astropy Table
    #ra dec in degrees, radius in arcsec
    def getGaiaCatalog(self,ra,dec,radius=60,save=False,output=""):
        print("Get GAIA Catalog: ",ra,dec)
        gaia_archive=GAIAAchiveCP()
        table=gaia_archive.getCatalog(ra, dec, radius)
        if save:
            if output == "":
                raise ValueError("output file path is required")
            else:
                self.saveCatalog(table,output)
        return table




    def generateIndex(self,catalog_path):

        print("Generate index files")
        file_name=os.path.basename(catalog_path)
        dir_path=os.path.dirname(catalog_path)+"/"
        astrometry_data = dir_path + "index"
        os.chdir(os.path.dirname(dir_path))
        Utils.createPath(astrometry_data)

        #astrometry_data=Config.getPath("astrometry-data")


        for i in range(-5,4,1):
            sign="P"
            if i < 0:
                sign = "M"

            index_output=file_name[:-5]+"_"+sign+str(abs(i))+".index"
            command = self.build_index.format(file_name,index_output,str(i))
            print(os.getcwd())
            runcommand=command.split(" ")
            r=Utils.runCommand(runcommand)
            if os.path.exists(dir_path+index_output):
                print("copy file")
                shutil.move(dir_path+index_output, astrometry_data+"/"+index_output)
                print("copy index to",dir_path+index_output, astrometry_data+"/"+index_output)
            print(r)

        config_cfg = "# In which directories should we search for indices?\nadd_path "+astrometry_data+"\n# Load any indices found in the directories listed above.\nautoindex"

        config_file = open(dir_path+"backend.cfg","w+")
        config_file.write(config_cfg)
        config_file.close()

        # set astromery
        # generate index
        """
                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_M5.index -P -5 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_M4.index -P -4 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_M3.index -P -3 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_M2.index -P -2 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_M1.index -P -1 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_0.index -P 0 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_p1.index -P 1 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_p2.index -P 2 -S MAG -E

                        build-astrometry-index -i SDSSJ133515.79+4330_GAIA2deg.fits -e1 -o SDSSJ133515.79+4330_GAIA2deg_p3.index -P 3 -S MAG -E
                        
                        solve-field --scale-units arcsecperpix --scale-low 0.20 --scale-high 0.30 SDSSJ133515.79+4330_bandH.fits
        """

    def fixHeaderError(self,file):
        with fits.open(file, mode='update') as hdul:
            hdr = hdul[0].header
            if "LIRLAMP1" in hdr:
                del hdr["LIRLAMP1"]
                hdul.flush()
            if "LIRLAMP2" in hdr:
                del hdr["LIRLAMP2"]
                hdul.flush()
            hdul.close()


    def doAstrometry(self,file_image,low_scale=0.20,high_scale=0.30):
        print("do astrometry")
        file_name = os.path.basename(file_image)
        dir_path = os.path.dirname(file_image) + "/"
        os.chdir(dir_path)

        conf_backend=dir_path+"gaia/backend.cfg"

        command=self.solved_file.format(low_scale,high_scale,file_name,conf_backend)
        r=Utils.runCommand(command.split(" "))
        print(r)

    def saveCatalog(self,astropy_table,output_path):
        dir_path=os.path.dirname(output_path)
        if not os.path.exists(dir_path):
            Utils.createPath(dir_path)
        format = "fits"

        t = Table()
        t["RA"] = list(astropy_table["ra"])
        t["DEC"] = list(astropy_table["dec"])
        t["MAG"] = list(astropy_table["phot_g_mean_mag"])
        t.write(output_path, format=format, overwrite=True)



if __name__ == "__main__":
    gc=GaiaAstrometry()
    ra = 214.87670833333337
    dec = 43.6916111111111


    file_path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/S1419+4341_J.fits"
    file_path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/S1419+4341_h.fits"
    file_path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/S1419+4341_ks.fits"

    file_path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/L2317+1234/L2317+1234_h.fits"
    file_path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/L2317+1234/L2317+1234_j.fits"
    file_path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/L2317+1234/L2317+1234_ks.fits"

    path = "/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/gaia/"
    gc.fixHeaderError(file_path)
    Utils.createPath(path)
    file = "S1419+4341_GAIA2deg.fits"
    output_path = path + file
    if not os.path.exists(output_path):
        print("Catalog already")
        catalog = gc.getGaiaCatalog(ra, dec,save=True,output=output_path)
        # format="fits"
        # t = Table()
        # t["RA"] = list(catalog["ra"])
        # t["DEC"] = list(catalog["dec"])
        # t["MAG"] = list(catalog["phot_g_mean_mag"])
        # t.write(output_path,format=format,overwrite=True)

    #gc.generateIndex(output_path)

    #gc.doAstrometry(file_path,low_scale=0.25,high_scale=0.31)

