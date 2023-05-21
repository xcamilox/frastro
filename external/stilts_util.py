import subprocess
from core.utils.utils import Utils
from core.utils.config import Config
import os
import os.path as path
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
import pandas as pd

class StiltsUtil():
    RUN_SWARP_CM = Config.getExtApp("stilts")
    MATCH_TWO_CM = RUN_SWARP_CM+' tmatch2 in1={0} in2={1} out={2} matcher=sky params="{3}"'
    MATCH_TWO_XY = RUN_SWARP_CM + ' tmatch2 in1={0} in2={1} out={2} matcher=2d params="{3}" suffix1=_a suffix2=_b'


    def __init__(self):
        pass

    def crossMatchTwoCatalogs(self,cat1_table,cat2_table,output,separation,list_str_params=[]):

        #save astropy tables in VOTables for stilts math
        cmd=str(self.MATCH_TWO_CM).format(cat1_table,cat2_table,output,separation)

        cmd=cmd.split(" ")
        cmd=cmd+list_str_params

        return self.runCommand(cmd)

    def crossMatchTwoCatalogsXY(self, cat1_table, cat2_table, output,separation, list_str_params=[]):
        # save astropy tables in VOTables for stilts math
        cmd = str(self.MATCH_TWO_XY).format(cat1_table, cat2_table, output, separation)

        cmd = cmd.split(" ")
        cmd = cmd + list_str_params
        print(" ".join(cmd))
        return self.runCommand(cmd)

    def runCommand(self,command):
       print("run command ",command)
       r = subprocess.run(command,stdout=subprocess.PIPE)

       return r.stdout



if __name__== "__main__":
    epoachs = ["2018-11-16",
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
    df=None
    for epoch in epoachs:
        refcat="in1=/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/coordiantes.csv"
        pathfile = f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/science/catalog_deep_model.fits'
        output=f'out=/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/science/match_model_catalog.fits'
        """
        Find all catalogs from the residual
        catalog=fits.open(pathfile,mode="update")

        if len(catalog)>2:
            del catalog[1]
            del catalog[1]
            del catalog[1]
            catalog.flush()
        catalog.close()
        """

        params=['java', '-jar', '/Applications/astro/stilts.jar', 'tmatch2',
         refcat,
         "in2="+pathfile,
         'ifmt1=csv','ifmt2=fits',
         output, 'matcher=sky', 'params="3"',"find=all",
         'values1=RA DEC', 'values2=ALPHA_J2000 DELTA_J2000', 'ofmt=fits-basic','progress=log']
        #stilts=StiltsUtil()
        #stilts.runCommand(params)
        t = Time(epoch, format='isot')
        table_file = f'/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/observations/{epoch}/science/match_model_catalog.fits'
        dat = Table.read(table_file, format='fits')
        del dat['MAG_APER']
        del dat['MAGERR_APER']
        del dat['FLUX_APER']
        del dat['FLUXERR_APER']
        del dat["RA"]
        del dat["DEC"]
        dat["id"] = dat["id"].astype(str)

        if df is None:
            df = dat.to_pandas()
            df["mjd"]=t.mjd
        else:
            df_temp = dat.to_pandas()
            df_temp["mjd"] = t.mjd
            df=pd.concat([df, df_temp], axis=0)
    full_sample="/Users/camilojimenez/PHD/thesis/data/chaper_2/S4TMJ1051+4439/full_sample.csv"
    df.to_csv(full_sample,index=False)

        #catalog = fits.open(pathfile, mode="update")
