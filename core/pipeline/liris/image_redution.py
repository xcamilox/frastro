from astropy.io import fits
import numpy as np
import shutil
import os
from frastro.core.utils.utils import Utils
from astropy import coordinates as coords
from astropy import units as u
import matplotlib.pyplot as plt

def setupImage(folder_path,filter="c"):

    row_data=folder_path
    images=Utils.getFilesList(row_data,".fit")

    for img in images:

        try:
            header = fits.getheader(row_data+img,0)
            print(img,img[0:1])
            if img[0:1] == filter:
                obtype = header["OBSTYPE"]
                imagetype = header["IMAGETYP"] #
                objName = header["CAT-NAME"] #object cat
                band = header["LIRF2NAM"] #filter name
                ra = header["RA"]
                dec = header["DEC"]
                date=header["DATE-OBS"]

                path=folder_path+objName+"/"+obtype+"/"+band+"/"
                Utils.createPath(path)
                print(path+img)
                if not os.path.exists(path+img):
                    shutil.copy2(row_data+img, path+img)
                    print(obtype,imagetype,objName,band,ra,dec,date,img)
                else:
                    print("already",path+img)
        except:
            print("error",row_data+img)
            log=open(folder_path+"error.log","w+")
            log.write(row_data+img + "\n")
            log.close()



def fixHeaderError(file):
    with fits.open(file, mode='update') as hdul:
        hdr = hdul[0].header
        if "LIRLAMP1" in hdr:
            del hdr["LIRLAMP1"]
            hdul.flush()
        if "LIRLAMP2" in hdr:
            del hdr["LIRLAMP2"]
            hdul.flush()
        hdul.close()


def checkCoordinates(folder_path):
    row_data = folder_path
    images = Utils.getFilesList(row_data, ".fits")
    list=[]
    steps="R1"
    toplot={"R1":[]}
    all=[]
    for img in images:

        header = fits.getheader(row_data + img, 0)

        if img[0:1] == "g":
            obtype = header["OBSTYPE"]
            imagetype = header["IMAGETYP"]  #
            objName = header["CAT-NAME"]  # object cat
            band = header["LIRF2NAM"]  # filter name
            ra = header["RA"]
            dec = header["DEC"]
            mjd= header["MJD-OBS"]
            object = header["OBJECT"]


            coordinate=coords.SkyCoord(ra=ra, dec=dec,unit=(u.hourangle, u.deg))
            list.append([coordinate.ra.deg,coordinate.dec.deg,mjd])

            print(coordinate)
            date = header["DATE-OBS"]
            #print(obtype,imagetype,objName,band,date,object)
            end = str(object).find("/")
            print(ra,dec)
            all.append([coordinate.ra.deg, coordinate.dec.deg, mjd,object,img])
            if object[0:end] != steps:
                steps = object[0:end]
                toplot[object[0:end]]=[[coordinate.ra.deg, coordinate.dec.deg, mjd,object,img]]

            else:
                toplot[object[0:end]].append([coordinate.ra.deg,coordinate.dec.deg,mjd,object,img])

    path = "/Users/cjimenez/Documents/PHD/data/liris_enero_2019/SDSSJ1419+4341/ks/files/"
    data = np.array(all)
    data = data[np.argsort(data[:, 2])]

    plt.plot(data[:, 0], data[:, 1], "o-")
    plt.show()

    """
    for key in toplot.keys():

        data=np.array(toplot[key])
        data = data[np.argsort(data[:, 2])]
        print(data[:,4])
        file=open(path+key+".txt","w+")
        file.write("\n".join(data[:,4]))
        file.close()
        plt.plot(data[:,0],data[:,1],"o-")
        plt.title(key+" "+str(len(toplot[key])))
        plt.show()
    """
def runswarp(file):
    pass

def fixObjectName(file):
    with fits.open(file, mode='update') as hdul:
        hdr = hdul[0].header
        if "OBJECT" in hdr:
            if ":" in hdr["OBJECT"]:
                hdr["OBJECT"]=str(hdr["OBJECT"])[str(hdr["OBJECT"]).find(":")+1:]
            hdul.flush()
        hdul.close()

if __name__ == "__main__":

    path = "/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/liris_2018_06_23/list_image.lst"
    pathcecilia ="/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/iris/cecilia/files.lst"
    pathraine = "/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/iris/raine/files.lst"
    path= "/Users/cjimenez/Documents/PHD/data/liris_enero_2019/"



    path="/Users/cjimenez/Documents/PHD/data/unlens_liris/SDSSJ1013+4650/SDSSJ1013+4650_Ks.new"
    path="/Users/cjimenez/Documents/PHD/data/unlens_liris/SDSSJ1013+4650/coadd.fits"
    path="/Users/cjimenez/Documents/PHD/data/wht_SDSSj085038.86+154917.8/"
    path="/Users/cjimenez/Documents/PHD/data/liris_enero_2019/SDSSJ1013+4650/stacking_Ks.new"

    path="/Users/cjimenez/Documents/PHD/data/liris_enero_2019/SDSSJ1157+0113/h/"
    path="/Users/cjimenez/Documents/PHD/data/liris_enero_2019/J003011/ks/"
    #checkCoordinates(path)
    #setupImage(path)
    #fixObjectName(path)
    path="/Users/cjimenez/Documents/PHD/data/liris_enero_2019/J003011/ks/stack/J003011_ks3.new"
    path="/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/liris/observations/"
    path="/Users/cjimenez/Documents/PHD/data/liris_june2019/data/"
    path="/Users/cjimenez/Documents/PHD/data/liris_june2019/20190614/"
    path="/Users/cjimenez/Documents/PHD/data/liris_june2019/20190619/"
    path="/Users/cjimenez/Documents/PHD/data/liris_june2019/20190615/"
    path = "/Users/cjimenez/Documents/PHD/data/liris_june2019/20190616/"
    path="/Users/cjimenez/Documents/PHD/data/liris_june2019/20190617/"
    path = "/Users/cjimenez/Documents/PHD/data/liris_june2019/20190618/"
    setupImage(path,"c")


    #astrometry find catalogs with 1deg radio
    # set astromery
        #generate index
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
"""
        #fixHeaderError() remove LIRLAMP1,LIRLAMP2 in the header
        #run astromery
            #solve-field --scale-units arcsecperpix --scale-low 0.20 --scale-high 0.30 SDSSJ133515.79+4330_bandH.fits









        #fix header wcs
            #cd /Applications/astro/sip_tpv-master/sip_tpv
            #python sip_to_pv.py /Users/cjimenez/Documents/PHD/data/unlens_liris/SDSSJ1013+4650/SDSSJ1013+4650_Ks.new /Users/cjimenez/Documents/PHD/data/unlens_liris/SDSSJ1013+4650/SDSSJ1013+4650_Ks_astro.fits --write_tan_si
        #run swarpt to fix SWC
            #swarp @files.lst -c default.swarp -COPY_KEYWORDS=INSTRUME,TELESCOP,OBSERVAT,PIXSCALE,FILTER,FILTER1,LIRF2NAM,OBJEC


