from astropy.io import fits
import numpy as np

def setupImage(image_lst_path):

    images = np.loadtxt(image_lst_path,"str")

    for img in images:
        header = fits.getheader(img,0)
        obtype = header["OBSTYPE"]
        imagetype = header["IMAGETYP"] #
        objName = header["CAT-NAME"] #object cat
        band = header["LIRF2NAM"] #filter name
        ra = header["RA"]
        dec = header["DEC"]
        date=header["DATE-OBS"]

        print(obtype,imagetype,objName,band,ra,dec,date,img)





if __name__ == "__main__":

    path = "/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/liris_2018_06_23/list_image.lst"
    pathcecilia ="/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/iris/cecilia/files.lst"
    pathraine = "/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/iris/raine/files.lst"

    #setupImage(pathcecilia)
    setupImage(pathcecilia)