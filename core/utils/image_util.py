import os
from io import BytesIO
from os import path

import io

import matplotlib.pyplot as plt
import base64
import matplotlib
from astropy import utils, io
import numpy as np
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.colors as mcolors
from astropy.visualization import (MinMaxInterval, SqrtStretch,ImageNormalize)
from astropy.io import fits
from astropy import units as u
from astropy import coordinates as coords
from astropy import wcs
from astropy.io import fits
from astropy.nddata import Cutout2D

from core.fits.FITSFIle import FITSFile

matplotlib.use('Agg')
class ImageUtils():
    @staticmethod
    def imgNpArrayToBase64(img_np_array,mark=""):
        buffer= BytesIO()
        plt.clf()
        fig1 = plt.gcf()

        w=img_np_array.shape[1]/96.
        h=img_np_array.shape[0]/96. #3.125 # 300pixel = 3.125inches (1in == 2.54cm) 1 inch = 96 pixel
        fig1.set_size_inches(w, h)

        plt.imshow(img_np_array,origin='lower', interpolation='none', cmap=plt.get_cmap("Greys"),norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=img_np_array.max()))

        plt.axis("off")
        radii=100
        if mark != "":
            plt.scatter(mark[0],mark[1],s=radii,facecolors='none', edgecolors='r')
        #plt.show()
        fig1.savefig(buffer,format='png')
        str="data:image/png;base64,{0}"
        base64str=base64.encodebytes(buffer.getvalue()).decode()
        return str.format(base64str)

    @staticmethod
    def saveImageFromNPArray(img_np_array,output_path, mark=""):

        plt.imshow(img_np_array, cmap='gray')
        plt.colorbar()
        plt.show()
        NBINS = 1000
        img_np_array=np.where(np.isnan(img_np_array),-99,img_np_array)
        print(img_np_array.max(), img_np_array.min())
        flatten=img_np_array.flatten()

        histogram = plt.hist(flatten, NBINS)
        plt.show()
        # plt.clf()
        # fig1 = plt.gcf()
        #
        # # Create an ImageNormalize object
        # norm = ImageNormalize(img_np_array, interval=MinMaxInterval(),
        #                       stretch=SqrtStretch())
        #
        # w = img_np_array.shape[1] / 96.
        # h = img_np_array.shape[0] / 96.  # 3.125 # 300pixel = 3.125inches (1in == 2.54cm) 1 inch = 96 pixel
        # fig1.set_size_inches(w, h)
        #
        # plt.imshow(img_np_array, origin='lower', vmin=100, vmax=6000,interpolation='none', norm=mcolors.PowerNorm(0.3),cmap=plt.get_cmap("Greys"))#norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),cmap=plt.get_cmap("Greys"))
        #
        # #plt.axis("off")
        # radii = 100
        # if mark != "":
        #     plt.scatter(mark[0], mark[1], s=radii, facecolors='none', edgecolors='r')
        # plt.show()
        # fig1.savefig(output_path, format='png')
        # return output_path

    @staticmethod
    def saveImageFromFits(fits_path,mark=""):
        data = ImageUtils.getImageFromFits(fits_path)

    def getFirstImageHDUIndex(file_path, beg=0):
        file = fits.open(file_path)

        header_idx = beg
        if len(file) > 1:
            for idx, h in enumerate(file):
                if idx >= beg:
                    header = h.header
                    if ("NAXIS" in header and header["NAXIS"] > 0) and (
                            "NAXIS1" in header and header["NAXIS1"] > 0 and "NAXIS2" in header and header["NAXIS2"] > 0):
                        header_idx = idx
                        break
        return header_idx

    @staticmethod
    def getFirstImageFromFitsFile(file_path):
        idx=ImageUtils.getFirstImageHDUIndex(file_path)
        file = fits.open(file_path)
        image=file[idx]
        return image


    @staticmethod
    def imcopy(image_path,output,origin=(0,0),crop=('full','full'),addcounts=0,n_ext=None):
        # Change these values to your desired filenames
        test_data = image_path

        ext_list=FITSFile.getImage(test_data)
        # Create iterable list of tuples to feed into Cutout2D,
        # seperate list for extensions with wcs, as feeding the wcs
        # back into the FITS file takes more work.
        if n_ext!=None:
            ext_list=[ext_list[n_ext]]

        for indx,ext in enumerate(ext_list):
            orig_wcs = wcs.WCS(ext.header)
            img_size = [ext.data.shape[0], ext.data.shape[1]]
            if crop[0]!="full":
                img_size[0]=crop[0]
            if crop[1]!="full":
                img_size[1]=crop[1]

            cutout = Cutout2D(ext.data, origin, img_size, wcs=orig_wcs)
            print(cutout.shape)
            ext.data = cutout.data
            ext.header.update(cutout.wcs.to_header())
            FITSFile.saveFile(output,cutout.data+addcounts,ext.header)
        return cutout.data+addcounts

    @staticmethod
    def getImageFromFits(url):
        try:
            idx=url.index("html")
            image = io.fits.getdata(utils.data.download_file(url, cache=True, show_progress=False, timeout=120))
        except ValueError:
            image = io.fits.getdata(url)
        return image

    @staticmethod
    def getBase64FromFitsURL(url,mark=""):
        base64 =ImageUtils.imgNpArrayToBase64(ImageUtils.getImageFromFits(url),mark=mark)
        return base64


    @staticmethod
    def getImageFromUrl(url):
        local_Path= utils.data.download_file(url, cache=True, show_progress=False, timeout=120)
        image = io.fits.getdata(local_Path)

    @staticmethod
    def getPixelFromRADEC():
        pass

    @staticmethod
    def getCoordinateFromPixel(fits_image_path, x,y):
        hdulist = fits.open(fits_image_path)
        data = hdulist[0].data
        w = WCS(hdulist[0].header)
        print(data.shape)
        ymax = data.shape[0]
        xmax = data.shape[1]

        world = w.all_pix2world(x, y, 1)
        skycoor = coords.SkyCoord(ra=world[0], dec=world[1], unit="deg")

        return skycoor

    @staticmethod
    def getCenterCoordinate(fits_image_path):
        hdulist = fits.open(fits_image_path)
        data = hdulist[0].data
        w = WCS(hdulist[0].header)

        ymax = data.shape[0]
        xmax = data.shape[1]

        skycoor = ImageUtils.getCoordinateFromPixel(fits_image_path,xmax / 2, ymax / 2)

        return skycoor

    @staticmethod
    def getImageRadiiSizeArcmin(fits_image_path):
        hdulist = fits.open(fits_image_path)
        data=hdulist[0].data
        w = WCS(hdulist[0].header)
        print(data.shape)
        ymax=data.shape[0]
        xmax = data.shape[1]


        pixcrd = np.array([[0, 0] ,[xmax/2,ymax/2],[xmax,ymax]] , np.float_)
        print(pixcrd)
        world = w.all_pix2world(pixcrd, 1)

        distance1 = np.sqrt(((world[0][0] - world[2][0]) ** 2) + ((world[0][1] - world[2][1]) ** 2))
        distance2 = np.sqrt(((world[0][0] - world[1][0]) ** 2) + ((world[0][1] - world[1][1]) ** 2))
        distance3 = np.sqrt(((world[1][0] - world[2][0]) ** 2) + ((world[1][1] - world[2][1]) ** 2))


        d1=distance1*u.deg
        d2 = distance2 * u.deg
        d3 = distance3 * u.deg
        rad=np.max([d2.to("arcmin").value,d3.to("arcmin").value])

        return rad



    @staticmethod
    def getPixelFromCoordinate(fits_image_path,ra,dec):
        hdulist = fits.open(fits_image_path)
        data = hdulist[0].data
        w = WCS(hdulist[0].header)
        skycoor = coords.SkyCoord(ra=ra, dec=dec, unit="deg")

        pixels=w.wcs_world2pix(np.array([[skycoor.ra.deg,skycoor.dec.deg]],np.float_),1)
        return pixels

    @staticmethod
    def getFilterImage(fits_image_path):
        hdul = fits.open(fits_image_path)
        header = hdul[0].header
        filters = []
        for idx, hdr in enumerate(header.keys()):
            if "FILTER" in hdr and "clear" not in str(header[hdr]).lower():
                filters.append([hdr, header[hdr]])
            elif "LIRF2NAM" in hdr and "clear" not in str(header[hdr]).lower():
                filters.append([hdr, header[hdr]])
        return filters


    @staticmethod
    def getImageThumbnail(file_path,ra,dec,image_index=0,size=(200,200)):
        ref_file_fits =fits.open(file_path)[image_index]
        skycoor = coords.SkyCoord(ra=ra, dec=dec, unit="deg")

        # ref image data
        orig_ref_wcs = wcs.WCS(ref_file_fits.header)

        pixel_ref = orig_ref_wcs.wcs_world2pix(np.array([[skycoor.ra.deg, skycoor.dec.deg]], np.float_), 1)
        origin_ref = (pixel_ref[0][0], pixel_ref[0][1])
        cutout_ref = Cutout2D(ref_file_fits.data, origin_ref, size, wcs=orig_ref_wcs)
        data_ref = cutout_ref.data
        return data_ref

    @staticmethod
    def removeBorders(file,output_name="",border_pix=300):

        image=FITSFile.open(file)
        data=image[0].data.shape
        width = data[0]-border_pix
        height = data[1] - border_pix
        base_path = os.path.dirname(file)
        saveFile = base_path +"/"+ output_name
        ImageUtils.imcopy(file, saveFile, origin=(data[0]/2,data[1]/2),crop=(width,height),addcounts=0)
        return saveFile

    @staticmethod
    def cropAllFiles(self,file,ra,dec,size=(100,100),addcounts=0):
        pixel = ImageUtils.getPixelFromCoordinate(file, ra, dec)
        ImageUtils.imcopy(file, file, (pixel[0][0], pixel[0][1]), size,addcounts=addcounts)

    @staticmethod
    def getBase64FromPyPlotFig(format="png"):

        pic_IObytes = io.BytesIO()
        plt.savefig(pic_IObytes, format=format)
        pic_IObytes.seek(0)
        pic_hash = base64.b64encode(pic_IObytes.read())
        return pic_hash

if __name__ == "__main__":
    fits_file = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/20181005_1/h_e_20181005_75_3_1_1.fits"
    fits_file = "/Users/cjimenez/Documents/PHD/data/tmp/349.4015416666667_12.582483333333334/panstarss/panstars_band_g.fits"
    fits_file = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ0143-1006/hst/SLACSJ0143-1006_12210_55568_F814W_1.fits"
    #rad=ImageUtils.getFirstImageFromFitsFile(fits_file)
    #print(rad)
    path_folder = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/science/"
    #onlyfiles = [f for f in os.listdir(path_folder) if path.isfile(path_folder + "/" + f)]
    #result=[]
    ra=52.55058333333333
    dec=-0.3477499999999987
    file="stacking_12_11_2018.fits"
    #pixel = ImageUtils.getPixelFromCoordinate(path_folder+file, ra, dec)


    #ImageUtils.imcopy(path_folder +file, path_folder + path.splitext(file)[0] + "_crop.fits",
    #                  (pixel[0][0], pixel[0][1]), (100, 100))



    #for file in onlyfiles:
    #    if path.splitext(file)[1] == ".fits" or path.splitext(file)[1] == ".fit":

    #        pixel=ImageUtils.getPixelFromCoordinate(path_folder+file,ra,dec)
    #        ImageUtils.imcopy(path_folder+file,path_folder+path.splitext(file)[0]+"_crop.fits",(pixel[0][0],pixel[0][1]),(100,100))
    #        print(path_folder+path.splitext(file)[0]+"_crop.fits")


    file1 = "/Users/cjimenez/Documents/PHD/data/adfs27/science/adfs_27_fors2_r_band.fits"
    file2 = "/Users/cjimenez/Documents/PHD/data/adfs27/science/adfs_27_fors2_z_band.fits"

    file1 = "/Users/cjimenez/Documents/PHD/data/LCO/coj1m011-fa12-20190514-0118-e91.fits.fz"
    file2 = "/Users/cjimenez/Documents/PHD/data/LCO/coj1m011-fa12-20190514-0119-e91.fits.fz"
    file3 = "/Users/cjimenez/Documents/PHD/data/LCO/coj1m011-fa12-20190514-0120-e91.fits.fz"

    files=[file1,file2,file3]
    for file in files:
        hdulist = fits.open(file, mode='update')
        data = hdulist[1].data
        shape = data.shape
        w = shape[0] - 109
        h = shape[1] - 109


        # telescope="vlt"
        #filter="R_SPECIAL"
        filter="SDSS-R"
        #instrument = "FORS2"
        #hdulist[0].header.set("INSTRUME", instrument)
        # hdulist[0].header.set("TELESCOP", telescope)
        hdulist[1].header.set("FILTER", filter)
        hdulist.flush()
        hdulist.close()

        header=hdulist[1].header
        FITSFile.saveFile(file[:-8]+"fix.fits",data,header)

        #ImageUtils.imcopy(file, file[:-5] + "_crop.fits", (shape[0] / 2, shape[1] / 2), (w, h))