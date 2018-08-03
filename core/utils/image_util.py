from io import BytesIO
import matplotlib.pyplot as plt
import base64
import matplotlib
from astropy import utils, io
import numpy as np
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.colors as mcolors
from astropy.visualization import (MinMaxInterval, SqrtStretch,ImageNormalize)


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
