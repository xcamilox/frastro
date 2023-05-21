from astropy.io import fits

from astropy import wcs
import numpy as np
from reproject import reproject_interp
import matplotlib.pyplot as plt
from astropy.wcs import WCS


class FITSFile():

    @staticmethod
    def open(filepath):
        hdul = fits.open(filepath)
        return hdul

    @staticmethod
    def header(filepath):
        file = FITSFile.open(filepath)
        headers = []
        for hd in file:
            headers.append(hd.header)
        return headers

    @staticmethod
    def getImage(filepath,hdul_idx=None):
        file=FITSFile.open(filepath)
        if hdul_idx!=None:
            return file[hdul_idx]
        else:
            images=[]
            if len(file) > 0:
                for hdl in file:
                    if "NAXIS" in hdl.header:
                        if hdl.header["NAXIS"] > 0:
                            images.append(hdl)
                        else:
                            continue
            return images



    @staticmethod
    def cropToFile(reference_file_fits,file_fits):
        reference = FITSFile.open(reference_file_fits)
        reference_header=reference[1].header+reference[0].header
        reference_data = reference[1].data

        size_reference = (reference_header["NAXIS1"],reference_header["NAXIS2"])

        file = FITSFile.open(file_fits)
        file_header = file[1].header+file[0].header
        file_data = file[1].data
        size_file = (file_header["NAXIS1"], file_header["NAXIS2"])
        print(size_file)

        wcs_reference= wcs.WCS(reference_header)
        wcs_file = wcs.WCS(file_header)

        reference_x1y1 = wcs_reference.wcs_pix2world(0,0,0)
        reference_x2y2 = wcs_reference.wcs_pix2world(size_reference[0], 0, 0)
        reference_x3y3 = wcs_reference.wcs_pix2world(0, size_reference[1], 0)
        reference_x4y4 = wcs_reference.wcs_pix2world(size_reference[0], size_reference[1], 0)



        crop_pix_x1y1 = wcs_file.wcs_world2pix(reference_x1y1[0],reference_x1y1[1],0)
        crop_pix_x2y2 = wcs_file.wcs_world2pix(reference_x2y2[0], reference_x2y2[1], 0)
        crop_pix_x3y3 = wcs_file.wcs_world2pix(reference_x3y3[0], reference_x3y3[1], 0)
        crop_pix_x4y4 = wcs_file.wcs_world2pix(reference_x4y4[0], reference_x4y4[1], 0)


        x_axes = np.array([crop_pix_x1y1[0],crop_pix_x2y2[0],crop_pix_x3y3[0],crop_pix_x4y4[0]])
        y_axes = np.array([crop_pix_x1y1[1], crop_pix_x2y2[1], crop_pix_x3y3[1], crop_pix_x4y4[1]])

        centerX = ((x_axes.max() - x_axes.min())/2)+x_axes.min()
        centerY = ((y_axes.max() - y_axes.min()) / 2) + y_axes.min()

        coordinate = wcs_file.wcs_pix2world(centerX,centerY,0)

        new_data=file_data[int(y_axes.min()):int(y_axes.max()),int(x_axes.min()):int(x_axes.max())]


        pixposition_x = int(x_axes.min())
        pixposition_y = int(y_axes.min())

        print(pixposition_x,pixposition_y)

        coordinate_loc = wcs_file.wcs_pix2world(pixposition_x, pixposition_y, 0)
        #print(float(coordinate[0])-float(coordinate_loc[0]))
        #print(float(coordinate[1]) - float(coordinate_loc[1]))

        file_header["CRVAL1"] = float(coordinate_loc[0])
        file_header["CRVAL2"] = float(coordinate_loc[1])
        file_header["CRPIX1"] = 0
        file_header["CRPIX2"] = 0

        file_header["NAXIS1"]=new_data.shape[0]
        file_header["NAXIS2"]=new_data.shape[1]
        path="/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/hotpants/crop.fits"
        FITSFile.saveFile(path,npa_data=new_data,header=file_header)


    @staticmethod
    def saveFile(path,npa_data=[],header=[],overwrite=True):
        hdu=fits.PrimaryHDU(data=npa_data,header=header)
        hdul = fits.HDUList([hdu])
        hdul.writeto(path,overwrite=overwrite)




    @staticmethod
    def scaleImage(path,outpath, min, max,deg=1,hdu=0):
        file=FITSFile.open(path)
        data=file[hdu].data
        image_min=data.min()
        image_max=data.max()
        x=[image_min,image_max]
        if min == "file":
            y = [image_min,max]
        else:
            y = [min, max]
        fit = np.poly1d(np.polyfit(x, y, deg))

        newdata=fit(data)


        header=file[hdu].header
        FITSFile.saveFile(outpath,npa_data=newdata,header=header)



    @staticmethod
    def trimImage(path,rec_size_pixel,ra_dec_refToCut):
        file = FITSFile.open(path)
        data = file[0].data
        header = file[0].header




    @staticmethod
    def reprojectionImage(reference_path,target_path):
        ref = fits.open(reference_path)
        target = fits.open(target_path)

        data=target[0]
        new_data, footprint = reproject_interp(data,ref[0].header)

        ax1 = plt.subplot(1, 2, 1, projection=WCS(ref[0].header))
        ax1.imshow(new_data, origin='lower', vmin=-2.e-4, vmax=5.e-4)
        ax1.coords.grid(color='white')
        ax1.coords['ra'].set_axislabel('Right Ascension')
        ax1.coords['dec'].set_axislabel('Declination')
        ax1.set_title('Reprojected MSX band E image')

        ax2 = plt.subplot(1, 2, 2, projection=WCS(ref[0].header))
        ax2.imshow(ref[0].data, origin='lower', vmin=-2.e-4, vmax=5.e-4)
        ax2.coords.grid(color='white')
        ax1.coords['ra'].set_axislabel('Right Ascension')
        ax1.coords['dec'].set_axislabel('Declination')
        ax2.coords['dec'].set_axislabel_position('r')
        ax2.coords['dec'].set_ticklabel_position('r')
        ax2.set_title('MSX band E image footprint')
        plt.show()

        fits.writeto('/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181109_26_4_1_1_rotate.fits', new_data, ref[0].header, clobber=True)



if __name__ == "__main__":

    #des 0.263 arcsecond/pixel scale
    #irac 1.213 arcsecond/pixel scale
    ref = "/Users/cjimenez/Documents/PHD/data/tmp/37.304500000000004_1.1664916666666683/des/des_band_g.fits"
    file = "/Users/cjimenez/Documents/PHD/data/tmp/37.304500000000004_1.1664916666666683/spitzer/swarp/irac2_37.3045_1.16649_swarp.fits"

    ref = "/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/hst/j9op03010_drz.fits"
    file = "/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/decals/image-decam-215899-N29-z.fits"

    ref ="/Users/cjimenez/Documents/PHD/data/bg1429_not_gtc/not_lyman_alpha_bg1429_astro_gaia_0p08.fits"
    file="/Users/cjimenez/Documents/PHD/data/bg1429_not_gtc/sdss1429_g_band_gtc_330sec_ipf_astro_gaia_dr1_using_decals_grz_0p1rms.fits"
    ref = "/Users/cjimenez/Documents/PHD/data/bg1429_not_gtc/gtctrim3.fits"

    #SLACSJ1538+5817
    #decals
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ1538+5817/decals/legacysurvey-2345p582-image-r.fits"

    #liverpool
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ1538+5817/h_e_20180717_47_1_1_1.fits"
    #newFile = FITSFile.scaleImage(ref, ref + "_rescale.fits", "file", 50000, 1, hdu=0)
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ1538+5817/h_e_20180801_27_1_1_1.fits"
    #newFile = FITSFile.scaleImage(ref, ref + "_rescale.fits", "file", 50000, 1, hdu=0)
    #hst
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ1538+5817/SLACSJ1538+5817_10587_53834_F814W_1.fits"
    #newFile = FITSFile.scaleImage(ref, ref + "_rescale.fits", "file", 50000, 1, hdu=1)
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/hst/SLACSJ0330-0020_10886_53994_F814W_4.fits"
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/decals/legacysurvey-0526m002-image-r.fits"
    #image=FITSFile.getImage(ref,hdul_idx=1)
    #FITSFile.saveFile(ref+"_1.fits",image.data,image.header)

    file="/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/liris/SDSSJ122040.72+084238.1/SDSSJ122040.72+084238.1_bandH.fits"
    file_out = "/Users/cjimenez/Documents/PHD/data/unlens_laes_rui/liris/SDSSJ122040.72+084238.1/SDSSJ122040.72+084238.1_bandH_scale.fits"

    #newFile = FITSFile.scaleImage(file, file_out, -0.9396, 50000, 1)

    #FITSFile.reprojectionImage(ref,file)

    file_header="/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_obs/2018-11-06/h_e_20181105_27_4_1_1.fits"

    file_data = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/pydia/nov06/d_h_e_20181105_27_4_1_1.fits"
    out_data = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/pydia/nov06/f_d_h_e_20181105_27_4_1_1.fits"


    """

    data_header=FITSFile.open(file_header)
    data_data = FITSFile.open(file_data)

    FITSFile.saveFile(out_data,data_data[0].data,data_header[0].header)
    """

    """
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180814_25_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180816_75_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180816_75_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180816_75_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180822_91_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180822_91_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180822_91_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180906_96_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180906_96_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180906_96_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180917_78_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180917_78_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180917_78_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181005_75_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181005_75_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181005_75_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181016_70_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181016_70_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181016_70_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181016_70_4_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181016_70_5_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_27_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_27_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_27_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_27_4_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_28_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_28_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_28_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_28_4_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181105_28_5_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181109_26_1_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181109_26_2_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181109_26_3_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181109_26_4_1_1.fits
    /Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181109_26_5_1_1.fits
    """


    ref="/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20180917_78_2_1_1.fits"
    file="/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/h_e_20181109_26_4_1_1.fits"
    file="/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1031+3026/decals/cutout_157.8375_30.4494_30.fits"

    image=FITSFile.open(file)
    data=image[0].data
    min=data.min()
    max=data.max()

    file_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1031+3026/pydia/2018-12-10/obs/2018-12-10_1.39_stack.fits"
    image_ref = FITSFile.open(file_ref)
    data_ref = image_ref[0].data
    min_ref = data_ref.min()
    max_ref = data_ref.max()


    file_out="/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1031+3026/decals/resize_decals_r.fits"

    #FITSFile.saveFile(file_out, data + 900, image[0].header)
    FITSFile.scaleImage(file, file_out, "file", max_ref, 0)

