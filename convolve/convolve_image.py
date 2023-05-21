from astropy.wcs import wcs
from scipy import signal
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.ndimage import gaussian_filter

from frastro import PSFexUtil, ImageUtils, FITSFile
from frastro.core.data.photometry.photometry_util import PhotometryUtil
from photutils import create_matching_kernel
from photutils import (HanningWindow, TukeyWindow, CosineBellWindow,
                       SplitCosineBellWindow, TopHatWindow)

from astropy.convolution import convolve

class ConvolveImage():
    def __init__(self):
        pass

    def convolveImage_scipy(self,np_data,np_kernel,mode="same",method="auto"):
        return signal.convolve(np_data, np_kernel, mode=mode, method=method)

    def convolveImage_astropy(self,np_data,np_kernel):
        return convolve(np_data,np_kernel)

    #Convolved Image using both methods, astropy convolve and scipy
    def convolveImage(self,np_data,np_kernel):
        scipy = self.convolveImage_scipy(np_data,np_kernel)
        astropy = self.convolveImage_astropy(np_data, np_kernel)
        return scipy, astropy


    #alpha and beta parameters only are required for turkey,cosinebell, splitcosinebell, and tophat functions
    def createKernel(self,windows_funtion="cosinebell",alpha=0.5,beta=0.3):
        wfuntion = windows_funtion.lower()
        if  wfuntion == "hanning":
            w = HanningWindow()
        if wfuntion == "tukey":
            w = TukeyWindow(alpha=alpha)
        if wfuntion == "cosinebell":
            w = CosineBellWindow(alpha=alpha)
        if wfuntion == "splitcosinebell":
            w = SplitCosineBellWindow(alpha=alpha, beta=beta)
        if wfuntion == "tophat":
            w = TopHatWindow(beta=beta)
        return w

    def convolve_img(image, kernel, padding="same", stride=1):
        image_h = image.shape[0]
        image_w = image.shape[1]

        kernel_h = kernel.shape[0]
        kernel_w = kernel.shape[1]

        h = kernel_h // 2
        w = kernel_h // 2

        if padding == "same":
            padding = (kernel_w - 1) / 2
            size = (int(image_h + padding), int(image_w + padding))

        image_ref = np.zeros(size)
        image_conv = np.zeros(image.shape)
        for i in range(image_h):
            for j in range(image_w):
                image_ref[i + h][j + h] = image[i][j]

        for i in range(image_h - h):
            for j in range(image_w - w):
                sum = 0
                for m in range(kernel_h):
                    for n in range(kernel_w):
                        sum = sum + kernel[m][n] * image_ref[i + m][j + n]

                image_conv[i][j] = sum

        return image_conv



def convolve_image():
    # SLACSJ0157-0056
    ra = 29.495583333333332
    dec = -0.9405833333333322
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-09-19/science/stack_best.fits"
    ref_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-09-19/science/stack_best_crop.fits"
    psf_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-09-19/science/psf_best.fits"
    cat_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-09-19/science/catalog_best_deep.fits"
    zp = 29.81058550612531
    mag = 19.549888610839844
    flux_ref = 12713.88
    # disaline
    # ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-10-17/science/stack_best.fits"
    # ref_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-10-17/science/stack_best_crop.fits"
    # psf_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-10-17/science/psf_best.fits"
    zp1 = 30.381520940420206
    mag1 = 19.457260131835938

    target = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2019-01-01/science/stack_best.fits"
    target_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2019-01-01/science/stack_best_crop.fits"
    psf_target = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2019-01-01/science/psf_best.fits"
    zp2 = 30.304233586279405
    mag2 = 19.44410514831543
    flux_target = 22082.65

    # BELLSJ1221+3806
    ra = 185.466325546679
    dec = 38.1029236161019
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/BELLSJ1221+3806/observations/2019-03-05/science/stack_best.fits"
    ref_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/BELLSJ1221+3806/observations/2019-03-05/science/stack_best_crop.fits"
    psf_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/BELLSJ1221+3806/observations/2019-03-05/science/psf_best.fits"
    zp = 30.2217960357666
    mag = 20.398475646972656
    flux_ref = 8498.216796875
    seeing_ref = 1.439271068572998

    target = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/BELLSJ1221+3806/observations/2019-02-01/science/stack_best.fits"
    target_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/BELLSJ1221+3806/observations/2019-02-01/science/stack_best_crop.fits"
    psf_target = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/BELLSJ1221+3806/observations/2019-02-01/science/psf_best.fits"
    zp2 = 30.256074905395508
    mag2 = 20.425241470336914
    flux_target = 8557.240234375
    seeing_target = 1.4405250549316406

    psf = PSFexUtil()
    ref_output_psf = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0157-0056/observations/2018-09-19/science/psf_best_img.fits"

    # S4TMJ1051+4439

    ra = 162.78920833333333
    dec = 44.65236388888888
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/astrometry/good/swarp_2018-11-20_1.fits"
    ref_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/astrometry/good/2018-11-20_1_crop.fits"
    psf_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2018-11-20/science/psf_best.fits"
    cat_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2018-11-20/science/catalog_best_deep.fits"

    zp = 30.010107046117383
    mag = 17.317863
    flux_ref = 119370.54
    seeing_ref = 1.3008779525756835

    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-04-23/science/stack_best.fits"
    ref_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-04-23/science/stack_best_crop.fits"
    psf_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-04-23/science/psf_best.fits"
    cat_ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-04-23/science/catalog_best_deep.fits"
    zp2 = 30.259170532226562
    mag2 = 17.600826
    flux_target = 115701.07
    seeing_target = 1.0545900106430053

    target = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-05-01/science/stack_best.fits"
    target_crop = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-05-01/science/stack_best_crop.fits"
    psf_target = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-05-01/science/psf_best.fits"

    zp2 = 30.20319175720215
    mag2 = 17.57571029663086
    flux_target = 112458.59375
    seeing_target = 0.8817330121994018

    scale = calcScale(cat_ref, ref, target)

    ref_file = fits.open(ref)[0].data
    ref_file = ref_file - ref_file.mean()
    sky_ref = ref_file.mean()
    sky_max = ref_file.max()
    sky_std_ref = ref_file.std()

    target_file = fits.open(target)[0].data
    target_file = target_file - target_file.mean()
    target_sky = target_file.mean()
    target_max = target_file.max()
    target_std_ref = target_file.std()

    print("ref", "sky", sky_ref, "max", sky_max, "std", sky_std_ref)
    print("target", "sky", target_sky, "max", target_max, "std", target_std_ref)

    psf = PSFexUtil()
    ref_output_psf = psf_ref[:-5] + "_img.fits"
    psf.saveSinglePSFFile(psf_ref, ref_output_psf)

    target_output_psf = psf_target[:-5] + "_img.fits"
    psf.saveSinglePSFFile(psf_target, target_output_psf)

    psf_ref_org = fits.open(ref_output_psf)[0].data
    psf_ref_img = psf_ref_org
    psf_target_img = fits.open(target_output_psf)[0].data

    fov = 80
    pixel = ImageUtils.getPixelFromCoordinate(ref, ra, dec)
    crop_ref = ImageUtils.imcopy(ref, ref_crop, (pixel[0][0], pixel[0][1]), (fov, fov))
    pixel_target = ImageUtils.getPixelFromCoordinate(target, ra, dec)
    crop_target = ImageUtils.imcopy(target, target_crop, (pixel_target[0][0], pixel_target[0][1]), (fov, fov))

    # crop_ref=crop_ref-crop_ref.mean()
    ref_total_flux = crop_ref.sum()

    mag_ref = -2.5 * np.log10(ref_total_flux / zp)

    crop_ref_mean = crop_ref.mean()
    crop_ref_std = crop_ref.std()

    # crop_target = crop_target-crop_target.mean()
    target_total_flux = crop_target.sum()

    mag_target = -2.5 * np.log10(target_total_flux / zp2)
    crop_target_mean = crop_target.mean()
    crop_target_std = crop_target.std()

    ref_img = fits.open(ref_crop)[0].data
    target_img = fits.open(target_crop)[0].data

    sky_refcrop = crop_ref.mean()
    sky_targetcrop = crop_target.mean()

    std_refcrop = crop_ref.std()
    std_targetcrop = crop_target.std()

    print("STD ref", std_refcrop)
    print("STD target", std_targetcrop)

    w1 = HanningWindow()
    w2 = TukeyWindow(alpha=0.5)
    w3 = CosineBellWindow(alpha=0.5)
    w4 = SplitCosineBellWindow(alpha=0.4, beta=0.3)
    w5 = TopHatWindow(beta=0.4)

    kernel_HanningWindow = create_matching_kernel(psf_ref_img, psf_target_img, window=w1)


    kernel_TukeyWindow = create_matching_kernel(psf_ref_img, psf_target_img, window=w2)


    kernel_CosineBellWindow = create_matching_kernel(psf_ref_img, psf_target_img, window=w3)


    kernel_SplitCosineBellWindow = create_matching_kernel(psf_ref_img, psf_target_img, window=w4)


    kernel_TopHatWindow = create_matching_kernel(psf_ref_img, psf_target_img, window=w5)


    diff_simple = crop_ref - crop_target

    sky_diff = sky_targetcrop - sky_refcrop
    # crop_ref += sky_diff


    psf_conv = signal.convolve(psf_ref_img, psf_target_img, mode='same', method="auto")

    gaussian_filter_diff = seeing_target - seeing_ref
    blur_conv = gaussian_filter(crop_ref, sigma=gaussian_filter_diff)

    img_conv = signal.convolve(crop_ref, psf_conv , mode='same', method="auto")
    target_conv = signal.convolve(crop_target, psf_conv, mode='same', method="auto")



    header_ref_crop = FITSFile.header(ref_crop)[0]
    header_target_crop = FITSFile.header(target_crop)[0]

    conv_ref_path = ref_crop[:-5] + "_conv.fits"
    conv_target_path = target_crop[:-5] + "_conv.fits"

    # save crop reference image convolved with target psf
    FITSFile.saveFile(conv_ref_path, img_conv, header_ref_crop)

    # save crop target image convolved with ref psf
    FITSFile.saveFile(conv_target_path, target_conv, header_target_crop)

    scale = calcScale(cat_ref, conv_ref_path, conv_target_path)
    print("scale",scale)

    w_ref = wcs.WCS(header_ref_crop)
    w_target = wcs.WCS(header_target_crop)
    diff=np.zeros(img_conv.shape)
    for i in range(img_conv.shape[0]):
        for j in range(img_conv.shape[1]):
            ra_dec = w_ref.wcs_pix2world([[i,j]],1)
            pixcrd2 = w_target.wcs_world2pix(ra_dec, 1)[0]
            #Dij=Iij-Mij
            #print(i,j, int(round(pixcrd2[0])),int(round(pixcrd2[1])))
            diff[i,j]=target_conv[int(round(pixcrd2[0])),int(round(pixcrd2[1]))]-(img_conv[i,j]*scale)

    diference = img_conv - target_conv
    diference_out = img_conv - crop_target




    conv_ref_path_diff = conv_ref_path[:-5] + "_diff.fits"
    # residuos from conv ref - conv target
    FITSFile.saveFile(conv_ref_path_diff, diff, header_ref_crop)


    psf_conv_path = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/astrometry/good/psf_conv.fits"

    header_con_psf = FITSFile.header(psf_target)[0]
    FITSFile.saveFile(psf_conv_path, psf_conv, header_con_psf)

    print(conv_target_path)
    print(conv_ref_path)
    print(conv_ref_path_diff)
    print(psf_conv_path)



#made plots
#def makePlot():

    print(conv_ref_path, conv_target_path, conv_ref_path_diff)

    plt.clf()
    plt.imshow(psf_ref_img)
    plt.title("psf ref")
    print("sum", psf_ref_img.sum())
    plt.colorbar()
    plt.show()
    plt.imshow(psf_target_img)
    plt.title("psf target")
    plt.colorbar()
    plt.show()
    # plt.imshow(psf_conv)
    # plt.title("psf convo ref-target")
    # plt.colorbar()
    # plt.show()

    plt.imshow(crop_ref)
    plt.title("Reference Image")
    plt.colorbar()
    plt.show()

    plt.imshow(blur_conv)
    plt.title("Image gausean")
    plt.colorbar()
    plt.show()

    plt.imshow(crop_target)
    plt.title("Target Image")
    plt.colorbar()
    plt.show()

    plt.imshow(img_conv)
    plt.title("Ref convol psf-target")
    plt.colorbar()
    plt.show()

    plt.imshow(diference_out)
    plt.title("Diference ref-conv-psf-target - target")
    plt.colorbar()
    plt.show()

    plt.imshow(diff)
    plt.title("Diference ref-conv-psf-target - target-conv-psf-ref")
    plt.colorbar()
    plt.show()

    """
    plt.imshow(diff_simple)
    plt.colorbar()
    plt.show()
    """


def readCatalog(cat_path):
    hdul = FITSFile.open(cat_path)
    catalog = hdul[1].data
    return catalog


def calcScale(cat_ref, img_ref, img_target):
    print(cat_ref)
    ref_cat = readCatalog(cat_ref)
    filter = ref_cat["CLASS_STAR"] > 0.8
    catalog_ref = ref_cat[filter]

    filter = catalog_ref["MAG_AUTO"] > 15.
    catalog_ref = catalog_ref[filter]

    filter = catalog_ref["BACKGROUND"] < 0.09
    catalog_ref = catalog_ref[filter]

    # image_ref=FITSFile.open(img_ref)
    ref_photometry = PhotometryUtil.doPhotometry(ra=catalog_ref["ALPHA_J2000"], dec=catalog_ref["DELTA_J2000"],
                                                 image=img_ref, r_arcsec_aperture=5)

    # image_target = FITSFile.open(img_target)
    target_photometry = PhotometryUtil.doPhotometry(ra=catalog_ref["ALPHA_J2000"], dec=catalog_ref["DELTA_J2000"],
                                                    image=img_target, r_arcsec_aperture=5)

    ra = 162.78920833333333
    dec = 44.65236388888888

    # image_ref=FITSFile.open(img_ref)
    source_ref_photometry = PhotometryUtil.doPhotometry(ra=[ra], dec=[dec],
                                                 image=img_ref, r_arcsec_aperture=5)

    # image_target = FITSFile.open(img_target)
    source_target_photometry = PhotometryUtil.doPhotometry(ra=[ra], dec=[dec],
                                                    image=img_target, r_arcsec_aperture=5)



    xcenter = ref_photometry['xcenter'] - target_photometry['xcenter']
    ycenter = ref_photometry['ycenter'] - target_photometry['ycenter']
    aperture_sum =  target_photometry['aperture_sum'] / ref_photometry['aperture_sum']

    print(aperture_sum)
    print(np.mean(aperture_sum),np.median(aperture_sum),np.std(aperture_sum))
    list =target_photometry['aperture_sum'] / ref_photometry['aperture_sum']


    scale = np.mean(list)
    scale_obj=source_target_photometry['aperture_sum'] / source_ref_photometry['aperture_sum']
    print(np.mean(xcenter))
    print(np.mean(ycenter))
    print(np.mean(ref_photometry['aperture_sum']), np.mean(target_photometry['aperture_sum']), np.mean(aperture_sum))
    return scale_obj





if __name__ == "__main__":
    convolve_image()
