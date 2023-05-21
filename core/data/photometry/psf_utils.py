from astropy.io import fits
from photutils import find_peaks
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from photutils.psf import extract_stars
from astropy.table import Table
from photutils import EPSFBuilder
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt

class PSFUtils():
    def __init__(self):
        pass
    def checkPSFOutput(self,file):
        hdu=fits.open(file)
        data=hdu[1].data[0][0]
        print(data)


    def getPSF(self,image_path):
        file = fits.open(image_path)
        data = file[0].data
        #find bright stars
        peaks_tbl = find_peaks(data, threshold=400.)
        #get sky value
        mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)
        #remove sky
        data -= median_val
        #create a ndarray , is a numpy narray with header iformation as wcs etc.
        nddata = NDData(data=data)
        #create table with stars position
        stars_tbl = Table()
        stars_tbl['x'] = peaks_tbl['x_peak']
        stars_tbl['y'] = peaks_tbl['y_peak']
        #extract stars from image
        stars = extract_stars(nddata, stars_tbl, size=25)

        nrows = 5
        ncols = 5
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20),squeeze = True)
        ax = ax.ravel()
        for i in range(nrows * ncols):
            norm = simple_norm(stars[i], 'log', percent=99.)
            ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')
        plt.show()


        epsf_builder = EPSFBuilder(oversampling=4, maxiters=3,progress_bar = False)
        epsf, fitted_stars = epsf_builder(stars)
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        plt.show()

if __name__=="__main__":
    psf=PSFUtils()
    file="/Users/cjimenez/Documents/PHD/data/liverpool_lens/SLACSJ0330-0020/liverpool_allobs/reduction/psfex/h_e_20181109_26_1_1_1/2018-11-13/catalog.psf"
    ref = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/S4TMJ1051+4439/observations/2019-04-23/science/stack_best.fits"
    psf.getPSF(ref)