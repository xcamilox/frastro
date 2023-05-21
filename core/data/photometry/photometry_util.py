from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from photutils import SkyCircularAperture
from photutils import aperture_photometry
from astropy import wcs

class PhotometryUtil():

    @staticmethod
    def doPhotometry(ra,dec,image,r_arcsec_aperture=5.):
        unit = "deg"
        positions = SkyCoord(ra=ra, dec=dec, unit=unit)
        apertures = SkyCircularAperture(positions, r=r_arcsec_aperture * u.arcsec)

        file=fits.open(image)
        data=file[0].data
        header=file[0].header
        wcs_img = wcs.WCS(header)
        phot_table = aperture_photometry(data, apertures,wcs=wcs_img)
        return phot_table
