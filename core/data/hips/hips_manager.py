from astropy.coordinates import SkyCoord
from hips import WCSGeometry
from frastro import HipsSkyMaps
from hips import make_sky_image
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
"""
    implementation base on hips api
    https://hips.readthedocs.io/en/latest/about.html
"""
class HipsManager():

    def __init__(self):
        #projections support http://docs.astropy.org/en/stable/wcs/#supported-projections

        ra = 214.641841
        dec = 51.2152

        geometry = WCSGeometry.create(
            skydir=SkyCoord(ra, dec, unit='deg'),
            width=100, height=100, fov="0.03 deg",
            coordsys='icrs', projection='SIN',
        )
        coor=SkyCoord(ra, dec, unit='deg', frame='galactic')
        print(coor.to_string('hmsdms'))
        # geometry = WCSGeometry.create(
        #     skydir=SkyCoord(280.4652, -32.8884, unit='deg', frame='galactic'),
        #     width=1600, height=1000, fov="2 deg",
        #     coordsys='icrs', projection='SIN',
        # )

#        hips_survey = HipsSkyMaps.getMap("CFHT")

        hips_survey = 'CDS/P/SPITZER/MIPS2'
        hips_survey = 'CDS/P/CFHTLS/W/r'
        #"P/CFHTLS/W/r) 14 18 42.36 +51 25 56.7 0.5 deg"
        #result = make_sky_image(geometry, hips_survey, 'jpeg')

        result = make_sky_image(geometry, hips_survey, 'fits',progress_bar=False)
        image=result.image
        plt.clf()
        fig1 = plt.gcf()
        w = h = 3.125  # 300pixel = 3.125inches (1in == 2.54cm) 1 inch = 96 pixel
        fig1.set_size_inches(w, h)

        plt.imshow(result.image, origin='lower', interpolation='none', cmap=plt.get_cmap("Greys"),
                   norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=result.image.max()))
        plt.axis("off")

        #plt.figure(figsize=(20, 10))

        result.plot()
        plt.show()