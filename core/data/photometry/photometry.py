from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from photutils import datasets
from frastro import FITSFile
from photutils import aperture_photometry
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils import SkyCircularAperture
from astropy import utils, io
from astropy import wcs
import math
from frastro import VOTableUtil,WaveLenghtCover
from astropy.table import Column, Table

class Photometry():

    def __init__(self):
        pass

    def sourceDetection(self,fits_file_path):


        data = io.fits.getdata(fits_file_path, 0)
        header = io.fits.getheader(fits_file_path,0)

        ra=17.689916666666665
        dec=-5.027574999999997

        # header= io.fits.
        w = wcs.WCS(header)

        mean, median, std = sigma_clipped_stats(data, sigma=10.0, iters=5)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std)
        sources = daofind(data - median)
        print((mean, median, std))

        radec= w.wcs_pix2world(sources["xcentroid"],sources["ycentroid"],0)

        sources.add_column(Column(name='ra', data=radec[0]))
        sources.add_column(Column(name='dec', data=radec[1]))

        VOTableUtil.saveFromTable(sources, fits_file_path+"votable_cat_10sigma.xml")
        """
        pixcrd2 = w.wcs_world2pix([[ra, dec]], 1)
        position = SkyCoord(ra, dec, unit='deg', frame='icrs')

        radii = [3., 4., 5.]
        apertures = [SkyCircularAperture(position, r=r * u.arcsec) for r in radii]
        phot_table = aperture_photometry(data, apertures,wcs=w)
        print(phot_table)
        print(pixcrd2)
        """

    def Magnitud_Absolute_UV(self,magnitud_ab,DL_pc,redshift,mu_magnification):
        MAV=magnitud_ab-5*math.log10(DL_pc)+5+2*math.log10(1+redshift)+2.5*math.log10(mu_magnification)
        return MAV

    def UVBeta_slope(self,ma1_ab,ma2_ab,lab_m1,lab_m2):
        BetaUV =-((ma1_ab-ma2_ab)/(2.5*math.log10(lab_m1/lab_m2)))-2
        return BetaUV

    def UVLuminosity(self,mag_ab,redshift,mu_magnification=1):
        mag_ab=float(mag_ab)
        redshift=float(redshift)
        DL_cm=self.luminosity_distance(redshift)["DL_cm"]
        const=9.5214*10**36
        L_uv= ((math.pi*DL_cm**2*const)/(1+redshift))*(1/mu_magnification)*10**-0.4*(mag_ab+48.6)
        return L_uv

    def Auv(self,UVbeta_slope):
        auv=4.43+1.99*UVbeta_slope
        return auv

    def SFRuv(self,uv_luminosity,a_uv):
        sfr=1.48*10**-28*uv_luminosity*10**(0.4*a_uv)
        return sfr


    def luminosity_distance(self,redshift, Ho=69.6, Omega_m=0.286, Omega_vac=0.714,verbose=0):
        #A Cosmology Calculator for the World Wide Web
        #Edward L. Wright (UCLA)
        #arXiv:astro-ph/0609593

        #input values = redshift, Ho=75, Omega_m, Omega_vac
        #ouput values = age at z, distance in Mpc, kpc/arcsec, apparent to abs mag conversion

        """
        # if no values, assume Benchmark Model, input is z
        if length == 2:
            if float(sys.argv[1 + verbose]) > 100:
                z = float(sys.argv[1 + verbose]) / 299790.  # velocity to redshift
            else:
                z = float(sys.argv[1 + verbose])  # redshift
            H0 = 75  # Hubble constant
            WM = 0.3  # Omega(matter)
            WV = 1.0 - WM - 0.4165 / (H0 * H0)  # Omega(vacuum) or lambda

        # if one value, assume Benchmark Model with given Ho
        elif length == 3:
            z = float(sys.argv[1 + verbose])  # redshift
            H0 = float(sys.argv[2 + verbose])  # Hubble constant
            WM = 0.3  # Omega(matter)
            WV = 1.0 - WM - 0.4165 / (H0 * H0)  # Omega(vacuum) or lambda

        # if Univ is Open, use Ho, Wm and set Wv to 0.
        elif length == 4:
            z = float(sys.argv[1 + verbose])  # redshift
            H0 = float(sys.argv[2 + verbose])  # Hubble constant
            WM = float(sys.argv[3 + verbose])  # Omega(matter)
            WV = 0.0  # Omega(vacuum) or lambda

        # if Univ is General, use Ho, Wm and given Wv
        elif length == 5:
            z = float(sys.argv[1 + verbose])  # redshift
            H0 = float(sys.argv[2 + verbose])  # Hubble constant
            WM = float(sys.argv[3 + verbose])  # Omega(matter)
            WV = float(sys.argv[4 + verbose])  # Omega(vacuum) or lambda

        # or else fail
        else:
            print
            'need some values or too many values'
            sys.exit()
        """
        # initialize constants

        z=float(redshift)
        H0=float(Ho)
        WM=float(Omega_m)
        WV=float(Omega_vac)



        WR = 0.  # Omega(radiation)
        WK = 0.  # Omega curvaturve = 1-Omega(total)
        c = 299792.458  # velocity of light in km/sec
        Tyr = 977.8  # coefficent for converting 1/H into Gyr
        DTT = 0.5  # time from z to now in units of 1/H0
        DTT_Gyr = 0.0  # value of DTT in Gyr
        age = 0.5  # age of Universe in units of 1/H0
        age_Gyr = 0.0  # value of age in Gyr
        zage = 0.1  # age of Universe at redshift z in units of 1/H0
        zage_Gyr = 0.0  # value of zage in Gyr
        DCMR = 0.0  # comoving radial distance in units of c/H0
        DCMR_Mpc = 0.0
        DCMR_Gyr = 0.0
        DA = 0.0  # angular size distance
        DA_Mpc = 0.0
        DA_Gyr = 0.0
        kpc_DA = 0.0
        DL = 0.0  # luminosity distance
        DL_Mpc = 0.0
        DL_Gyr = 0.0  # DL in units of billions of light years
        V_Gpc = 0.0
        a = 1.0  # 1/(1+z), the scale factor of the Universe
        az = 0.5  # 1/(1+z(object))

        h = H0 / 100.
        WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species, T0 = 2.72528
        WK = 1 - WM - WR - WV
        az = 1.0 / (1 + 1.0 * z)
        age = 0.
        n = 1000  # number of points in integrals
        for i in range(n):
            a = az * (i + 0.5) / n
            adot = math.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
            age = age + 1. / adot

        zage = az * age / n
        zage_Gyr = (Tyr / H0) * zage
        DTT = 0.0
        DCMR = 0.0

        # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
        for i in range(n):
            a = az + (1 - az) * (i + 0.5) / n
            adot = math.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
            DTT = DTT + 1. / adot
            DCMR = DCMR + 1. / (a * adot)

        DTT = (1. - az) * DTT / n
        DCMR = (1. - az) * DCMR / n
        age = DTT + zage
        age_Gyr = age * (Tyr / H0)
        DTT_Gyr = (Tyr / H0) * DTT
        DCMR_Gyr = (Tyr / H0) * DCMR
        DCMR_Mpc = (c / H0) * DCMR

        # tangential comoving distance

        ratio = 1.00
        x = math.sqrt(abs(WK)) * DCMR
        if x > 0.1:
            if WK > 0:
                ratio = 0.5 * (math.exp(x) - math.exp(-x)) / x
            else:
                ratio = math.sin(x) / x
        else:
            y = x * x
            if WK < 0: y = -y
            ratio = 1. + y / 6. + y * y / 120.
        DCMT = ratio * DCMR
        DA = az * DCMT
        DA_Mpc = (c / H0) * DA
        kpc_DA = DA_Mpc / 206.264806
        DA_Gyr = (Tyr / H0) * DA
        DL = DA / (az * az)
        DL_Mpc = (c / H0) * DL
        DL_Gyr = (Tyr / H0) * DL

        # comoving volume computation

        ratio = 1.00
        x = math.sqrt(abs(WK)) * DCMR
        if x > 0.1:
            if WK > 0:
                ratio = (0.125 * (math.exp(2. * x) - math.exp(-2. * x)) - x / 2.) / (x * x * x / 3.)
            else:
                ratio = (x / 2. - math.sin(2. * x) / 4.) / (x * x * x / 3.)
        else:
            y = x * x
            if WK < 0: y = -y
            ratio = 1. + y / 5. + (2. / 105.) * y * y
        VCM = ratio * DCMR * DCMR * DCMR / 3.
        V_Gpc = 4. * math.pi * ((0.001 * c / H0) ** 3) * VCM

        if verbose == 1:
            print('For H_o = ' + '%1.1f' % H0 + ', Omega_M = ' + '%1.2f' % WM + ', Omega_vac = ')
            print('%1.2f' % WV + ', z = ' + '%1.3f' % z)
            print('It is now ' + '%1.1f' % age_Gyr + ' Gyr since the Big Bang.')
            print('The age at redshift z was ' + '%1.1f' % zage_Gyr + ' Gyr.')
            print('The light travel time was ' + '%1.1f' % DTT_Gyr + ' Gyr.')
            print('The comoving radial distance, which goes into Hubbles law, is')
            print('%1.1f' % DCMR_Mpc + ' Mpc or ' + '%1.1f' % DCMR_Gyr + ' Gly.')
            print('The comoving volume within redshift z is ' + '%1.1f' % V_Gpc + ' Gpc^3.')
            print('The angular size distance D_A is ' + '%1.1f' % DA_Mpc + ' Mpc or')
            print('%1.1f' % DA_Gyr + ' Gly.')
            print('This gives a scale of ' + '%.2f' % kpc_DA + ' kpc/".')
            print('The luminosity distance D_L is ' + '%1.1f' % DL_Mpc + ' Mpc or ' + '%1.1f' % DL_Gyr + ' Gly.')
            print('The distance modulus, m-M, is ' + '%1.2f' % (5 * math.log10(DL_Mpc * 1e6) - 5))
        """
        else:
            print('%1.2f' % zage_Gyr)
            print('%1.2f' % DCMR_Mpc)
            print('%1.2f' % kpc_DA)
            print('%1.2f' % (5 * math.log10(DL_Mpc * 1e6) - 5))
        """
        DL_cm=DL_Mpc*(3.08568*1024)
        #1 Mpc = 1,000,000 parsecs = 3.08568*1024 cm, or 3,261,566 light years.
        return {"DL_cm":DL_cm,"DL_Mpc":DL_Mpc,"DA_Gyr":DA_Gyr}

"""
if __name__ == "__main__":
    p = Photometry()
    image = "/Users/cjimenez/Documents/PHD/data/tmp/17.689916666666665_-5.027574999999997/panstarss/panstars_band_g.fits"
    image = "/Users/cjimenez/Documents/PHD/data/tmp/37.304500000000004_1.1664916666666683/spitzer/swarp/irac2_37.3045_1.16649_swarp.fits"
    p.sourceDetection(image)
"""
