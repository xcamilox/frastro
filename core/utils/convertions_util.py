import math
from astropy import units as u

class Convertion():
    @staticmethod
    def janskyToAb(fluxJansky):
        AB = 2.5 * (23 - math.log10(fluxJansky)) - 48.6
        return AB

    @staticmethod
    def abToJansky(magAB):
        Flux_jansky = 10(23 - (magAB + 48.6) / 2.5)
        return Flux_jansky

    @staticmethod
    def degreeToArcsec(degrees):
        d=degrees*u.degrees
        return d.to("arcsec")