import math
from astropy import units as u
from astropy.coordinates import Distance
import numpy as np
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
        d=degrees*u.deg
        return d.to("arcsec")

    @staticmethod
    def degreeToHMS(degrees):
        d = degrees * u.deg
        return d.to("hmsdms")

    @staticmethod
    def aparentToAbsoluteMagnitud(abmagnitud, z="",distance_pc=""):

        if z=="" and distance_pc =="":
            raise Exception('Distance in pc or redshift is necessary')

        if z!="":

            distance =  Distance(z=z,allow_negative=True).to(u.parsec).value
        if distance_pc != "":
            distance = distance_pc


        M = abmagnitud - 5 * np.log10(distance) + 5

        return M