from astropy import coordinates as coords
from astropy import units as u
import numpy as np

class CoordinateParser():
    @staticmethod
    def splitCoordinates(coordinates_str):
        radec = coordinates_str.split(",")
        ra = dec = 0
        if len(radec) > 1:
            ra = radec[0]
            dec = radec[1]
        elif len(coordinates_str.split(" ")) > 1:
            radec = coordinates_str.split(" ")
            ra = radec[0]
            dec = radec[1]
        return ra, dec

    @staticmethod
    def validateCoordinates(coordinates_str):
        if coordinates_str is type(coords):
            return coordinates_str

        ra, dec = CoordinateParser.splitCoordinates(coordinates_str)
        if "s" in coordinates_str or "h" in coordinates_str:  # validate with hou
            skycoor = coords.SkyCoord(ra=ra, dec=dec)
        else:
            unit = "deg"
            skycoor = coords.SkyCoord(ra=ra, dec=dec, unit=unit)
        return skycoor

    @staticmethod
    def getMinToDegree(arcmin):
        min=arcmin*u.arcmin
        return min.to(u.degree).value

    @staticmethod
    def getSecToDegree(arcsec):
        sec=arcsec*u.arcsec
        return sec.to(u.degree).value

    @staticmethod
    def getDegreeToMin(degree):
        deg=degree*u.degree
        return  deg.to(u.arcmin).value

    @staticmethod
    def getDegreeToSec(degree):
        deg = degree * u.degree
        return deg.to(u.arcsec).value

    @staticmethod
    def getNearPositionIndex(ra,dec,ra_list,dec_list):
        loss = np.sqrt((ra - ra_list) ** 2 + (dec - dec_list) ** 2)
        index = np.where(loss.min() == loss)[0][0]
        return index,loss.min()

