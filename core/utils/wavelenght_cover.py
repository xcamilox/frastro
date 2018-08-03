import astropy.units as u
import numpy as np
import pandas as pd

class WaveLenghtCover():
    @staticmethod
    def sdss():
        return {"u": (356.18 * u.nm).value, "g": (471.89 * u.nm).value, "r": (618.52 * u.nm).value, "i": (749.97 * u.nm).value, "z": (896.15 * u.nm).value, "cover_u":	(98.0* u.nm).value, "cover_g":	(176.6* u.nm).value, "cover_r":	(157.4* u.nm).value, "cover_i":	(170.0* u.nm).value, "cover_z":	(287.3* u.nm).value}

    @staticmethod
    def decals():
        return {"g": (475 * u.nm).value, "r": (635 * u.nm).value, "z": (925 * u.nm).value,"cover_g":(162.5*u.nm).value, "cover_r":(173.2* u.nm).value, "cover_z":(181.8* u.nm).value}

    @staticmethod
    def cfht():
        return {"u": (335 * u.nm).value, "g": (475 * u.nm).value, "r": (640 * u.nm).value, "i": (776 * u.nm).value, "z": (925 * u.nm).value,"cover_u":	(190.3* u.nm).value, "cover_g":	(182.2* u.nm).value, "cover_r":	(163.0* u.nm).value, "cover_i":	(189.8* u.nm).value, "cover_z":	(221.7* u.nm).value}

    @staticmethod
    def des():
        return  {"g": (475 * u.nm).value, "r": (635 * u.nm).value, "i": (775 * u.nm).value, "z": (925 * u.nm).value, "y": (1000 * u.nm).value,"cover_g":(162.5*u.nm).value, "cover_r":(173.2* u.nm).value, "cover_i":(169.6* u.nm).value, "cover_z":(181.8* u.nm).value, "cover_y":(137.5* u.nm).value}

    @staticmethod
    def panstars():
        return {"g": (486 * u.nm).value, "r": (621 * u.nm).value, "i": (754 * u.nm).value, "z": (867 * u.nm).value, "y": (963 * u.nm).value,"cover_g": (165.0*u.nm).value,"cover_r": (165.0*u.nm).value,"cover_w": (437.8*u.nm).value,"cover_i": (152.6*u.nm).value,"cover_z": (131.8*u.nm).value,"cover_y": (173.8*u.nm).value}

    @staticmethod
    def spitzer():
        return {"irac1": (3557.26 * u.nm).value, "irac2": (4504.93 * u.nm).value, "irac3": (5738.57 * u.nm).value,
                        "irac4": (7927.37 * u.nm).value, "mips24": (23843.31 * u.nm).value, "mips70": (72555.53 * u.nm).value,
                        "mips160": (156962.71 * u.nm).value,"cover_irac1":	(831.8*u.nm).value,"cover_irac2":	(1138.8*u.nm).value,"cover_irac3":	(1610.6*u.nm).value,"cover_irac4":	(3288.2*u.nm).value,"cover_mips24":	(11049.4*u.nm).value,"cover_mips70":(59052.7*u.nm).value,"cover_mips160":(78452.7*u.nm).value}

    @staticmethod
    def wise():
        #return {"w1": (3352.60 * u.nm).value, "w2": (4602.80 * u.nm).value, "w3": (11560.80 * u.nm).value, "w4": (22088.30 * u.nm).value,"cover_w1": (1118.3*u.nm).value,"cover_w2": (1378.1*u.nm).value,"cover_w3": (9818.3*u.nm).value,"cover_w4": (8390.6*u.nm).value,"j": (1251.02 * u.nm).value, "h": (1637.71 * u.nm).value, "k": (2208.28 * u.nm).value,"cover_j":(204.0*u.nm).value,"cover_h":(371.0*u.nm).value,"cover_k":(474.5*u.nm).value}

        return {"w1": (3352.60 * u.nm).value, "w2": (4602.80 * u.nm).value, "w3": (11560.80 * u.nm).value,
            "w4": (22088.30 * u.nm).value, "cover_w1": (1118.3 * u.nm).value, "cover_w2": (1378.1 * u.nm).value,
            "cover_w3": (9818.3 * u.nm).value, "cover_w4": (8390.6 * u.nm).value}

    @staticmethod
    def ukidss():
        return {"z": (883.07 * u.nm).value, "y": (1031.87 * u.nm).value, "j": (1251.02 * u.nm).value, "h": (1637.71 * u.nm).value, "k": (2208.28 * u.nm).value,"cover_z":(122.5*u.nm).value,"cover_y":(139.0*u.nm).value,"cover_j":(204.0*u.nm).value,"cover_h":(371.0*u.nm).value,"cover_k":(474.5*u.nm).value}

    @staticmethod
    def vhs():
        return {"z": (878.05 * u.nm).value,"nb980": (978.13 * u.nm).value,"nb990": (990.87 * u.nm).value, "y": (1021.12 * u.nm).value,"nb118": (1190.76 * u.nm).value, "j": (1254.09 * u.nm).value,
                "h": (1646.37 * u.nm).value, "ks": (2148.77 * u.nm).value,"cover_z": (124.3*u.nm).value,"cover_nb980": (27.4*u.nm).value,"cover_nb990": (28.6*u.nm).value,"cover_y": (155.0*u.nm).value,"cover_nb118": (35.4*u.nm).value,"cover_j": (233.2*u.nm).value,"cover_h": (381.8*u.nm).value,"cover_ks": (434.1*u.nm).value}


    @staticmethod
    def waveLenghtCoverByFilter(lambda_nm, redshift=0):

        sdss = WaveLenghtCover.sdss()
        decals = WaveLenghtCover.decals()
        cfht = WaveLenghtCover.cfht()
        des = WaveLenghtCover.des()
        panstars = WaveLenghtCover.panstars()
        spitzer = WaveLenghtCover.spitzer()
        wise = WaveLenghtCover.wise()
        ukidss = WaveLenghtCover.ukidss()
        vhs = WaveLenghtCover.vhs()

        if type(lambda_nm) is not list and type(lambda_nm) is not np.array:
            lambda_nm = np.array([lambda_nm])
        elif type(lambda_nm) is list:
            lambda_nm = np.array(lambda_nm)

        if type(redshift) is list or type(redshift) is np.array:
            lambda_nm_tmp = {}
            for z in redshift:
                lambda_nm_tmp[str(z)]=lambda_nm * (1 + z)
            lambda_nm = lambda_nm_tmp
        else:

            lambda_nm = {str(redshift):lambda_nm * (1+redshift)}
            redshift=[redshift]

        filters = [sdss,decals,cfht,des,panstars,spitzer,wise,ukidss,vhs]
        archive = ["SDSS", "decals", "CFHT", "DES", "PanSTARRS", "spitzer", "wise_unwise", "Ukidss", "vhs"]
        utils_filters=[]
        for z in redshift:
            lambda_nm_val = lambda_nm[str(z)]
            for idx, filter in enumerate(filters):
                for filter_range in filter:
                    if "cover" not in filter_range:
                        min_wev = filter[filter_range] - (filter["cover_" + filter_range] / 2)
                        max_wev = filter[filter_range] + (filter["cover_" + filter_range] / 2)
                        r = np.where(np.logical_and(lambda_nm_val >= min_wev, lambda_nm_val <= max_wev))[0]
                        for result in r:
                            utils_filters.append([lambda_nm_val[result],str(z),archive[idx],filter_range,(filter[filter_range]* u.nm),(filter["cover_"+filter_range]* u.nm)])



        df=pd.DataFrame(data=utils_filters,columns=["WaveLenght","redshift","archive","band","lambda_mean","lambda_cover"])

        return df

if __name__ == "__main__":

    #get filter thant cover the rest wavelenght, for diferents redshifts
    redshift=[0.0,1.0,2.7,3.0,4.0]
    wavelenghts_nm=[430.0,460.0]
    filters=WaveLenghtCover.waveLenghtCoverByFilter(wavelenghts_nm,redshift)
    #example to conver units like nm to angstrom
    print(filters["archive"],filters["band"],filters["lambda_mean"].apply(lambda x:x.to(u.AA)))