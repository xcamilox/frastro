from astropy.table import Column

from frastro import ContentProvider
from astroquery.ukidss import Ukidss
import astropy.units as u
from frastro import CoordinateParser
from frastro import AstroSource
from frastro import ImageSource
from frastro import SpectraSource
from frastro import TAPManager
import numpy as np

from frastro import CatalogSource
from frastro import VOTableUtil
from frastro import Utils
from astropy.table import Table
from frastro import WaveLenghtCover

class UkidssArchiveCP(ContentProvider):
    __service_provider = {
        "tap_url": "http://tapvizier.u-strasbg.fr/TAPVizieR/tap"
    }

    __simple_rec_query = 'SELECT * FROM "{0}" as w WHERE ' + "1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {1},{2}, {3} ))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/ukidss/"

    __wavelenght = {"z": 883.07 * u.nm, "y": 1031.87 * u.nm, "j": 1251.02 * u.nm, "h": 1637.71 * u.nm, "k": 2208.28 * u.nm}

    __key_band = {"mag_z": "mag_z", "mag_y": "mag_y", "mag_j": "mag_j", "mag_h": "mag_h","mag_ks": "mag_k" }

    def __init__(self,catalog_provider='Ukidss'):
        self.__catalog_provider=catalog_provider
        self.onCreate()


    def getBand(self, band):

        key_band = ""
        if band in self.__key_band:
            key_band = self.__key_band[band]
        return key_band

    def onCreate(self):

        all_databases = ("UKIDSSDR9PLUS", "UKIDSSDR8PLUS", "UKIDSSDR7PLUS",
                         "UKIDSSDR6PLUS", "UKIDSSDR5PLUS", "UKIDSSDR4PLUS",
                         "UKIDSSDR3PLUS", "UKIDSSDR2PLUS", "UKIDSSDR1PLUS",
                         "UKIDSSDR1", "UKIDSSEDRPLUS", "UKIDSSEDR", "UKIDSSSV",
                         "WFCAMCAL08B", "U09B8v20120403", "U09B8v20100414","UKIDSSDR10PLUS","UHSDR1","UKIDSSDR11PLUSUDSONLY")

        ukidss_programmes_short = {'ALL':'all',
                                   'LAS': 101,
                                   'GPS': 102,
                                   'GCS': 103,
                                   'DXS': 104,
                                   'UDS': 105,
                                   'UHS': 107}

        ukidss_programmes_long = {
                                'ALL Survey':'all',
                                'Large Area Survey': 101,
                                'Galactic Plane Survey': 102,
                                'Galactic Clusters Survey': 103,
                                'Deep Extragalactic Survey': 104,
                                'Ultra Deep Survey': 105,
                                'UKIRT Hemisphere Survey, UHS':107}


        Ukidss.all_databases = all_databases
        Ukidss.programmes_long = ukidss_programmes_long
        Ukidss.programmes_short = ukidss_programmes_short
        Ukidss.ukidss_programmes_long=ukidss_programmes_long
        Ukidss.ukidss_programmes_short=ukidss_programmes_short


    def query(self, **kwargs):

        # check for data release to use
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"]) + "," + str(kwargs["dec"])
        else:
            raise ValueError("Not valid coordinates found. Used coordinates key or ra dec keys")

        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)

        program_id = "ALL"


        respond = self.getCatalog(self.__coordinates.ra.degree, self.__coordinates.dec.degree,radius)
        result = AstroSource(self.__coordinates)

        program_id = self.program_id

        if len(respond)>0:


            Utils.createPath(self.__save_path)
            catalog = CatalogSource("ukidss", "Search in radius " + str(radius))
            VOTableUtil.saveFromTable(respond, self.__save_path + "catalog.xml")
            catalog.addFile("ukidss", self.__save_path + "catalog.xml", "vo")

            result.addCatalog(catalog)

            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          respond["RAJ2000"], respond["DEJ2000"])

            ukidss = WaveLenghtCover.ukidss()
            result.addSummaryParams("ra", respond["RAJ2000"][index])
            result.addSummaryParams("dec", respond["DEJ2000"][index])
            result.addSummaryParams("distance", loss)

            magnitudes={"Y":-99,"J":-99,"H":-99,"K":-99,"Z":-99}

            try:
                data={"lambda": ukidss["h"], "ab": str(respond["mag_h"][index]), "err": str(respond["e_Hmag"][index])}
                result.addSummaryParams("mag_h", data)
                magnitudes["H"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_h found")

            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["mag_j"][index]), "err": str(respond["e_Jmag"][index])}
                result.addSummaryParams("mag_j", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_j found")

            """
            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["Jmag1"][index]), "err": str(respond["e_Jmag1"][index])}
                result.addSummaryParams("mag_j1", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_j1 found")

            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["Jmag2"][index]),"err": str(respond["e_Jmag2"][index])}
                result.addSummaryParams("mag_j2", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_j2 found")
            
            try:
                data = {"lambda": ukidss["j"], "ab": str(respond["pJmag"][index]),"err": str(respond["e_pJmag"][index])}
                result.addSummaryParams("pmag_j", data)
                magnitudes["J"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No pmag_j found")
            
            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["pKmag"][index]),"err": str(respond["e_pKmag"][index])}
                result.addSummaryParams("pmag_k", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No pmag_k found")
            """
            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["mag_k"][index]),"err": str(respond["e_Kmag"][index])}
                result.addSummaryParams("mag_k", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_k found")
            """
            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["Kmag1"][index]), "err": str(respond["e_Kmag1"][index])}
                result.addSummaryParams("mag_k1", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_k1 found")
            try:
                data = {"lambda": ukidss["k"], "ab": str(respond["Kmag2"][index]),"err": str(respond["e_Kmag2"][index])}
                result.addSummaryParams("mag_k2", data)
                magnitudes["K"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_k2 found")
            """
            try:
                data = {"lambda": ukidss["y"], "ab": str(respond["mag_y"][index]),"err": str(respond["e_Ymag"][index])}
                result.addSummaryParams("mag_y", data)
                magnitudes["Y"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_y found")
            try:
                data = {"lambda": ukidss["z"], "ab": str(respond["mag_z"][index]), "err": str(respond["e_Zmag"][index])}
                result.addSummaryParams("mag_z", data)
                magnitudes["Z"] = data["ab"] + "+/-" + data["err"]
            except KeyError:
                print("No mag_z found")

            images = Ukidss.get_image_list(self.__coordinates,waveband='all', frame_type='stack',
                                            image_width=radius * u.arcmin,
                                            database="UKIDSSDR10PLUS",programme_id=program_id)

            display_urls = [link.replace("getFImage", "getImage")for link in images]
            image_urls = [link.replace("getFImage", "getJImage") for link in images]

            image = ImageSource(result.getSummary()["id"], self.__catalog_provider)


            for index in range(len(images)):
                img= image_urls[index]
                band = str(img[img.find("band")+5:img.find("band")+6]).replace(" ","")
                mag=magnitudes[band]
                image.addCutout(image_urls[index],300, band, magnitude=mag,link=display_urls[index])

                local_path = self.__save_path + "ukidss_band_" + band + ".fits"
                image.addFile(name="band-" + band, url=images[index], download=True,
                              local_path=local_path, uncompress=False, thumbnail=True,external=False)

                #image.addFile("band-" + band,images[index] , "fits")

            result.addImage(image)

        return result

    def getCatalog(self,ra,dec,radius):
        self.__coordinates = str(ra) + "," + str(dec)
        radius = radius  # arcmin
        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        large_suvey = "II/319/las9"
        galactic_cluster_suvey = "II/319/gcs9"
        deep_suvey = "II/319/dxs9"
        galatic_plane_survey = "II/316/gps6"

        radius_degree = CoordinateParser.getMinToDegree(radius)


        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        query = self.__simple_rec_query.format(large_suvey, self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                               radius_degree)
        respond = self.getTapRequest(query=query)


        if len(respond) <= 0:
            query = self.__simple_rec_query.format(galactic_cluster_suvey, self.__coordinates.ra.degree,
                                                   self.__coordinates.dec.degree, radius_degree)
            respond = self.getTapRequest(query=query)
            self.program_id = "GCS"

        if len(respond) <= 0:
            query = self.__simple_rec_query.format(deep_suvey, self.__coordinates.ra.degree,
                                                   self.__coordinates.dec.degree, radius_degree)
            respond = self.getTapRequest(query=query)
            self.program_id = "DXS"

        if len(respond) <= 0:
            query = self.__simple_rec_query.format(galatic_plane_survey, self.__coordinates.ra.degree,
                                                   self.__coordinates.dec.degree, radius_degree)
            respond = self.getTapRequest(query=query)
            self.program_id = "GPS"


        #covert all magnitudes in AB system
        """
        The UKIRT Infrared Deep Sky Survey ZY JHK
        Photometric System: Passbands and Synthetic Colours
        P. C. Hewett / 2006
        
        Table 7. Conversion to AB mag.
        band λeff S flux density AB
              µm  Jy Wm−2µ−1    offset
        u 0.3546 1545 3.66 × 10−8 0.927
        g 0.4670 3991 5.41 × 10−8 -0.103
        r 0.6156 3174 2.50 × 10−8 0.146
        i 0.7471 2593 1.39 × 10−8 0.366
        z 0.8918 2222 8.32 × 10−9 0.533
        Z 0.8817 2232 8.59 × 10−9 0.528
        Y 1.0305 2026 5.71 × 10−9 0.634
        J 1.2483 1530 2.94 × 10−9 0.938
        H 1.6313 1019 1.14 × 10−9 1.379
        K 2.2010 631 3.89 × 10−10 1.900
        
        """

        if len(respond)>0:

            h_ab_offset = 1.379
            y_ab_offset = 0.634
            j_ab_offset = 0.938
            k_ab_offset = 1.900
            hmag = respond.field("Hmag")
            jmag = respond.field("Jmag1")
            ymag = respond.field("Ymag")
            kmag = respond.field("Kmag")

            mag_h = hmag + h_ab_offset
            mag_y = ymag + y_ab_offset
            mag_j = jmag + j_ab_offset
            mag_k = kmag + k_ab_offset

            col_h = Column(name='mag_h', data=mag_h.data,description="magnitud h ab, using ab offset +1.379 P.C.Hewett/2006")
            col_y = Column(name='mag_y', data=mag_y.data,description="magnitud y ab, using ab offset +0.634 P.C.Hewett/2006")
            col_j = Column(name='mag_j', data=mag_j.data,description="magnitud j ab, using ab offset +0.938 P.C.Hewett/2006")
            col_k = Column(name='mag_k', data=mag_k.data,description="magnitud k ab, using ab offset +1.900 P.C.Hewett/2006")
            respond.add_columns([col_h,col_y,col_j,col_k])


        #http://www.ukidss.org/technical/photom/hewett-ukidss.pdf

        return respond

    def getCatalogAstroquery(self,ra,dec,radius=1):
        # check for data release to use
        self.__coordinates = str(ra) + "," + str(dec)
        radius = radius  # arcmin

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        ukidss_server=Ukidss
        ukidss_server.programme_id="UHS"
        ukidss_server.database="UHSDR1"

        result=ukidss_server.query_region(self.__coordinates,radius=radius * u.arcmin)

        if len(result) <=0:

            ukidss_server.programme_id = "ALL"
            ukidss_server.database = "UKIDSSDR10PLUS"
            result = ukidss_server.query_region(self.__coordinates, radius=radius * u.arcmin)



        return result




    def getTapRequest(self, query=""):

        tap = TAPManager()
        tap_url = self.__service_provider["tap_url"]
        tap.connect(url=tap_url)

        result=tap.sync_query(query=query)
        r=result.get_results()
        return r

    def contentUrl(self):
        pass
    def delete(self):
        pass
    def getType(self):
        pass
    def insert(self):
        pass
    def update(self):
        pass


if __name__ == "__main__":
    ukisscp=UkidssArchiveCP()
    file_path_lens_liverpool = "/Users/cjimenez/Documents/PHD/data/liverpool_lens/sources_all_ra_dec.csv"
    file_path_laes="/Users/cjimenez/Documents/PHD/data/liris_enero_2019/laes.csv"


    file = np.loadtxt(file_path_laes, delimiter=',')

    for target in file:
        r = ukisscp.getCatalogAstroquery(target[0],target[1], 10)
        print(target,len(r))



