from frastro import ContentProvider
from frastro import TAPManager
from frastro import AstroSource
from frastro import CoordinateParser
import astropy.units as u
from frastro import ImageSource
import numpy as np
from frastro import NOAOArchiveCP
from frastro import CatalogSource
from frastro import VOTableUtil
from frastro import Utils
from astropy.table import Table
from frastro import WaveLenghtCover
import urllib.parse

class VHSAchiveCP(ContentProvider):
    __catalog_provider = "vhs"

    __service_provider = {
        "vhs": {"url":"http://wfaudata.roe.ac.uk/vhsDR3-dsa/TAP","table":"vhsSource"},
        "viking":{"url":"http://wfaudata.roe.ac.uk/vikingDR3-dsa/TAP","table":"vikingSource"},
        "video":{"url":"http://wfaudata.roe.ac.uk/videoDR3-dsa/TAP","table":"videoSource"},
        "image_url":"http://horus.roe.ac.uk:8080/vdfs/GetImage?"
    }

    __simple_rec_query = "SELECT * FROM {0} WHERE((ra between {1} and {2}) and (dec between {3} and {4}))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/vhs/"

    __wavelenght = {"g": 475 * u.nm, "r": 635 * u.nm, "i": 775 * u.nm, "z": 925 * u.nm, "y": 1000 * u.nm}

    def __init__(self):
        pass

    def query(self, **kwargs):
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"]) + "," + str(kwargs["dec"])
        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        self.__save_path=self.__save_path.format(str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))



        degrees_range=CoordinateParser.getMinToDegree(radius)

        ra_min = self.__coordinates.ra.degree - (degrees_range/2)
        ra_max = self.__coordinates.ra.degree + (degrees_range/2)
        dec_min = self.__coordinates.dec.degree - (degrees_range/2)
        dec_max = self.__coordinates.dec.degree + (degrees_range/2)


        table = self.__service_provider["video"]["table"]
        query = self.__simple_rec_query.format(table, ra_min, ra_max, dec_min, dec_max)
        url = self.__service_provider["video"]["url"]
        respond=self.getTapRequest(url=url,query=query)
        program="video"
        if len(respond) <= 0:
            table = self.__service_provider["viking"]["table"]
            query = self.__simple_rec_query.format(table, ra_min, ra_max, dec_min, dec_max)
            url = self.__service_provider["viking"]["url"]
            respond = self.getTapRequest(url=url, query=query)
            program = "viking"

        if len(respond) <= 0:
            table = self.__service_provider["vhs"]["table"]
            query = self.__simple_rec_query.format(table, ra_min, ra_max, dec_min, dec_max)
            url = self.__service_provider["vhs"]["url"]
            respond = self.getTapRequest(url=url, query=query)
            program = "vhs"


        result = AstroSource(self.__coordinates)
        if len(respond) > 0:
            vhs = WaveLenghtCover.vhs()
            Utils.createPath(self.__save_path)
            catalog = CatalogSource("vhs", "Search in radius " + str(radius))
            VOTableUtil.saveFromTable(respond, self.__save_path + "catalog.xml")
            catalog.addFile("vhs", self.__save_path + "catalog.xml", "vo")

            result.addCatalog(catalog)

            index, loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree,
                                                                self.__coordinates.dec.degree,
                                                                respond["ra"], respond["dec"])

            result.addSummaryParams("ra", respond["ra"][index])
            result.addSummaryParams("dec", respond["dec"][index])
            result.addSummaryParams("distance", loss)
            mag = {"Y": str(respond["yAperMag3"][index])+"+/-"+str(respond["yAperMag3Err"][index]), "J": str(respond["jAperMag3"][index])+"+/-"+str(respond["jAperMag3Err"][index]),
                   "H": str(respond["hAperMag3"][index])+"+/-"+str(respond["hAperMag3Err"][index]), "Ks": str(respond["ksAperMag3"][index])+"+/-"+str(respond["ksAperMag3Err"][index])}
            if program=="video" or program=="viking":
                result.addSummaryParams("magpetro_z", {"lambda": vhs["z"], "ab": str(respond["zPetroMag"][index]),"err":str(respond["zPetroMagErr"][index])})
                result.addSummaryParams("magAper_z", {"lambda": vhs["z"], "ab": str(respond["zAperMag3"][index]),"err":str(respond["zAperMag3Err"][index])})
                mag["Z"] = str(respond["zAperMag3"][index])+"+/-"+str(respond["zAperMag3Err"][index])

            result.addSummaryParams("magpetro_y", {"lambda": vhs["y"], "ab": str(respond["yPetroMag"][index]),"err":str(respond["yPetroMagErr"][index])})
            result.addSummaryParams("magpetro_j", {"lambda": vhs["j"], "ab": str(respond["jPetroMag"][index]),"err":str(respond["jPetroMagErr"][index])})
            result.addSummaryParams("magpetro_h", {"lambda": vhs["h"], "ab": str(respond["hPetroMag"][index]),"err":str(respond["hPetroMagErr"][index])})
            result.addSummaryParams("magpetro_ks", {"lambda": vhs["ks"], "ab": str(respond["ksPetroMag"][index]),"err":str(respond["ksPetroMagErr"][index])})


            result.addSummaryParams("magAper_y", {"lambda": vhs["y"], "ab": str(respond["yAperMag3"][index]),"err":str(respond["yAperMag3Err"][index])})
            result.addSummaryParams("magAper_j", {"lambda": vhs["j"], "ab": str(respond["jAperMag3"][index]),"err":str(respond["jAperMag3Err"][index])})
            result.addSummaryParams("magAper_h", {"lambda": vhs["h"], "ab": str(respond["hAperMag3"][index]),"err":str(respond["hAperMag3Err"][index])})
            result.addSummaryParams("magAper_ks", {"lambda": vhs["ks"], "ab": str(respond["ksAperMag3"][index]),"err":str(respond["ksAperMag3Err"][index])})

            result.addSummaryParams("type", str(respond["mergedClass"][index]))

            images= self.getImages(program)

            image = ImageSource(result.getSummary()["id"], self.__catalog_provider)



            saved_path={}

            for idx in range(len(images)):
                img = images[idx]
                band = img["band"]
                mag_val = mag[band]
                image.addCutout(img["image"], 256, band, mag_val)
                index_image=""
                if band in saved_path:
                    saved_path[band]+=1
                    index_image=saved_path[band]
                else:
                    saved_path[band] = 0

                local_path = self.__save_path + "vhs_band_" + band +str(index_image)+ ".fits"
                image.addFile(band, img["fits"], "fits", download=True, local_path=local_path,
                              uncompress=True, thumbnail=False, external=False)

            result.addImage(image)

        return result

    def getTapRequest(self,url, query=""):

        tap = TAPManager()
        tap_url = url
        tap.connect(url=tap_url)

        result = tap.sync_query(query=query)
        r = result.get_results()
        print("results",len(r))
        return r

    def getImages(self,program):

        if program=="vhs":
            programmeID = 110
            database = 'VHSDR5'
        elif program=="viking":
            programmeID = 140
            database = 'VIKINGDR4'
        else:

            programmeID = 150
            database = 'VIDEODR5'


        params={
            'archive' : 'VSA',
            'programmeID' : programmeID,
            'database' : database,
            'ra' : self.__coordinates.ra.degree,
            'dec' : self.__coordinates.dec.degree,
            'sys' : 'J',
            'filterID' : 'all',
            'xsize' : 15,
            'ysize' : 15,
            'obsType' : 'object',
            'frameType' : 'tilestack',
            'mfid' : '',
            'fsid' : ''
        }
        url = self.__service_provider["image_url"]

        url_request=url + urllib.parse.urlencode(params)

        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        browser = webdriver.Chrome('/Applications/chromedriver')
        browser.set_window_size(1120, 550)
        browser.get(url_request)
        filter = (By.XPATH, "//a")
        element = WebDriverWait(browser, 1).until(
            EC.presence_of_element_located(filter)
        )

        data = element.find_elements_by_xpath("//a[@href]")
        list_img=[]
        for item in data:
            uri = item.get_attribute("href")
            uri_parse=urllib.parse.urlparse(uri)
            data=urllib.parse.parse_qs(uri_parse.query)
            query=uri_parse.query.replace("/","%2F").replace("&",";")
            jpeg_uri="http://horus.roe.ac.uk/wsa/cgi-bin/getJImage.cgi?"+query
            fits_uri = "http://horus.roe.ac.uk/wsa/cgi-bin/getFImage.cgi?" + query

            band = data["band"][0]
            list_img.append({"fits": fits_uri, "band": band, 'image':jpeg_uri })

        browser.quit()
        return list_img

    def onCreate(self):
        pass

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