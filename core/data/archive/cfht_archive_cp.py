from frastro import ContentProvider
from frastro import TAPManager
from frastro import AstroSource
from frastro import CoordinateParser
import astropy.units as u
from frastro import ImageSource
from frastro import CatalogSource
import numpy as np
from frastro.core.utils.html_parse import HtmlParser
import requests
import re
from bs4.element import Tag
from frastro import ImageUtils
from frastro import VOTableUtil
from frastro import Utils
from frastro import WaveLenghtCover

class CFHTAchiveCP(ContentProvider):
    __catalog_provider = "CFHT"

    __service_provider = {
        "tap_url": "http://tapvizier.u-strasbg.fr/TAPVizieR/tap",
        "width_cutout":"http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/megapipe/imc.pl" #"?lang=en&object=&ra={0}&dec={1}&size=256"
    }

    __simple_rec_query = 'SELECT * FROM "II/317/cfhtls_{0}" as w WHERE '+"1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {1},{2}, {3} ))"
    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/cfht/"
    #all unit in micras
    __wavelenght={"u":335*u.nm,"g":475*u.nm,"r":640*u.nm,"i":776*u.nm,"z":925*u.nm}


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
        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        width_survey = "w"
        deep_survey = "d"

        radius_degree = CoordinateParser.getMinToDegree(radius)

        query = self.__simple_rec_query.format(deep_survey, self.__coordinates.ra.degree, self.__coordinates.dec.degree, radius_degree)
        respond = self.getTapRequest(query=query)
        self.wide=False

        if len(respond) <= 0:
            query = self.__simple_rec_query.format(width_survey, self.__coordinates.ra.degree,
                                                   self.__coordinates.dec.degree, radius_degree)
            respond = self.getTapRequest(query=query)
            self.wide = True

        result = AstroSource(self.__coordinates)
        if len(respond)>0:
            Utils.createPath(self.__save_path)

            catalog = CatalogSource("cfht", "Search in radius " + str(radius))
            VOTableUtil.saveFromTable(respond, self.__save_path + "catalog.xml")
            catalog.addFile("allwise", self.__save_path + "catalog.xml", "vo")

            index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                          respond["RAJ2000"], respond["DEJ2000"])

            result.addSummaryParams("ra", respond["RAJ2000"][index])
            result.addSummaryParams("dec", respond["DEJ2000"][index])

            result.addSummaryParams("distance", loss)
            result.addSummaryParams("survey","wide" if self.wide else "deep")

            cfht = WaveLenghtCover.cfht()
            data = {"lambda": cfht["u"], "ab": str(respond["umag"][index]), "err": str(respond["e_umag"][index])}

            result.addSummaryParams("mag_u", data)

            data = {"lambda": cfht["g"], "ab": str(respond["gmag"][index]), "err": str(respond["e_gmag"][index])}
            result.addSummaryParams("mag_g", data)

            data = {"lambda": cfht["r"], "ab": str(respond["rmag"][index]), "err": str(respond["e_rmag"][index])}
            result.addSummaryParams("mag_r", data)

            data = {"lambda": cfht["i"], "ab": str(respond["imag"][index]), "err": str(respond["e_imag"][index])}
            result.addSummaryParams("mag_i", data)

            data = {"lambda": cfht["z"], "ab": str(respond["zmag"][index]), "err": str(respond["e_zmag"][index])}
            result.addSummaryParams("mag_z", data)

            result.addSummaryParams("flux_u", str(respond["uFpk"][index]))
            result.addSummaryParams("flux_g", str(respond["gFpk"][index]))
            result.addSummaryParams("flux_r", str(respond["rFpk"][index]))
            result.addSummaryParams("flux_i", str(respond["iFpk"][index]))
            result.addSummaryParams("flux_z", str(respond["zFpk"][index]))

            result.addCatalog(catalog)

            if self.wide:

                #img=self.getWideImages()
                url=self.__service_provider["width_cutout"]+"?lang=en&object=&ra={0}&dec={1}&size=256"
                url=url.format(self.__coordinates.ra.degree, self.__coordinates.dec.degree)
                images = self.delayRequest(url)
                image = ImageSource(result.getSummary()["id"], self.__catalog_provider)
                mag={"U":str(respond["umag"][index])+"+/-"+str(respond["e_umag"][index]),"G":str(respond["gmag"][index])+"+/-"+str(respond["e_gmag"][index]),"R":str(respond["rmag"][index])+"+/-"+str(respond["e_rmag"][index]),"I":str(respond["imag"][index])+"+/-"+str(respond["e_imag"][index]),"Z":str(respond["zmag"][index])+"+/-"+str(respond["e_zmag"][index]),"2":-99}
                for idx in range(len(images)):
                    img=images[idx]
                    band = img["band"]
                    mag_val=mag[band]
                    image.addCutout(img["image"], 256, band,mag_val)

                    local_path = self.__save_path + "cfht_band_" + band + ".fits"
                    image.addFile(band, img["fits"], "fits", download=True, local_path=local_path,
                                  uncompress=False,thumbnail=True,external=False)

                result.addImage(image)


        return result


    def getWideImages(self):

        #   parameters "?lang=en&object=&ra={0}&dec={1}&size=256"

        img_size=256
        params = {
            'lang': 'en',
            'object': '',
            'ra': self.__coordinates.ra.degree,
            'dec': self.__coordinates.dec.degree,
            'size': img_size}

        r = requests.get(self.__service_provider["width_cutout"], params=params)


        try:
            error = HtmlParser.getById(r.text, "url0")[0].contents[0]
            sr = re.search("Unable to determine", error.string).group(0)
            success = 300
        except AttributeError:
            success = 200

        tables = HtmlParser.findAllElements(r.text, "tr")
        pos = str(params["ra"]) + "," + str(params["deg"])
        image = ImageSource(pos, self.__catalog_provider)


        for index in range(len(tables[0].contents)):
            imgDescript = tables[0].contents[index]
            if type(imgDescript) == Tag:
                img_source = tables[1].contents[index].contents[0]
                links = []

                # image.addCutout("url":url,"size":size,"band":band,"magnitude")

                # image.addFile({"url": url, "name": name, "type": type})
                display_link = ""
                for tags in imgDescript.contents:
                    if type(tags) == Tag and tags.name == "a":
                        if "Display" == tags.contents[0]:
                            display_link = tags.attrs["href"]
                            if "http:" not in display_link:
                                display_link = "http:" + display_link
                        elif "FITS" == tags.contents[0]:
                            url_file = tags.attrs["href"]
                            if "http://ps1images.stsci.edu/" not in url_file:
                                url_file = "http://ps1images.stsci.edu/" + url_file

                            name_file = tags.attrs["title"]
                            type_file = "FITS"
                            # links.append({"title":tags.attrs["title"],"href":tags.attrs["href"]})

                            band = re.search('(y/i/g|g|r|i|z|y)', imgDescript.contents[0]).group(0)
                            name = imgDescript.contents[0]
                            # img={'name':re.search('(y/i/g|g|r|i|z|y)',imgDescript.contents[0]).group(0),
                            #     'link':links,
                            #     'img':img_source.attrs["src"]}
                            url_img = "http:" + img_source.attrs["src"] if "http:" not in img_source.attrs["src"] else \
                            img_source.attrs["src"]
                            local_path=self.__save_path+"cfht_band_"+band+".fits"
                            image.addCutout(url_img, img_size, band, -99, display_link)
                            image.addFile(name_file, url_file, type_file,download=True,local_path=local_path, uncompress=False)

        return image


    def delayRequest(self,link):
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        browser = webdriver.Chrome('/Applications/chromedriver')
        browser.set_window_size(1120, 550)
        browser.get(link)
        element = WebDriverWait(browser, 1).until(
            EC.presence_of_element_located((By.ID, "url0"))
        )

        data = element.find_elements_by_xpath("//tr/td/div/a")
        list_img=[]
        for item in data:
            url=item.get_attribute("href")

            position=url.find('.fits')
            band=url[position-1:position]
            list_img.append({"fits":url,"band":band,'image':ImageUtils.getBase64FromFitsURL(url)})

        browser.quit()
        return list_img

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
    def onCreate(self):
        pass
    def update(self):
        pass
