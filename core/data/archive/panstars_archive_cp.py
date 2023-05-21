import requests
import re
from frastro import ContentProvider
from astropy.io.votable import parse_single_table
from frastro.core.utils.html_parse import HtmlParser
from bs4.element import Tag
from frastro import ImageSource
from frastro import AstroSource
from frastro import CoordinateParser
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd

from frastro import CatalogSource
from frastro import VOTableUtil
from frastro import Utils
from astropy.table import Table
from frastro import WaveLenghtCover

class PanSTARRSArchiveCP(ContentProvider):
    """
    This code extract the html form result from ps1images.stsci.edu
    example request: http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=20.48371+0.4223&filter=g&filter=z&filter=y&filetypes=stack&auxiliary=data&size=40&output_size=512&verbose=0&autoscale=99.500000&catlist=
    if the html output change the code will't work due this is extract from html page code and not from web/rest services
    """

    __url_img_cutout = "http://ps1images.stsci.edu/cgi-bin/ps1cutouts"
    __server_catalog = 'https://archive.stsci.edu/panstarrs/search.php'
    __catalog_provider = "Panstarss"

    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/panstarss/"

    __wavelenght = {"g": 486 * u.nm, "r": 621 * u.nm, "i": 754 * u.nm, "z": 867 * u.nm, "y": 963 * u.nm}
    __key_band = {"mag_g": "gMeanApMag", "mag_r": "rMeanApMag", "mag_i": "iMeanApMag", "mag_z": "zMeanApMag","mag_y": "yMeanApMag",}

    #def query_region(self,):

    def query(self,**kwargs):
        # check for data release to use
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"]) + "," + str(kwargs["dec"])
        radius = 1 # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]
        #ra = 0, dec = 0, rad = 3, mindet = 1,
        #maxsources = 101,



        """
        Query Pan-STARRS DR1 @ MAST
        parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                              radius in degrees
                    mindet: minimum number of detection (optional)
                    maxsources: maximum number of sources
                    server: servername
        returns: astropy.table object
        """



        #result = self.panstars_request(server,
        #                 params={'RA': ra_deg, 'DEC': dec_deg,
        #                         'SR': rad_deg, 'max_records': maxsources,
        #                         'outputformat': 'JSON',
        #                         'ndetections': ('>%d' % mindet)})

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)

        self.__save_path = self.__save_path.format(
            str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        result = AstroSource(self.__coordinates)
        #cover limit -30 degrees
        if self.__coordinates.dec.degree>=-30:
            table,cat=self.getCatalog(self.__coordinates.ra.degree, self.__coordinates.dec.degree,radius)
            if len(cat)>0:
                #conver json data to pandas for pick up ra dec columns, each column was converted to degree values
                df = pd.DataFrame(cat.json())
                coord_list = [SkyCoord(df["raMean"][index].replace(" ",":")+" "+df["decMean"][index].replace(" ",":"), unit=(u.hour, u.deg), frame='icrs')for index in range(df["raMean"].size)]
                ra_list=[ra.ra.degree for ra in coord_list]
                dec_list = [dec.dec.degree for dec in coord_list]

                #table = Table.from_pandas(df)
                Utils.createPath(self.__save_path)
                catalog = CatalogSource("panstars", "Search in radius " + str(radius))
                VOTableUtil.saveFromTable(table, self.__save_path + "catalog.xml")
                catalog.addFile("panstars", self.__save_path + "catalog.xml", "vo")

                result.addCatalog(catalog)
                panstars = WaveLenghtCover.panstars()

                index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree,
                                                              self.__coordinates.dec.degree, np.array(ra_list),
                                                              np.array(dec_list))
                respond=cat[index]
                result.addSummaryParams("ra", ra_list[index])
                result.addSummaryParams("dec", dec_list[index])
                result.addSummaryParams("distance", loss)


                result.addSummaryParams("mag_g", {"lambda": panstars["g"], "ab": str(respond["gMeanApMag"]),"err":str(respond["gMeanApMagErr"])})
                result.addSummaryParams("mag_r", {"lambda": panstars["r"], "ab": str(respond["rMeanApMag"]),"err":str(respond["rMeanApMagErr"])})
                result.addSummaryParams("mag_i", {"lambda": panstars["i"], "ab": str(respond["iMeanApMag"]),"err": str(respond["iMeanApMagErr"])})
                result.addSummaryParams("mag_z", {"lambda": panstars["z"], "ab": str(respond["zMeanApMag"]),"err": str(respond["zMeanApMagErr"])})
                result.addSummaryParams("mag_y", {"lambda": panstars["y"], "ab": str(respond["yMeanApMag"]),"err": str(respond["yMeanApMagErr"])})

                mags={"g":str(respond["gMeanApMag"])+"+/-"+str(respond["gMeanApMagErr"]),"r":str(respond["rMeanApMag"])+"+/-"+str(respond["rMeanApMagErr"]),"i":str(respond["iMeanApMag"])+"+/-"+str(respond["iMeanApMagErr"]),"z":str(respond["zMeanApMag"])+"+/-"+str(respond["zMeanApMagErr"]),"y":str(respond["yMeanApMag"])+"+/-"+str(respond["yMeanApMagErr"])}

                """
                result.addSummaryParams("mag_g", str(respond["gMeanApMag"]))
                result.addSummaryParams("mag_r", str(respond["rMeanApMag"]))
                result.addSummaryParams("mag_i", str(respond["iMeanApMag"]))
                result.addSummaryParams("mag_z", str(respond["zMeanApMag"]))
                result.addSummaryParams("mag_y", str(respond["yMeanApMag"]))

                result.addSummaryParams("magerr_g", str(respond["gMeanApMagErr"]))
                result.addSummaryParams("magerr_r", str(respond["rMeanApMagErr"]))
                result.addSummaryParams("magerr_i", str(respond["iMeanApMagErr"]))
                result.addSummaryParams("magerr_z", str(respond["zMeanApMagErr"]))
                result.addSummaryParams("magerr_y", str(respond["yMeanApMagErr"]))
                """

            img=self.cutout(self.__coordinates.ra.degree, self.__coordinates.dec.degree,radius,mags)
            if type(img) is ImageSource:
                result.addImage(img)
            else:
                result.addSummaryParams("error","Any results was found in Panstarss Archive")
        return result

    def getBand(self, band):

        key_band = ""
        if band in self.__key_band:
            key_band = self.__key_band[band]
        return key_band

    def getCatalog(self,ra_deg, dec_deg, rad_arcmin,filter=['color','g','r','i','z','y']):

        """

        :param ra_deg:
        :param dec_deg:
        :param rad_arcmin:
        :param filter:
        :return:
        """

        """"
        target: 
        resolver: Resolve
        radius: 3.0
        ra: 8.049666
        dec: -29.872637
        equinox: J2000
        nDetections: > 1
        selectedColumnsCsv: objname,objid,ramean,decmean,rameanerr,decmeanerr,ndetections,randomid,projectionid,skycellid,objinfoflag,qualityflag,rastack,decstack,rastackerr,decstackerr,epochmean,nstackdetections,ng,nr,ni,nz,ny,gqfperfect,gmeanpsfmag,gmeanpsfmagerr,gmeankronmag,gmeankronmagerr,gmeanapmag,gmeanapmagerr,gflags,rqfperfect,rmeanpsfmag,rmeanpsfmagerr,rmeankronmag,rmeankronmagerr,rmeanapmag,rmeanapmagerr,rflags,iqfperfect,imeanpsfmag,imeanpsfmagerr,imeankronmag,imeankronmagerr,imeanapmag,imeanapmagerr,iflags,zqfperfect,zmeanpsfmag,zmeanpsfmagerr,zmeankronmag,zmeankronmagerr,zmeanapmag,zmeanapmagerr,zflags,yqfperfect,ymeanpsfmag,ymeanpsfmagerr,ymeankronmag,ymeankronmagerr,ymeanapmag,ymeanapmagerr,yflags,ang_sep
        availableColumns: objname
        ordercolumn1: ang_sep
        ordercolumn2: objid
        ordercolumn3: 
        coordformat: sex
        outputformat: JSON
        max_records: 50001
        max_rpp: 5000
        action: Search
        """
        params = {
            'target':'',
            'resolver': 'Resolve',
            'radius': rad_arcmin,
            'ra': ra_deg,
            'dec': dec_deg,
            'equinox': 'J2000',
            'nDetections': '> 1',
            'selectedColumnsCsv':'objname,objid,ramean,decmean,rameanerr,decmeanerr,ndetections,randomid,projectionid,skycellid,objinfoflag,qualityflag,rastack,decstack,rastackerr,decstackerr,epochmean,nstackdetections,ng,nr,ni,nz,ny,gqfperfect,gmeanpsfmag,gmeanpsfmagerr,gmeankronmag,gmeankronmagerr,gmeanapmag,gmeanapmagerr,gflags,rqfperfect,rmeanpsfmag,rmeanpsfmagerr,rmeankronmag,rmeankronmagerr,rmeanapmag,rmeanapmagerr,rflags,iqfperfect,imeanpsfmag,imeanpsfmagerr,imeankronmag,imeankronmagerr,imeanapmag,imeanapmagerr,iflags,zqfperfect,zmeanpsfmag,zmeanpsfmagerr,zmeankronmag,zmeankronmagerr,zmeanapmag,zmeanapmagerr,zflags,yqfperfect,ymeanpsfmag,ymeanpsfmagerr,ymeankronmag,ymeankronmagerr,ymeanapmag,ymeanapmagerr,yflags,ang_sep',
            'availableColumns': 'objname',
            'ordercolumn1': 'ang_sep',
            'ordercolumn2': 'objid',
            'ordercolumn3':'',
            'coordformat': 'sex',
            'outputformat': 'JSON',
            'max_records': 50001,
            'max_rpp': 5000,
            'action': 'Search'
        }

        r = self.panstars_request(self.__server_catalog, params=params)
        print("panstars",r)
        table=[]
        json=r.json()
        if r.status_code ==200:
            df = pd.DataFrame(r.json())
            df["gMeanApMag"]=pd.to_numeric(df["gMeanApMag"], errors='coerce')
            df["rMeanApMag"] = pd.to_numeric(df["rMeanApMag"], errors='coerce')
            df["iMeanApMag"] = pd.to_numeric(df["iMeanApMag"], errors='coerce')
            df["zMeanApMag"] = pd.to_numeric(df["zMeanApMag"], errors='coerce')
            df["yMeanApMag"] = pd.to_numeric(df["yMeanApMag"], errors='coerce')
            table = Table.from_pandas(df)

        return table,r



    def cutout(self,ra_deg, dec_deg, rad_arcmin,mags,filter=['color','g','r','i','z','y'],img_size=256):
        """
        Parse the hrml respond, find html table output using BeautifulSoup
        the first row in table have links to download fits and display image in interative image display(http://ps1images.stsci.edu/cgi-bin/display)
        the second row have cutout image, the link have a no clear pather to replecate with out request
        this is the reason for parse all html output and take the content from scratch

        :param ra_deg: (float) Right ascension in degree
        :param dec_deg: (float) Declination in degree
        :param rad_arcsec:(int) radio of search from the source
        :param mags:(dic) dictionary with magnitudes
        :param filter: (list) list of filters to search ['color','g','r','i','z','y']
        :param img_size: pizel size of image output are avalible '256,512,1024' pixels
        :return: ImageSource of objects with {name: name of image, link:list of links fits-files complete fits file
        and cutout fits file, some case are link to display navigation image, img: is the link to JPG cutout}
        """


        pos = str(ra_deg)+","+str(dec_deg)
        size = (rad_arcmin*60)*15/3.75

        params={
            'pos':pos,
            'filter':filter,
            'filetypes':'stack',
            'auxiliary':'data',
            'size':int(size),
            'output_size':img_size,#'256,512,1024'
            'verbose':'0',
            'autoscale':'99.500000'}

        r = self.panstars_request(self.__url_img_cutout,params=params)



        try:
            error = HtmlParser.findAllElements(r.text, "h2")[0].contents[0]
            sr = re.search("Unable to determine",error.string).group(0)
            success = 300
        except AttributeError:
            success = 200



        tables = HtmlParser.findAllElements(r.text,"tr")


        image = ImageSource(str(ra_deg)+"_"+str(dec_deg), self.__catalog_provider)

        for index in range(len(tables[0].contents)):
            imgDescript = tables[0].contents[index]
            if type(imgDescript) == Tag:
                img_source = tables[1].contents[index].contents[0]
                links=[]

                #image.addCutout("url":url,"size":size,"band":band,"magnitude")

                #image.addFile({"url": url, "name": name, "type": type})
                display_link = ""
                for tags in imgDescript.contents:
                    if type(tags) == Tag and tags.name == "a" :
                        if "Display" == tags.contents[0]:
                            display_link = tags.attrs["href"]
                            if "http:" not in display_link:
                                display_link = "http:"+display_link
                        elif "FITS" == tags.contents[0]:
                            url_file = tags.attrs["href"]
                            if "http://ps1images.stsci.edu/" not in url_file:
                                url_file="http://ps1images.stsci.edu/"+url_file

                            name_file = tags.attrs["title"]
                            type_file = "FITS"
                        #links.append({"title":tags.attrs["title"],"href":tags.attrs["href"]})

                            band = re.search('(y/i/g|g|r|i|z|y)',imgDescript.contents[0]).group(0).replace(" ","")
                            name = imgDescript.contents[0]
                #img={'name':re.search('(y/i/g|g|r|i|z|y)',imgDescript.contents[0]).group(0),
                #     'link':links,
                #     'img':img_source.attrs["src"]}
                            url_img =  "http:"+img_source.attrs["src"] if "http:" not in img_source.attrs["src"] else img_source.attrs["src"]

                            image.addCutout(url_img, img_size, band,mags[band],display_link)
                            local_path = self.__save_path + "panstars_band_" + band + ".fits"
                            image.addFile(name="band-" + band, url=url_file, download=True,
                                          local_path=local_path, uncompress=False,thumbnail=False,external=True)
                            #image.addFile(name_file, url_file, type_file)

        return image

    def panstars_request(self,url,params):

        try:
            r = requests.get(url,params=params)
        except ConnectionError:
            r=[]
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

