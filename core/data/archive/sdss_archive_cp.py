from frastro.frastro.core.data.interface.content_provider import ContentProvider
from astroquery.sdss import SDSS
from frastro.frastro.core.utils.coordinates_parser import CoordinateParser
from frastro.frastro.core.data.astrosource.astro_source import AstroSource
from frastro.frastro.core.data.astrosource.image_source import ImageSource
from frastro.frastro.core.data.astrosource.catalog_source import CatalogSource
from frastro.frastro.core.data.astrosource.spectra_source import SpectraSource
import astropy.units as u
from frastro.frastro.core.utils.votable_util import VOTableUtil
from frastro.frastro.core.utils.utils import Utils
from frastro.frastro.core.utils.wavelenght_cover import WaveLenghtCover
import json
from astropy import coordinates
from collections import OrderedDict
import requests
from astropy.coordinates import SkyCoord
from astropy import units as u



import pandas as pd

class SDSSArchiveCP(ContentProvider):

    __service_provider = {
        "cutout": "http://skyserver.sdss.org/dr{0}/SkyServerWS/ImgCutout/getjpeg?ra=ravalue&dec=decvalue&scale=scaleval&width=widthval&height=heightval&opt=optvalue&query=queryvalue",
        "spect": "http://skyserver.sdss.org/dr{0}/en/get/SpecById.ashx?id=specobjid",
        "imageformat":"https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/rerun/run/camcol/frame-g-00run-camcol-field.fits.bz2",
        "image": "https://dr12.sdss.org/sas/dr{0}/boss/photoObj/frames/{1}/{2}/{3}/frame-{4}-{5}-{6}-{7}.{8}",
        "display":"http://skyserver.sdss.org/dr14/en/tools/explore/summary.aspx?ra={0}&dec={1}",
        "tapprovider":"http://skyserver.sdss.org/dr14/en/tools/search/x_results.aspx?searchtool=SQL&TaskName=Skyserver.Search.SQL&syntax=NoSyntax&cmd={0}",
        "sql_server":"http://skyserver.sdss.org/dr{0}/en/tools/search/x_results.aspx"

    }

    __save_path = "/Users/cjimenez/Documents/PHD/data/tmp/{0}/sdss/"


    __wavelenght = {"u":354 *u.nm ,"g": 477 * u.nm, "r": 623 * u.nm, "i": 762 * u.nm, "z": 913 * u.nm}

    __key_band = {"mag_u":"u","mag_g":"g","mag_r":"r","mag_i":"i","mag_z":"z"}

    __coordinates=000

    __data_release=16
    __image_release=12

    def __init__(self,catalog_provider='SDSS'):
        self.__catalog_provider=catalog_provider

        self.onCreate()

    def onCreate(self):
        pass

    def getCatalog(self,ra,dec,radius):
        self.__coordinates = str(ra) + "," + str(dec)
        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)
        respond = SDSS.query_region(self.__coordinates, spectro=False,
                                    photoobj_fields=['ra', 'dec', 'objid', 'run', 'type', 'rerun', 'camcol', 'field',
                                                     'u', 'g',
                                                     'r', 'i', 'z', 'err_u', 'err_g', 'err_r', 'err_i', 'err_z'],
                                    data_release=self.__data_release,
                                    radius=radius * u.arcmin)


        return respond


    def query(self, **kwargs):


        # check for data release to use
        if "coordinates" in kwargs:
            self.__coordinates = kwargs["coordinates"]
        elif "ra" in kwargs and "dec" in kwargs:
            self.__coordinates = str(kwargs["ra"])+","+str(kwargs["dec"])
        else:
            raise ValueError("Not valid coordinates found. Used coordinates key or ra dec keys")

        self.__coordinates = CoordinateParser.validateCoordinates(self.__coordinates)

        self.__save_path=self.__save_path.format(str(self.__coordinates.ra.degree) + "_" + str(self.__coordinates.dec.degree))

        Utils.createPath(self.__save_path)

        radius = 1  # arcmin
        if "radius" in kwargs:
            radius = kwargs["radius"]

        #check for data release to use
        if "release" in kwargs:
            self.__data_release = kwargs["release"]


        respond = SDSS.query_region(self.__coordinates, spectro=True,
                                    photoobj_fields=['ra', 'dec', 'objid', 'run','type', 'rerun', 'camcol', 'field', 'u', 'g',
                                                     'r', 'i', 'z','err_u','err_g','err_r','err_i','err_z'],
                                    specobj_fields=['z','zErr','zWarning', 'plate', 'mjd', 'fiberID', 'specobjid', 'run2d',
                                                    'instrument', 'targetObjID'], data_release=self.__data_release,radius=radius*u.arcmin)


        #sp=SDSS.get_spectra(matches=respond,data_release=self.__data_release)
        #im = SDSS.get_images(matches=respond,data_release=self.__data_release)

        # print(respond)
        # for col in respond.columns:
        #     print(respond[col])





        result = AstroSource(self.__coordinates)

        sdss = WaveLenghtCover.sdss()
        index,loss = CoordinateParser.getNearPositionIndex(self.__coordinates.ra.degree, self.__coordinates.dec.degree,
                                                      respond["ra"], respond["dec"])
        if respond != None:
            Utils.createPath(self.__save_path)
            VOTableUtil.saveFromTable(respond, self.__save_path + "catalog.xml")
            catalog = CatalogSource("sdss", "Search in radius(arcmin) " + str(radius))
            catalog.addFile("sdss", self.__save_path + "catalog.xml", "vo")





            result.addSummaryParams("ra", respond["ra"][index])
            result.addSummaryParams("dec", respond["dec"][index])

            result.addSummaryParams("distance", loss)

            data={"lambda": sdss["u"], "ab": str(respond["u"][index]), "err": str(respond["err_u"][index])}
            result.addSummaryParams("mag_u", data)

            data = {"lambda": sdss["g"], "ab": str(respond["g"][index]), "err": str(respond["err_g"][index])}
            result.addSummaryParams("mag_g", data)

            data = {"lambda": sdss["r"], "ab": str(respond["r"][index]), "err": str(respond["err_r"][index])}
            result.addSummaryParams("mag_r", data)

            data = {"lambda": sdss["i"], "ab": str(respond["i"][index]), "err": str(respond["err_i"][index])}
            result.addSummaryParams("mag_i", data)

            data = {"lambda": sdss["z"], "ab": str(respond["z"][index]), "err": str(respond["err_z"][index])}
            result.addSummaryParams("mag_z", data)

            result.addSummaryParams("type", str(respond["type"][index]))
            result.addSummaryParams("z", respond["z1"][index])
            result.addSummaryParams("zErr", str(respond["zErr"][index]))
            result.addSummaryParams("zWarning", str(respond["zWarning"][index]))

            image = ImageSource(result.getSummary()["id"],self.__catalog_provider)

            bands = ["u","g","r","i","z"]
            displayUrl = self.__service_provider["display"]

            imagesList = self.downloadImages(respond["rerun"][index],respond["run"][index],respond["camcol"][index],respond["field"][index],bands)

            result.addCatalog(catalog)

            for img in imagesList:
                band=str(img["band"]).replace(" ","")
                image.addFile(name="band-"+img["band"],url=img["url"],download=True,local_path=img["local_path"],uncompress=True,thumbnail=True,external=False)
                image.addCutout(self.getImageCutout(respond["ra"][index],respond["dec"][index],query="G"+str(band).upper()),300,band,magnitude=str(respond[band][index])+"+/-"+str(respond["err_"+band][index]),link=displayUrl.format(respond["ra"][index],respond["dec"][index]))


            result.addImage(image)

            if respond["specobjid"][index]!= None:
                spec = self.getSpecImage(respond["specobjid"][index])
                spectra= SpectraSource(result.getSummary()["id"],self.__catalog_provider)
                spectra.addCutout(spec,700,"visible")
                result.addSpectra(spectra)
        else:
            result.addSummaryParams("error","Any results was found in SDSS Archive")

        return result

    def getBand(self,band):

        key_band=""
        if band in self.__key_band:
            key_band=self.__key_band[band]
        return key_band

    def getImageCutout(self,ra,dec,scale=0.396127,width=300,height=300,opt="G",query="GA",low_mag=0,high_mag=30):


        """
        this implementation are based on http://skyserver.sdss.org/dr14/en/tools/chart/chartinfo.aspx.
        This application is based on an underlying web service, ImgCutout.asmx which can be called in many different ways, using the SOAP protocol,
        or just using the standard HTTP GET and PUT interfaces. The formal description is contained in the WSDL, Web Service Description Language document.
        The getjpeg service can be directly called from any web page through the HTTP GET protocol. In order to build a dynamic cutout into your own web page,
        insert the following example. Naturally, replace the parameter values with your own.
        <IMG SRC="http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra=179.689293428354&dec=-0.454379056007667&scale=0.79224&width=400&height=400&opt=GST&query=SR(10,20)">

        :param scale(float): arcsec/pixel (the natural scale of SDSS is 0.396127)
        :param width(int): image width in pixels, limited to [64..2048]
        :param height(int): image height in pixels, limited to [64..2048]
        :param opt(str): options string, a set of upper-case characters, like 'GPST'.
                Drawing options:
                    The characters present will select the corresponding option
                    from the list below. Characters not in the list are ignored.
                    G	Grid	Draw a N-S E-W grid through the center
                    L	Label	Draw the name, scale, ra, and dec on image
                    P	PhotoObj	Draw a small cicle around each primary photoObj
                    S	SpecObj	Draw a small square around each specObj
                    O	Outline	Draw the outline of each photoObj
                    B	Bounding Box	Draw the bounding box of each photoObj
                    F	Fields	Draw the outline of each field
                    M	Masks	Draw the outline of each mask considered to be important
                    Q	Plates	Draw the outline of each plate
                    I	Invert	Invert the image (B on W)
        :param query(str): This option will draw a triangle on top of objects selected by a marking string.
                    Objects must be inside the field of view of the image to be displayed.
                    The format of the string can be from the following choices:
                        ObjType:	S | G | P
                        marks Stars, Galaxies or PhotoPrimary objects.
                        Band:	U | G | R | I | Z | A
                        restricts marks to objects with Band BETWEEN low_mag AND high_mag
                        Band 'A' will mark all objects within the specified magnitude range in any band (ORs composition).
                        Examples:	S
                        SR(0.0, 23.5)
                        GA(20, 30)

        :param low_mag(float): restricts marks to objects with Band BETWEEN low_mag AND high_mag
        :param high_mag(float): restricts marks to objects with Band BETWEEN low_mag AND high_mag
        :return imageUrl(str): url composition for cutout image
        """



        imageUrl=self.__service_provider["cutout"]
        imageUrl=imageUrl.format(self.__data_release)


        imageUrl=imageUrl.replace("ravalue",str(ra)[0:8])
        imageUrl = imageUrl.replace("decvalue", str(dec)[0:8])

        imageUrl = imageUrl.replace("scaleval", str(scale))
        imageUrl = imageUrl.replace("widthval", str(width))
        imageUrl = imageUrl.replace("heightval", str(height))
        imageUrl = imageUrl.replace("optvalue", str(opt))

        queryOpt=query+"("+str(low_mag)+","+str(high_mag)+")"
        imageUrl = imageUrl.replace("queryvalue", queryOpt)

        return imageUrl


    def getSpecImage(self,specobjid):

        specImageUrl = self.__service_provider["spect"]
        specImageUrl=specImageUrl.format(self.__data_release)
        imageUrl = specImageUrl.replace("specobjid", str(specobjid))
        return imageUrl

    def getSpectra(self,rerun,run,camcol,field,bands):
        url='http://data.sdss3.org/sas/dr14/boss/spectro/redux/v5_10_0/spectra/6370/spec-6370-56238-0791.fits'
        url='http://dr14.sdss.org/optical/spectrum/view/data/format=fits?plateid=6370&mjd=56238&fiberid=791&reduction2d=v5_7_0'


    def getImage(self,rerun,run,camcol,field,bands):
        """
        :param rerun:
        :param run:
        :param camcol:
        :param field:
        :return:
        """

        """
        run 6 digits
        camcol 1 digit
        field 4 digits
        """

        format_file="fits.bz2"

        imageUrl=self.__service_provider["image"]


        run_format="0"*(6-len(str(run)))+str(run)
        field_format="0"*(4-len(str(field)))+str(field)

        if type(bands) is list:
            images=[]
            for band in bands:
                if band=="irg":
                    format_file="jpg"
                else:
                    format_file = "fits.bz2"

                imageUrl = imageUrl.format(self.__image_release, rerun, run, camcol, bands, run_format, camcol,
                                           field_format,format_file)
                images.append({"band":band,"url":imageUrl})
            return images

        else:
            if bands == "irg":
                format_file = "jpg"

            imageUrl=imageUrl.format(self.__image_release, rerun, run, camcol,bands,run_format,camcol,field_format,format_file)
            return imageUrl


    #def getCatalog(self,coordinates):

    def downloadImages(self, rerun, run, camcol,field,bands):

        imagePath="http://dr14.sdss.org/sas/dr14/eboss/photoObj/frames/{0}/{1}/{2}/frame-{3}-{4}-{5}-{6}.fits.bz2"
        paths=[]
        for band in bands:

            run2 = run if len(str(run))>=5 else str("0"*(6-len(str(run))))+str(run)
            field2=field if len(str(field))>=4 else str("0"*(4-len(str(field))))+str(field)
            url=imagePath.format(rerun,run,camcol,band,run2,camcol,field2)

            #if Utils.validURL(url):
            path_tmp = self.__save_path+"sdss_band_"+band+".fits"
            paths.append({"band": band, "url": url,"local_path":path_tmp})
            #local_path=Utils.saveLocalFile(url,path_tmp,True)
            #print(local_path)

        return paths

    def getRedshift(self,ra,dec,radio=15):
        radio=radio*u.arcsec
        radio=radio.to("deg").value
        co = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))



        ra_min =co.ra.degree - (radio / 2)
        ra_max = co.ra.degree + (radio / 2)
        dec_min = co.dec.degree - (radio / 2)
        dec_max = co.dec.degree + (radio / 2)
        query="""SELECT p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,p.err_u, p.err_g, p.err_r, p.err_i, p.err_z, p.run, p.rerun, p.camcol,p.type as photo_type ,p.field,s.specobjid,s.sourceType as s_type, s.class as spec_type, s.z as spec_z,s.zErr as spec_zerr, s.plate, s.mjd, s.fiberid,pz.z as photo_z, pz.zErr as photo_zerr
        FROM PhotoObjAll AS p
           left outer join SpecObj AS s ON p.objid = s.bestobjid
           left outer join Photoz as pz on p.objid = pz.objid
        WHERE 
           p.ra BETWEEN {0} AND {1}
           AND p.dec BETWEEN {2} AND {3}"""

        json_result = self.SQLRequest(query.format(ra_min,ra_max,dec_min,dec_max))

        dp_result = pd.DataFrame(json_result[0]["Rows"])

        if len(json_result[0]["Rows"])>0:

            ra_ref = dp_result["ra"].tolist()
            dec_ref = dp_result["dec"].tolist()
            cref = SkyCoord([ra], [dec], frame='icrs', unit='deg')
            c1 = SkyCoord(ra_ref, dec_ref, frame='icrs', unit='deg')
            desi_distance = cref.separation(c1).arcsec
            dp_result["separation"] = desi_distance

        return dp_result



    def SQLRequest(self,command):
        url_service=self.__service_provider["sql_server"].format(self.__data_release)
        params= {"cmd":command,"format":"json","searchtool":"SQL","TaskName":"Skyserver.Search.SQL","syntax":"NoSyntax"}
        json_file = OrderedDict(params)
        r=requests.get(url_service,params=json_file)
        return r.json()

    @staticmethod
    def getSpecType(id):

        types = ["GALAXY", "QSO", "STAR_BHB", "STAR_CARBON", "STAR_BROWN_DWARF", "STAR_SUB_DWARF", "STAR_CATY_VAR", "STAR_RED_DWARF", "STAR_WHITE_DWARF", "REDDEN_STD","SPECTROPHOTO_STD","HOT_STD","ROSAT_A","ROSAT_B","ROSAT_C","ROSAT_D","SERENDIPITY_BLUE","SERENDIPITY_FIRST","SERENDIPITY_RED","SERENDIPITY_DISTANT","SERENDIPITY_MANUAL","QA","SKY","NA","STAR_PN"]
        if type(id) is int and id < len(types):
            return types[id]
        else:
            return "ANYTYPE"


    @staticmethod
    def getPhotoType(id):


        types=["UNKNOWN","COSMIC_RAY","DEFECT","GALAXY","GHOST","KNOWNOBJ","STAR","TRAIL","SKY","NOTATYPE"]
        if type(id) is int and id < len(types):
            return types[id]
        else:
            return "ANYTYPE"

    def contentUrl(self):
        pass
    def delete(self, **kwargs):
        pass
    def getType(self):
        pass
    def insert(self, **kwargs):
        pass
    def update(self, **kwargs):
        pass