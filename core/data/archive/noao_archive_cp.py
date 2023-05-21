import numpy as np
import pylab as plt
import matplotlib
from astropy import utils, io
from getpass import getpass
from astropy.visualization import make_lupton_rgb
from pyvo.dal import sia
from frastro import CoordinateParser, ImageUtils
from dl import authClient as ac, queryClient as qc
from dl import storeClient as sc, helpers
import json
from astropy import units as u
from frastro import ImageSource

class NOAOArchiveCP():

    # Datalab and related imports
    # You'll need at least these for authenticating and for issuing database queries

    # You'll want storeClient if you plan to use virtual storage or myDB
    # Get helpers for various convenience function


    # To get image cutouts, you'll need the VO-based SIA package, and define which SIA service
    # to use
    __TOKEN_SESSION=""
    DEF_ACCESS_URL = "http://datalab.noao.edu/sia/des_dr1"  # DES SIA service URL


    def __init__(self):
        self._svc = sia.SIAService(self.DEF_ACCESS_URL)


    def login(self,user="",password=""):
        if user != "" and password !="":
            token = ac.login(user,password)
        else:
            token = ac.login('anonymous')

        self.__TOKEN_SESSION=token
        return token

    def logout(self,token=""):
        if token=="":
            token=self.__TOKEN_SESSION

        ac.logout(token)


    def getAvaliableDataset(self):
        # set default profile
        qc.set_profile('default')

        # these schemas are not astronomical datasets
        _remove = set(['ivao', 'ivao_smash', 'tap_schema', 'schema'])

        # get all schemas from DB
        _schemas = set(qc.query(self.__TOKEN_SESSION, sql='SELECT schema FROM tbl_stat').split())

        # remove non-astro schemas
        _alldatasets = sorted(list(_schemas - _remove))
        print("Datasets available in Data Lab (with current profile):\n", _alldatasets)

    def getDesCatalog(self):
        # The pre-release DES DR1 profile
        try:
            qc.set_profile('des-proto')
        except Exception as e:
            print(e)

        try:
            schema=qc.schema('des_dr1', format='json', profile='des-proto')
            if type(schema) is not json:
                tmp=schema.decode('utf8').replace("'", '"')
                data = json.loads(tmp)
                s = json.dumps(data, indent=4, sort_keys=True)
                print(s)
            print(schema)
        except Exception as e:
            print(e)

    def getToken(self):
        if self.__TOKEN_SESSION=="":
            self.login()
        return self.__TOKEN_SESSION

    def desQuery(self,ra,dec,radius_arcmin=1,columns="*"):

        radius_degree = CoordinateParser.getMinToDegree(radius_arcmin)

        query_template = "SELECT {0:s} FROM des_dr1.main WHERE q3c_radial_query(ra,dec,{1:f},{2:f},{3:f})"
        query = query_template.format(columns, ra, dec, radius_degree)
        df = None
        try:
            result = qc.query(self.getToken(), sql=query)  # by default the result is a CSV formatted string
            result = result.decode('utf8').replace("'", '"')
            df = helpers.convert(result, 'pandas')
        except Exception as e:
            print(e)
        return df

    def download_deepest_image(self,ra, dec, radius_arcmin=1, path=''):
        fov=CoordinateParser.getMinToDegree(radius_arcmin)
        size=fov / np.cos(np.array(dec) * np.pi / 180)
        imgTable = self._svc.search((ra, dec),  (size, fov), verbosity=2)
        imgTable=imgTable.votable.to_table()
        print("The full image list contains", len(imgTable), "entries")
        #sel0 = imgTable['obs_bandpass'].astype(str) == band
        #sel0 &
        sel =  ((imgTable['proctype'].astype(str) == 'Stack') & (
                    imgTable['prodtype'].astype(str) == 'image'))  # basic selection
        Table = imgTable[sel]  # select


        imageSource = None
        if (len(Table) > 0):
            imageSource = ImageSource(str(ra)+"_"+str(dec), "DES")

            for row in Table:
                band=row['obs_bandpass'].decode()

                if band=='':
                    band="stackbands"
                url = row['access_url'].decode()
                #image = io.fits.getdata(utils.data.download_file(url, cache=True, show_progress=False, timeout=120))
                #base64 = ImageUtils.imgNpArrayToBase64(image)
                local_path = path + "des_band_" + band + ".fits"
                imageSource.addFile(band, url, "fits", download=True, local_path=local_path,uncompress=False, thumbnail=True,external=False)

                #imageSource.addFile("band-" + band, url, "fits")
                thumbnail=imageSource.getFiles()
                pos =len(thumbnail)-1
                base64img=thumbnail[pos]["thumbnail"]
                imageSource.addCutout(base64img, 300, band)


            # if (len(Table) > 0):
            #     row = Table[np.argmax(Table['exptime'].data.data.astype('float'))]  # pick image with longest exposure time
            #     url = row['access_url'].decode()  # get the download URL
            #     print('downloading deepest stacked image...')
            #     image = io.fits.getdata(utils.data.download_file(url, cache=True, show_progress=False, timeout=120))
            #
            # else:
            #     print('No image available.')
            #     image = None

        return imageSource





    # multi panel image plotter
    def plot_images(self,images, geo=None, panelsize=4, bands=list('gri'), cmap=matplotlib.cm.gray_r):
        n = len(images)
        if geo is None: geo = (n, 1)

        fig = plt.figure(figsize=(geo[0] * panelsize, geo[1] * panelsize))
        for j, img in enumerate(images):
            ax = fig.add_subplot(geo[1], geo[0], j + 1)
            if img is not None:
                print(img.min(), img.max())
                ax.imshow(img, origin='lower', interpolation='none', cmap=cmap,
                          norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=img.max()))
                ax.set_title('%s band' % bands[j])
                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)
