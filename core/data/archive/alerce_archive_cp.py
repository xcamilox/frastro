
import json
from astropy.time import Time

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import requests
from PIL import Image
from io import BytesIO
import logging

import pandas as pd
from astropy.table import Table, Column, vstack

class AlerceArchive():

    query_url = "http://ztf.alerce.online/query"
    cross_match_url="https://catshtm.alerce.online/crossmatch_all"


    get_stapms = "http://avro.alerce.online/get_stamp"
    get_url_detections = "http://ztf.alerce.online/get_detections"
    get_url_nodetections = "http://ztf.alerce.online/get_non_detections"
    get_url_probabilities = "http://ztf.alerce.online/get_probabilities"
    get_url_features = "http://ztf.alerce.online/get_features"
    get_url_stats = "http://ztf.alerce.online/get_stats"

    result=pd.DataFrame()

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    @classmethod
    def getStamps(self,ZTFObjectId,candid,format=["png"],type=["science", "template", "difference"]):
        params={
            "oid": ZTFObjectId,
            "candid": candid,
            "format": format, # ["png", "fits"]
            "type": type # ["science", "template", "difference"]
        }
        r = requests.get(self.get_stapms,params=params)
        return r

    @classmethod
    def getLightCurve(self,ZTFObjectId):

        params = {
            "oid": ZTFObjectId
        }
        r = requests.post(self.get_url_detections,json=params)
        return r.json()

    @classmethod
    def getStats(self,ZTFObjectId):

        params = {
            "oid": ZTFObjectId
        }
        r = requests.post(self.get_url_stats,json=params)
        return r.json()

    @classmethod
    def getFeatures(self, ZTFObjectId):

        params = {
            "oid": ZTFObjectId
        }
        r = requests.post(self.get_url_features, json=params)
        return r.json()

    #ra ,dec in degrees and radius in arcsec
    @classmethod
    def crossMatch(self, ra, dec, radius):
        parameters={
            "ra":ra,
            "dec":dec,
            "radius":radius
        }

        r = requests.get(self.cross_match_url,params=parameters)
        if r.status_code == 500:
            return {"error":r.reason}
        else:
            return r.json()

    @classmethod
    def getProbabilities(self,ZTFObjectId):
        params = {
            "oid": ZTFObjectId
        }
        r = requests.post(self.get_url_probabilities, json=params)
        return r.json()

    @classmethod
    def showStamps(self,ZTFObjectId,observationId):
        imgs = self.getStamps(ZTFObjectId, observationId, type=["science"])
        im = Image.open(BytesIO(imgs.content))

        plt.subplot(131)
        plt.imshow(im)

        imgs = self.getStamps(ZTFObjectId, observationId, type=["template"])
        im = Image.open(BytesIO(imgs.content))

        plt.subplot(132)
        plt.imshow(im)

        imgs = self.getStamps(ZTFObjectId, observationId, type=["difference"])
        im = Image.open(BytesIO(imgs.content))

        plt.subplot(133)
        plt.imshow(im)
        plt.show()

    @classmethod
    def getCandidates(self, days_ago=15,page=1):
#        self.logger.info('getting Alerce detections in last {0} days'.format(str(days_ago)))
        total = 1000000
        records_per_page = 5000

        sortBy = "min_magap_g"
        nobsmin = 1
        nobsmax = 40
        classearly = [19]
        lateclass=[10,11,12,13,14]
        pclassrf = 0.8
        d = datetime.today() - timedelta(days=days_ago)

        t = Time(d.strftime("%Y-%m-%dT%H:%M:%S"), format='isot')
        print(t.mjd,t)
        params = {
        "total": total,
        "records_per_pages": records_per_page,
        "page": page,
        "sortBy": sortBy,
        "query_parameters": {
            "dates": {
                "firstmjd": {
                    "min": int(t.mjd)
                }
            },
            "filters": {
                "nobs": {
                    "min": 3

                },
                # "min_magap_g":{
                #     "min": 17.0
                # },
                # "min_magap_r": {
                #     "min": 17.0
                # },
                "classearly":classearly,
                #"classrf": lateclass,
                #"pclassearly": pclassrf
                }
            }
        }

        query_results = requests.post(self.query_url, json=params)




        data = query_results.json()

        # with open('alerce.json', 'w') as outfile:
        #     json.dump(data, outfile)

        alerceGoodCandidates = [data["result"][target] for target in data["result"]]
        pd_result=pd.DataFrame(alerceGoodCandidates)

        print("total",len(self.result),len(pd_result))
        if data["page"] > 1:
            self.result = pd.concat([self.result, pd_result])
        else:
            self.result = pd_result
        if len(pd_result) >= 5000:
            page += 1
            return self.getCandidates(days_ago=days_ago, page=page)
        else:
            return self.result






#getLast15days()