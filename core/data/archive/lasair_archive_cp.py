import requests
import logging
import json
import pandas as pd
from astropy.table import Table, Column, vstack
import numpy as np
from scipy import stats

from frastro.frastro.core.database.mongodb.mongodb_manager import MongodbManager
from frastro.frastro.core.utils.convertions_util import Convertion

class LasairArchive():

    url_lasair = "https://lasair.roe.ac.uk/object/"
    #min  num of candidate > 2
    #min 0.6 days between last detections
    filterByPassDetections={
       "selected":"""objects.objectId as oid, objects.ramean as meanra, objects.decmean as meandec, objects.jdmin - 2400000.5 AS mjdmin, objects.jdmax - 2400000.5 AS mjdmax, objects.jdmax - 2400000.5 AS lastmjd, objects.maggmin, objects.magrmax, objects.magrmin, objects.maggmax, objects.maggmin, latestrmag, latestgmag, comments.content as comment, sherlock_classifications.classification, IF(objects.distpsnr1 < 2 AND objects.sgscore1 > 0.49,"Within 2arcsec of PS1 star","Not Near PS1 star") score, objects.ncand as detections, sherlock_classifications.classification as classearl""",
       "tables":"objects,sherlock_classifications,comments",
       "put":'on',
       "conditions":"""objects.objectId = comments.objectId AND objects.primaryId = sherlock_classifications.transient_object_id AND objects.objectId = objects.objectId AND objects.primaryId = sherlock_classifications.transient_object_id AND objects.objectId = objects.objectId AND comments.objectId = objects.objectId AND objects.primaryId = sherlock_classifications.transient_object_id AND objects.jdmin > JDNOW() - {0} AND objects.jdmax - objects.jdmin > 0.6 AND sherlock_classifications.classification NOT IN ("VS" , "AGN", "CV", "BS") AND objects.ncand > 2 AND objects.ncandgp > 1 ORDER BY score , mjdmin DESC""",
       "check_days_ago":"on","days_ago":"{0}","page":"{0}","json":"on"
    }

    result=pd.DataFrame()

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    @classmethod
    def getLastDetectionsCandidates(self):
        detections = self.getLastDetections(20)
        result = iter([row for index, row in detections])
        return result


    @classmethod
    def getMinMaxDetections(self,measures):
        # peak is min value
        mjds=measures[:,0]
        gdetections = [str(item)[0:str(item).find(".")] for item in mjds]
        gdetectionsdays = np.unique(np.array(gdetections))
        numObs = len(gdetectionsdays)

        peak_idx = np.where(measures[:, 1] == np.min(measures[:, 1]))
        peak = measures[peak_idx][0]
        mindetections_idx = np.where(measures[:, 0] < peak[0])
        if mindetections_idx[0].shape[0]<=0:
            mindetections_idx = np.where(measures[:, 0] <= peak[0])
        rise_peak = measures[mindetections_idx]
        mindetections = np.where(rise_peak[:, 1] == np.max(rise_peak[:, 1]))
        min_peak = rise_peak[mindetections][0]
        mjd_peak = peak[0]
        fall_curve = np.where(measures[:,0]>mjd_peak,1,0)
        isfall = fall_curve.sum(axis=0)
        status="fall" if isfall>0 else "rise"
        if numObs<=1:
            status="point"

        return {"peak":peak,"mindet":min_peak,"peakidx":peak_idx[0].tolist()[0],"midx":mindetections[0].tolist()[0],"status":status,"numdetections":numObs }


    @classmethod
    def getPeakLightCurve(self, candidates):

        data=[]

        for index,detections in enumerate(candidates):

            if "drb" in detections.keys() or "rb" in detections.keys(): #good detections

                mjd=detections["mjd"]
                magpsf=detections["magpsf"]
                sigmapsf = detections["sigmapsf"]
                fid = detections["fid"]
                data.append([mjd,magpsf,sigmapsf,fid])
        data = np.array(data)

        #Filter ID (1=g; 2=r; 3=i)
        gfilter = np.where(data[:, 3] == 1)
        gmeasures = data[gfilter]

        rfilter = np.where(data[:,3]==2)
        rmeasures = data[rfilter]

        peak_data ={"peak":{},"min_det":{},"stats":{},"lightcurve":{},"status":{"r":"anydet","g":"anydet"}}

        if gmeasures.shape[0]>0:
            gpeak=self.getMinMaxDetections(gmeasures)
            peak_data["peak"]["g"]=gpeak["peak"].tolist()
            peak_data["min_det"]["g"] = gpeak["mindet"].tolist()
            peak_data["peak"]["g"].append(gpeak["midx"])
            peak_data["status"]["g"]=gpeak["status"]
            peak_data["min_det"]["g"].append(gpeak["peakidx"])
            x_g = [0,(peak_data["peak"]["g"][0] - peak_data["min_det"]["g"][0])]
            y_g = [peak_data["min_det"]["g"][1],peak_data["peak"]["g"][1]]
            slope_g, intercept_g, r_value_g, p_value_g, std_err_g = stats.linregress(x_g, y_g)
            angle_g = np.arctan(slope_g)+360
            #gmeasures.sort(axis=0)
            peak_data["stats"]["g"]={"days_peak":(peak_data["peak"]["g"][0] - peak_data["min_det"]["g"][0]),
                                     "slope":slope_g,"angle":angle_g,"x":x_g,"y":y_g,
                                     "deltax":(peak_data["peak"]["g"][0] - peak_data["min_det"]["g"][0]),"deltay":(peak_data["min_det"]["g"][1]-peak_data["peak"]["g"][1]),"peakmag":min(y_g),"peakmjd":peak_data["peak"]["g"][0]}
            peak_data["lightcurve"]["g"] = {"days":(gmeasures[:,0]-gmeasures[:,0].min()).tolist(),"mjd": (gmeasures[:, 0]).tolist(), "mag": gmeasures[:, 1].tolist(), "magerr": gmeasures[:, 2].tolist()}
            peak_data["lightcurve"]["g"]["detections"]=gpeak["numdetections"]

        if rmeasures.shape[0]>0:
            rpeak=self.getMinMaxDetections(rmeasures)
            peak_data["peak"]["r"] = rpeak["peak"].tolist()
            peak_data["min_det"]["r"] = rpeak["mindet"].tolist()
            peak_data["peak"]["r"].append(rpeak["midx"])
            peak_data["status"]["r"] = rpeak["status"]
            peak_data["min_det"]["r"].append(rpeak["peakidx"])
            x_r = [0, (peak_data["peak"]["r"][0] - peak_data["min_det"]["r"][0])]
            y_r = [peak_data["min_det"]["r"][1], peak_data["peak"]["r"][1]]
            slope_r, intercept_r, r_value_r, p_value_r, std_err_r = stats.linregress(x_r, y_r)
            angle_r = np.arctan(slope_r)+360

            peak_data["stats"]["r"] = {"days_peak": (peak_data["peak"]["r"][0] - peak_data["min_det"]["r"][0]),
                                       "slope": slope_r, "angle": angle_r, "x": x_r, "y": y_r,
                                       "deltax":(peak_data["peak"]["r"][0] - peak_data["min_det"]["r"][0]),"deltay":(peak_data["min_det"]["r"][1]-peak_data["peak"]["r"][1]),"peakmag":min(y_r),"peakmjd":peak_data["peak"]["r"][0]}
            peak_data["lightcurve"]["r"]={"days":(rmeasures[:,0]-rmeasures[:,0].min()).tolist(),"mjd":(rmeasures[:,0]).tolist(),"mag":(rmeasures[:,1]).tolist(),"magerr":(rmeasures[:,2]).tolist()}
            peak_data["lightcurve"]["r"]["detections"]=rpeak["numdetections"]

        if rmeasures.shape[0]>0 and gmeasures.shape[0]>0:
            delta_g_r_magpeak= peak_data["stats"]["g"]["peakmag"] - peak_data["stats"]["r"]["peakmag"]
            peak_data["stats"]["delta_g_r_magpeak"] = delta_g_r_magpeak
            delta_g_r_mjdpeak = peak_data["stats"]["g"]["peakmjd"] - peak_data["stats"]["r"]["peakmjd"]
            peak_data["stats"]["delta_g_r_mjdpeak"] = delta_g_r_mjdpeak

        return peak_data

    @classmethod
    def getObjectInfo(self, idobject):
        r = requests.get(self.url_lasair + idobject + "/json")
        return r.json()

    @classmethod
    def generateLightCurveImage(self,lightcurve_json):
        pass

    # get NULL,UNCLEAR, ORPHAN, and SN
    @classmethod
    def getLastDetections(self, days_ago=14,page=0):
        print('getting lasair last detections in last {0} days'.format(days_ago))
        url_lasair = "https://lasair.roe.ac.uk/objlist/"
        DAYS_AGO = str(days_ago)
        self.filterByPassDetections["conditions"]=self.filterByPassDetections["conditions"].format(DAYS_AGO)
        self.filterByPassDetections["days_ago"]=self.filterByPassDetections["days_ago"].format(DAYS_AGO)
        self.filterByPassDetections["page"] = self.filterByPassDetections["page"].format(page)
        r = requests.post(url_lasair, data=self.filterByPassDetections)

        pd_result = pd.DataFrame(r.json())
        print(len(self.result))
        if len(self.result)>0:
            self.result = pd.concat([self.result, pd_result])
        else:
            self.result = pd_result
        if len(pd_result)>=1000:
            page+=1
            return self.getLastDetections(days_ago=days_ago, page=page)
        else:
            return self.result



if __name__ == "__main__":
    db = MongodbManager()

    db.setDatabase("iacbroker")
    db.setCollection("tnssn")
    data_sn=db.getData(filter={"lightcurve":{"$exists":True}},projection={"lightcurve":1,"Redshift":1,"Name":1})
    for idenx,row in enumerate(data_sn):
        name = row["Name"]
        print(name)
        ls = LasairArchive()
        pick = ls.getPickLightCurve(row)

        db.update({"Name":{"$eq":name}},{ "$addToSet": {"pick": pick}})