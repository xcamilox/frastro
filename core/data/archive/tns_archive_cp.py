from django import forms
import requests
import csv
import pandas as pd
from astropy.table import Table
import json
from collections import OrderedDict

from frastro.frastro.core.data.archive.sussex_archive_cp import SussexArchive
from frastro.frastro.core.data.archive.lasair_archive_cp import LasairArchive
from frastro.frastro.core.database.mongodb.mongodb_manager import MongodbManager

tns_broker_url = "https://wis-tns.weizmann.ac.il/api/get/search"
search_sn = "https://wis-tns.weizmann.ac.il/search?&isTNS_AT=all&classified_sne=1&num_page=500&format=csv&page={0}"

############################# PARAMETERS #############################
# API key for Bot                                                    #
api_key = "5443cf3bf665514d3059e991f05fec790b78ee91"
api_key="b2fbc4c3e3fa2d4b3b7b9ac12b24617ec3a8598e" #13 february
# list that represents json file for search obj                      #
search_obj=[("ra",""), ("dec",""), ("radius",""), ("units",""),      #
            ("objname",""), ("internal_name","")]                    #
# list that represents json file for get obj                         #
get_obj=[("objname",""), ("photometry","0"), ("spectra","1")]        #
######################################################################

#############################    URL-s   #############################
# url of TNS and TNS-sandbox api                                     #
url_tns_api="https://wis-tns.weizmann.ac.il/api/get"                 #
url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get"     #
######################################################################


class TNSServices():

    downloadeddata=None
    results=[]

    @classmethod
    def downloadCSV(self,currentPage):

        print("current page downloades",currentPage)
        URLPAGE = search_sn.format(currentPage)
        print(URLPAGE)
        with requests.Session() as s:
            download = s.get(URLPAGE)

            decoded_content = download.content.decode('utf-8')

            cr = csv.reader(decoded_content.splitlines(), delimiter=',')
            data_list = list(cr)

            header = data_list.pop(0)

            self.results += data_list

            if len(data_list) <= 0:


                self.downloadeddata = pd.DataFrame(self.results,columns=header)
                base_path = "/Users/camiloj/projects/broker_data/all/"
                self.downloadeddata.to_csv(base_path+'tns_sndetections.csv',sep=",")
                return self.downloadeddata
            else:
                next_page = currentPage+1
                return self.downloadCSV(next_page)


    @classmethod
    def getDownloadSNDetections(self):

        current_page = 0

        data_table = self.downloadCSV(current_page)
        return data_table

    @classmethod
    def get_data_object(self,objectId):
        try:
            search_url = url_tns_api + '/object'
            ######################################################################
            # api key for your Bot
            get_obj = [("objname", objectId), ("photometry", "1"), ("spectra", "1")]
            json_file = OrderedDict(get_obj)
            get_data = [('api_key', (None, api_key)),  #
                        ('data', (None, json.dumps(json_file)))]  #
            # get obj using request module
            response = requests.post(search_url, files=get_data)

            # parsed = json.loads(response.text, object_pairs_hook=OrderedDict)  #
            # result = json.dumps(parsed, indent=4)  #
            return response.json()

        except Exception as e:  #
            return [None, 'Error message : \n' + str(e)]




    @classmethod
    def getLightCurves(self, file_path):
        db = MongodbManager()


        db.setDatabase("iacbroker")
        db.setCollection("tnssn")
        df_tns_reports = pd.read_csv(file_path)
        df_tns_reports=df_tns_reports.fillna(" ")
        lasair = LasairArchive()

        for index, row in df_tns_reports.iterrows():
            if index<4748:
                print("already")
                continue
            try:
                list_data = json.loads(row.to_json())
                print(index, row["objindex"], row["DiscInternalName"])
                if type(row["DiscInternalName"])==str and "ZTF" in row["DiscInternalName"]:
                    print("get light curve", row["DiscInternalName"])
                    curve=lasair.getObjectInfo(row["DiscInternalName"])
                    list_data["lightcurve"]=curve

                print("save db", row["DiscInternalName"])
                db.saveData(list_data)
            except Exception as e:
                print("error with ",row["objindex"], row["DiscInternalName"], e.args)

        print("saved desi table")

    @classmethod
    def machTNSMassivegalaxiesDESI(self,file_path,xml_path_to_save):
        df_tns_reports = pd.read_csv(file_path)
        table_tns_report = Table.from_pandas(df_tns_reports)
        dataarchive = SussexArchive()

        desi_targetsvo, desi_targetstable = dataarchive.getTNSdetectionOnDesiMassiveGalaxies(table_tns_report)
        desi_targetsvo.to_xml(xml_path_to_save)
        print("saved desi table")

    @classmethod
    def searchOBJ(self, ra,dec,radio=5):

        params = {"ra": ra,
                  "dec": dec,
                  "radius": radio,
                  "units": "arcsec",
                  "objname": "",
                  "internal_name": ""}
        try:  #
            # url for search obj                                             #
            search_url = url_tns_api + '/search'  #
            # change json_list to json format                                #
            json_file = OrderedDict(params)  #
            # construct the list of (key,value) pairs                        #
            search_data = [('api_key', (None, api_key)),  #
                           ('data', (None, json.dumps(json_file)))]  #
            # search obj using request module                                #
            response = requests.post(search_url, files=search_data)  #
            # return response                                                #
            return response.json()  #
        except Exception as e:  #
            return [None, 'Error message : \n' + str(e)]  #


if __name__ == "__main__":



    tns=TNSServices()
    #tns.getDownloadSNDetections()
    base_path = "/Users/camiloj/projects/broker_data/all/"
    #tns.getLightCurves(base_path+"tns_sndetections_ok.csv")
    tns.machTNSMassivegalaxiesDESI(base_path+"tns_sndetections_todesi.csv",base_path+"tns_vs_desi1000.xml")


    # data=tns.searchOBJ(2.55475, -19.69235)
    # #data=tns.searchOBJ(45.990417,43.401)
    # for items in data['data']['reply']:
    #     objname=items["objname"]
    #     content_data=tns.get_data_object(objname)
    #     print(content_data)
