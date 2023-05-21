import sncosmo
import pandas as pd
import requests
from astropy.table import Table
import numpy as np
import json
import math
from pandas.io.json import json_normalize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.font_manager import FontProperties

base_data_path = "/Users/camiloj/projects/broker_data/"
file_path_report = "/Users/camiloj/projects/broker_data/reports/"

class SNCosmosFit():

    url_lasair = "https://lasair.roe.ac.uk/object/"

    light_curve=""
    def __init__(self):
        self.avaliable_models=sncosmo.builtins._SOURCES.get_loaders_metadata()
        plt.rcParams.update({'font.size': 8})

    def getAvaliableModels(self):
        return self.avaliable_models

    def loadTnsSn(self,file_path_csv):
        df_sn_list=pd.read_csv(file_path_csv)
        num_ztf=0
        for index,row in df_sn_list.iterrows():

            if type(row["DiscInternalName"])==str and "ZTF" in row["DiscInternalName"]:
                num_ztf+=1
                self.fitSNByZTFId(row["DiscInternalName"], save_on=file_path_report,
                                           reference={"redshift": row["Redshift"], "type": 0, "label": str(row["Name"]).replace("SN ","")+row["ObjType"],"type_ob":row["ObjType"]})


    def generatePlots(self,data_path,mark="",plot_pdf=""):
        data=pd.read_csv(data_path)
        fig, ax = plt.subplots(2,1,constrained_layout=True)
        fields=data.groupby("typeidx")


        for groups_type,groups_data in fields:
            ax[0].scatter(groups_data["redshift"].tolist(), groups_data["chi2"].tolist(), alpha=0.5,label=groups_data["type"].tolist()[0])


        if mark!="":
            ax[0].scatter([mark["redshift"]],[0.0],color="black",alpha=1,marker="+",label=mark["label"]+" "+mark["type_ob"])
        ax[0].legend(loc="upper right",title="SN type")
        ax[0].grid(True)
        ax[0].set_xlabel("Redshift",fontsize=9)
        ax[0].set_ylabel("Chi2",fontsize=9)
        fontP = FontProperties()
        fontP.set_size('small')

        ax[0].set_title(mark["title"],fontsize=8)
        ax[0].tick_params(which='minor', width=0.75, length=2.5, labelsize=5)


        if "best" in mark:
            detections_points=mark["best"]["data_table"]
            ax[0].scatter([mark["best"]["redshift"]], [mark["best"]["chi2"]], color="black", alpha=0.7, marker="x", label="Best: "+mark["best"]["model"]+" "+mark["best"]["type"])
            color = detections_points["band"]
            color_index_green = np.where(color == "ztfg",True,False)
            color_index_red = np.where(color == "ztfr", True, False)
            green_color="#31FB3B80"
            red_color = "#FC202480"
            ax[1].errorbar(detections_points["time"][color_index_green],detections_points["mag"][color_index_green], yerr=detections_points["mag_err"][color_index_green], fmt='o',color=green_color)
            ax[1].errorbar(detections_points["time"][color_index_red], detections_points["mag"][color_index_red],
                           yerr=detections_points["mag_err"][color_index_red], fmt='o', color=red_color)
            ax[1].invert_yaxis()
            ax[1].grid(True)

            sncosmo.plot_lc(mark["best"]["data_table"], model=mark["best"]["fitted_model"], errors=mark["best"]["errors"], figtext="BEST FIT\n"+mark["best"]["fit_title"], fname=plot_pdf,
                            format='pdf', figtextsize=2.)
        ax[0].legend(prop=fontP, markerscale=0.7, loc='center left', bbox_to_anchor=(1, 0.5))

        if plot_pdf != "":
            plot_pdf.savefig(fig)
            plot_pdf._file.pageList.insert(0,plot_pdf._file.pageList[-1])
            plot_pdf._file.pageList.insert(1, plot_pdf._file.pageList[-2])
            plot_pdf._file.pageList.pop()
            plot_pdf._file.pageList.pop()

            print("save in pdf")
        else:
            plt.show()







    def fitSNByZTFId(self,ZTFID,save_on="",reference=""):

        obj_info = self.getLasairObInfo(ZTFID)
        table_info = self.getTableFromJson(obj_info)
        pp=None
        if save_on!="":

            pp = PdfPages(file_path_report+ZTFID+'.pdf')
        if "redshift" in reference:
            report,best=self.fitLightCurve(table_info, title=ZTFID,report=pp,redshift_range=(reference["redshift"]-0.1,reference["redshift"]+0.1))
        else:
            report, best = self.fitLightCurve(table_info, title=ZTFID, report=pp)
        pd.DataFrame(report).to_csv(file_path_report+ZTFID+'.csv')



        if pp!=None:
            txt_result = "best fit:" + str(best["chi2"]) + " model: " + best["model"] + " version:" + best[
                "version"] + " type: " + best["type"]
            if reference != "":
                reference["title"] = txt_result
                reference["best"] = best
                self.generatePlots(file_path_report + ZTFID + '.csv', mark=reference, plot_pdf=pp)



            pp.attach_note(txt_result)
            pp.close()


    def fitLightCurve(self,data_table,title="",redshift_range=(0.005, 0.2),report=None,usemodels=[]):

        models=self.getAvaliableModels()

        results=[]
        best=[]
        min_result={"model":"xx","chi2":999,"type":""}
        types=[]
        for index,model_ref in enumerate(models):
            plt.clf()
            name_model = model_ref["name"]
            if len(usemodels)>0 and name_model not in usemodels:
                continue
            supernovea_type = model_ref["type"]
            model_version= model_ref["version"]
            zp_mean = np.mean(data_table["zp"])
            if supernovea_type not in types:
                types.append(supernovea_type)

            type_indx=types.index(supernovea_type)


            try:
                model = sncosmo.Model(source=name_model)  # model.set(z=0.064)
                res, fitted_model = sncosmo.fit_lc(data_table, model, model.param_names, bounds={'z': redshift_range},guess_z=True)
                if res.ndof > 0:
                    redchi2_allpts = res.chisq / res.ndof  # reduced chi2
                else:
                    redchi2_allpts = res.chisq  # reduced chi2

                if res.success:
                    title_fig=title+'\nmodel:'+name_model+'\ntype:'+supernovea_type+'\nchisq:' + str(round(redchi2_allpts, 2))
                    if report != None:
                        sncosmo.plot_lc(data_table, model=fitted_model, errors=res.errors, figtext=title_fig, fname=report, format='pdf',figtextsize=2.,zp=zp_mean)
                    else:
                        sncosmo.plot_lc(data_table, model=fitted_model, errors=res.errors,
                                    figtext=title_fig,figtextsize=2.,zp=zp_mean) #redchi2_allpts

                        plt.show()

                else:
                    print("any solution found: "+title, name_model, supernovea_type,res.chisq,res.errors)




                results.append({"model":name_model,"type":supernovea_type,"version":model_version,"chi2":redchi2_allpts,"redshift":fitted_model.parameters[0],"typeidx":type_indx,"min" :abs(redchi2_allpts-1)})
                best.append(redchi2_allpts)
                #print("Fiting model", name_model, model_version, supernovea_type)
                #print(abs(redchi2_allpts-1) < min_result["chi2"], min_result["chi2"],abs(redchi2_allpts-1),redchi2_allpts)
                #print("model fit:", redchi2_allpts, " redshift ", fitted_model.parameters[0])
                if abs(redchi2_allpts-1) < abs(min_result["chi2"]-1):
                    min_result["chi2"] = redchi2_allpts
                    min_result["model"] = name_model
                    min_result["type"] = supernovea_type
                    min_result["version"]=model_version
                    min_result["redshift"]= fitted_model.parameters[0]
                    min_result["typeidx"]=type_indx
                    #min_result["data_table"]=data_table
                    #min_result["fitted_model"]=fitted_model
                    #min_result["errors"] = res.errors
                    #min_result["fit_title"]=title_fig



            except Exception as error_msg:
                print("error with model ",name_model,supernovea_type,model_version,error_msg)


        return results,min_result

    def getTableFromJson(self,data_info):

        # ************************************************************
        # get fields of interest from json and produce input table for sncosmo
        # ************************************************************


        table_json = json_normalize(data_info)

        data = Table.from_pandas(table_json)
        # define fields for sncosmo imput table:

        time = []  # np.array([])
        band = []  # np.str([])
        mag = []  # np.array([])
        mag_err = []  # np.array([])
        zp = []  # np.array([])
        zpsys = []  # np.str([])
        flux = []  # np.array([])
        # define fields from data to be scanned

        for ind,row in enumerate(data):
            if row["candid"]>0:
                time.append(row["mjd"])
                mag.append(row["magpsf"])
                mag_err.append(row["sigmapsf"])
                zp.append(row["magzpsci"])

                zpsys.append('ab')
                # fid 1=g; 2=r; 3=i
                if row['fid'] == 1:

                    band.append('ztfg')
                else:
                    band.append('ztfr')


        # convert photometry to flux inputs with errors for sncosmo:

        ab = sncosmo.get_magsystem('ab')
        n = len(mag)
        index = np.arange(n)

        for i in index:
            flux.append(ab.band_mag_to_flux(mag[index[i]], band[index[i]]))

        flux_err = np.multiply(np.multiply(np.multiply(flux, 0.4), mag_err), np.log(10))

        # produce table:
        data_draft = Table([time, band, flux, flux_err, zp, zpsys,mag,mag_err],
                           names=('time', 'band', 'flux', 'flux_err', 'zp', 'zpsys','mag','mag_err'))
        return data_draft


    def getLasairObInfo(self,idobject,save=True,force=True,path=""):

            path_file=base_data_path+"lightcurve/"+idobject+".json"
            if not os.path.exists(path_file) or force:
                r = requests.get(self.url_lasair + idobject + "/json")
                jsondata_file = r.json()
                if save:
                    with open(path_file, 'w') as json_file:
                        json.dump(jsondata_file,json_file)
                        json_file.close()
            else:

                # read file
                with open(path_file, 'r') as json_file:
                    data = json_file.read()


                    # parse file
                    jsondata_file = json.loads(data)
                    json_file.close()


            self.light_curve = jsondata_file
            return jsondata_file


if __name__ == "__main__":
    file_path="/Users/camiloj/projects/broker_data/tns_sndetections.csv"
    file_path_report = "/Users/camiloj/projects/broker_data/reports/"
    sn_classifier = SNCosmosFit()

    ztf_id="ZTF20aahbfmf"
    redshift_ob=0.123
    type_ob="SN Ia"
    tns_name="Fred"
    sn_classifier.fitSNByZTFId(ztf_id, save_on=file_path_report,
                               reference={"redshift": 0.123 , "type": 0, "label": tns_name, "type_ob": type_ob})
    # sn_classifier.loadTnsSn(file_path)
    # ztf_id="ZTF20aadvbni"
    # redshift_ob=0.106
    # type_ob="SN Ia"
    # tns_name="2020uo"
    # sn_classifier.fitSNByZTFId(ztf_id,save_on=file_path_report,reference={"redshift":0.106,"type":0,"label":tns_name,"type_ob":type_ob})
    # sn_classifier.generatePlots(file_path_report+"ZTF20aadvbni.csv",mark={"redshift":0.106,"type":0,"label":tns_name})
