import os
import os.path as path

from astropy.io import fits
from astropy.table import Table

from frastro import MongodbManager, Config, Utils
from frastro.core.utils.report_util import ReportUtil
from frastro.core.utils.image_util import ImageUtils
import matplotlib.pyplot as plt
from pylab import figure, cm
from matplotlib.gridspec import GridSpec


class Generatereports():
    @staticmethod
    def galfitReport(filter={'epoach': {'$exists': 'true'}}):
        mgdb = MongodbManager()
        report = ReportUtil()
        mgdb.setCollection("liverpool_sn_v3")
        all_data = mgdb.getData(filter)
        report_output=""
        for target in all_data:
            try:
                html_text = ""
                ra = target["ra"]
                dec = target["dec"]
                html_text+=report.addTitle(target["id"])
                html_text_content=""
                for epoach in target["epoach"]:
                    stack = target["epoach"][epoach]["stacking"]["best"]
                    if "model" in stack:

                        model = stack["model"]
                        psf = stack["psf_path"]
                        model_path = os.path.dirname(model)
                        model_data = fits.open(model)
                        model_image = model_data[1].data
                        mean_img = model_image.mean()
                        std_img = model_image.std()

                        model_model = model_data[2].data
                        mean_model = model_model.mean()
                        std_model = model_model.std()

                        model_residuos = model_data[3].data
                        mean_residuos = model_residuos.mean()
                        std_residuos = model_residuos.std()

                        fig = plt.figure()

                        # show reference

                        gs1 = GridSpec(1, 3)
                        gs1.update(wspace=0.05)
                        ax1 = plt.subplot(1, 3, 1)
                        plt.imshow(model_image, cmap=cm.Oranges, vmin=mean_img - std_img, vmax=mean_img + std_img,
                                   origin="lower")
                        plt.axis('off')

                        ax2 = plt.subplot(1, 3, 2)
                        plt.subplots_adjust(wspace=0.5)
                        plt.imshow(model_model, cmap=cm.Oranges, vmin=mean_model - std_model, vmax=mean_model + std_model,
                                   origin="lower")
                        plt.axis('off')

                        ax3 = plt.subplot(1, 3, 3)
                        plt.subplots_adjust(wspace=0.01)
                        plt.imshow(model_residuos, cmap=cm.Oranges, vmin=mean_residuos - std_residuos, vmax=mean_residuos + std_residuos,
                                   origin="lower")

                        plt.axis('off')
                        image_output = model_path + "/model.png"
                        plt.savefig(image_output)
                        html_text+=report.addSubTitle(epoach)
                        html_text+=report.addImage(image_output,width="1500")
                        html_text+=report.addText(model)
                        html_text += report.addText(psf)
                        html_text_content+=report.wrapContent(html_text)
                        html_text=""
                        fig.clf()

                path = os.path.dirname(target["stacking"]["best"]["file_path"])
                report.saveReport(path+"/model_report.html",html_text_content)
                report_output+=html_text_content
            except:
                print("Error",target["id"])

        path = Config.getPath("liverpool_sn")
        report.saveReport(path + "/models_report.html", report_output)
        report.createPDF(path + "/models_report.pdf", report_output)

    @staticmethod
    def dataReports(filter={'epoach': {'$exists': 'true'}}):
        mgdb = MongodbManager()
        report = ReportUtil()
        mgdb.setCollection("liverpool_sn_v3")
        all_data = mgdb.getData(filter)
        path = Config.getPath("liverpool_sn")
        report_output=""
        table = report.createTable(["id","epoach","group_id","mjd","airmass","seeing","exptime","source"])
        report_out=""
        img_paths=path+"reportimg/"

        ids=[]
        ras=[]
        decs = []
        group_ids=[]
        mjds=[]
        airmasss=[]
        seeings=[]
        exptimes=[]
        epoachs=[]
        for target in all_data:
            try:
                html_text = ""
                id = target["id"]
                ra = target["ra"]
                dec = target["dec"]
                for epoach in target["epoach"]:
                    if "stacking" in target["epoach"][epoach]:
                        if "error" not in target["epoach"][epoach]["stacking"]:
                            group_id = target["epoach"][epoach]["obs"][0]["group_id"]
                            mjd = target["epoach"][epoach]["obs"][0]["mjd"]
                            airmass = target["epoach"][epoach]["obs"][0]["airmass"]
                            seeing = target["epoach"][epoach]["stacking"]["best"]["seeing_cal"]
                            exptime = target["epoach"][epoach]["stacking"]["best"]["exp_time"]
                            source= target["epoach"][epoach]["stacking"]["best"]["cat_deep_maglimit"]["n_sources"]
                            path = target["epoach"][epoach]["stacking"]["all"]["file_path"]

                            # image = ImageUtils.getImageThumbnail(path,ra,dec,0,(500,500))
                            # vmin= image.mean()-image.std()
                            # vmax = image.mean() + image.std()
                            #
                            # plt.imshow(image, cmap=cm.Oranges, vmin=vmin,vmax=vmax,origin="lower")
                            # plt.axis('off')
                            # image_output = img_paths + id+"_"+epoach+".png"
                            # plt.savefig(image_output)

                            #print(id,epoach,group_id,mjd,airmass,seeing)
                            report_out+=report.addToTable([id,epoach,group_id,mjd,airmass,str(seeing)[0:5],exptime,source] ) #,img=image_output)

                        else:
                            group_id = target["epoach"][epoach]["obs"][0]["group_id"]
                            mjd = target["epoach"][epoach]["obs"][0]["mjd"]
                            airmass = target["epoach"][epoach]["obs"][0]["airmass"]
                            seeing = target["epoach"][epoach]["obs"][0]["seeing_cal"]
                            exptime = target["epoach"][epoach]["obs"][0]["exptime"] * len(target["epoach"][epoach]["obs"])
                            source = 0
                            report_out += report.addToTable([id, epoach, group_id, mjd, airmass, str(seeing)[0:5], exptime, source])

                        ids.append(id)
                        ras.append(ra)
                        decs.append(dec)
                        epoachs.append(epoach)
                        group_ids.append(group_id)
                        mjds.append(mjd)
                        airmasss.append(airmass)
                        seeings.append(seeing)
                        exptimes.append(exptime)

                    else:
                        print("not stacking",id,epoach)

            except Exception as inst:
                print("error",id,epoach)
                print(type(inst))  # the exception instance
                print(inst.args)  # arguments stored in .args
                print(inst)

        path = Config.getPath("liverpool_sn")
        datatable=Table([ids,ras,decs,epoachs,group_ids,mjds,airmasss,seeings,exptimes],names=("id","ra","dec","epoch","group_id","mjd","airmass","seeing","exptime"))
        datatable.write(path + "data_report_list.fits", format='fits')

        #table_out=report.saveTable(table+report_out)
        #report.saveReport(path + "data_report_list.html", table_out)
        #report.createPDF(path + "data_report_list.pdf", table_out)



    @staticmethod
    def checkDataFolder(folder_path):
        paths = []
        onlyfiles = []
        for folder in os.listdir(folder_path):

            if path.isdir(folder_path+folder):

                files = [folder_path+folder + "/" +f for f in os.listdir(folder_path+folder) if path.isfile(folder_path+folder + "/" + f)]
                onlyfiles.extend(files)

            elif path.isfile(folder_path+folder):
                 onlyfiles.append(folder_path+folder)

        for file in onlyfiles:
            if path.splitext(file)[1] == ".fits" or path.splitext(file)[1] == ".fit":
                paths.append(file)
        return paths

if __name__ == "__main__":
    Generatereports.galfitReport()
    # Generatereports.dataReports()
    # path_files="/Users/cjimenez/Documents/PHD/data/liverpool_lens/observations/"
    # obs="/Users/cjimenez/Documents/PHD/data/liverpool_lens/alldatareduction"
    # fileslist=Generatereports.checkDataFolder(path_files)
    # filesobs=Utils.getFilesList(obs)
    # for files in fileslist:
    #     if files[files.rfind("/")+1:] not in filesobs:
    #         print(files)
    #         fi=fits.open(files)
    #         header=fi[0].header
    #         #print(header["DATE"])
    #         #print(header["OBJECT"])
