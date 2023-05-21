import subprocess
from frastro import Utils, Config, LOGUtil, ImageUtils
from frastro.external.psfex_util import PSFexUtil
import os
import os.path as path
import numpy as np
import shutil
from astropy.io import fits



class SwarpUtil():

    RUN_SWARP_CM = Config.getExtApp("swarp")
    DEFAULT_PARAMS_CM= [RUN_SWARP_CM,"-d"]


    def __init__(self):
       pass

    def createDefaultParamsfiles(self,path):

        default=self.runCommand(self.DEFAULT_PARAMS_CM)
        f = open(path + "default.swarp", "w+")
        f.write(default.decode("utf-8"))
        f.close()
        return path + "default.swarp"

    def runCommand(self,command):
       r = subprocess.run(command,stdout=subprocess.PIPE)
       return r.stdout

    def createFileListFromBestPSF(self,files_path,output_path,file_type=".fits",max_error=0.3):
        LOGUtil.log("============Create file List from best files to Swarp===========", output_path)

        #Utils.createFileList(files_path, output_file=output_path + "files.lst")
        folder_path = files_path if files_path[-1] == "/" else files_path + "/"
        only_fits_files = [folder_path + f for f in os.listdir(folder_path) if
                           path.isfile(folder_path + f) and path.splitext(f)[1] == file_type]

        psf = PSFexUtil()
        seeing_lst=[]
        psf_items=[]
        for file in only_fits_files:

            fwhm = psf.calcSeeing(file)
            seeing_lst.append(float(fwhm[0]))
            psf_items.append(fwhm[1])

        seeing_lst = np.array(seeing_lst,dtype=float)
        min=seeing_lst.min()

        if min < 999:
            #bad files observations

            index_best=np.where(seeing_lst==min)[0][0]

            LOGUtil.log("Best Seeing: "+str(min), output_path)
            LOGUtil.log("Best Seeing file: " +  only_fits_files[index_best], output_path)
            LOGUtil.log("Best Seeing PSF file: " + psf_items[index_best], output_path)
            LOGUtil.log("Good Files: ", output_path)

            std=max_error
            index=np.where(seeing_lst< min+(min*std),True,False)
            good_files=np.array(only_fits_files)[index]

            str_log="=========== Files with seeing < "+str(min+(min*std))+"\n"
            for idx,value in enumerate(only_fits_files):
                state= "accepted" if index[idx] else "rejected"
                str_log += "seeing: "+str(seeing_lst[idx]) + "  " +state +  " File: " +only_fits_files[idx] + " PSFfile: "+ psf_items[idx]+"\n"

            LOGUtil.log(str_log, output_path)

            with open(output_path, 'w+') as f:
                for file in good_files:
                    data=str(file)+"\n"
                    f.write(data)
                f.close()
                print("Save file:", output_path)

            print(good_files)
            return output_path
        else:
            return None



    def run(self,file_list,output_path,**kargs):
        if not path.exists(output_path):
            os.makedirs(output_path)

        commands = []
        current_path = self.runCommand(["pwd"])
        commands.append(self.RUN_SWARP_CM)

        os.chdir(output_path)
        new_path_list = output_path+"files.lst"
        shutil.copy2(file_list,new_path_list)

        commands.append("@files.lst")
        commands.append("-c")
        default_path = self.createDefaultParamsfiles(output_path)
        commands.append(default_path)
        kargs["-COPY_KEYWORDS"] = "INSTRUME,TELESCOP,OBSERVAT,PIXSCALE,FILTER,FILTER1,OBJECT"
        for arg in kargs:
            if arg != "ra" and arg != "dec" and arg != "size":
                commands.append(arg)
                commands.append(kargs[arg])

        #kargs["-RESAMPLING_TYPE"] = "NEAREST"  # No iterpolate
        #kargs["-IMAGEOUT_NAME"] = "coadd_nearest_sapling.fits"  # No iterpolate
        #kargs["-WEIGHTOUT_NAME"] = "weightout_coadd_nearest_sapling.fits"  # No iterpolate

        r = self.runCommand(commands)
        file_name ="coadd.fits"
        if "-IMAGEOUT_NAME" in kargs:
            file_name=kargs["-IMAGEOUT_NAME"]
        print("swarpt ")
        print(" ".join(commands))

        hdl=fits.open(output_path+file_name)
        data=hdl[0].data
        shape = data.shape
        outpath=output_path+"crop_"+file_name
        coordinate = (shape[0]/2,shape[1]/2)
        size = (shape[0] - 150, shape[1] - 150)
        if "ra" in kargs and "dec" in kargs:
            coordinate=ImageUtils.getPixelFromCoordinate(output_path + file_name,kargs["ra"],kargs["dec"])[0]

        if "size" in kargs:
            size = kargs["size"]




        ImageUtils.imcopy(output_path+file_name, outpath, coordinate, size, addcounts=0)

        return outpath

    #run swarp from a files path. get all fits files from a path
    #if quality is best, will be calculate the psf for each file and take the best fwhm and reject file with over hwfm+(hwfm*0.3)
    #this value could change by karg[max_error]
    def runSwarp(self,files_path,quality="all",**kargs):

        file_path = ""
        if not os.path.isdir(files_path):
            file_path=files_path
            files_path=os.path.dirname(file_path)

        LOGUtil.log("============Run swarp===========", files_path)
        commands=[]
        current_path = self.runCommand(["pwd"])

        output_path = files_path

        files_path = files_path if files_path[-1] == "/" else files_path + "/"

        if not path.exists(files_path+"reduction/swarp"):
            os.makedirs(files_path+"reduction/swarp")

        if "/reduction/swarp" not in files_path:
            output_path=files_path+"reduction/swarp/"



        commands.append(self.RUN_SWARP_CM)

        if file_path!="":
            filelist = file_path
        else:
            if quality == "best":
                file_list_path=self.createFileListFromBestPSF(files_path,output_path+"best_files.lst")
                if file_list_path == None:
                    return None
                filelist="@best_files.lst"
            else:
                LOGUtil.log("============Create file List for all files to Swarp===========", output_path)
                Utils.createFileList(files_path,output_file=output_path+"files_all.lst")
                filelist = "@files_all.lst"
        os.chdir(output_path)
        print(os.getcwd())
        commands.append(filelist)
        commands.append("-c")
        default_path = self.createDefaultParamsfiles(output_path)
        commands.append(default_path)
        kargs["-COPY_KEYWORDS"]="INSTRUME,TELESCOP,OBSERVAT,PIXSCALE,FILTER,FILTER1,LIRF2NAM,OBJECT,AIRMASS"
        for arg in kargs:
            commands.append(arg)
            commands.append(kargs[arg])




        LOGUtil.log("========Run Swarp Command===========")
        LOGUtil.log(" ".join(commands), output_path)
        r = self.runCommand(commands)
        LOGUtil.log(str(r), output_path)

        file_name = "coadd.fits"
        if "-IMAGEOUT_NAME" in kargs:
            file_name = kargs["-IMAGEOUT_NAME"]

        kargs["-RESAMPLING_TYPE"]="NEAREST" #No iterpolate
        kargs["-IMAGEOUT_NAME"] = file_name[0:-5]+"_nearst.fits"  # No iterpolate
        kargs["-WEIGHTOUT_NAME"] = "weight_coadd_nearest_sapling.fits"  # No iterpolate



        r = self.runCommand(commands)

        return output_path+file_name

        #os.chdir(current_path.decode("utf-8").replace("\n",""))



if __name__ == "__main__":
    sw=SwarpUtil()
    path_files="/Users/cjimenez/Documents/PHD/data/unlens_liris/S1419+4341/ks/"
    kargs={"-IMAGEOUT_NAME":"S1419+4341_ks_swarp.fits"}
    kargs = {}
    sw.runSwarp(path_files,**kargs)