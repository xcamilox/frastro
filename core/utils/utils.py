from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from astropy import utils, io
import urllib.request
import requests
import bz2
import gzip
from astropy import wcs
import os.path as path
import subprocess
import os
from subprocess import Popen,PIPE,STDOUT,call

from core.utils.image_util import ImageUtils


class Utils():
    @staticmethod
    def validURL(str_url):


        state=None
        try:
            req = requests.get(str_url)
            if req.status_code<=200:
                state = True
            else:
                state = False
        except HTTPError as e:
            print('The server couldn\'t fulfill the request.')
            print('Error code: ', e.code)

        except URLError as e:
            print('We failed to reach a server.')
            print('Reason: ', e.reason)
        except ConnectionError as e:
            print('The server couldn\'t fulfill the request.')
        else:
            if state == None:
                state = True
        return state

    @staticmethod
    def saveLocalFile(url,output_path,uncompress=False,external=True):
        print("downloading ",url)

        #used an external process, create a parallel process unlinked
        if external:

            Utils.createPath(output_path)
            path = '"{1}"'
            index = output_path.rindex("/")
            file = output_path[index + 1:-1]
            path = path.format(output_path[0:index], url, file)
            cmd = ["cd", output_path[0:index] + ";", "curl", path.format(url), "-o", file] #
            #subprocess.run(cmd,shell=True)
            cmd="cd {0}; curl '{1}' -o {2} &"
            cmd=cmd.format(output_path[0:index],url,file)
            p = subprocess.Popen(cmd, stdout=None, stderr=None, shell=True)
            #output = p.communicate()
            #print(output)

        else:
            r=requests.get(url, verify=False,timeout=60)
            data=r.content
            if uncompress:
                data=Utils.unzipDataFiles(r.content)

            if data is not None:
                Utils.createPath(output_path)

                with open(output_path,'wb') as f:
                    f.write(data)
                    f.close()
                    print("Save file:",output_path)

                return output_path
        return None

    @staticmethod
    def createPath(path):
        #if not os.path.isdir(path):
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    @staticmethod
    def runCommand(command):
        print("commnad to run:"," ".join(command))
        try:
            r = subprocess.run(command, stdout=subprocess.PIPE)
            return r.stdout
        except:
            p = subprocess.Popen(command, stdout=None, stderr=None, shell=True)
            output = p.communicate()
            return output




    @staticmethod
    def createFileList(folder_path,output_file="files.lst",file_type=".fits"):
        folder_path = folder_path if folder_path[-1]=="/" else folder_path+"/"
        only_fits_files = [folder_path + f for f in os.listdir(folder_path) if
                           path.isfile(folder_path + f) and path.splitext(f)[1] == file_type]

        if path.exists(output_file):
            os.remove(output_file)

        f = open(output_file, "w+")
        for file in only_fits_files:
            f.write(file+"\n")
        f.close()

        return folder_path+output_file

    @staticmethod
    def unzipDataFiles(data):
        try:
            bites=bz2.BZ2Decompressor().decompress(data)
        except OSError:
            print("no bz files")
            bites=Utils.gzipDataFile(data)
        return bites

    @staticmethod
    def gzipDataFile(data):
        try:
            bites = gzip.decompress(data)
        except OSError:
            print("no gz files")

        return bites


    @staticmethod
    def getPixelPositionFromRaDec(ra,dec,filepath):

        header = io.fits.getheader(filepath)

        # header= io.fits.
        w = wcs.WCS(header)
        try:
            pixcrd2 = w.wcs_world2pix([[ra, dec]], 1)
        except ValueError:
            return ""
        return pixcrd2[0]


    @staticmethod
    def delayParserRequest(link,xpath_filter="a",xpath_condition=""):
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        browser = webdriver.Chrome('/Applications/chromedriver')
        browser.set_window_size(1120, 550)
        browser.get(link)
        filter=(By.XPATH, xpath_condition)
        element = WebDriverWait(browser, 1).until(
            EC.presence_of_element_located(filter)
        )

        data = element.find_elements_by_xpath(xpath_filter)
        list_img=[]
        for item in data:
            url=item.get_attribute("href")

            position=url.find('.fits')
            band=url[position-1:position]
            list_img.append({"fits":url,"band":band,'image':ImageUtils.getBase64FromFitsURL(url)})
        browser.quit()
        return data

    @staticmethod
    def getFilesList(folder_path,file_extension=""):
        onlyfiles = [f for f in os.listdir(folder_path) if path.isfile(folder_path + "/" + f)]
        files=[]
        if file_extension != "":
            for file in onlyfiles:
                if path.splitext(file)[1] == file_extension:
                    files.append(file)
            return files
        else:
            return onlyfiles


if __name__ == "__main__":
    file_path="/Users/cjimenez/Documents/PHD/data/liverpool_lens/lensed_sn_hst"
    files=Utils.getFilesList(file_path)
    print(files)