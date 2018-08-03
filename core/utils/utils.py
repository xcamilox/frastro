from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from astropy import utils, io
import os
import urllib.request
import requests
import bz2
import gzip
from astropy import wcs

import subprocess
import os
from subprocess import Popen,PIPE,STDOUT,call
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
            index = path.rindex("/")
            main_path = path[0:index]
            paths=main_path.split("/")
            check_path="/"
            for path_index in paths:
                check_path+=path_index
                if not os.path.exists(check_path):
                    os.mkdir(check_path)
                if len(check_path)>1:
                    check_path += "/"





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
        for item in data:
            url=item.get_attribute("href")

            position=url.find('.fits')
            band=url[position-1:position]
            list_img.append({"fits":url,"band":band,'image':ImageUtils.getBase64FromFitsURL(url)})
        browser.quit()
        return data