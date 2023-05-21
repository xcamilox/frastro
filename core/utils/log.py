from core.utils.config import Config
from core.utils.singleton import Singleton
import os
from datetime import datetime


class LOGUtil(Singleton):
    log_path = "./"
    log_error = []
    logs=[]
    verbose = True
    file = None
    current_path=""

    def openLog(self,log_path=""):

        if log_path=="":
            log_path = Config.getPath("log")
        else:
            log_path = os.path.dirname(log_path)
            log_path = log_path + "/log.log"

        try:

            if os.path.exists(log_path):
                self.file = open(log_path, "a")
            else:
                output_path = os.path.dirname(log_path)
                os.makedirs(output_path)
                self.file = open(log_path, "w+")
        except PermissionError:
            print("error print log")

    @staticmethod
    def log(str,path="",force_close=True):

        """
        lo=LOGUtil()
        try:
            lo.openLog(path)
            now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            lo.file.write(now+": "+str + "\n")
            if lo.verbose:
                print(str)
                lo.logs.append(str)
            lo.file.close()
        except ValueError:
            print("error generate log")
        """
        print(str)


