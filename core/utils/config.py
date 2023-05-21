import json
from core.utils.singleton import Singleton
import pathlib

class Config(Singleton):

    def __init__(self,path="config.json"):

        filepath = pathlib.Path(pathlib.Path(pathlib.Path(__file__).parent).parent).parent
        self.path = str(filepath)+"/"+path
        self.__loadFile()

    def __loadFile(self):
        with open(self.path) as json_data_file:
            self.data = json.load(json_data_file)
        json_data_file.close()

    @staticmethod
    def getPath(str):
        self = Config()
        if str in self.data["paths"].keys():
            return self.data["paths"][str]
        else:
            return False

    @staticmethod
    def getExtApp(str):
        self = Config()
        if str in self.data["external"].keys():
            return self.data["external"][str]
        else:
            return False


    @staticmethod
    def getQuery(str):
        self = Config()
        if str in self.data["tapqueries"].keys():
            return self.data["tapqueries"][str]
        else:
            return False

    @staticmethod
    def getServices(str):
        self = Config()
        if str in self.data["services"].keys():
            return self.data["services"][str]
        else:
            return False

    @staticmethod
    def getCommand(str):
        self = Config()
        if str in self.data["commands"].keys():
            return self.data["commands"][str]
        else:
            return False

    @staticmethod
    def getDatabase(name):
        self = Config()
        if name in self.data["database"].keys():
            return self.data["database"][name]
        else:
            return False

