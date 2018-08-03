from frastro import AstroSource
from astropy.time import Time
from frastro import MongodbManager


class LiverPoolObsservation():
    __target = AstroSource
    __group = ""
    __texp = 0.0
    __seeing= 0.0
    __date=0
    __mjd=0
    __collection="liverpool_sn"
    __mgdb = MongodbManager()


    def __init__(self,**kargs):
        if "source" in kargs:
            self.__target=kargs["source"]

        if "group" in kargs:
            self.__group=kargs["group"]

        if "texp" in kargs:
            self.__texp==kargs["texp"]

        if "seeing" in kargs:
            self.__seeing==kargs["seeing"]

        if "date" in kargs:
            self.__date = kargs["date"]
            date = Time(self.__date, format='isot', scale='utc')
            self.__mjd = date.mjd

        self.__mgdb.setCollection(self.__collection)

    def getObservations(self,id=None):
        if id==None:
            list = self.__mgdb.getData()

        return list




    def saveObservation(self):

        if self.__collection != "":

            self.__mgdb.saveData(self.__dict__())


    def __dict__(self):
        return {"target":"","group":self.__group,"texp":self.__texp,"seeing":self.__seeing,"date":self.__date,"mjd":self.__mjd}







