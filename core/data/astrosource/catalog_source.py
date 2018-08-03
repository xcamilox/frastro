import json
from astropy.table import Table
from astroquery.utils import TableList
class CatalogSource(object):
    __name = ""
    __description = ""
    __summary={}
    __files = []
    __provider = ""
    __cache_date = False
    __results=[]

    def __init__(self,provider,description):
        self.__provider=provider
        self.__description=description


    def setCatalaog(self,result):
        if result != None:
            if type(result) is TableList:
                list=[]
                for table in TableList:
                    dic = table.to_pandas().to_json(orient='records')
                    dict_with_ints = dict((k, str(v)) for k, v in dic)
                    list.append(dict_with_ints)
                self.__results.append(list)
            elif type(result) is dict:
                self.__results.append(result)
            else:
                dic_result = json.loads(result.to_pandas().to_json(orient='records'))
                list_items = []
                for dic in dic_result:
                    dict_with_ints = dict((k, str(v)) for k, v in dic.items())
                    list_items.append(dict_with_ints)
                self.__results.append(list_items)

    def addFile(self,name,url,type="fits"):
        self.__files.append({"url": url, "name": name, "type": type})



    def __dict__(self):
        catalog = {"name":self.__name,"description":self.__description,"provider":self.__provider,"cache":self.__cache_date,"summary":self.__summary,"files":[]}
        for file in self.__files:
            catalog["files"].append(file)
        return catalog