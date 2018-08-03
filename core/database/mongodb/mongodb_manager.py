from pymongo import MongoClient
class MongodbManager():
    __client=None
    __host = 'localhost'
    __port = 27017
    __main_database = "frastro"
    __main_collection = "v3"
    __db=None

    def __init__(self):
        self.getClient()

    def getClient(self):
        if self.__client == None:
            self.connect()
        return  self.__client

    def connect(self):
        self.__client = MongoClient(self.__host,self.__port)

    def getDB(self,id=None):
        if self.__db == None:
            client = self.getClient()
            if id is not None:
                self.__db = client[id]
            else:
                self.__db = client[self.__main_database]
        return self.__db

    def getCollection(self,collection=""):
        collection = collection if collection != "" else self.__main_collection
        return self.getDB()[collection]

    def setCollection(self,collection):
        if type(collection) == str and collection !="":
            self.__main_collection=collection

    def insertAstrosource(self,astrosource,collection=""):
        if astrosource != None:
            indb= self.getAstrosourceByID(astrosource["id"])
            if indb is None:
                self.getCollection(collection).insert_one(astrosource)
                return True
            else:
                return self.updateAstrosource(astrosource)

    def saveData(self,data,collection=""):
        if collection=="":
            collection= self.__main_collection
        self.getCollection(collection).insert_one(data)
        return True

    def updateAstrosource(self,astrosource):
        source = astrosource
        collection = self.getCollection()
        return collection.update_one({"id":source['id']},{ "$addToSet": {"archives": source['archives']}})

    def update(self,id,query):
        collection = self.getCollection()

        return collection.update_one({"id": id}, query)

    def getAstrosourceByID(self,id):
        return self.getCollection().find_one({"id":id})

    def query(self,filter=None,projection={}):
        r=self.getCollection().find(filter, projection)
        return r

    def getData(self,filter=None,projection={}):
        result = []
        if filter == None:
            for obj in self.getCollection().find():
                obj.pop("_id")
                result.append(obj)
        else:
            for obj in self.getCollection().find(filter, projection):
                obj.pop("_id")
                result.append(obj)
        return result

    def getAstroSources(self,filter=None,projection={}):
        return self.getData(filter=filter,projection=projection)

