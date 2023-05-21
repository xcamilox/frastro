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

    def getDatabase(self,database=""):
        db = database if database != "" else self.__main_database
        return self.getDB()[db]

    def setDatabase(self,db):
        if type(db) == str and db !="":
            self.__main_database=db

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
            collection = self.__main_collection
        id = self.getCollection(collection).insert_one(data)
        return id.inserted_id

    def updateAstrosource(self,astrosource):
        source = astrosource
        collection = self.getCollection()
        return collection.update_one({"id":source['id']},{ "$addToSet": {"archives": source['archives']}})

    def delete(self,query,collection=""):
        if collection!="":
            self.setCollection(collection)
            collection = self.getCollection()
        else:
            collection = self.getCollection()

        return collection.delete_many(filter=query)

    def update(self,filter,query,collection=""):
        if collection!="":
            self.setCollection(collection)
            collection = self.getCollection()
        else:
            collection = self.getCollection()
        if type(filter) == dict:
            return collection.update_one(filter, query)
        else:
            return collection.update_one({"id": filter}, query)


    def updateMany(self,filter,query,collection=""):
        if collection!="":
            self.setCollection(collection)
            collection = self.getCollection()
        else:
            collection = self.getCollection()
        if type(filter) == dict:
            return collection.update_many(filter, query)
        else:
            return collection.update_many({"id": filter}, query)


    def getAstrosourceByID(self,id):
        return self.getCollection().find_one({"id":id})


    def query(self,filter=None,projection=None,collection=""):
        if collection !="":
            self.setCollection(collection)
        r=self.getCollection().find(filter, projection)

        return r

    def getData(self,filter={},projection={},skip=0):
        result = []
        if filter!={} and projection !={}:
            r=self.getCollection().find(filter, projection)
        elif filter!={}:
            r = self.getCollection().find(filter)
        elif projection!={}:
            r = self.getCollection().find(projection=projection)
        else:
            r=self.getCollection().find()

        r.skip(skip)
        for obj in r:
            obj.pop("_id")
            result.append(obj)

        return result

    def getAstroSources(self,filter=None,projection={}):
        return self.getData(filter=filter,projection=projection)


    def command(self,collection="",pipeline={}):
        if collection !="":
               self.setCollection(collection)
        r = self.getCollection().aggregate(pipeline=pipeline)
        result = []
        for obj in r:
            try:
                obj.pop("_id")
                result.append(obj)
            except Exception as err:
                print(err)
        return result
