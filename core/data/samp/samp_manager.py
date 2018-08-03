import time

from astropy.samp import SAMPIntegratedClient
from astropy.samp import SAMPHubServer


class SAMPManager():

    __client=None
    __received = False
    __params = {}
    __hub =None
    __server_name = "frastro server"
    __server_description = "frastro Web client"
    __server_description_html = "<h1>frastro Web client</h1>"
    __author = "Camilo E. Jimenez-Angel"
    __instritution = "IAC"
    __email_contact = "camilo.jimenez@iac.es"
    __url_icon = "http://localhost:8000/static/img/frastro_icon.jpg"


    def __init__(self,server_name=""):
        # Instantiate the client and connect to the hub
        self.__server_name = server_name if server_name != "" else self.__server_name
        self.connect()

    def starHubServer(self,web_profile = True):
        if self.__hub is None:
            self.__hub = SAMPHubServer(web_profile=web_profile,label=self.__server_name)
        self.__hub.start()

    def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
        self.__params = params
        self.__received = True
        self.__client.reply(msg_id, {"samp.status": "samp.ok", "samp.result": {}})

    def receive_notification(self, private_key, sender_id, mtype, params, extra):
        self.__params = params
        self.__received = True


    def bind_to_server(self):
        # Listen for any instructions to load a table
        self.__client.bind_receive_call("table.load.votable", self.receive_call)
        self.__client.bind_receive_notification("table.load.votable", self.receive_notification)

    def sampMetadata(self):
        meta={
            "samp.name":self.__server_name,
            "samp.description.text":self.__server_description,
            "samp.description.html":self.__server_description_html,
            "samp.icon.url":self.__url_icon,
            "author.affiliation":self.__instritution,
            "author.name":self.__author,
            "author.email":self.__email_contact
        }
        return meta

    def getRegisteredClients(self):
        clients=self.__client.get_registered_clients()
        run_clients = {}
        for client in clients:
            client_data = self.getMetadata(client)
            run_clients[client] = client_data
        return run_clients

    def getMetadata(self,client_id):
        return self.__client.get_metadata(client_id)

    def sendTable(self,message,recipient="all"):
        if recipient == "all":
            self.__client.notify_all(message)
        else:
            self.__client.notify(recipient,message)



    def connect(self):
        self.__client = SAMPIntegratedClient(name=self.__server_name,metadata=self.sampMetadata())
        self.__client.connect()
        self.bind_to_server()

    def sendImage(self,params,id=""):
        message = {}
        message["samp.mtype"] = "image.load.fits"
        message["samp.params"] = params
        if id=="" or id=="all":
            self.__client.notify_all(message)
        else:
            self.__client.notify(recipient_id=id,message=message)

    def sendMessage(self,params,id=""):
        message = {}
        message["samp.mtype"] = "table.load.votable"
        message["samp.params"] = params
        if id=="" or id=="all":
            self.__client.notify_all(message)
        else:
            self.__client.notify(recipient_id=id,message=message)



    def disconnect(self):
        self.__client.disconnect()

    def __del__(self):
        self.disconnect()

