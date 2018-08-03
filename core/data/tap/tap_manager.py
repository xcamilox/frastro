from astroquery.utils.tap.core import TapPlus


class TAPManager():
    __tap_connector=None

    def __init__(self):
        pass


    def connect(self,url=""):
        if self.__tap_connector==None:
            if url=="":
                raise Exception("provide URL TAP Server")
            else:
                self.__tap_connector = TapPlus(url=url)
        else:
            if url != "":
                self.__tap_connector = TapPlus(url=url)
            else:
                return self.__tap_connector

    def getAvalibleTables(self):
        conector=self.connect()
        tables=conector.load_tables()

        return tables

    """
    There is a limit of 2000 rows. If you need more than that, you must use asynchronous queries method (async_query).        
    """

    def sync_query(self,query,dump_to_file=False):
        if query=="" or type(query) is not str:
            raise Exception("you need a query ADQL format")
        conector = self.connect()
        try:
            job=conector.launch_job(query,dump_to_file=dump_to_file)
            print("rows: ",len(job.get_results()))
        except ValueError as e:
            print(e)
            job = JOBError()
        return job

    def async_query(self,query,dump_to_file=False):
        if query=="" or type(query) is not str:
            raise Exception("you need a query ADQL format")
        conector = self.connect()
        job = conector.launch_job_async(query,dump_to_file=dump_to_file)
        return job



class JOBError():
    def get_results(self):
        print("Tap query request error service")
        return []