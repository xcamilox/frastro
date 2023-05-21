from antares_client.search import search
from antares_client.search import download


class AntaresArchive():
    name = 'Alerce'
    ztfid=""

    query_single_item={
          "query": {
            "wildcard": {
              "properties.ztf_object_id.keyword": ztfid
            }
        }
    }

    def getItem(self,ztfId):
        self.ztfid = ztfId
        self.fetch_alerts(self.query_single_item)


    @classmethod
    def fetch_alerts(self, query):
        result_set = search(query,progress_callback=self.progress_callback)
        result_set = download(query, "results.json.gz", output_format="json", decompress=True)
        return result_set

    @classmethod
    def progress_callback(self,status, **kwargs):
        print("Search status: {}; Received kwargs: {}".format(status, kwargs))

