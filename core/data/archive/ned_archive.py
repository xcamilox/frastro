import requests
import json
import pyvo as vo
from astropy import coordinates as coords
from astropy import units as u
import math
from astropy.coordinates import SkyCoord
import numpy as np
import time
from astropy.table import QTable
import pandas as pd
from astropy.table import Table, join, vstack, QTable

class NEDArchive():

    @staticmethod
    def getCatalog(ra, dec,radio=5):
        ra = str(ra)
        dec = str(dec)
        url = "https://ned.ipac.caltech.edu/tap/"

        radio=radio*u.arcsec
        radio_degree = str(radio.to("degree").value)


        query = "SELECT * FROM NEDTAP.objdir WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', " + ra + ", " + dec + ", "+radio_degree+"))"

        service = vo.dal.TAPService(url)
        match=service.search(query)

        # job = service.submit_job(query=query)
        # job.run()
        # job_url = job.url
        # job_result = vo.dal.tap.AsyncTAPJob(job_url)
        # start_time = time.time()
        #
        # job_result.wait(phases=["COMPLETED","ABORTED", "ERROR"])
        # # wait = 5.
        # # while job.phase == 'EXECUTING':
        # #     print(job.phase)
        # #     time.sleep(wait)  # wait and try again
        # #     # wait *= 2
        #
        # desi_photozvo, desi_photoztable = self.getDesiResult(job_result)

        # print('Job {} running after {} seconds with candidates{}'.format(job.phase, round(time.time() - start_time)),len(desi_photoztable))
        table_result = match.to_table()


        table_result["prefname"] = table_result["prefname"].astype(str)

        ra_ref = table_result["ra"].tolist()
        dec_ref = table_result["dec"].tolist()
        cref = SkyCoord([ra], [dec], frame='icrs', unit='deg')
        c1 = SkyCoord(ra_ref, dec_ref, frame='icrs', unit='deg')
        desi_distance = cref.separation(c1).arcsec
        table_result["separation"] = desi_distance
        pd_result = Table(table_result, masked=False).to_pandas()

        return pd_result







