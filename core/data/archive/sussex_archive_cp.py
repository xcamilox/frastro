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
#data from
#http://batc.bao.ac.cn/~zouhu/online-data/desi_photoz
class SussexArchive():
    def getSDSS(self,ra, dec, radii=10 * u.arcsec):
        # pos = coords.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
        # xid = SDSS.query_region(pos, spectro=True,data_release=15,photoobj_fields=["ra","dec","objid","run","rerun","camcol","field","u","g","r","i","z","err_u","err_g","err_r","err_i","err_z"],specobj_fields=["plate","fiberID","specObjID","z","zErr","zWarning"],radius=3*u.arcsec )

        # ra, dec are in degrees, r is in arc minutes.
        query = ''' SELECT TOP 1000
                                p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,
                                p.run, p.rerun, p.camcol, p.field,
                                f.z as photoz, f.zErr as photoZerr, GN.distance,
                                s.specobjid, s.class, s.z as redshift, s.zErr as errredshift, s.plate, s.mjd, s.fiberid
                                FROM PhotoObj AS p
                                JOIN dbo.fGetNearbyObjEq({0},{1}, {2}) AS GN on p.objid = GN.objid 
                                join Photoz as f ON f.objid = GN.objid 
                                LEFT OUTER JOIN SpecObj AS s ON s.bestobjid = GN.objid'''

        params = {"searchtool": 'SQL',
                  "TaskName": 'Skyserver.Search.SQL',
                  "syntax": 'NoSyntax',
                  "ReturnHtml": 'true',
                  "cmd": query.format(str(ra), str(dec), str(radii.to("arcmin").value)),

                  "format": 'json',
                  "TableName": ''
                  }
        url = "http://skyserver.sdss.org/dr12/en/tools/search/x_results.aspx"
        r = requests.get(url, params=params)
        return r.json()

    def getDESI(self,ra, dec,radio=5):
        ra = str(ra)
        dec = str(dec)
        url = "https://herschel-vos.phys.sussex.ac.uk/__system__/tap/run/tap"

        radio=radio*u.arcsec
        radio_degree = str(radio.to("degree").value)


        query = "SELECT * FROM desi_photoz.main WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', " + ra + ", " + dec + ", "+radio_degree+"))"

        service = vo.dal.TAPService(
            "https://herschel-vos.phys.sussex.ac.uk/__system__/tap/run/tap"
        )
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

        return match.to_table()




    def getDesiPhotoZfromTable(self,table):
        print("Get desi crossmatch")
        cross_query = """
        SELECT
            tc.id as id,
            db.ID as desiid,
            db.RA as desira,
            db.DEC as desidec,
            db.field,
            db.photo_z,
            db.photo_zerr,
            db.MASS_BEST,
            MASS_INF,
            MASS_SUP,
            SFR_BEST,
            SFR_INF,
            SFR_SUP,
            AGE_BEST,
            AGE_INF,
            AGE_SUP,
            tc.ramean,
            tc.decmean
        FROM desi_photoz.main AS db
        JOIN TAP_UPLOAD.t1 AS tc
        ON 1=CONTAINS(POINT('ICRS', db.ra, db.dec),
        CIRCLE('ICRS', tc.ramean, tc.decmean, 5.0/3600.))
        """
        #SQRT(POWER((db.RA-tc.ramean),2)+POWER((db.DEC-tc.decmean),2)) as distance
        # construct a service; I’ve taken the URL from TOPCAT’s
        # TAP service browser # ("Selected TAP Service" near the
        # foot of the dialog)
        service = vo.dal.TAPService(
            "https://herschel-vos.phys.sussex.ac.uk/__system__/tap/run/tap"
        )
        job = service.submit_job(query=cross_query,
                                 uploads={'t1': table})
        job.run()
        job_url = job.url
        job_result = vo.dal.tap.AsyncTAPJob(job_url)
        start_time = time.time()
        job_result.wait(phases=["COMPLETED"])
        desi_photozvo,desi_photoztable = self.getDesiResult(job_result)

        print('Job {} running after {} seconds with candidates'.format(job.phase, round(time.time() - start_time)))

        return desi_photozvo,desi_photoztable



    def getTNSdetectionOnDesiMassiveGalaxies(self,table):
        print("Get desi crossmatch")
        cross_query = """
        SELECT
            tc.id as id,
		    tc.Name,
		    tc.objtype,
		    tc.Redshift,
		    tc.DiscInternalName,
		    tc.HostRedshift,
            db.ID as desiid,
            db.RA as desira,
            db.DEC as desidec,
            db.field,
            db.photo_z,
            db.photo_zerr,
            db.MASS_BEST,
            db.MASS_INF,
            db.MASS_SUP,
            db.SFR_BEST,
            db.SFR_INF,
            db.SFR_SUP,
            db.AGE_BEST,
            db.AGE_INF,
            db.AGE_SUP,
            tc.ra,
            tc.dec
        FROM desi_photoz.main AS db
        JOIN TAP_UPLOAD.t1 AS tc
        ON 1=CONTAINS(POINT('ICRS', db.ra, db.dec),
        CIRCLE('ICRS', tc.ra, tc.dec, 5.0/3600.))
        """
        #SQRT(POWER((db.RA-tc.ramean),2)+POWER((db.DEC-tc.decmean),2)) as distance
        # construct a service; I’ve taken the URL from TOPCAT’s
        # TAP service browser # ("Selected TAP Service" near the
        # foot of the dialog)
        service = vo.dal.TAPService(
            "https://herschel-vos.phys.sussex.ac.uk/__system__/tap/run/tap"
        )
        job = service.submit_job(query=cross_query,
                                 uploads={'t1': table})
        job.run()
        job_url = job.url
        job_result = vo.dal.tap.AsyncTAPJob(job_url)
        start_time = time.time()
#        print("result",job_result.phases)

        job_result.wait(phases=["COMPLETED"])
        desi_photozvo,desi_photoztable = self.getDesiResult(job_result)

        print('Job {} running after {} seconds with candidates'.format(job.phase, round(time.time() - start_time)))

        return desi_photozvo,desi_photoztable

    def getDesiResult(self,job):

        try:
            result = job.fetch_result()
            desi_photozvo = result.votable
            desi_photoztable = result.table
        except AttributeError as err:
            print("error call again afeter 1 sec",err)
            time.sleep(1)
            desi_photoz = self.getDesiResult(job)


        return desi_photozvo,desi_photoztable



    def getDESITargets(self, ra, dec, radio=10):

        ra = ra*u.deg
        dec = dec*u.deg
        radio = radio*u.arcsec

        ra_min = ra-radio
        ra_max = ra + radio

        dec_min = dec - radio
        dec_max = dec + radio
        #22.4173&rahi=22.4707&declo=19.7097&dechi=19.7239
        url = "http://legacysurvey.org/viewer/targets-dr8/1/cat.json"
        params={"ralo":ra_min.value,"rahi":ra_max.value,"declo":dec_min.value,"dechi":dec_max.value}
        # params = {"ralo": 22.4173, "rahi": 22.4707, "declo": 19.7097, "dechi": 19.7239}
        r = requests.get(url, params=params)
        result=r.json()
        radec = result["rd"]
        positions = np.array(radec)
        catalog = SkyCoord(positions[:,0], positions[:,1], frame='icrs', unit='deg')
        cref = SkyCoord(ra, dec, frame='icrs', unit='deg')



        desi_distance_arc = cref.separation(catalog).arcsec
        indx = np.where(desi_distance_arc == desi_distance_arc.min())[0][0]

        #result["mag_ab"] = self.fluxToMag(result["fluxes"])
        result["distance"] = desi_distance_arc
        result["closed_idx"] = indx


        return result

    def fluxToMag(self,desiTargetFlux):
        #desi flux in maggies 1 x 10^-9
        desiTargetFlux = np.array(desiTargetFlux)
        return -2.5 * (np.log(desiTargetFlux)/np.log(10.) - 9)


if __name__ == "__main__":
    ra = 152.13029797142855
    dec = 9.239761214285716
    dataArchive = SussexArchive()
    r = dataArchive.getDESI(ra,dec,5)
    print(len(r),r)
    print("ok")