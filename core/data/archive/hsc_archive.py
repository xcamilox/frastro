import json
import urllib.request, urllib.error, urllib.parse
from astropy.coordinates import SkyCoord
import pandas as pd

class HSCArchive():
    user = "camiloja"
    passw = "dvhnWS8fz6/RMrIOJDr4ykxdr4nCqmbNVqv3Wj10"
    version = 20190514.1

    wide_query="""
-- Get object_id, ra, dec, CModel mags and mags_err on g,r,i,z,y bands
-- which are within 30 arcsec radius from (RA, Dec) = (150.93, 1.93) in degrees
-- and with good CModel magnitude measurement in g and r
-- skip ohter band due no in ztf

SELECT
        object_id
      , ra
      , dec
      , g_cmodel_mag
      , g_cmodel_magsigma
      , r_cmodel_mag
      , r_cmodel_magsigma
      , i_cmodel_mag
      , i_cmodel_magsigma
      , z_cmodel_mag
      , z_cmodel_magsigma
      , y_cmodel_mag
      , y_cmodel_magsigma
    , specz_name 
    , specz_redshift
    , specz_redshift_err
    , specz_flag_homogeneous
    , specz_original_flag
    , photoz_best
    , photoz_std_best
    , stellar_mass
    , sfr
    , prob_gal
    , prob_qso
    , prob_star
    FROM
  pdr2_wide.forced
    LEFT JOIN pdr2_wide.forced2 USING (object_id)
    LEFT JOIN pdr2_wide.forced3 USING (object_id)
    LEFT JOIN pdr2_wide.forced4 USING (object_id)
    LEFT JOIN pdr2_wide.forced5 USING (object_id)
    LEFT JOIN pdr2_wide.specz USING (object_id)
  LEFT JOIN pdr2_wide.photoz_mizuki USING (object_id)
    WHERE
    isprimary
    AND conesearch(coord, {0}, {1}, {2})
    AND NOT g_cmodel_flag
    AND NOT r_cmodel_flag
;
    """
    deep_query = """
-- Get object_id, ra, dec, CModel mags and mags_err on g,r,i,z,y bands
-- which are within 30 arcsec radius from (RA, Dec) = (150.93, 1.93) in degrees
-- and with good CModel magnitude measurement in r and i band
-- and have (r - i) >= 2.0 in CModel magnitudes.

SELECT
        object_id
      , ra
      , dec
      , g_cmodel_mag
      , g_cmodel_magsigma
      , r_cmodel_mag
      , r_cmodel_magsigma
      , i_cmodel_mag
      , i_cmodel_magsigma
      , z_cmodel_mag
      , z_cmodel_magsigma
      , y_cmodel_mag
      , y_cmodel_magsigma
      ,n387_cmodel_mag
      ,n816_cmodel_mag
      ,n921_cmodel_mag
      ,n387_cmodel_magsigma
      ,n816_cmodel_magsigma
    ,n921_cmodel_magsigma
    , specz_name
    , specz_redshift
    , specz_redshift_err
    , specz_flag_homogeneous
    , specz_original_flag
    , photoz_best
    , photoz_std_best
    , stellar_mass
    , sfr
    , prob_gal
    , prob_qso
    , prob_star
    FROM
        pdr2_dud.forced
    LEFT JOIN pdr2_dud.forced2 USING (object_id)
    LEFT JOIN pdr2_dud.forced3 USING (object_id)
    LEFT JOIN pdr2_dud.forced4 USING (object_id)
    LEFT JOIN pdr2_dud.forced5 USING (object_id)
    LEFT JOIN pdr2_dud.specz USING (object_id)
  LEFT JOIN pdr2_dud.photoz_mizuki USING (object_id)
    WHERE
    isprimary
    AND conesearch(coord, {0}, {1}, {2})
    AND NOT g_cmodel_flag
    AND NOT r_cmodel_flag
    
;
    """


    args = None

    def search(self,ra,dec,radio=5):
        query_str = self.deep_query.format(str(ra),str(dec),str(radio))
        credential = {'account_name': self.user, 'password': self.passw}
        results = self.preview(credential, query_str)
        catalog="deep"
        if len(results["rows"]) <= 0:
            query_str = self.wide_query.format(str(ra), str(dec), str(radio))
            results = self.preview(credential, query_str)
            catalog="wide"

        r=pd.DataFrame(results["rows"])
        if len(results["rows"]) > 0:

            dp_result = pd.DataFrame(results["rows"],columns=results["fields"])

            ra_ref = dp_result["ra"].tolist()
            dec_ref = dp_result["dec"].tolist()
            cref = SkyCoord([ra], [dec], frame='icrs', unit='deg')
            c1 = SkyCoord(ra_ref, dec_ref, frame='icrs', unit='deg')
            hsc_distance = cref.separation(c1).arcsec

            dp_result["separation"] = hsc_distance
            dp_result["catalog"] = catalog
            r = dp_result



        return r


    def httpJsonPost(self,url, data):
        data['clientVersion'] = self.version
        postData = json.dumps(data)
        return self.httpPost(url, postData, {'Content-type': 'application/json'})

    def httpPost(self,url, postData, headers):
        req = urllib.request.Request(url, postData.encode('utf-8'), headers)
        res = urllib.request.urlopen(req)
        return res

    def preview(self,credential, sql):
        url_api = "https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/"
        url = url_api + 'preview'
        catalog_job = {
            'sql': sql,
            'release_version': "pdr2",
        }
        postData = {'credential': credential, 'catalog_job': catalog_job}
        res = self.httpJsonPost(url, postData)
        result = json.load(res)

        print("HST found ",len(result['result']['rows']))
        if result['result']['count'] > len(result['result']['rows']):
            raise self.QueryError('only top %d records are displayed !' % len(result['result']['rows']))
        return result['result']


class QueryError(Exception):
    pass