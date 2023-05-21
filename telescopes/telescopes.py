from core.fits.FITSFIle import FITSFile
from telescopes.liverpool import LivepoolTelescope
from telescopes.decals import DecalsTelescope
from telescopes.panstarrs import PanstarrsTelescope
from telescopes.wht import WHTTelescope
from telescopes.hst import HSTelescope
from telescopes.las_cumbres import LasCumbres
from telescopes.vlt import VLTelescope
from telescopes.gtc import GTCTelescope


class Telescopes():

    telescopes = {
        "liverpool":LivepoolTelescope(),
        "decals":DecalsTelescope(),
        "panstarrs":PanstarrsTelescope(),
        "wht":WHTTelescope(),
        "0m4-03":LasCumbres(),
        "0m4-06": LasCumbres(),
        "hst":HSTelescope(),
        "vlt": VLTelescope(),
        "1m0-11":LasCumbres(),
        "gtc":GTCTelescope()

    }

    keys={
        "ctio4.0-mtelescope":"decals",
        "ctio": "decals",
        "ps1":"panstarrs",
        "liverpooltelescope":"liverpool",
        "wht":"wht",
        "0m4-03":"0m4-03",
        "0m4-06": "0m4-06",
        "hst":"hst",
        "vlt":"vlt",
        "1m0-11":"1m0-11",
        "gtc":"gtc"

    }

    @staticmethod
    def getTelescopeData(filepath,header_idx=0):
        file = FITSFile.open(filepath)
        telescopes = Telescopes()
        name = ""
        instrument=""
        if len(file)>1:
            for idx,h in enumerate(file):
                header = h.header
                if "NAXIS" in header and header["NAXIS"] > 0:
                    header_idx=idx
                    break

        if header_idx != None:
            file[header_idx].header
            header = file[header_idx].header
            name = telescopes.getTelescopeName(header)
            instrument = telescopes.getInstrument(header)
        else:
            for hd in file:
                header = file.header
                name = telescopes.getTelescopeName(header)
                instrument = telescopes.getInstrument(header)

                if name!="":
                    break

        print(name,instrument)
        if name != "":
            tele = telescopes.getTelescopeById(name)
        if instrument!="":
            params = tele.getParams(instrument)
        else:
            params = tele.getParams()

        params["instrument"] = instrument
        params["instance"] = tele


        return params

    def getTelescopeById(self,id):
        id=id.lower().replace(" ","")
        if id in self.keys:
            name=self.keys[id]
            tel = self.telescopes[name]
            return tel

    def getInstrument(self,header):
        instrument_name = ""
        if "INSTRUME" in header:
            instrument_name  = header["INSTRUME"]


        if instrument_name=="":
            list_keys = header.keys()
            string_keys = " ".join(list_keys)
            index = string_keys.find("INSTRUME")
            index_in = string_keys.find(" ", index - 10)
            index_out = string_keys.find(" ", index_in + 1)
            key = string_keys[index_in:index_out]
            instrument_name = header[key]

        return str(instrument_name).replace(" ","").upper()

    def getTelescopeName(self,header):
        telescope_name = ""



        if "TELESCOP" in header:
            telescope_name = header["TELESCOP"]
        elif "OBSERVAT" in header:
            telescope_name = header["OBSERVAT"]

        if telescope_name=="":
            list_keys = header.keys()
            string_keys = " ".join(list_keys)
            index = string_keys.find("TELESCOP")
            index_in = string_keys.find(" ", index - 10)
            index_out = string_keys.find(" ", index_in + 1)
            key = string_keys[index_in:index_out]
            telescope_name = header[key]




        return str(telescope_name).replace(" ","").upper()

    @staticmethod
    def getPixelScale(file):
        header = FITSFile.header(file)

        try:
            pixscale = header["PIXSCALE"]
        except Exception:
            tele = Telescopes.getTelescopeData(file)
            pixscale = tele["pixelscale"]
        return pixscale
