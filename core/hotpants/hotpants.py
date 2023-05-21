from astropy.io import fits


def generateFiles(**kargs):
    files={}
    if "tmplim" in kargs.keys():
        template = fits.open(kargs["tmplim"])
        files["template"]=template
        data_tmp = template[1].data
        header_tmp = template[1].header
        newTmp=fits.PrimaryHDU(data=data_tmp,header=header_tmp)
        hdul = fits.HDUList([newTmp])
        hdul.writeto('/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/hotpants/template.fits',overwrite=True)

    if "inim" in kargs.keys():
        comparation = fits.open(kargs["inim"])
        files["comparation"] = comparation
        data_comp = comparation[1].data
        header_comp = comparation[1].header
        newTmp = fits.PrimaryHDU(data=data_comp, header=header_comp)
        hdul = fits.HDUList([newTmp])
        hdul.writeto('/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/hotpants/compare.fits',overwrite=True)

    print(files)




if __name__=="__main__":

    data={}
    data["tmplim"] = "/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/hst/j9op03010_drz.fits"
    data["inim"] = "/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/decals/image-decam-215899-N29-z.fits"
    data["tmplim"] = "/Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/decals/image-decam-510449-N31-r.fits"


    generateFiles(tmplim=data["tmplim"],inim=data["inim"])




#    "hotpants -inim /Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/decals/image-decam-510449-N31-r.fits.gz -tmplim /Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/hst/j9op03010_drz.fits -outim /Users/cjimenez/Documents/PHD/data/lens_SN/SLACSJ0157-0056/hotpants/out.fits"