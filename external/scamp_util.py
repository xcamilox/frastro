import subprocess
from frastro import Utils, Config, LOGUtil, ImageUtils
from frastro.external.psfex_util import PSFexUtil
import os
import os.path as path
import numpy as np
import shutil
from astropy.io import fits



class ScampUtil():
    RUN_SWARP_CM = Config.getExtApp("scamp")
    DEFAULT_PARAMS_CM = [RUN_SWARP_CM, "-d"]

