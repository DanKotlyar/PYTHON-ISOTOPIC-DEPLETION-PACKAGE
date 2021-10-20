from pyIsoDep.functions import *
from pyIsoDep.data import *
from pyIsoDep.tests import *

import os

_ROOT = os.path.abspath(os.path.dirname(__file__))
def setDataPath(path):
    return os.path.join(_ROOT, 'data', path)
