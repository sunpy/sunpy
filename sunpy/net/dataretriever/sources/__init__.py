from __future__ import absolute_import, division, print_function
from ..client import GenericClient

__all__ = ['EVEClient','GOESClient','LYRAClient','NOAAIndicesClient',
           'NOAAPredictClient','NoRHClient','RHESSIClient','ERNEClient',
           'SEPTClient','HETClient','SITClient','PLASTICClient',
           'MAGClient','LETClient']

from .eve import EVEClient
from .goes import GOESClient
from .lyra import LYRAClient
from .noaa import NOAAIndicesClient,NOAAPredictClient
from .norh import NoRHClient
from .rhessi import RHESSIClient
from .soho import ERNEClient
from .stereo import SEPTClient,HETClient,SITClient,PLASTICClient,MAGClient,LETClient


