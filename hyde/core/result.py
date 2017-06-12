""" Class for processing HyDe output. """

from __future__ import print_function
import numpy as np
import pandas as pd

class HydeResult:
    """
    A class for reading in and working with results from a HyDe analysis
    using the C++ interface.
    """
    def __init__(self, infile):
        self.infile = infile
        
