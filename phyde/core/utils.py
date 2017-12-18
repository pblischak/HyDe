"""Utilities for running HyDe analyses"""

from __future__ import print_function

def expand_prefix(prefix):
    """
    Check for '~', '.', and '..' in file path names
    and expand to full path. Modeled code from Deren
    Eaton's pyRAD (pyRAD.py, lines 155-164).
    """
    if prefix[0] == "~":
        home = os.path.expanduser("~")
        return home + prefix[1:]
    elif prefix[0:2] == "..":
        dirname = os.path.abspath("..")
        return dirname + prefix[2:]
    elif prefix[0] == "." and prefix[1] != ".":
        dirname = os.path.abspath(".")
        return dirname + prefix[1:]
    else:
        return prefix
