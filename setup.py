from __future__ import print_function
from setuptools import setup, find_packages
import sys

missing_modules = []
INSTALL_ERROR   = False

# Try importing necessary modules to test
# and see if they are installed.
try:
    from Cython.Build import cythonize
except ImportError:
    missing_modules.append('cython')

try:
    import numpy
except ImportError:
    missing_modules.append('numpy')

try:
    import scipy
except ImportError:
    missing_modules.append('scipy')

try:
    import pandas
except ImportError:
    missing_modules.append('pandas')

try:
    import matplotlib
except ImportError:
    missing_modules.append('matplotlib')

try:
    import seaborn
except ImportError:
    missing_modules.append('seaborn')

if len(missing_modules) > 0:
    INSTALL_ERROR = True
    print("\nERROR:")
    print("  You are missing the following required modules:")
    for m in missing_modules:
        print("  \t", m)
    print("\n")

if INSTALL_ERROR:
    print("\nERROR:")
    print("  Unable to install hyde.")
    print("  Please see the documentation at http://hybridization-detection.rtfd.io/.\n")
    sys.exit()
else:
    setup(name="phyde",
        version="0.1.0",
        description="Hybridization detection using phylogenetic invariants",
        long_description=open('README.rst').read(),
        url="https://github.com/pblischak/HyDe",
        author="Paul Blischak & Laura Kubatko",
        author_email="blischak.4@osu.edu",
        license="GPL",
        packages=find_packages(),
        ext_modules=cythonize("phyde/**/*.pyx", compiler_directives={'--cplus': True}),
        include_dirs=[numpy.get_include()],
        scripts=['scripts/run_hyde.py', 'src/hyde_cpp'],
        zip_safe=False
    )
