from setuptools import setup, find_packages
from setuptools.extension import Extension

from Cython.Build import cythonize
import numpy


setup(
    packages=find_packages(),
    ext_modules=cythonize(
        [
            Extension(
                "phyde.data",
                ["phyde/data.pyx"],
                include_dirs=[numpy.get_include()],
                language="c++",
            ),
        ]
    ),
    scripts=[
        "scripts/run_hyde.py",
        "scripts/run_hyde_mp.py",
        "scripts/individual_hyde.py",
        "scripts/individual_hyde_mp.py",
        "scripts/bootstrap_hyde.py",
        "scripts/bootstrap_hyde_mp.py",
        "scripts/hyde_gui.py",
    ],
    zip_safe=False,
)
