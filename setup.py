from setuptools import setup, find_packages
from setuptools.extension import Extension

from Cython.Build import cythonize
import numpy


setup(
    name="phyde",
    version="1.0.1",
    description="Hybridization detection using phylogenetic invariants",
    long_description=open("README.rst").read(),
    url="https://github.com/pblischak/HyDe",
    author="Paul Blischak & Laura Kubatko",
    author_email="paul.blischak@gmail.com",
    packages=find_packages(),
    install_requires=[
        "Cython",
        "multiprocess",
        "numpy",
    ],
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
    license="MIT",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
    ],
    zip_safe=False,
)
