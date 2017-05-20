from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

setup(name="HyDe",
      version="0.1.0",
      description="Hybrid detection using phylogenetic invariants",
      long_description=open('README.rst').read(),
      url="https://github.com/pblischak/HyDe",
      author="Paul Blischak & Laura Kubatko",
      author_email="blischak.4@osu.edu",
      license="GPL",
      packages=find_packages(),
      ext_modules=cythonize("hyde/**/*.pyx"),
      include_dirs=[numpy.get_include()],
      scripts=['scripts/run_hyde.py', 'src/hyde_cpp'],
      install_requires = ['numpy', 'pandas', 'Cython', 'matplotlib'],
      zip_safe=False
)
