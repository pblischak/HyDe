from setuptools import setup

setup(name="HyDe",
      version="0.1.0a",
      description="Hybrid detection using phylogenetic invariants",
      long_description=open('README.md').read(),
      url="https://github.com/pblischak/HyDe",
      author="Paul Blischak & Laura Kubatko",
      author_email="blischak.4@osu.edu",
      license="GPL",
      packages=['hyde'],
      scripts=['src/hyde_cpp', 'scripts/run_hyde.py'],#, 'scripts/hyde'],
      install_requires = ['numpy', 'pandas'],
      zip_safe=False
)
