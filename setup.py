from setuptools import setup

setup(name="HyDe",
      version="0.1.0a",
      description="Hybrid detection using phylogenetic invariants",
      long_description=open('README.rst').read(),
      url="https://github.com/pblischak/HyDe",
      author="Paul Blischak & Laura Kubatko",
      author_email="blischak.4@osu.edu",
      license="GPL",
      packages=['hyde'],
      scripts=['scripts/run_hyde.py'],
      data_files=[('bin', ['src/hyde_cpp'])],
      install_requires = ['numpy', 'pandas'],
      zip_safe=False
)
