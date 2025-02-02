from setuptools import setup
import sys
sys.path.insert(0, ".")
from PyOAE import __version__

setup(
    name='PyOAE',
    version=__version__,
    author='Greg Pelletier',
    py_modules=['PyOAE'], 
    install_requires=['numpy','scipy','PyCO2SYS'],
)