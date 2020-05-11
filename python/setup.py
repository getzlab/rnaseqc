import os
import re
from setuptools import setup, find_packages
with open("rnaseqc/__init__.py") as reader:
    __version__ = re.search(
        r'__version__ ?= ?[\'\"]([\w.]+)[\'\"]',
        reader.read()
    ).group(1)
with open(os.path.join(os.path.dirname(__file__), 'README.md')) as r:
    long_description = r.read()

# Setup information
setup(
    name = 'rnaseqc',
    version = __version__,
    packages = find_packages(),
    description = 'Multi-sample visualization of metrics from RNA-SeQC',
    long_description = long_description,
    long_description_content_type='text/markdown',
    install_requires = [
        'numpy',
        'pandas',
        'matplotlib',
        'seaborn',
        'qtl',
        'agutil',
        'nbformat'
    ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
