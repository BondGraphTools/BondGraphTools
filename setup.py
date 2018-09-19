import setuptools
from BondGraphTools.config import VERSION

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BondGraphTools",
    version=VERSION,
    author="Pete Cudmore",
    author_email="peter.cudmore@uqconnect.edu.au",
    description="Bond Graph Modelling Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/peter-cudmore/BondGraph",
    classifiers=(
        "Intended Audience :: Science/Research",
        "Intended Audience :: End Users/Desktop",
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ),
    keywords="modelling control engineering",
    packages=['BondGraphTools'],
    package_dir={'BondGraphTools': 'BondGraphTools'},
    package_data={'BondGraphTools': ['components/*.json']}
)
