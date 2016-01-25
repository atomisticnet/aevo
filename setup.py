from setuptools import setup, find_packages
from numpy import get_include

with open("README.rst") as f:
    long_desc = f.read()
    ind = long_desc.find("\n")
    long_desc = long_desc[ind + 1:]

setup(
    name="pygs",
    packages=find_packages(),
    version="0.1.0a1",
    install_requires=["numpy>=1.5"],
    author="Alexander Urban",
    author_email="aurban@lbl.gov",
    maintainer="Alexander Urban",
    license="MIT",
    description="A collection of tools for the pymatgen library.",
    long_description=long_desc,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    include_dirs=[get_include()]
)
