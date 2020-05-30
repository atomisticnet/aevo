from setuptools import setup, find_packages
from codecs import open
import os
import glob

__author__ = "Alexander Urban, Nongnuch Artrith"
__email__ = "aurban@atomistic.net, nartrith@atomistic.net"

here = os.path.abspath(os.path.dirname(__file__))
package_name = 'aevo'
package_description = 'Atomistic evolution'

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as fp:
    long_description = fp.read()

# Get version number from the VERSION file
with open(os.path.join(here, package_name, 'VERSION')) as fp:
    version = fp.read().strip()

setup(
    name=package_name,
    packages=find_packages(),
    version=version,
    install_requires=["numpy>=1.5"],
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    license="MPL2",
    description=package_description,
    long_description=long_description,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
    scripts=glob.glob(os.path.join("scripts", "*.py"))
)
