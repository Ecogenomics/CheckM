#!/usr/bin/env python3

import os
from setuptools import setup


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'checkm', 'VERSION'))
    return versionFile.readline().strip()


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='checkm-genome',
    version=version(),
    author='Donovan Parks, Michael Imelfort, Connor Skennerton',
    author_email='donovan.parks@gmail.com',
    packages=['checkm', 'checkm.plot', 'checkm.test', 'checkm.util'],
    scripts=['bin/checkm'],
    package_data={'checkm': ['VERSION', 'DATA_CONFIG']},
    include_package_data=True,
    url='http://pypi.python.org/pypi/checkm-genome/',
    license='GPL3',
    description='Assess the quality of putative genome bins.',
    long_description=readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=[
        "numpy >= 1.21.3",
        "scipy >= 1.7.3",
        "matplotlib >= 3.5.1",
        "pysam >= 0.19.0",
        "dendropy >= 4.5.2",
        "setuptools"],
    zip_safe=False
)
