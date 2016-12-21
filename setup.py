#!/usr/bin/env python

from setuptools import setup

from distutils.core import Extension

import glob

install_requires=[ "networkx >= 1.7", "pysam >= 0.8.4" ]

scripts = glob.glob("src/py_scripts/*.py")

setup(name='falcon_unzip',
      version='0.4.0',
      description='Falcon unzip',
      author='Jason Chin',
      author_email='jchin@pacificbiosciences.com',
      maintainer='Christopher Dunn',
      maintainer_email='pb.cdunn@gmail.com',
      packages=['falcon_unzip'],
      package_dir={'falcon_unzip':'falcon_unzip/'},
      scripts = scripts,
      zip_safe = False,
      install_requires=install_requires
     )

