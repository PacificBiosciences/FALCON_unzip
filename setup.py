#!/usr/bin/env python

from setuptools import setup

from distutils.core import Extension

import glob

install_requires=[ "networkx >= 1.7" ]

scripts = glob.glob("src/py_scripts/*.py")

setup(name='falcon_unzip',
      version='0.1.0',
      description='Falcon unzip',
      author='Jason Chin',
      author_email='jchin@pacificbiosciences.com',
      packages=['falcon_unzip'],
      package_dir={'falcon_unzip':'src/py/'},
      scripts = scripts,
      zip_safe = False,
      install_requires=install_requires
     )

