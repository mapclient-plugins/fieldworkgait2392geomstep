"""
MAP Client, a program to generate detailed musculoskeletal models for OpenSim.
    Copyright (C) 2012  University of Auckland
    
This file is part of MAP Client. (http://launchpad.net/mapclient)

    MAP Client is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MAP Client is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MAP Client.  If not, see <http://www.gnu.org/licenses/>..
"""

from setuptools import setup, find_packages
import sys, os, io

def readfile(filename, split=False):
    with io.open(filename, encoding="utf-8") as stream:
        if split:
            return stream.read().split("\n")
        return stream.read()

readme = readfile("README.md", split=True)[3:]  # skip title
requires = readfile("requirements.txt", split=True)
license = readfile("LICENSE")
package_data = {
  'mapclientplugins.fieldworkgait2392geomstep': [
    'data/*',
  ],
}


setup(name=u'mapclientplugins.fieldworkgait2392geomstep',
    version='0.1.2',
    description='',
    long_description='\n'.join(readme) + license,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
    ],
    author=u'Ju Zhang',
    author_email='',
    url='https://github.com/mapclient-plugins/fieldworkgait2392geomstep',
    license='APACHE',
    packages=find_packages(exclude=['ez_setup',]),
    namespace_packages=['mapclientplugins'],
    include_package_data=True,
    package_data=package_data,
    zip_safe=False,
    install_requires=requires,
    )
