from setuptools import setup, find_packages
import sys, os

setup(name=u'mapclientplugins.fieldworkgait2392geomstep',
      version='0.1',
      description='',
      long_description="",
      classifiers=[],
      author=u'Ju Zhang',
      author_email='',
      url='https://github.com/mapclient-plugins/fieldworkgait2392geomstep',
      license='GPL',
      packages=find_packages(exclude=['ez_setup',]),
      namespace_packages=['mapclientplugins'],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'numpy',
          'scipy',
          'transforms3d',
      ],
      )
