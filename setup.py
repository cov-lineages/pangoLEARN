from setuptools import setup, find_packages
import glob
import os
import pkg_resources
# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from pangoLEARN import __version__, _program

setup(name='pangoLEARN',
      version=__version__,
      packages=find_packages(),
      scripts=[],
      package_data={'pangoLEARN':['data/*']},
      description='trained lineage data',
      url='https://github.com/cov-lineages/pangoLEARN',
      author='cov-lineages group',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = pangoLEARN.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)