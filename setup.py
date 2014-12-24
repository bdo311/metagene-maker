#!/usr/bin/env python

from setuptools import setup

setup(name='metagene-maker',
      version='0.1',
      description='Obtain average profiles of NGS datasets over your favorite regions',
      
	  requires = ['numpy (>=1.7)'],
	  
	  packages=['metagene_maker'],
	  package_dir={'metagene_maker': 'src/'},
	  package_data={'': ['README*'], 'metagene_maker': ['src/*.r']},
	  
	  entry_points = {'console_scripts': ['metagene_maker = metagene_maker.metagene_maker:main']},

	  author='Brian Do',
      author_email='bdo@stanford.edu',
      url='https://www.github.com/bdo311/metagene-maker',
     )