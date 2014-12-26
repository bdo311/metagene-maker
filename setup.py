#!/usr/bin/env python

from setuptools import setup

setup(name='metagenemaker',
      version='0.2',
      description='Obtain average profiles of NGS datasets over your favorite regions',
      
	  install_requires = ['numpy>=1.7'],
	  
	  packages=['metagene_maker', 'knownGenes'],
	  package_dir={'metagene_maker': 'metagene_maker', 'knownGenes': 'knownGenes'},
	  package_data={'metagene_maker': ['rscript/*.r']},
	  
	  entry_points = {'console_scripts': ['metagene_maker = metagene_maker.metagene_maker:main', 'knownGenes = knownGenes.knowngenes_to_transcript_regions:main']},
	  
	  author='Brian Do',
      author_email='bdo@stanford.edu',
      url='https://www.github.com/bdo311/metagene-maker',
	  license = "GPL2",
     )
