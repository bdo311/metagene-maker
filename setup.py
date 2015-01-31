#!/usr/bin/python

from setuptools import setup
import sys, site, os

setup(name='metagenemaker',
      version='0.6.0a1',
      description='Obtain average profiles of NGS datasets over your favorite regions',
      
	  install_requires = ['numpy>=1.7', 'pandas>=0.14'],
	  
	  packages=['metagene_maker', 'extractTranscriptRegions'],
	  package_dir={'metagene_maker': 'metagene_maker', 'extractTranscriptRegions': 'extractTranscriptRegions'},
	  
	  entry_points = {'console_scripts': ['metagene_maker = metagene_maker.metagene_maker:main', 'extractTranscriptRegions = extractTranscriptRegions.extract_transcript_regions:main',
	  'metagene_subsets = metagene_maker:metagene_subsets:main']},
	  
	  author='Brian Do',
      author_email='bdo@stanford.edu',
      url='https://www.github.com/bdo311/metagene-maker',
	  license = "GPL2",
     )
