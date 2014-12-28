#!/usr/bin/env python

from setuptools import setup
import sys, site, os

#sys.executable = '/usr/bin/env%s' % " python"

	
setup(name='metagenemaker',
      version='0.2',
      description='Obtain average profiles of NGS datasets over your favorite regions',
      
	  install_requires = ['numpy>=1.7'],
	  
	  packages=['metagene_maker', 'knownGenes'],
	  package_dir={'metagene_maker': 'metagene_maker', 'knownGenes': 'knownGenes'},
	  package_data={'metagene_maker': ['rscript/*.r']},
	  
	  entry_points = {'console_scripts': ['metagene_maker = metagene_maker.metagene_maker:main', 'knownGenes = knownGenes.knowngenes_to_transcript_regions:main']},
	  #options = {'build_scripts': {'script_text': "#!/usr/bin/env python", 'executable': '/usr/bin/env'}},
	  
	  author='Brian Do',
      author_email='bdo@stanford.edu',
      url='https://www.github.com/bdo311/metagene-maker',
	  license = "GPL2",
     )

print ("Fixing hashbangs")
print site.getuserbase()
exit()
def fix_hashbang(script_list):
    sed = "sed -i 's/#!\/usr\/bin\/python/#!\/usr\/bin\/env python/'"
    cmd = ' '.join([sed, ' '.join(script_list)])
    os.system(cmd)
	
parentDir = sys.prefix + '/bin/'
list_of_scripts = [parentDir + x for x in ["knownGenes", "metagene_maker"]]
fix_hashbang(list_of_scripts)