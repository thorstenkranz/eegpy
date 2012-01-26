#!/usr/bin/env python

from distutils.core import setup
import sys

import eegpy 


# For some commands, use setuptools
if len(set(('develop', 'sdist', 'release', 'bdist_egg', 'bdist_rpm',
           'bdist', 'bdist_dumb', 'bdist_wininst', 'install_egg_info',
           'build_sphinx', 'egg_info', 'easy_install',
            )).intersection(sys.argv)) > 0:
    from setupegg import extra_setuptools_args

# extra_setuptools_args is injected by the setupegg.py script, for
# running the setup with setuptools.
if not 'extra_setuptools_args' in globals():
    extra_setuptools_args = dict()


setup(name='eegpy',
      version=eegpy.__version__,
      summary='Library for the analysis of EEG and MEG.',
      author='Thorsten Kranz',
      author_email='thorstenkranz@gmail.com',
      url='http://www.thorstenkranz.de',
      description="""
      Library for the analysis of EEG and MEG
""",
      long_description=file('README').read(),
      license='BSD',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Education',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Topic :: Utilities',
      ],
      platforms='any',
      package_data={
          'eegpy.ui.viewer': ['images/*_32.png'],
          'eegpy.plot' : ['sample_locs/*.ced','sample_locs/*,locs'] 
                    },
      packages=['eegpy',
                'eegpy.analysis',
                'eegpy.analysis.entropy',
                'eegpy.analysis.mvpa',
                'eegpy.cnn',
                'eegpy.computing',
                'eegpy.filter',
                'eegpy.formats',
                'eegpy.models',
                'eegpy.plot',
                'eegpy.stats',
                'eegpy.study',
                'eegpy.tools',
                'eegpy.ui',
                'eegpy.ui.icon',
                'eegpy.ui.tools',
                'eegpy.ui.viewer',
                'eegpy.ui.widgets',
                'eegpy.ui',
                ],
      #scripts=['bin/pylocator'],
      **extra_setuptools_args)

