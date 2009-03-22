#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

# TODO: add tests for required packages
# numpy, pyrex, matplotlib and scipy

# build the pure c++ part
libcadda = Extension('libcadda', 
                     ['src/cadda_parameters.cpp',
                      'src/cadda_io.cpp',
                      'src/cadda_optimise.cpp',
                      'src/cadda_index.cpp', 
                      'src/cadda_convert.cpp',
                      'src/cadda_mst.cpp',])
                               
# build the interface
cadda = Extension(
    "cadda",                   # name of extension
    [ "src/cadda.pyx",        # filename of our Pyrex/Cython source
        'src/cadda_parameters.cpp',
        'src/cadda_io.cpp',
        'src/cadda_optimise.cpp',
        'src/cadda_index.cpp', 
        'src/cadda_convert.cpp',
        'src/cadda_mst.cpp',],
    library_dirs=[],
    libraries=[],              
    language="c++",               # this causes Pyrex/Cython to create C++ source
    )

Components = Extension(
     "Components",                   # name of extension
     [ "src/Components.pyx",        # filename of our Pyrex/Cython source
       'src/connected_components.cpp',],
     library_dirs=[],
     libraries=[],
     language="c++",               # this causes Pyrex/Cython to create C++ source
     )

setup(name='Adda',
      version='1.0',
      description='ADDA - automatic domain detection algorithm',
      author='Andreas Heger',
      author_email='andreas.heger@helsinki.fi',
      url='http://wwwfgu.anat.ox.ac.uk/~andreas',
      packages = ["Adda",],
      package_dir = {'Adda': 'lib'},
      package_data = { "Adda" : [ 'adda_blast_parser.pl',] },
      scripts=['scripts/adda.py' ],
      ext_modules=[cadda, Components ],
      cmdclass = {'build_ext': build_ext}
     )
    
