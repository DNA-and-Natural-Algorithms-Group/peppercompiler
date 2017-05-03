#!/usr/bin/env python

from setuptools import setup, find_packages
from distutils.command.build import build
from setuptools.command.develop import develop

class build_with_spurious(build):
    def run(self):
        import os
        print "Compiling spuriousSSM"
        os.system("cc -Wall -O3 SpuriousDesign/spuriousSSM.c -o pepper/_spuriousSSM -lm")
        
        build.run(self)

class develop_with_spurious(develop):
    def run(self):
        import os
        print "Compiling spuriousSSM"
        os.system("cc -Wall -O3 SpuriousDesign/spuriousSSM.c -o pepper/_spuriousSSM -lm")
        
        develop.run(self)

setup(
    name = "pepper",
    version = "v0.0.1dev",
    packages = ['pepper'],

    install_requires = ["pyparsing","xdg","six"],

    include_package_data=True,
    package_data={'pepper': ['_spuriousSSM']},
    
    cmdclass={'build': build_with_spurious, 'develop': develop_with_spurious},
    
    entry_points={ 'console_scripts': [
        'pepper-compiler = pepper.compiler:main',
        'pepper-design-spurious = pepper.design.spurious_design:main',
        'pepper-finish = pepper.finish:main',
        'pepper-config = pepper.config:main',
        'spuriousSSM = pepper._spuriousSSM_wrapper:main']},

    author = "Constantine Evans et al (this version)",
    author_email = "cgevans@evans.foundation",
    description = "PepperCompiler in a pythonic form"
)
