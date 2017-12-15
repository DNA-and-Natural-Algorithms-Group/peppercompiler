#!/usr/bin/env python

from setuptools import setup, find_packages
from distutils.command.build import build
from setuptools.command.develop import develop

class build_with_spurious(build):
    def run(self):
        import os
        if "CC" in os.environ:
            cc = os.environ['CC']
        else:
            cc = "cc"
        os.system("{} -Wall -O3 peppercompiler/SpuriousDesign/spuriousSSM.c -o peppercompiler/_spuriousSSM -lm".format(cc))

        build.run(self)

class develop_with_spurious(develop):
    def run(self):
        import os
        os.system("cc -Wall -O3 peppercompiler/SpuriousDesign/spuriousSSM.c -o peppercompiler/_spuriousSSM -lm")

        develop.run(self)

setup(
    name = "peppercompiler",
    version = "0.1.0",
    packages = ['peppercompiler', 'peppercompiler.design'],

    install_requires = ["pyparsing","six"],

    include_package_data=True,
    package_data={'peppercompiler': ['_spuriousSSM', 'SpuriousDesign/spuriousSSM.c']},

    test_suite='peppercompiler.tests',

    cmdclass={'build': build_with_spurious, 'develop': develop_with_spurious},

    entry_points={ 'console_scripts': [
        'pepper-compiler = peppercompiler.compiler:main',
        'pepper-design-spurious = peppercompiler.design.spurious_design:main',
        'pepper-finish = peppercompiler.finish:main',
        'spuriousSSM = peppercompiler._spuriousSSM_wrapper:main']},

    author = "Constantine Evans et al (this version)",
    author_email = "cge@dna.caltech.edu",
    description = "PepperCompiler in a pythonic form"
)
