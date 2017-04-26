#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = "pepper",
    version = "v0.0.1dev",
    packages = ['pepper'],

    install_requires = [],

    include_package_data=True,

    #entry_points={ 'console_scripts': [] },

    author = "Constantine Evans et al (this version)",
    author_email = "cgevans@evans.foundation",
    description = "PepperCompiler in a pythonic form"
}
