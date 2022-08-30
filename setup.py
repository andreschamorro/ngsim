# setup file for building with distutils.
#
# Eli Bendersky [http://eli.thegreenplace.net]
# This code is in the public domain.
from distutils.core import setup, Extension

module1 = Extension('ngsim', sources = ['ngsimmodule.c'], libraries = ['z',  'm'])

setup (name = 'NGsim',
        version = '0.1',
        description = 'This is a demo package',
        ext_modules = [module1])
