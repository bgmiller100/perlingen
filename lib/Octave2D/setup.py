# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 18:31:11 2019

@author: Ben
"""

from distutils.core import setup, Extension 

setup(name = "Octave2D", ext_modules = [Extension('Octave2D',sources=['Octave2D.c'])])
    