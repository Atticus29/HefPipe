#!/usr/bin/env python

from distutils.core import setup

setup(name='HeFPipe',
	version='1.0',
	description='Run heterozygosity-fitness correlation analyses',
	author='Mark A. Fisher',
	author_email='mark.aaron.fisher@gmail.com',
	url='https://github.com/Atticus29/HefPipe_repos',
	py_modules=['HefPipe_modules','pyper'], # copy to site-packages
	)

