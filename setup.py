from setuptools import setup

import glob
import os

folder = os.path.dirname(os.path.realpath(__file__))

silent_exe = []
for file in glob.glob(folder + "/silent*"):
   file = os.path.basename(file)
   if '_' in file:
      continue
   if '.' in file:
      continue
   silent_exe.append(file)

__version__ = "1.0.0"

setup(
   name='silent_tools',
   version='1.0.0',
   description='A bunch of shell utilities for dealing with silent files',
   license='MIT',
   author='Brian Coventry',
   author_email='bcoventry77@gmail.com',
   url="https://github.com/bcov77/silent_tools",
   download_url="https://github.com/bcov77/silent_tools/archive/refs/tags/v1.0.0.tar.gz",
   packages=['_helpers_silent'],  #same as name
   py_modules=['silent_tools'],
   install_requires=['numpy'], #external packages as dependencies
   scripts=silent_exe + ['_extract_flags.sh', 'pyjd2', 'pyextract_pdbs']
)