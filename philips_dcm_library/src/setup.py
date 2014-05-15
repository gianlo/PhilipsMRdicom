'''
Created on 15 May 2014

@author: glf12
'''

from distutils.core import setup

setup(name='philips_dcm',
      version='1.0',
      description='A python interface for Philips MR multiframe dicom files',
      author='Gianlorenzo Fagiolo',
      url='https://github.com/gianlo/PhilipsMRdicom',
      license='License :: OSI Approved :: MIT License',
      requires=['numpy', 'dicom'],
      packages=['philips_dcm'])
