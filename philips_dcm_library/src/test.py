'''
Created on 22 Jan 2013

@author: gfagiolo
'''

import philips_dcm as PD
import os
import numpy as np
from pynii import Nifti1Data
##import nibabel as nb


TEST_DATA_PATH = r'D:\var\common\mrs'

if __name__ == '__main__':
    fname = os.path.join(TEST_DATA_PATH,'MRS.dcm')
    pdc = PD.PhilipsMultiFrameDcm(fname)
    #print pdc.get_minimal_header().get_nifti_affine()
    nii = Nifti1Data()
    nii.setAffine(pdc.get_minimal_header().get_nifti_affine())
    nii.setData(100*np.ones((1,1,1), dtype=np.int16))
    Nifti1Data.save(nii, fname + '.nii.gz')