'''
Created on 22 Jan 2013

@author: gfagiolo
'''

import dicom
# import logging

#===============================================================================
# CONSTANTS
#===============================================================================

#SOPClassUID = (0008, 0016)
TAG_SOP_CLASS_UID = dicom.tag.Tag(dicom.datadict.NameDict['SOPClassUID'])
#DimensionIndexes = (0020, 9222)
TAG_DIMENSION_INDEX_SEQUENCE = dicom.tag.Tag(dicom.datadict.NameDict['DimensionIndexes'])
#PerframeFunctionalGroups = (5200, 9230) 
TAG_PER_FRAME_FUNCTIONAL_GROUPS_SEQUENCE = dicom.tag.Tag(dicom.datadict.NameDict['PerframeFunctionalGroups'])
UID_ENHANCED_MR_IMAGE_STORAGE = '1.2.840.10008.5.1.4.1.1.4.1'
UID_MR_SPECTROSCOPY_STORAGE = '1.2.840.10008.5.1.4.1.1.4.2'
#DTI related content
DTIDIRNO = 'NONE'
DTIDIRYES = 'DIRECTIONAL'
DTIDIRISO = 'ISOTROPIC'

METADATA_DESCRIPTION ="""#DICOM metadata in flat txt-format
#'key' = 'value'
#where key is a hierarchical set of  ('.'-separated) dicom tags and value is the corresponding value
#Parameters common to all frames can be found as either direct parameters (i.e. 'BodyPartExamined')
#or under the sequence 'SharedFunctionalGroupsSequence' (i.e. 'RepetitionTime')
#frame specific parameters can be found under the sequence 'PerFrameFunctionalGroupsSequence', 
#note that the number 'n' in 'PerFrameFunctionalGroupsSequence[n]' refers to the frame described
#Frame specific tags are for example 'EffectiveEchoTime'
"""

ANONYMISE_FIELDS = set(
('PatientName',
'PatientID',
'PatientBirthDate', 
'PatientSex',
'OperatorsName',
'RequestingPhysician',
'ScheduledPerformingPhysicianName',
'ReferringPhysicianName',))

#===============================================================================
# FUNCTIONS
#===============================================================================
    
def inspect_dicom(fname):
    fp = open(fname, 'rb')
    _ = fp.read(0x80)
    magic = fp.read(4)
    if magic != "DICM":
        #logging.debug(fname,'not a dicom file')
        raise dicom.filereader.InvalidDicomError
    meta = dicom.filereader._read_file_meta_info(fp)
    fp.close()
    return meta

def is_storagesopclassuid(df, uid):
    if isinstance(df, dicom.dataset.FileDataset):
        return df[TAG_SOP_CLASS_UID].value == uid
    elif isinstance(df, str):
        try:
            return inspect_dicom(df).MediaStorageSOPClassUID == uid
        except  dicom.filereader.InvalidDicomError:
            #this is not a dicom file
            return False
    else:
        return False

def is_accepted_dicom(df):
    if isinstance(df, dicom.dataset.FileDataset):
        return df[TAG_SOP_CLASS_UID].value in [UID_ENHANCED_MR_IMAGE_STORAGE, UID_MR_SPECTROSCOPY_STORAGE]
    elif isinstance(df, str):
        try:
            return inspect_dicom(df).MediaStorageSOPClassUID in [UID_ENHANCED_MR_IMAGE_STORAGE, UID_MR_SPECTROSCOPY_STORAGE]
        except  dicom.filereader.InvalidDicomError:
            #this is not a dicom file
            return False
    else:
        return False
    
    
def is_multiframe(df):
    return is_storagesopclassuid(df, UID_ENHANCED_MR_IMAGE_STORAGE)

def is_mrspectroscopystorage(df):
    return is_storagesopclassuid(df, UID_MR_SPECTROSCOPY_STORAGE)

def get_a_frame(df, frame_number=0):
    return df[TAG_PER_FRAME_FUNCTIONAL_GROUPS_SEQUENCE][frame_number]

def get_frame_content(frame):
    return frame.FrameContentSequence[0]

def get_shared_functional_group_sequence(df):
    return df.SharedFunctionalGroupsSequence[0]

def get_shared_functional_group_sequence_repetion_time(SharedFunctionalGroupsSequence):
#    logging.debug("SharedFunctionalGroupsSequence %d items"%len(df.SharedFunctionalGroupsSequence))
    return SharedFunctionalGroupsSequence.MRTimingAndRelatedParametersSequence[0].RepetitionTime

def dicomobj_to_str(seq, level=0, prefix=None, only_non_private=True, anonymise=True):
    def tag_to_name(atag):
        try:
            return dicom.datadict.DicomDictionary[atag][-1]
        except KeyError:
            return str(atag).replace(' ','')
    def my_repr(self):
        repVal = self.repval
        s = '%s = %s'%(tag_to_name(self.tag), repVal)         
        return s
    strings = []
    for data_element in seq:
        name = tag_to_name(data_element.tag)
        if anonymise and name in ANONYMISE_FIELDS:
            continue
        if data_element.VR == "SQ":
            #this is a sequence, use its name as the start of the fields contained
            if only_non_private and data_element.tag in dicom.datadict.DicomDictionary:   # a sequence
                #skip private fields
                continue
            #name = tag_to_name(data_element.tag)
            for indx, dataset in enumerate(data_element.value):
                if prefix is None:
                    nprefix = name + '[%d].'%(indx+1)
                else:
                    nprefix = prefix + name + '[%d].'%(indx+1)    
                strings.append(dicomobj_to_str(dataset, level + 1, nprefix, only_non_private))
        else:
            if only_non_private and data_element.tag in dicom.datadict.DicomDictionary:
                #only non-private fields
                if prefix is None:
                    strings.append(my_repr(data_element))
                else:
                    strings.append(prefix + my_repr(data_element))
            elif not only_non_private:
                if prefix is None:
                    strings.append(my_repr(data_element))
                else:
                    strings.append(prefix + my_repr(data_element))

    return "\n".join(strings)

def get_repetion_time(df):
#    logging.debug("SharedFunctionalGroupsSequence %d items"%len(df.SharedFunctionalGroupsSequence))
    return df.SharedFunctionalGroupsSequence[0].MRTimingAndRelatedParametersSequence[0].RepetitionTime

def get_frame_repetion_time(frame):
    return frame.MRTimingAndRelatedParametersSequence[0].RepetitionTime

def define_DimensionIndexValues(df):
    desc = []
    for el in df[TAG_DIMENSION_INDEX_SEQUENCE]:
        desc.append((el.DimensionIndexPointer, el.FunctionalGroupPointer))
    return desc

def DimensionIndexes_to_tagnames(df):
    def find_tag_name(t):
        if t in dicom.datadict.DicomDictionary:
            tagname = dicom.datadict.DicomDictionary[t][2]
        else:
            tagname = str(t)
        return tagname
    desc = define_DimensionIndexValues(df)
    return [(find_tag_name(e[0]), find_tag_name(e[1])) for e in desc]

def get_frame_index(frame, idx):
    return [frame[e.FunctionalGroupPointer][0][e.DimensionIndexPointer] for e in idx]
        
def get_dimension_index_description_position(idx, desc):
    return [x[1] for x in filter(lambda x:x[0].DimensionDescriptionLabel == desc, zip(idx, range(idx.VM)))][0]

class DiffusionFrameInfo(object):

    FRAME_NO = None
    DIFF_TYPE = None
    DIFF_BVALUE = None
    DIFF_BVEC = None
    
    def __init__(self, frame, frame_no):
        self.set_frame_no(frame_no)
        dti_info = diffusionInfo(frame)
        self.set_diff_type(dti_info[0])
        self.set_diff_bvalue(dti_info[1])
        self.set_diff_bvec(dti_info[2])
        
    def get_frame_no(self):
        return self.__FRAME_NO

    def get_diff_type(self):
        return self.__DIFF_TYPE

    def get_diff_bvalue(self):
        return self.__DIFF_BVALUE

    def get_diff_bvec(self):
        return self.__DIFF_BVEC

    def set_frame_no(self, value):
        self.__FRAME_NO = value

    def set_diff_type(self, value):
        self.__DIFF_TYPE = value

    def set_diff_bvalue(self, value):
        self.__DIFF_BVALUE = value

    def set_diff_bvec(self, value):
        self.__DIFF_BVEC = value

    def isDirectional(self):
        return self.get_diff_type() == DTIDIRYES

    def isUnweighted(self):
        return self.get_diff_type() == DTIDIRNO

    def toList(self):
        return [self.get_diff_type(), self.get_diff_bvalue(), self.get_diff_bvec(), self.get_frame_no()]

    def __str__(self):
        return str(self.toList())

    def __repr__(self):
        return str(self)
    
    FRAME_NO = property(get_frame_no, set_frame_no, None, "FRAME_NO's docstring")
    DIFF_TYPE = property(get_diff_type, set_diff_type, None, "DIFF_TYPE's docstring")
    DIFF_BVALUE = property(get_diff_bvalue, set_diff_bvalue, None, "DIFF_BVALUE's docstring")
    DIFF_BVEC = property(get_diff_bvec, set_diff_bvec, None, "DIFF_BVEC's docstring")

def diffusionInfo(frame):
    out = [frame.MRDiffusionSequence[0].DiffusionDirectionality,
           None,
           None]
    if 'DiffusionBValue' in frame.MRDiffusionSequence[0]:
        out[1] = frame.MRDiffusionSequence[0].DiffusionBValue
    if 'DiffusionGradientDirectionSequence' in frame.MRDiffusionSequence[0]:
        out[2] = list(
            frame.MRDiffusionSequence[0].DiffusionGradientDirectionSequence[0].DiffusionGradientOrientation)
    return out
