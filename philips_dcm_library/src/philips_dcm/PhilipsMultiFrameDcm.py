'''
This library provides an interface for Philips MR Multiframe dicoms 
through the PhilipsMultiFrameDcm class.

Created on 22 Jan 2013
@date:$Date$
@author: gfagiolo
'''
#python standard libs
from operator import itemgetter
import logging

#libs from the internet
import numpy
import dicom
logging.debug("pydicom %s"%(str(dicom.__version__)))
#print 'pydicom version',dicom.__version__

#custom libs
from philips_dcm.helpers import (is_multiframe,
    TAG_PER_FRAME_FUNCTIONAL_GROUPS_SEQUENCE, get_frame_content,
    TAG_DIMENSION_INDEX_SEQUENCE, DimensionIndexes_to_tagnames, 
    DiffusionFrameInfo,# get_frame_repetion_time, get_repetion_time,
    get_shared_functional_group_sequence,
    get_shared_functional_group_sequence_repetion_time,
    is_mrspectroscopystorage, get_a_frame, is_accepted_dicom)
#from philips_dcm import __version__

FRAME_INDEX_LABEL = 'PhilipsMultiFrameDcm:Frame Index'

class FrameOrientation(object):
    
    """
    
    This object represents the Frame orientation, position and pixel and slice measurements (pixel spacing and slice thickness)
    It can compute an affine transformation from pixel coordinates to Patient Coordinates
    
    """
    
    def __init__(self, a_frame, frame_no, spacing_between_slices=1.):
        self.set_frame_no(frame_no)
        self.__initialise(a_frame, spacing_between_slices)

    def get_last_position(self):
        return self.__LAST_POSITION

    def set_last_position(self, frame):
        self.__LAST_POSITION = numpy.array(frame.PlanePositionSequence[0].ImagePositionPatient,
                             dtype=numpy.float64)

    def get_slice_thickness(self):
        return self.__SLICE_THICKNESS

    def set_slice_thickness(self, value):
        self.__SLICE_THICKNESS = value

    def get_pixel_spacing(self):
        return self.__PIXEL_SPACING

    def set_pixel_spacing(self, value):
        self.__PIXEL_SPACING = value

    def get_depth(self):
        return self.__DEPTH

    def set_depth(self, value):
        self.__DEPTH = value

    def get_vec_row(self):
        return self.__VEC_ROW

    def get_vec_col(self):
        return self.__VEC_COL

    def get_position(self):
        return self.__POSITION

    def set_vec_row(self, value):
        self.__VEC_ROW = value

    def set_vec_col(self, value):
        self.__VEC_COL = value

    def set_position(self, value):
        self.__POSITION = value

    def get_frame_no(self):
        return self.__FRAME_NO

    def set_frame_no(self, value):
        self.__FRAME_NO = value
        
    def __initialise(self, frame, spacing_between_slices):
        orient = numpy.array(frame.PlaneOrientationSequence[0].ImageOrientationPatient,
                             dtype=numpy.float64)
        self.set_vec_row(orient[:3])
        self.set_vec_col(orient[3:])
        depth = numpy.cross(orient[:3], orient[3:])
        depth /= numpy.linalg.norm(depth)
        self.set_depth(depth)
        positi = numpy.array(frame.PlanePositionSequence[0].ImagePositionPatient,
                             dtype=numpy.float64)
        self.set_position(positi)
        self.set_pixel_spacing(numpy.array(frame.PixelMeasuresSequence[0].PixelSpacing, 
                                           dtype=numpy.float64))
        try:
            self.set_slice_thickness(numpy.float64(frame.PixelMeasuresSequence[0].SliceThickness))
        except AttributeError as e:
            self.set_slice_thickness(spacing_between_slices)
            #TODO: log to warning system
#            print 'WARNING: ' + str(RuntimeWarning(str(e) + " Using 'Spacing Between Slices' = %f"%spacing_between_slices))
            msg = str(RuntimeWarning(str(e) + " Using 'Spacing Between Slices' = %f"%spacing_between_slices))
            logging.warning(msg)
            
    def get_position_matrix(self):
        mat = numpy.eye(4)
        mat[:3,3] = self.get_position()
        return mat
    
    def get_rotation_matrix(self):
        mat = numpy.eye(4)
        mat[:3, 0] = self.get_vec_row()
        mat[:3, 1] = self.get_vec_col()
        mat[:3, 2] = self.get_slice_direction()
        return mat
    
    def get_resolution_matrix(self):
        return numpy.diag([self.get_pixel_spacing()[0],
                           self.get_pixel_spacing()[1],
                           self.get_slice_thickness(), 1.])
    
    def get_frame_affine(self):
        out = numpy.dot(self.get_rotation_matrix(), self.get_resolution_matrix())
        out[:3,3] = self.get_position()
        return out
    
    def get_slice_direction(self):
        lp = self.get_last_position()
        if lp is not None:
            sd = lp-self.get_position()
            return sd/numpy.linalg.norm(sd)
        else:
            return self.get_depth()
    
    FRAME_NO = property(get_frame_no, set_frame_no, None, "FRAME_NO's docstring")
    VEC_ROW = property(get_vec_row, set_vec_row, None, "VEC_ROW's docstring")
    VEC_COL = property(get_vec_col, set_vec_col, None, "VEC_COL's docstring")
    POSITION = property(get_position, set_position, None, "POSITION's docstring")
    DEPTH = property(get_depth, set_depth, None, "DEPTH's docstring")
    PIXEL_SPACING = property(get_pixel_spacing, set_pixel_spacing, None, "PIXEL_SPACING's docstring")
    SLICE_THICKNESS = property(get_slice_thickness, set_slice_thickness, None, "SLICE_THICKNESS's docstring")
    LAST_POSITION = property(get_last_position, set_last_position, None, "LAST_POSITION's docstring")

class FrameIndexList(list):
    """
        A list of DimensionIndexes where the last element of each item is the frame number
    """
    
    def __init__(self, *arg):
        super(FrameIndexList, self).__init__(*arg)
        #copy also labels (if arg[0] is a FrameIndexList)
        if len(arg)>0 and isinstance(arg[0], FrameIndexList):
            self.set_indices_labels(arg[0].get_indices_labels())        
    
    def get_frames_indices_only(self):
        """
        returns only the slice number in the current order as a list
        """
        return map(itemgetter(-1), self)
    
    def sorted_slices(self, order):
        """
            order is a tuple/list specifying the indexes to sort  
            returns the sorted FrameIndexList 
        """
        out = FrameIndexList(sorted(self, key=itemgetter(*order)))
        out.set_indices_labels(self.get_indices_labels())
        return out
    
    def select_slices(self, idx, value):
        """
        select only slices for which the field identified by idx is equal to value 
        """
        if isinstance(idx, int):
            out = FrameIndexList(filter(lambda x:x[idx] == value, self))
            out.set_indices_labels(self.get_indices_labels())
            return out
        elif isinstance(idx, str):
            try:
                cidx = self.get_indices_labels().index(idx)
                out = FrameIndexList(filter(lambda x:x[cidx] == value, self))
                out.set_indices_labels(self.get_indices_labels())
                return out
            except ValueError:
                #idx is not in a valid label
                pass
        #if everything failed return empty FrameIndexList
        return FrameIndexList()
    
    def get_slices_excluding_frame_list(self, excluded_frame_list):
        """
        Returns a FrameIndexList excluding excluded_frame_list (a list of frame numbers)
        """
        out = FrameIndexList([x for x in self if x[-1] not in excluded_frame_list])
        out.set_indices_labels(self.get_indices_labels())
        return out

    def get_indices_labels(self):
        return self.__INDICES_LABELS

    def set_indices_labels(self, value):
        self.__INDICES_LABELS = value

    def __str__(self):
        outs = str([])
        if self:
            outs = ''
            #set standard template
            outt = '# '.join(['{%d:5}'%n for n in range(len(self[0]))])
            if self.get_indices_labels():
                #modify template to fit labels
                outt = '# '.join(['{%d:%d}'%(n, max([len(x), 5])) for n, x in enumerate(self.get_indices_labels())])
                outs += outt.format(*self.get_indices_labels()) + '\n'
            for e in self:
                outs += outt.format(*map(str, e)) + '\n'
        return outs 
    
    def toCSV(self):
        outs = ''
        if self:
            outs = ''
            if self.get_indices_labels():
                outs += ','.join(['"%s"'%l for l in self.get_indices_labels()]) + '\n'
            for e in self:
                outs += ','.join(map(str, e)) + '\n'
        return outs 
        
    INDICES_LABELS = property(get_indices_labels, set_indices_labels, None, "Dimension Index Labels")

class IntensityRescaleInfo(object):
    
    """
        Given a frame from PerFrameFunctionalGroup, this class extracts intensity rescale slope and intercept
    """
    
    def __init__(self, frame):
        self.set_slope(numpy.float32(frame.PixelValueTransformationSequence[0].RescaleSlope))
        self.set_intercept(numpy.float32(frame.PixelValueTransformationSequence[0].RescaleIntercept))

    def get_slope(self):
        return self.__SLOPE

    def get_intercept(self):
        return self.__INTERCEPT

    def set_slope(self, value):
        self.__SLOPE = value

    def set_intercept(self, value):
        self.__INTERCEPT = value
    
    SLOPE = property(get_slope, set_slope, None, "Intensity Rescaling Slope")
    INTERCEPT = property(get_intercept, set_intercept, None, "Intensity Rescaling Intercept")
       
    
class MinimalHeaderInfo(object):
        
    def __init__(self, dicom_object):
        self.__initialise(dicom_object)

    def get_shared_functional_group_sequence(self):
        return self.__SHARED_FUNCTIONAL_GROUP_SEQUENCE


    def set_shared_functional_group_sequence(self, value):
        self.__SHARED_FUNCTIONAL_GROUP_SEQUENCE = value


    def get_repetition_time(self):
        return self.__REPETITION_TIME


    def set_repetition_time(self, value):
        self.__REPETITION_TIME = value


    def get_dti_info_list(self):
        return self.__DTI_INFO_LIST

    def set_dti_info_list(self, value):
        self.__DTI_INFO_LIST = value

    def get_dimension_indices_specs(self):
        return self.__DIMENSION_INDICES_SPECS

    def set_dimension_indices_specs(self, value):
        self.__DIMENSION_INDICES_SPECS = value

    def get_intensity_scaling(self):
        return self.__STACK_INTENSITY_SCALING

    def set_intensity_scaling(self, value):
        self.__STACK_INTENSITY_SCALING = value

    def get_dimension_indices_counts(self):
        return self.__DIMENSION_INDICES_COUNTS

    def set_dimension_indices_counts(self, value):
        self.__DIMENSION_INDICES_COUNTS = value

    def get_dimension_indices_labels(self):
        return self.__DIMENSION_INDICES_LABELS

    def set_dimension_indices_labels(self, value):
        self.__DIMENSION_INDICES_LABELS = value

    def get_indices_list(self):
        return self.__INDICES_LIST

    def set_indices_list(self, value):
        self.__INDICES_LIST = value

    def get_no_rows(self):
        return self.__NO_ROWS

    def get_no_cols(self):
        return self.__NO_COLS

    def get_spacing_between_slices(self):
        return self.__SPACING_BETWEEN_SLICES

    def set_no_rows(self, value):
        self.__NO_ROWS = value

    def set_no_cols(self, value):
        self.__NO_COLS = value

    def set_spacing_between_slices(self, value):
        self.__SPACING_BETWEEN_SLICES = value

    def get_index_stack_id(self):
        return self.__INDEX_STACK_ID

    def get_index_in_stack_position_number(self):
        return self.__INDEX_IN_STACK_POSITION_NUMBER

    def set_index_stack_id(self, value):
        self.__INDEX_STACK_ID = value

    def set_index_in_stack_position_number(self, value):
        self.__INDEX_IN_STACK_POSITION_NUMBER = value

    def get_dimension_indices_info(self):
        return self.__DIMENSION_INDICES_INFO

    def set_dimension_indices_info(self, value):
        self.__DIMENSION_INDICES_INFO = value

    def get_stack_orientation(self):
        return self.__STACK_ORIENTATION

    def set_stack_orientation(self, value):
        self.__STACK_ORIENTATION = value
        
    STACK_ORIENTATION = property(get_stack_orientation, set_stack_orientation, None, "STACK_ORIENTATION a dictionary with stack-id as key and Frame orientation as value")
    INTENSITY_SCALING = property(get_intensity_scaling, set_intensity_scaling, None, "INTENSITY_SCALING is a dictionary with a single key None for standard images or MRImageType keys (e.g. magnitude and phase, useful for field maps etc)")
    DIMENSION_INDICES_INFO = property(get_dimension_indices_info, set_dimension_indices_info, None, "DIMENSION_INDICES_INFO's docstring")
    DIMENSION_INDICES_LABELS = property(get_dimension_indices_labels, set_dimension_indices_labels, None, "DIMENSION_INDICES_LABELS's docstring")
    DIMENSION_INDICES_COUNTS = property(get_dimension_indices_counts, set_dimension_indices_counts, None, "DIMENSION_INDICES_COUNTS's docstring")
    DIMENSION_INDICES_SPECS = property(get_dimension_indices_specs, set_dimension_indices_specs, None, "DIMENSION_INDICES_SPECS is a list of length #of dimension indexes, each element is a dictionary whose keys are the unique index values and the relative value is a list containing frame numbers belonging to this particular index value")
    INDEX_STACK_ID = property(get_index_stack_id, set_index_stack_id, None, "INDEX_STACK_ID's docstring")
    INDEX_IN_STACK_POSITION_NUMBER = property(get_index_in_stack_position_number, set_index_in_stack_position_number, None, "INDEX_IN_STACK_POSITION_NUMBER's docstring")
    NO_ROWS = property(get_no_rows, set_no_rows, None, "NO_ROWS's docstring")
    NO_COLS = property(get_no_cols, set_no_cols, None, "NO_COLS's docstring")
    SPACING_BETWEEN_SLICES = property(get_spacing_between_slices, set_spacing_between_slices, None, "SPACING_BETWEEN_SLICES's docstring")
    INDICES_LIST = property(get_indices_list, set_indices_list, None, "INDICES_LIST's docstring")
    DTI_INFO_LIST = property(get_dti_info_list, set_dti_info_list, None, "DTI_INFO_LIST's docstring")
    REPETITION_TIME = property(get_repetition_time, set_repetition_time, None, "REPETITION_TIME's docstring")
    SHARED_FUNCTIONAL_GROUP_SEQUENCE = property(get_shared_functional_group_sequence, set_shared_functional_group_sequence, None, "SHARED_FUNCTIONAL_GROUP_SEQUENCE's docstring")

       
    def __initialise(self, df):
        """
        Loops over Frames to extract relevant information
        @type df: dicom.dataset.DataSet 
        """
        
        self.set_shared_functional_group_sequence(get_shared_functional_group_sequence(df))
        try:
            self.set_repetition_time(df.RepetitionTime)
        except Exception:
            self.set_repetition_time(None)
        if self.get_repetition_time() is None:
            try:
                self.set_repetition_time(get_shared_functional_group_sequence_repetion_time(self.get_shared_functional_group_sequence()))
            except Exception:
                pass
        self.set_no_rows(df.Rows)
        self.set_no_cols(df.Columns)
        self.set_spacing_between_slices(numpy.float64(df.SpacingBetweenSlices))
        idx = df[TAG_DIMENSION_INDEX_SEQUENCE]
        self.set_dimension_indices_info(idx)
        #the next line doesn't seem to work on all multiframes, hence the try statement
        try:
            dim_index_lab = [x.DimensionDescriptionLabel for x in idx]            
        except AttributeError:
            dim_index_lab = map(itemgetter(0), DimensionIndexes_to_tagnames(df))             
        self.set_dimension_indices_labels(dim_index_lab)
        all_frames = FrameIndexList()
        all_frames.set_indices_labels(dim_index_lab + [FRAME_INDEX_LABEL])
        self.set_stack_orientation(dict())
        imagetyppos = self.get_index_image_type_mr()
        self.set_intensity_scaling(dict())
        stackidpos = self.find_index('Stack ID')
        self.set_index_stack_id(stackidpos)
        instackpos = self.find_index('In-Stack Position Number')
        self.set_index_in_stack_position_number(instackpos)
        #count unique indexes per Dimension Index
        conts = list()
        ns = range(len(self.get_dimension_indices_labels()))
        for _ in ns:
            conts.append(dict())
        first_frame_in_stacks = dict()
        last_frame_in_stacks = dict()
        is_dti = self.is_dti()
        self.set_dti_info_list(DiffusionWeightedImageInfo())        
        for fn, frame in enumerate(df[TAG_PER_FRAME_FUNCTIONAL_GROUPS_SEQUENCE]):
            content = get_frame_content(frame)            
            indices = content.DimensionIndexValues
#            if self.get_repetition_time() is None:
#                try:
#                    self.set_repetition_time(get_frame_repetion_time(frame))
#                except Exception:
#                    pass
            if is_dti:
                dti_info = DiffusionFrameInfo(frame, fn)
                self.get_dti_info_list().append(dti_info)
            #get intensity scaling information
            if not imagetyppos is None:
                if not indices[imagetyppos] in self.get_intensity_scaling():
                    self.get_intensity_scaling()[indices[imagetyppos]] = IntensityRescaleInfo(frame)
            elif None not in self.get_intensity_scaling():
                self.get_intensity_scaling()[None] = IntensityRescaleInfo(frame)
            #update dimension index counters
            for i in ns:
                key = indices[i]
                if key in conts[i]:
                    conts[i][key].append(fn)
                else:
                    conts[i][key] = [fn]
            cstack = indices[stackidpos]
            #update stack orientation information
            if not cstack in first_frame_in_stacks:
                first_frame_in_stacks[cstack] = 1000000
                last_frame_in_stacks[cstack] = -1
            #ADD HEADER FOR FIRST FRAME IN STACK
            if indices[instackpos] < first_frame_in_stacks[cstack]:
                #current stack position is the smallest, update first_frame for this stack
                first_frame_in_stacks[cstack] = indices[instackpos]
                #update frame header for this stack
                self.get_stack_orientation()[cstack] = FrameOrientation(frame, fn, numpy.float64(df.SpacingBetweenSlices))
            #UPDATE Patient Position for last frame in stack
            if indices[instackpos] > last_frame_in_stacks[cstack]:
                #current stack position is the largest, update last_frame for this stack
                last_frame_in_stacks[cstack] = indices[instackpos]
                #update last_position for this stack                
                self.get_stack_orientation()[cstack].set_last_position(frame)
            all_frames.append(indices + [fn])
        self.set_indices_list(all_frames)
        self.set_dimension_indices_counts([len(c.keys()) for c in conts])
        self.set_dimension_indices_specs(conts)
    
    def __to_nifti_coord_matrix(self):
        """
        returns an affine matrix to transform dicom orientation axis into nifti
        """
        mat = numpy.eye(4)
        mat[0, 0] = -1.
        mat[1, 1] = -1.
#         mat[0, 3] = self.get_no_rows()-1
#         mat[1, 3] = self.get_no_cols()-1
        return mat
    
## AXIS PERMUTATIONS: ALGEBRA EXPLANATION
##Let's assume there's a permutation from I to I' such that I' = P*I and P is a 4x4 matrix composed of
## a 3x3 matrix with the axis permutations (det = +/-1) and a 3x1 vector representing a shift of
## Na-1 pixels in the a-direction where there's a sign change (e.g. in 1D for a pixel reversal i = N-1-i', since i goes from 0 to N-1)
## since the location pointed by both indexes must be the same:
##X = A*I = A'*I'
##we get that:
## A*I =  A'*P*I
## and that implies that
## A*inv(P) = A'*P*inv(P) = A'
## A' = A*inv(P)
    
    def get_nifti_affine_stack(self, stack_no):
        aff = self.get_affine_stack(stack_no)
#        return numpy.dot(aff, numpy.linalg.inv(self.__to_nifti_coord_matrix()))
#        return numpy.dot(numpy.linalg.inv(self.__to_nifti_coord_matrix()), aff)
        return numpy.dot(self.__to_nifti_coord_matrix(), aff)

    def get_affine_stack(self, stack_no):
        fo = self.get_stack_orientation()[stack_no]
        rot_mat = fo.get_rotation_matrix()
        res_mat = numpy.diag(numpy.hstack((fo.get_pixel_spacing(), 
                                self.get_spacing_between_slices(), numpy.float64(1.))))
        out = numpy.dot(rot_mat, res_mat)
        out[:3,3] = fo.get_position()
        return out
    
    def default_order(self, by_slice_last=False):
        """Analysis of the Dimension Indexes and decision on what order is needed to sort slices
            Stack ID: identifies a stack
            In-Stack Position Number: indentifies a slice in a stack (geometrical)
            Temporal Position Index: identifies a temporal instance
            Diffusion b-value: identifies the current b-value weighting
            Diffusion Gradient Orientation: the current gradient direction
            Effective Echo Time: identifies an echo time for multi-echo acquisitions
            Nominal Cardiac Trigger Delay Time: (possibly) identifies a phase in the cardiac cycle
            Private ImageTypeMR: identifies whether it is a magnitude, real part, imaginary part or a phase image (2005, 1011)
                NB this drives different Intensity Scaling slope and intercept! needs to be handled correctly.
                
            Private Image Plane Number: used in surveys and CSI.2D.T2.dyna.REFERENCE [IGNORED] 
                (not clear how it works, but assigns a progressive number to each slice, even if they belong to different stacks)
            (2001, 100a): same as 'Private Image Plane Number' [IGNORED] found in cardiac.dcm and in dti-16a 
                on http://incenter.medical.philips.com/doclib/getdoc.aspx?func=ll&objid=8162223&objaction=open it states that it is: 
                Slice Number MR'
        """
        #=======================================================================
        # Each frames has got a multidimensional index [dim1, dim2, ..., dimN] describing its location in a multi-dimensional loop.
        # 
        # It would be thus tempting to think of the data order as multidimensional matrix.
        # Something like, rows x columns x size(dim1) x size(dim2) x size(dim3) etc.
        # However that is not generally the case since dim1 item 1 could spawn only 1 item for dim2 
        # and item 2 from dim1 could spawn n items of dim2.
        # This is what happens for instance in dti sequences, where dim1 is bvalue and dim2 is gradient orientation
        # when dim1 = 0 then there's only 1 gradient orientation (i.e. an arbitrary [0,0,0], which is not even a direction)
        # whereas any dim1 != 0 would spawn N gradient orientation (common values for N are 7, 16, 32 or 64).
        # In order to represent the data as a matrix the maximum # of items per dim should be computed and a bigger than
        # the data matrix should be allocated
        #=======================================================================
        corder = range(len(self.get_dimension_indices_labels()))
        ispid = self.get_index_in_stack_position_number()
        #has the image multiple types (i.e. magnitude phase etc)
        imtid = self.get_index_image_type_mr()
        if imtid is None:
            order = list([self.get_index_stack_id(), ispid])
        else:
            #sort first by image type
            order = list([imtid, self.get_index_stack_id(), ispid])
            corder.remove(imtid)
        corder.remove(self.get_index_stack_id())
        corder.remove(ispid)
        #temporal position index
        ti = self.get_index_temporal()
        if not ti is None:
            #sort by time then by space             
            order.append(ti)
            corder.remove(ti)
        #cardiac temporal order
        tc = self.get_index_cardiac_trigger_delay_time()
        if not tc is None:
            #sort by cardiac trigger time then by space             
            order.append(tc)
            corder.remove(tc)                    
        #diffusion information
        
        #remove Private indexes from remaining order index list
        for li, lab in enumerate(self.get_dimension_indices_labels()):
            if lab.lower().find('image plane number')>-1:
                #ignore image plane number
                corder.remove(li)
            elif lab.lower().find('scanning sequence')>-1:
                #ignoring 'Private Scanning Sequence'
                corder.remove(li)
            elif lab.lower().find('(2001, 100a)')>-1:
                #ignoring '(2001, 100a)'
                corder.remove(li)            
        #what is not known needs to come last
        if by_slice_last:
            #remove in stack position from order and add it last to corder
            #for nifti libraries such as pynifti or pynii this gives the correct data ordering for 4d images
            order.remove(ispid)
            corder.append(ispid)
        return order + corder
    
    def sort_slices(self, order=None):
        if order is None:
            order = self.default_order()
        return self.get_indices_list().sorted_slices(order)
    
    def find_index(self, name):
        if name in self.get_dimension_indices_labels():
            return self.get_dimension_indices_labels().index(name)
        else:
            return None
    
    def get_number_of_slices(self):
        """
        returns the number of volumetric slices (not frames)
        """
        return self.get_dimension_indices_counts()[self.get_index_in_stack_position_number()]
    
    def get_echo_numbers(self):
        echid = self.get_index_effective_echo_time()
        if not echid is None: 
            return sorted(self.get_dimension_indices_specs()[echid].keys())
        else:
            return [] 
        
    def get_index_cardiac_trigger_delay_time(self):
        return self.find_index('Nominal Cardiac Trigger Delay Time')
    
    def get_index_temporal(self):
        return self.find_index('Temporal Position Index')

    def get_index_diffusion_bvalue(self):
        if not self.find_index('Diffusion b-Value') is None:
            return self.find_index('Diffusion b-Value')
        else:
            return self.find_index('Diffusion b-value')

    def get_index_diffusion_grad_orient(self):
        return self.find_index('Diffusion Gradient Orientation')

    def get_index_effective_echo_time(self):
        return self.find_index('Effective Echo Time')
    
    def get_index_image_type_mr(self):
        idx = self.find_index('Private ImageTypeMR')
        if idx is None:
            idx = self.find_index('(2005, 1011)')
        return idx
    
    def is_multi_image_type(self):
        return not self.get_index_image_type_mr() is None

    def is_multi_echo(self):
        return not self.get_index_effective_echo_time() is None
    
    def is_dti(self):
        return not self.get_index_diffusion_grad_orient() is None
    
    def format_as_4D(self):
        """
        boolean, True if data is to be reformatted in 4D (fmri/time repeated, cardiac, dti)
        """
        #4d if
        #there is a temporal index
        test = not self.get_index_temporal() is None
        #or a cardiac trigger index
        test = test or not self.get_index_cardiac_trigger_delay_time() is None
        #or it is a dti acquisition with different b-values orientation
        test = test or not self.get_index_diffusion_bvalue() is None
        #or it is a dti acquisition with different gradient orientation
        test = test or not self.get_index_diffusion_grad_orient() is None
        
        return  test
        

class MRSFrameOrientation(FrameOrientation):
    
    def __init__(self, shared_functional_groups, per_frame_functional_groups, slab_thickness):
        orient = numpy.array(shared_functional_groups.PlaneOrientationSequence[0].ImageOrientationPatient,
                             dtype=numpy.float64)
        self.set_vec_row(orient[:3])
        self.set_vec_col(orient[3:])
        depth = numpy.cross(orient[:3], orient[3:])
        depth /= numpy.linalg.norm(depth)
        self.set_depth(depth)
        positi = numpy.array(shared_functional_groups.PlanePositionSequence[0].ImagePositionPatient,
                             dtype=numpy.float64)
        self.set_position(positi)
        self.set_pixel_spacing(numpy.array(per_frame_functional_groups.PixelMeasuresSequence[0].PixelSpacing, 
                                           dtype=numpy.float64))
        self.set_slice_thickness(slab_thickness)
        
    def get_slice_direction(self):
        return self.get_depth()


class MRSpectroscopyMinimalHeaderInfo(MinimalHeaderInfo):
    
    #this only works with single voxel acquisitions
        
    def __init__(self, df):
        #produce 1 frame and its orientation information
        self.set_shared_functional_group_sequence(get_shared_functional_group_sequence(df))

        self.set_no_rows(df.Rows)
        self.set_no_cols(df.Columns)
        #voxel depth is equal to slab thicknes
        self.set_spacing_between_slices(numpy.float64(df.VolumeLocalizationSequence[0].SlabThickness))
        
        self.set_stack_orientation({0:MRSFrameOrientation(self.get_shared_functional_group_sequence(), get_a_frame(df, 0) , self.get_spacing_between_slices())})

    def get_nifti_affine(self):
        return self.get_nifti_affine_stack(0)
    
class PhilipsMultiFrameDcm(object):
    '''
    Wraps the pydicom dataset object and adds Philips EnhancedMRImageStorage header decoding useful for conversion purposes 
    '''
    def __init__(self, filename):
        '''
        filename is a string
        '''
        if is_accepted_dicom(filename):            
            self.__initialise(filename)
        else:
            raise ValueError('%s not a multiframe dicom file'%filename)

    def get_minimal_header(self):
        return self.__MINIMAL_HEADER

    def set_minimal_header(self, value):
        self.__MINIMAL_HEADER = value
        
    def get_pixel_array(self, dtype=numpy.int16):
        "Returns data in int16 format"
        if not dtype is None:
            return self.get_dicom_obj().pixel_array.astype(dtype)
        else:
            return self.get_dicom_obj().pixel_array

    def get_filename(self):
        return self.__FILENAME

    def set_filename(self, value):
        self.__FILENAME = value
        
    def get_dicom_obj(self):
        return self.__DICOM_OBJ

    def set_dicom_obj(self, value):
        self.__DICOM_OBJ = value

    DICOM_OBJ = property(get_dicom_obj, set_dicom_obj, None, "DICOM_OBJ is the dicom represented as a pydicom DataSet object")
    FILENAME = property(get_filename, set_filename, None, "FILENAME is the filename")
    MINIMAL_HEADER = property(get_minimal_header, set_minimal_header, None, "MINIMAL_HEADER is a MinimalHeaderInfo obj containing orientation and other info")    
        
    def __initialise(self, filename):
        df = dicom.read_file(filename)
        if not TAG_PER_FRAME_FUNCTIONAL_GROUPS_SEQUENCE in df:
            raise ValueError('%s is a multiframe dicom file that this software does not support [missing PerFrameFunctionalGroups (5200,9230) dicom tag]'%filename)
        self.set_filename(filename)
        self.set_dicom_obj(df)
        if is_mrspectroscopystorage(df):
            self.set_minimal_header(MRSpectroscopyMinimalHeaderInfo(df))
        else:
            self.set_minimal_header(MinimalHeaderInfo(df))
    
    def get_data(self):
        return self.get_dicom_obj().pixel_array
    
    def get_ordered_data(self, order=None):
        orsl = self.get_minimal_header().sort_slices(order).get_frames_indices_only()
        return self.get_data()[orsl, :, :]
    
    def get_sequence_custom_description_txt(self):
        try:
            pn = self.get_dicom_obj().ProtocolName
        except Exception:
            pn = 'UNKNOWN'
        sn = self.get_dicom_obj().SeriesNumber 
        return '_'.join([str(sn), pn]) 
    
    def is_mrs(self):
        return isinstance(self.get_minimal_header(), MRSpectroscopyMinimalHeaderInfo)
            
class DiffusionWeightedImageInfo(list):
    
    def get_only_acquired_data(self):
        """only select frames that were acquired (i.e. not processed/reconstructed data such as apparent diffusion image)"""
        return DiffusionWeightedImageInfo(filter(lambda x:x.isDirectional() or x.isUnweighted(), self))
    
    def get_processed_data(self):
        """only select frames that were acquired (i.e. not processed/reconstructed data such as apparent diffusion image)"""
        return DiffusionWeightedImageInfo(filter(lambda x:not x.isDirectional() and not x.isUnweighted(), self))
    
    def get_bvals(self):
        return numpy.array(map(lambda x:x.get_diff_bvalue(), self.get_only_acquired_data()))
    
    def get_bvecs(self):
        return numpy.array(map(lambda x:x.get_diff_bvec(), self.get_only_acquired_data()))
    
    def save_FSL_bval_bvecs(self, selected_frames, dti_idx, bvalname, bvecsname):
        """
            Saves bvals and bvecs in FSL format (in the correct volume order) 
            selected_frames is a FrameIndexList object
            dti_idx is (MinimalHeaderInfo().get_index_diffusion_bvalue(), MinimalHeaderInfo().get_index_diffusion_grad_orient())            
        """
        def nifti_transform(v):
            return [-v[0], -v[1], v[2]]
        bvals = list()
        bvecs = list()
        uniqSet = set()
        for x in self.get_only_acquired_data():
            
            try:
                idxs = filter(lambda z:z[-1]==x.get_frame_no(), selected_frames)[0]
            except IndexError:
                #current frame from get_only_acquired_data not in selected frames 
                continue
            cidx = '%d %d'%(idxs[dti_idx[0]], idxs[dti_idx[1]])
            if not cidx in uniqSet:
                #only get unique bvalues/bvecs
                bvals.append(x.get_diff_bvalue())
                grad = x.get_diff_bvec() if not x.get_diff_bvec() is None else [0., 0., 0.]
                bvecs.append(nifti_transform(grad))
                uniqSet.add(cidx)
        numpy.savetxt(bvalname, numpy.array(bvals).reshape(len(bvals), 1).T,fmt='%g')
        numpy.savetxt(bvecsname, numpy.array(bvecs).T,fmt='%g')
