
import scipy.io as sio
import numpy as np
import os
import re
import struct

### This file will contain the test and setup portion
## setup algoParams
gamma_Hz_per_Tesla = 42.577481e6;
#species(1).name = 'water';
#species(1).frequency = 4.70;
#species(1).relAmps = 1;
#species(2).name = 'fat';
#species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29]; # 9-peak model
#species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37]; # Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
decoupled_estimation = True; # flag for decoupled R2 estimation
Fibonacci_search = True; # flag for Fibonacci search
B0_smooth_in_stack_direction = False; # flag for B0 smooth in stack direction
multigrid = True; # flag for multi-level resolution pyramid
estimate_R2 = True; # flag to estimate R2star
verbose = True# flag for verbose status messages (default false)
process_in_3D = True; # flag to process in 3D (default True)
R2star_calibration = True; # flag to perform R2* calibration (default False)
ICM_iterations = 2; # ICM iterations
F = 100; # number of discretized B0 values
mu = 10; # regularization parameter
R2_stepsize = 1; # R2 stepsize in s^-1
max_R2 = 120; # maximum R2 in s^-1
max_label_change = 0.1; #
fine_R2_stepsize = 1.0; # Fine stepsize for discretization of R2(*) [s^-1] (used in decoupled fine-tuning step)
coarse_R2_stepsize = 10.0; # Coarse stepsize for discretization of R2(*) [s^-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only
water_R2 = 0.0; # Water R2 [sec-1]
fat_R2s = np.zeros((1,9)); # fat peak R2s [s^-1]
R2star_calibration_max = 800; # max R2* to use for calibration [s^-1] (default 800)
R2star_calibration_cdf_threshold = 0.98; ## threshold for R2* calibration cumulative density function [0,1] (default 0.98)

i = 0
path = '/Users/joseferpaz/Documents/Vanderbilt Internship/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases' #HARDCODED (needs to be replaced eventually)
file_list = os.listdir(path)

edit_path = path + '/' + file_list[i]
mat = sio.loadmat(edit_path)

num_array = mat['imDataParams']['images']
echo_times = mat['imDataParams']['TE']
Field_Strength = mat['imDataParams']['FieldStrength']
#F_Strength = np.squeeze(Field_strength)

#TEDIOUS and varies for some reason
#images
image_imDataParams = num_array[0,0] #HARDCODED replace this as well.. somehow

#echo times (TE)
echo_times = echo_times[0,0]
TE_imDataParams = np.squeeze(echo_times)

#Field Strength
#FieldStrength_imDataParams = np.squeeze(Field_Strength)

        ##getting dimensions
nx, ny, nz, nc, ne = np.shape(image_imDataParams)
print "The dimensions are:" , np.shape(image_imDataParams) #not sure if necessary since there is a loop down there



validParams = 0

# if 'myVar' in globals():  ???? USefule maybe? :D

#algoparams - verbose
if 'verbose' not in locals():
    verbose = False
else:
    if (verbose == 1 or verbose == True):
        verbose = True
    else:
        verbose = False

print "Verbose is :", verbose

#imDataParams - image
#if the images are not set as a variable, shoot a warningn
#otherwise if it isnt five dimensions, ask for them, but if there is 5 then we gucci

if not image_imDataParams.any():
    print 'Warning: (fw_i3cm1i_3pluspoint_berglund_QPBO: input must contain field ''images'' with size [nx ny nz ncoils nTE])'
elif len(np.shape(image_imDataParams)) != 5:
    print 'Warning: (fw_i3cm1i_3pluspoint_berglund_QPBO: imDataParams.images should have dimensions [nx ny nz ncoils nTE])' # also produces error for nTE=1 <-- Maybe not in python
else:
    [nx, ny, nz, ncoils, nTE] = np.shape(image_imDataParams)
    print nx, ny, nz, ncoils, nTE

#imDataParams - TE
#This one has to do with Tesla n stuffVC
#yell at it if there is not imdataparams.TE
#if there is, make sure the numer of TE matches the last dimension (nTE)

if not TE_imDataParams.any():
    print 'Warning( fw_i3cm1i_3pluspoint_berglund_QPBO: input must contain field ''TE'' with number of elements nTE)'
elif len(TE_imDataParams) != nTE:
    print 'Warning: fw_i3cm1i_3pluspoint_berglund_QPBO; TE_imDataParams has length %d (expected %d)' % (len(TE_imDataParams), ne)
print "TE parameters has length %d (expected %d)" % (len(TE_imDataParams), ne)

#imDataParams - FieldStrength
if 'FieldStrength_imDataParams' not in locals():
    FieldStrength_imDataParams = 1.5
    print "Warning: defaulting FieldStrength_imDataParams to 1.5"
else:
    print "It's in here."





#elif 'verbose' in globals():
#    print ' nein hier!'
#
#elif 'verbose' not in locals():
#    print ' the var isnt here'

 #   if not there:
