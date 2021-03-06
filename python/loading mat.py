import scipy.io as sio
import numpy as np
import os
import re
#import AlgoParams

############### Dixon app is really just a huge function ###########################

#########################################################################################
###########################  Defining variables
#########################################################################################
i = 0
path = '/Users/joseferpaz/Documents/Vanderbilt Internship/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases'
file_list = os.listdir(path)
edit_path = path + '/' + file_list[i]  #This is a list of all of the .mat files
vtk_max = 2^14-1

decoupled_estimation = true # flag for decoupled R2 estimation
Fibonacci_search = true # flag for Fibonacci search
B0_smooth_in_stack_direction = false # flag for B0 smooth in stack direction
multigrid = true # flag for multi-level resolution pyramid
estimate_R2 = true # flag to estimate R2star
verbose = true # flag for verbose status messages (default false)
R2_stepsize = 1
mu = 10
max_R2 = 120; # maximum R2 in s^-1
max_label_change = 0.1; #
fine_R2_stepsize = 1.0; # Fine stepsize for discretization of R2(*) [s^-1] (used in decoupled fine-tuning step)
coarse_R2_stepsize = 10.0; # Coarse stepsize for discretization of R2(*) [s^-1] (used for joint estimation step, value 0.0 --> Decoupled estimation only
water_R2 = 0.0; # Water R2 [sec-1]

#sub_num1 = num_array[0]
#image_array = num_array[0][0] # <---- This is kinda hardcoded so in the future find an alternative to it.
#sub_num3 = num_array[0][0][0]
#sub_num4 = num_array[0][0][0][0]
#sub_num5 = num_array[0][0][0][0][0]

DixonApp_parent_folder = os.path.dirname(os.path.abspath(__file__))
DixonApp_vtk_folder = DixonApp_parent_folder + "/vtk"
DixonApp_vtk_prefix_input  = 'DixonApp_Input'
DixonApp_vtk_prefix_output = 'DixonApp_Output'


#########################################################################################
################### Functions!
#########################################################################################

def write_VTK(vol,vtkfile):

    #where vol is the 3D matrix... but is it really?
    #it isn't but we ignore the rest because we don't need it
    sz = np.shape(vol)
    X = sz[1]
    Y = sz[2]
    Z = 1

    if (len(sz) == 3):
        Z = sz[3]

        # where s is a string


    fid = open(vtkfile,'a+')

    #fid = open(vtkfile....)???????!!!!!####,'w','b')
    ## Well, apparently it just tells the system that you'll be writting the strings
    ## to come in binary

    header = """
    # vtk DataFile Version 3.0);
    created by writeVTK (Matlab implementation by Erik Vidholm)
    BINARY
    DATASET STRUCTURED_POINTS
    """

    header += '%s%d%c%d%c%d\n' % ('DIMENSIONS ', X, ' ', Y, ' ', Z);
    header += '%s%f%c%f%c%f\n' % ('ORIGIN ', 0.0, ' ', 0.0, ' ', 0.0);
    header += '%s%f%c%f%c%f\n' % ('SPACING ', 1.0, ' ', 1.0, ' ', 1.0);
    header += '%s%d\n' % ('POINT_DATA ', X*Y*Z);

    vol_dataType = vol.dtype

    if vol_dataType == uint8 :
        header += 'SCALARS image_data unsigned_char\n'
    elif vol_dataType == int8 :
        header += 'SCALARS image_data char\n'
    elif vol_dataType == uint16 :
        header += 'SCALARS image_data unsigned_short\n'
    elif vol_dataType == int16 :
        header += 'SCALARS image_data short\n'
    elif vol_dataType == uint32 :
        header += 'SCALARS image_data unsigned_int\n'
    elif vol_dataType == int32 :
      header += 'SCALARS image_data int\n'
    elif vol_dataType == float32 :
      header += 'SCALARS image_data float\n'
    elif vol_dataType == float64 :
      header += 'SCALARS image_data double\n'

    header += 'LOOKUP_TABLE default\n' #Make sure to print 'header' to check

    header.encode('utf-16BE')
    fid.write(header)
    fid.close()

def read_VTK(vtkfile):

    fid = open(vtkfile, "r")

    fid.readline()  # vtk DataFile Version x.x
    fid.readline() # comments
    fid.readline() # BINARY
    fid.readline() # DATASET STRUCTURED_POINTS

    z = fid.readline() # DIMENSIONS NX NY NZ   #5

    #sz = sscanf(s, '%*s%d%d%d')
    sz = map(int, re.findall(r'\d+', z))    #5

    ###THE EXAMPLE TEST WAS READING LINES JUST Fine
    ### KEEP AN EYE ON THE BINARY stuffVC


    fid.readline() # ORIGIN OX OY OZ
    fid.readline() # SPACING SX SY SZ
    fid.readline() # POINT_DATA NXNYNZ


    s = fid.readline() # SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
    svstr = sscanf(s, '%s', 1);
    dtstr = sscanf(s, '%*s%*s%s');





#files in the directory
print file_list     #Just double checking...

#THIS IS READY TO DO EVERY FILE INDIVIDUALLY
#for filename in file_list:

#########################################################################################
############## Calling the information from the mat file
#########################################################################################

mat = sio.loadmat(edit_path)

num_array = mat['imDataParams']['images']
echo_times = mat['imDataParams']['TE']

image_array = num_array[0][0]

        ##getting dimensions
nx, ny, nz, nc, ne = np.shape(image_array)
print "The dimensions are:" , np.shape(image_array)

        ##get maximum values
max_real = np.amax(abs(np.real(image_array)))
imag_real = np.amax(abs(np.imag(image_array)))

if max_real > imag_real:
    images_max = max_real
else:
    images_max = imag_real

print "This is the MAX of the real # : " , max_real
print "This is the MAX of the imaginary # : " , imag_real
print "This is the MAX out of both # : " , images_max
print "-" * 100
print "These are the echo times: " , echo_times

#Attempting to get the values that are ultimately written in the file ----NOT FINISHED
for idx_echo in range(1,ne):

    image_re = (np.real(image_array[:,:,:,:,idx_echo])*(vtk_max/images_max))
    image_im = (np.imag(image_array[:,:,:,:,idx_echo])*(vtk_max/images_max))

    #vtk file names that have to be part of a loop
    vtk_filename_re = '%s/%s_echo%04d_re.vtk' % (DixonApp_vtk_folder, DixonApp_vtk_prefix_input, idx_echo)
    vtk_filename_im = '%s/%s_echo%04d_im.vtk' % (DixonApp_vtk_folder, DixonApp_vtk_prefix_input, idx_echo)

    write_VTK(image_re,vtk_filename_re)
    write_VTK(image_im,vtk_filename_im)

    #Now as part of the loop, write the vtk files!

    ## SO (one for the real and one for the imaginary ) for every iteration


#########################################################################################
########################### Making THE string
#########################################################################################

### Dixon stuff is being defined where the rest of the variables are

###This will have to be in a case statement because it depends on the operating system
path_to_DixonApp_executable = DixonApp_parent_folder + "/MAC/DixonApp_MAC.exe"
path_to_DixonApp_libraries = DixonApp_parent_folder + "/MAC/lib"
system_command_string_prefix = 'export DYLD_LIBRARY_PATH="%s"; "%s"' % (path_to_DixonApp_libraries, path_to_DixonApp_executable)

#from here down in the making string section it doesn't work yet
system_command_string_prefix = '%s "%s"' % (system_command_string_prefix, DixonApp_vtk_folder)
system_command_string_prefix = '%s "%s"' % (system_command_string_prefix, DixonApp_vtk_prefix_input)
system_command_string_prefix = '%s "%s"' % (system_command_string_prefix, DixonApp_vtk_prefix_output)

print system_command_string_prefix

system_command_string = system_command_string_prefix;

# flags (6)
system_command_string = '%s "%d"' % (system_command_string, verbose);
system_command_string = '%s "%d"' % (system_command_string, estimate_R2);
system_command_string = '%s "%d"' % (system_command_string, decoupled_estimation);
system_command_string = '%s "%d"' % (system_command_string, Fibonacci_search);
system_command_string = '%s "%d"' % (system_command_string, multigrid);
system_command_string = '%s "%d"' % (system_command_string, B0_smooth_in_stack_direction);

# floats (11)

F = mat['imDataParams']['FieldStrength']; FS = F[0][0]; ## temporary fix
system_command_string = '%s "%.5e"' % (system_command_string, FS);

#system_command_string = sprintf('%s "%.5e"', system_command_string, input.te_used);
#system_command_string = sprintf('%s "%.5e"', system_command_string, input.dte_used);
#system_command_string = sprintf('%s "%.5e"', system_command_string, input.water_chemical_shift);
system_command_string = sprintf('%s "%.5e"', system_command_string, mu);
system_command_string = sprintf('%s "%.5e"', system_command_string, R2_stepsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, max_R2);
system_command_string = sprintf('%s "%.5e"', system_command_string, max_label_change);
#ystem_command_string = sprintf('%s "%.5e"', system_command_string, input.InplaneOverThroughplaneVoxelsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, fine_R2_stepsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, coarse_R2_stepsize);
system_command_string = sprintf('%s "%.5e"', system_command_string, water_R2);

#DixonQPBO_input.verbose = algoParams.verbose;
#DixonQPBO_input.decoupled_estimation = algoParams.decoupled_estimation;
#DixonQPBO_input.Fibonacci_search = algoParams.Fibonacci_search;
#DixonQPBO_input.B0_smooth_in_stack_direction = algoParams.B0_smooth_in_stack_direction;
#DixonQPBO_input.multigrid = algoParams.multigrid;

#if ( (algoParams.estimate_R2) && (imDataParams.nTEeven>=4) ),
#    DixonQPBO_input.estimate_R2 = true;
#else
#    DixonQPBO_input.estimate_R2 = false;
#


#for i in range(n): ###DONT NEED A LOOP







### stuff #############
    #print mat['imDataParams']['TE']
    ##space()
    #print mat['imDataParams']['FieldStrength']
    #space()
    #print mat['imDataParams']['PrecessionIsClockwise']
    #space()
#################################
