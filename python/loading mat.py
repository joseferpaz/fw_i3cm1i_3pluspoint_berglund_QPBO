import scipy.io as sio
import numpy as np
import os
import AlgoParams

#########################################################################################
###########################  Defining variables
#########################################################################################

i = 0
path = '/Users/joseferpaz/Documents/Vanderbilt Internship/fw_i3cm1i_3pluspoint_berglund_QPBO/test_cases'
file_list = os.listdir(path)    #This is a list of all of the .mat files
vtk_max = 2^14-1

sub_num1 = num_array[0]
sub_num2 = num_array[0][0]
sub_num3 = num_array[0][0][0]
sub_num4 = num_array[0][0][0][0]
sub_num5 = num_array[0][0][0][0][0]


#files in the directory
print file_list     #Just double checking...

#THIS IS READY TO DO EVERY FILE INDIVIDUALLY
#for filename in file_list:
edit_path = path + '/' + file_list[i]


#########################################################################################
############## Calling the information from the mat file
#########################################################################################

mat = sio.loadmat(edit_path)

num_array = mat['imDataParams']['images']
echo_times = mat['imDataParams']['TE']

        ##getting dimensions
nx, ny, nz, ne = np.shape(sub_num2)
print "The dimensions are:" , np.shape(sub_num2)

        ##get maximum values
max_real = np.amax(abs(np.real(sub_num2)))
imag_real = np.amax(abs(np.imag(sub_num2)))

if max_real > imag_real:
    return images_max = max_real
else:
    return images_max = imag_real

print "This is the MAX of the real # : " , max_real
print "This is the MAX of the imaginary # : " , imag_real
print "This is the MAX out of both # : " , images_max
print "-" * 100
print "These are the echo times: " , echotimes

#Attempting to get the values that are ultimately written in the file ----NOT FINISHED
for i in range(1,ne):
    image_re = (np.real(sub_num2[:,:,:,i])*(vtk_max/images_max)
    image_im = (np.imag(sub_num2[:,:,:,i])*(vtk_max/images_max)

    #vtk file names that have to be part of a loop
    vtk_filename_re = '%s/%s_echo%04d_re.vtk' % (DixonApp_vtk_folder, DixonApp_vtk_prefix_input, e);
    vtk_filename_im = '%s/%s_echo%04d_im.vtk' % (DixonApp_vtk_folder, DixonApp_vtk_prefix_input, e);

    #Now as part of the loop, write the vtk files,
    ## SO (one for the real and one for the imaginary ) for every iteration


#########################################################################################
########################### Making THE string
#########################################################################################

DixonApp_parent_folder = os.path.dirname(os.path.abspath(__file__))
DixonApp_vtk_folder = DixonApp_parent_folder + "/vtk"
DixonApp_vtk_prefix_input  = 'DixonApp_Input'
DixonApp_vtk_prefix_output = 'DixonApp_Output'

###This will have to be in a case statement because it depends on the operating system
path_to_DixonApp_executable = DixonApp_parent_folder + "/MAC/DixonApp_MAC.exe"
path_to_DixonApp_libraries = DixonApp_parent_folder + "/MAC/lib"
system_command_string_prefix = 'export DYLD_LIBRARY_PATH="%s"; "%s"' % (path_to_DixonApp_libraries, path_to_DixonApp_executable)

system_command_string_prefix = '%s "%s"' % (system_command_string_prefix, DixonApp_vtk_folder)
system_command_string_prefix = '%s "%s"' % (system_command_string_prefix, DixonApp_vtk_prefix_input)
system_command_string_prefix = '%s "%s"' % (system_command_string_prefix, DixonApp_vtk_prefix_output)


system_command_string = '%s "%d"' % (system_command_string, AlgoParams.verbose)
system_command_string = '%s "%d"' % (system_command_string, AlgoParams.estimate_R2)
system_command_string = '%s "%d"' % (system_command_string, AlgoParams.decoupled_estimation)
system_command_string = '%s "%d"' % (system_command_string, AlgoParams.Fibonacci_search)
system_command_string = '%s "%d"' % (system_command_string, AlgoParams.multigrid)
system_command_string = '%s "%d"' % (system_command_string, AlgoParams.B0_smooth_in_stack_direction);

#for i in range(n): ###DONT NEED A LOOP







### stuff #############
    #print mat['imDataParams']['TE']
    ##space()
    #print mat['imDataParams']['FieldStrength']
    #space()
    #print mat['imDataParams']['PrecessionIsClockwise']
    #space()
#################################



def write_VTK(vol,vtkfile):

    #where vol is the 3D matrix... but is it really?
    sz = np.shape(vol)
    X = sz[1]
    Y = sz[2]
    Z = 1

    if (sz.len == 3):
        return Z = sz[3]

    fid = open(vtkfile....)???????!!!!!####,'w','b')

    fprintf(fid, '%s\n' % ('# vtk DataFile Version 3.0');
    fprintf(fid, '%s\n' % ('created by writeVTK (Matlab implementation by Erik Vidholm)');
    fprintf(fid, '%s\n' % ('BINARY');
    fprintf(fid, '%s\n' % ('DATASET STRUCTURED_POINTS');
    fprintf(fid, '%s%d%c%d%c%d\n' % ('DIMENSIONS ', X, ' ', Y, ' ', Z);
    fprintf(fid, '%s%f%c%f%c%f\n' % ('ORIGIN ', 0.0, ' ', 0.0, ' ', 0.0);
    fprintf(fid, '%s%f%c%f%c%f\n' % ('SPACING ', 1.0, ' ', 1.0, ' ', 1.0);
    fprintf(fid, '%s%d\n' % ('POINT_DATA ', X*Y*Z);
