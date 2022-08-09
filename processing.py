"""
Pre- and/or post processing function for in-silico pacemapping

This module provides methods to compute generate .csv paced file from CARP LF simulations, interpolate correlation maps from point stimuli to mesh of interest,
generate vtk files, generate output names for correlation files

"""

import numpy as np
import os

# Function to generate output file for correlation map
def output_name(folder, type_map, type_sig,flag=1):

    """
    Generate name for output file

    :param folder:          (str) name of output filder
    :param type_map:        (str) type of map (conventional, reference-less)
    :param type_sig:        (str) type of signal (ecg, egm)
    :param flag:            (int) flag to compute: (0) mean correlation of ALL leads; (1) mean correlation of highest (ALL leads - 2) ---> default: 1

    :return:                (str) name of output file
    """

    if flag==0:
        extra = 'mean_12_12'
    elif flag==1:
        extra = 'mean_10_12'
    elif flag==2:
        extra = 'mean_2'

    output = folder + type_map + '_' + type_sig + '_' + extra
    print('Output file (no extension): {}\n'.format(output))

    return output


# Function to interpolate correlation map from point stimuli to entire mesh + convert it to vtk
# THIS FUNCTION CAN ONLY BE EXEXCUTED IF MESHTOOL and VTKGlConvert ARE IN THE MAIN PATH
def interpolate_map(idat,odat,mesh,pts, flag=0):

    """

    :param idat:        (str) filename of correlation map .dat to interpolate (with .dat extension)
    :param odat:        (str) filename of output correlation map (with .dat extension)
    :param mesh:        (str) filename of mesh
    :param pts:         (str) filename of points to interpolate from
    :param flag:        (bool) open vtk file in paraview (1)

    """

    print('This function can only be executed if meshtool and VTKGlConvert are in main path !!')

    command_interpolate = 'meshtool interpolate clouddata -omsh={} -idat={} -odat={} -pts={}'.format(mesh,idat,odat,pts)
    print(command_interpolate)
    os.system(command_interpolate)

    command_convert_vtk = 'GlVTKConvert -m {} -n {} -o {} -F "bin"'.format(mesh,odat,odat[:-4])
    print(command_convert_vtk)
    os.system(command_convert_vtk)

