import matlab.engine
import pandas as pd
from ..pacemapping import PaceMaps
from ..processing import output_name, interpolate_map
import numpy as np

"""
Test script to compute conventional and reference-less pace-maps between reference VT signals and paced signals
"""

# INITIALISE FILES
ref_file = '/media/sm18/Seagate Backup Plus Drive/PhD/SBRT_Patient1/post_processed/ecg/VT_PVC2_VT_MULTI_BEAT_ecg.csv'
pts_pace_file = '/media/sm18/Seagate Backup Plus Drive/PhD/SBRT_Patient1/pacing/final_clinical_pacing.pts'
pace_file = '/media/sm18/Seagate Backup Plus Drive/PhD/SBRT_Patient1/pacing/COMBINED_ECG_paces.csv'
type_sig = 'ecg'
folder_out = './output_files/'
mesh_file = '/media/sm18/Seagate\ Backup\ Plus\ Drive/PhD/SBRT_Patient1/Mesh/myo_final'
flag = 1 # mean of 10/12
cycle = 40

########################################################################################################################
## TESTING MATLAB IMPLEMENTATION OF ALIGN AND CORRELATE
########################################################################################################################

# Computing conventional pace-map
outfile_conv = output_name(folder=folder_out, type_map = 'conventional', type_sig = type_sig, flag = flag)
eng = matlab.engine.start_matlab()
eng.addpath('../matlab/')
corr = eng.clinical_pacemapping(ref_file,pace_file,flag,outfile_conv + '.dat',type_sig,cycle)
interpolate_map(idat = outfile_conv + '.dat',odat = outfile_conv + '_myo.dat',mesh = mesh_file,pts = pts_pace_file.replace(' ','\ '), flag=flag)

# Computing reference-less pace-map
outfile_less = output_name(folder=folder_out, type_map = 'reference_less', type_sig = type_sig, flag = flag)
corr_less = eng.clinical_gradient_pacemapping(pace_file,pts_pace_file,outfile_less + '.dat',type_sig,cycle)
new_pts_pace_file = outfile_less + '.pts'
interpolate_map(idat = outfile_less + '.dat',odat = outfile_less + '_myo.dat',mesh = mesh_file,pts = new_pts_pace_file.replace(' ','\ '), flag=flag)

########################################################################################################################
## TESTING PYTHON IMPLEMENTATION OF ALIGN AND CORRELATE
########################################################################################################################
"""
# Read files
ref = pd.read_csv(ref_file, header=None).values
pace = pd.read_csv(pace_file, header=None).values
pts_pace = pd.read_csv(pts_pace_file, header=None, skiprows=1, delimiter=' ').values

# Downsample reference file to match timestep of simulated pace file (usually clinical ECGs/EGMs are acquired every millisecond, whereas simulations every 10 milliseconds)
ref = ref[:,::10]

# Check shape signals
print('Shape pace: ({},{})\nShape ref: ({},{})\nShape pace pts: ({},{})'.format(pace.shape[0],pace.shape[1],ref.shape[0],ref.shape[1],pts_pace.shape[0],pts_pace.shape[1]))

# Creating a PaceMap object
test_pace = PaceMaps(ref = ref,pace = pace,type_sig = type_sig, pts = pts_pace, flag = flag, cycle = cycle)

#Test aling and correlate for just one pace
mean_corr = test_pace.aling_and_correlate(ref,pace,flag = flag, cycle = cycle)
print('Mean correlation is: {} '.format(mean_corr*100))

# Computing conventional correlation map
corr, stim = test_pace.conventional()
outfile_conv = output_name(folder=folder_out, type_map = 'conventional', type_sig = type_sig, flag = flag)

print('Printing max and min correlation: {} and {}'.format(max(corr),min(corr)))

print('Writing out {}.dat pacemap'.format(outfile_conv))
np.savetxt(outfile_conv + '.dat',corr, fmt='%.2f')

print('Writing out {}.csv new stim file'.format(outfile_conv))
with open(outfile_conv + '.pts', 'w') as f:
    f.write('{}\n'.format(stim.shape[0]))
    np.savetxt(f,stim, delimiter=' ', fmt='%.3f %.3f %.3f')

# Interpolating pace-maps into mesh file of interest and convert to vtk
interpolate_map(idat = outfile_conv + '.dat',odat = outfile_conv + '_myo.dat',mesh = mesh_file,pts = outfile_conv + '.pts', flag=flag)

# Computing reference-less correlation map
corr_less,new_stim = test_pace.referenceless()
outfile_less = output_name(folder=folder_out, type_map = 'reference-less', type_sig = type_sig, flag = flag)

print('Printing max and min reference-less correlation: {} and {}'.format(max(corr_less),min(corr_less)))

print('Writing out {}.dat pacemap'.format(outfile_less))
np.savetxt(outfile_less+'.dat',corr_less, fmt='%.2f')

print('Writing out {}.csv new stim file'.format(outfile_less))
with open(outfile_less + '.pts', 'w') as f:
    f.write('{}\n'.format(new_stim.shape[0]))
    np.savetxt(f,new_stim, delimiter=' ', fmt='%.3f %.3f %.3f')

# Interpolating pace-maps into mesh file of interest and convert to vtk
interpolate_map(idat = outfile_less + '.dat',odat = outfile_less + '_myo.dat',mesh = mesh_file,pts = outfile_less + '.pts', flag=flag)
"""
