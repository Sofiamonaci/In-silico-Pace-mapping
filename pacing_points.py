#!/usr/bin/env python
# This function computes equally distanced pacing points in specific tags of a mesh, and in specific AHA segments --> can be used for pace-mapping or cnn pacing training/testing
#
#
# Sofia Monaci 

import os
import sys
import argparse
import numpy as np
import pandas as pd
import logging
import open3d as o3d
import matlab.engine

# Importing logger
logging.basicConfig(filename = 'output.log', format='%[levelname]s -> %(lineno)s: %(message)s', level = logging.DEBUG)


def main(args):

    # Read mesh points
    print('Reading %s.pts ...'%(args.mesh))
    pts = pd.read_csv(args.mesh + '.pts', header=None, skiprows=1, delimiter=' ').values

    # Read elem points
    print('Reading %s.elem ...'%(args.mesh))
    elem = pd.read_csv(args.mesh + '.elem', header=None, skiprows=1, delimiter=' ', usecols=[1,2,3,4,5]).values

    # Read segment-model file
    print('Reading %s ...'%(args.aha))
    aha = pd.read_csv(args.aha, header=None, delimiter=' ').values[:,0]

    # Check for consistency in pts/elem/aha dimensions
    logging.debug('Checking sizes of pts/elem/aha arrays ...')
    if not np.max(elem[:,:-2].flatten())==pts.shape[0]-1:
        logging.error('Max index in elem does not match total n of pts!')
    if not len(aha)==pts.shape[0]:
        logging.error('Aha length does not match total n of pts!')

    # Find nodes equal to tag of interest
    try:
        nodes = elem[elem[:,-1]==args.tag,:-2].flatten()
        new_pts = pts[np.unique(nodes),:]
    except IndexError:
        logging.error('There are no elements equal to {0}!'.format(args.tag))


    # Initialise output points
    out_pts = []
    for i in args.segs:
        # Find nodes equal to segment
        pts_seg = pts[aha==i,:]
        # Create point cloud from pts_seg
        ptsCloud = o3d.geometry.PointCloud()
        ptsCloud.points = o3d.utility.Vector3dVector(pts_seg)
        # Downsample according to distance provided
        out_ptCloud = ptsCloud.voxel_down_sample(voxel_size=args.dist)
        # Append to output points
        out_pts.append(np.asarray(out_ptCloud.points))

    # Find closest mesh points to the grid points (out_pts)
    try:
        # Write out grid points in format compatible for meshtool query idxlist
        print('Writing out cloud points ...')
        out_pts = np.concatenate((out_pts))
        with open(args.out, 'w') as f:
            f.write('%d\n'%(out_pts.shape[0]))
            for i in range(out_pts.shape[0]):
                f.write('%.3f %.3f %.3f %d\n'%(out_pts[i,0],out_pts[i,1],out_pts[i,2],args.res))
        # Use meshtool query idxlist to find mesh points
        os.system('meshtool query idxlist -msh='+args.mesh+' -coord='+args.out)
        # Read just found indices
        print('Reading indices')
        ind_out = pd.read_csv(args.out + '.out.txt',header=None, skiprows=1).values[:,0]
        final_pts = pts[ind_out,:]
        # Print out mesh points in .pts format
        print('Printing out final pacing points in %s ...\n'%(args.out))
        with open(args.out, 'w') as f:
            f.write('%d\n'%(out_pts.shape[0]))
            for i in range(out_pts.shape[0]):
                f.write('%.3f %.3f %.3f\n'%(out_pts[i,0],out_pts[i,1],out_pts[i,2]))
    except IndexError:
        logging.error('Empty output points')

# Defining arguments
if __name__ == '__main__':

    parser = argparse.ArgumentParser('Generating equally distanced pacing points for Pace-Mapping (or CNN pacing training)')
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--mesh', type = str, 
            default = "/home/sm18/syncdir/CT_data/Torso-Female-Segmentation/myo_final", 
            help = 'Path to mesh points/elem')
    
    parser.add_argument('--aha' , type = str, 
            default = "/home/sm18/syncdir/CT_data/Torso-Female-Segmentation/geom_17segs_myo.dat", 
            help = '17-segment aha file')
    
    parser.add_argument('--dist', type = int, 
            default = 800, 
            help = 'Distance in micrometers between pacing points')
    
    parser.add_argument('--out' , type = str, 
            default = "pacing_points.pts", 
            help = 'Path to Output file')
    
    parser.add_argument('--tag' , type = int, 
            default = 13, 
            help = "Tag of healthy/unhealthy myocardium to compute pacing points on")
    
    parser.add_argument('--segs', type = list, 
            default = [1,2,3], 
            help = "List of AHA segments to compute pacing points on")
    
    parser.add_argument('--res' , type = int, 
            default = 100, 
            help = "Resolution for mesh point serch (meshtool query idxlist")
    
    args = parser.parse_args()

    main(args)

