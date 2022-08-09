"""
Generate reference-less or conventional ECG/EGM-based pace-maps

This module provides methods to compute conventional pace-maps given a reference VT signal (ecgs or egm vectors) and paced beat traces,
and/or reference-less pace-maps between paced beat traces (within a certain neighbourhood).

"""

import numpy as np
from dataclasses import dataclass
from scipy import signal
import matlab.engine

@dataclass
class PaceMaps:
    """
    Class for keeping track of conventional and/or reference-less pace-maps, and additional grid points needed for the latter

    Args:
        ref (np.ndarray):   (Nleads,TimeSteps) array of the reference VT ecg/egm
        paced (np.ndarray): (Nleads*Npaces,TimeSteps) array of the paced beat traces
        type_sig (str):     'ecg' or 'egm' (lower or upper)
        flag (int):         flag to compute: (0) mean correlation of ALL leads; (1) mean correlation of highest (ALL leads - 2) ---> default: 1
        cycle (int):        basic cycle length (BCL) of paced beat ---> default: 300 ms
        pts (np.ndarray):   (N,3) array of pacing points ---> default: empty
        thres (float):      radius of points to consider to perform reference-less pace-mapping (in mm) ---> default: 20000 microm
        units (float):      units in mm ---> default: 100 (mm). Meshes are usually in micrometers

    """

    ref: np.ndarray
    pace:  np.ndarray
    type_sig: str
    flag: np.int = 1
    cycle: np.int = 30
    pts: np.ndarray = np.zeros([], dtype=float)
    thres: float = 20000
    units: float = 100

    def __post_init__(self):

        # Initialise number of leads according to type of signal
        if self.type_sig.lower() == 'ecg':
            self.Nleads = 12
        elif self.type_sig.lower() == 'egm':
            self.Nleads = 8

        # Initialise number of pacing sites
        self.Nsites = int(round(self.pace.shape[0]/self.Nleads))

        """ Setting up pacemap object """
        print('\n-------- Setting up object --------')
        print('Type signal: {}\tCycle paced beats: {}\tNumber of leads: {}\tNumber of pacing sites: {}\t'.format(self.type_sig,self.cycle,self.Nleads,self.Nsites))
        print('-----------------------------------\n')

    # Function to check and/or remove nan values (for referenceless pace-mapping)
    def remove_nan(self, corr, stim):

        """
            Check/Remove nan elements from correlation array, and modify point stimuli accordingly

            Args:
                corr (np.ndarray):  correlation array
                stim (np.ndarray):  point stimuli array

            Returns:
                new_corr (np.ndarray):  new correlation array
                new_stim (np.ndarray):  new point stimuli array

            """


        print('Checking for naan value ...\n')
        # Find corr values that are not nan
        not_nan_ind = np.where(~np.isnan(corr))[0]
        new_corr = corr[not_nan_ind]
        new_stim = stim[not_nan_ind,:]

        return new_corr, new_stim

    # Function to check and/or remove inf values (for referenceless pace-mapping)
    def remove_inf(self, corr, stim):

        """
            Check/Remove inf elements from correlation array, and modify point stimuli accordingly

            Args:
                corr (np.ndarray):  correlation array
                stim (np.ndarray):  point stimuli array

            Returns:
                new_corr (np.ndarray):  new correlation array
                new_stim (np.ndarray):  new point stimuli array

            """

        print('Checking for inf value ...\n')
        # Find corr values that are not nan
        not_nan_ind = np.where(~np.isinf(corr))[0]
        new_corr = corr[not_nan_ind]
        new_stim = stim[not_nan_ind,:]

        return new_corr, new_stim


    # Function to align QRSs and compute correlation
    def aling_and_correlate(self,x,y,flag=1,cycle=30):

        """
            Align two signals according to max cross correlation, crop qrs and compute correlation

            Args:
                x (np.ndarray):  reference signal Nleads x Nsteps
                y (np.ndarray):  pace signal Nleads x Nsteps
                flag (int):      same as self.flag. (0) mean correlation of ALL leads; (1) mean correlation of highest (ALL leads - 2) ---> default: 1
                cycle (int):     basic cycle length (BCL) of paced beat ---> default: 30 timesteps (300 ms)

            Returns:
                mean_corr (float):  final correlation value (according to flag)

            """

        # Initialise correlation
        corr = np.zeros([x.shape[0],])

        # Loop over each lead
        for i in range(x.shape[0]):

            """ Align signals according to cross-correlation """
            # Compute cross correlation signals
            crosscorr = signal.correlate(x[i,:],y[i,:])
            # Index of best correlation
            best_correlate_ind = np.argmax(crosscorr)
            # Shift positions
            shift_positions = np.arange(-len(x[i,:]) + 1,len(x[i,:])-1,1)
            # Shift amount
            shift = shift_positions[best_correlate_ind]
            shift_amount = int(round(shift))

            # Shift signals according to neg or pos value of shift
            if shift>0:
                xa = x[i,shift_amount:]
                ya = y[i,:]
            else:
                ya = y[i,abs(shift_amount):]
                xa = x[i, :]

            # Crop signals so they have same length
            #xa,ya = self.check_shape(xa,ya)

            """ Find peak(s) in reference signal and crop QRSs"""
            ref_peaks,_ = signal.find_peaks(xa**2, distance=cycle)
            pace_peaks,_ = signal.find_peaks(ya**2, distance=cycle)

            diff_peaks = abs(pace_peaks - ref_peaks[:, np.newaxis])

            """
            # Uncomment following to plot align signals and indices of max peaks
            #min_dist = np.argmin(diff_peaks)
            min_dist = np.unravel_index(np.argmin(diff_peaks, axis = None), diff_peaks.shape)
            print(min_dist)
            plt.plot(xa, '-k')
            plt.plot(ya, '--r')
            plt.plot(ref_peaks[min_dist[0]], xa[ref_peaks[min_dist[0]]], 'or')
            plt.plot(pace_peaks[min_dist[1]], ya[pace_peaks[min_dist[1]]], 'ok')
            plt.show()
            """


            # If more than 1 peaks
            if len(ref_peaks)>1:

                # Find best aligned peaks
                ind_max = np.unravel_index(np.argmin(diff_peaks, axis = None), diff_peaks.shape)
                j_max = pace_peaks[ind_max[1]]

                # Initialise cropping values
                ind_1 = j_max - cycle // 2
                ind_2 = j_max + cycle // 2

                """
                if ind_1<0:
                    ind_1 = 0
                if ind_2>len(ya):
                    ind_2 = len(ya)

                # Find first intersection
                for j in range(j_max,0,-1):
                    if ya[j]*ya[j_max]<=0.006:
                        ind_1 = j
                        break

                # Find second intersection
                for j in range(j_max,len(ya)):
                    if ya[j]*ya[j_max]<=0.006:
                        ind_2 = j
                        break
                """

                # Crop signals
                xa = xa[ind_1:ind_2]
                ya = ya[ind_1:ind_2]

                """
                # Uncomment to plot aligned and cropped signals
                plt.plot(xa, '-k')
                plt.plot(ya, '--r')
                plt.show()
                """


            """ Compute correlation for each pair of aligned signals """
            # Compute correlation
            R = np.corrcoef(xa, ya)
            corr[i] = R[0, 1]

            # Uncomment to plot correlation in the title of the figure
            # plt.title('Corr: {}'.format(corr[i]*100))


        # Flag 0 returns mean of all correlation, Flag 1 return mean of best 10
        """ Compute final correlation value across leads according to desired mode mean of 12/12 or mean of top 10/12 """
        if flag == 0:
            mean_corr = np.mean(corr)
        else:
            sorted_corr = sorted(corr, reverse=True)
            mean_corr = np.mean(sorted_corr[:-2])

        return mean_corr

    # Function to check whether signals have same length
    def check_shape(self,x,y):

        """
            Check size of signals and crop them so that they have same length

                Args:
                    x (np.ndarray):  reference signal Nleads x Nsteps
                    y (np.ndarray):  pace signal Nleads x Nsteps

                Returns:
                    x (np.ndarray):  cropped reference signal Nleads x Nsteps
                    y (np.ndarray):  cropped pace signal Nleads x Nsteps

                    """

        # Shortening 2D signals
        if x.ndim>1:
            m = 1
            if x.shape[m]>y.shape[m]:
                x = x[:,:y.shape[m]]
            elif x.shape[m]<y.shape[m]:
                y = y[:,:x.shape[m]]
        # Shortening 1D signals
        elif x.ndim==1:
            m = 0
            if x.shape[m]>y.shape[m]:
                x = x[:y.shape[m]]
            elif x.shape[m]<y.shape[m]:
                y = y[:x.shape[m]]
        return x,y


    # Function to compute conventional in-silico pace-maps
    def conventional(self):

        """
            Generate conventional in-silico pace-maps by looping over all pace signals and align/correlate with reference

            Args:
                 no arguments, see PaceMaps class

            Returns:
                corr (np.ndarray):      correlation array Nsites x 1
                stim (np.ndarray):      point stimuli Nsites x 3
        """

        print('Computing conventional pace-maps for {} ...'.format(self.type_sig))
        # Initialise counter
        l = 0
        # Initialise correlation array
        pacemap = np.zeros([self.Nsites,])
        # Loop over each pace signal
        for i in range(self.Nsites):
            # 12-lead or 8-vector pace ECG/EGM
            new_pace = self.pace[l : self.Nleads + l,:]
            # Align pace to reference and compute correlation
            pacemap[i] = self.aling_and_correlate(self.ref, new_pace, self.flag, self.cycle)
            # Increment counter
            l += self.Nleads

        # Check for nan
        new_corr, new_stim = self.remove_nan(pacemap, self.pts)
        # Check for inf
        new_corr, new_stim = self.remove_inf(new_corr, new_stim)

        return new_corr*100, new_stim

    # Function to compute referenceless in-silico pace-maps
    def referenceless(self):

        """
            Generate reference-less in-silico pace-maps by looping over all pace signals and align/correlate with reference

                Args:
                    no arguments, see PaceMaps class

                Returns:
                    new_corr (np.ndarray):      correlation array
                    new_stim (np.ndarray):      new interpolated point stimuli
                """

        print('Computing referenceless pace-maps for {} ...'.format(self.type_sig))
        # Initialise counter
        m = 0
        corr = []
        new_stim = []

        for j in range(self.Nsites):

            l = 0
            # Take one pacing signal and corresponding coordinates
            pace_1 = self.pace[m : self.Nleads + m,:]
            pts_1 = self.pts[j,:]

            # Compute correlation with another pace if not computed already (i>j) and within a certain radius (self.thres)
            for i in range(self.Nsites):

                # Take another pacing signal and corresponding coordinates
                pace_2 = self.pace[l : self.Nleads + l,:]
                pts_2 = self.pts[i,:]
                # Distance between the two paces
                dist = np.linalg.norm(pts_2 - pts_1)

                if i>j and dist<self.thres:
                    # Compute distance between the points (according to units)
                    dist /= self.units
                    # Compute correlation between the two points and divide by distance
                    corr.append(abs(self.aling_and_correlate(pace_1,pace_2,self.flag,self.cycle)*100 - 100)/dist)
                    # Compute mid point between the two
                    new_stim.append((pts_2 + pts_1)/2)

                l += self.Nleads

            m += self.Nleads

        # Check for nan
        new_corr, new_stim = self.remove_nan(np.asarray(corr), np.asarray(new_stim))
        # Check for inf
        new_corr, new_stim = self.remove_inf(new_corr, new_stim)

        return new_corr, new_stim





