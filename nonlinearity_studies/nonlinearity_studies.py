import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.signal import find_peaks as scipy_find_peaks
import math

from pathlib import Path
from glob import glob

#---------------- ANALYSIS FUNCTIONS ----------------------------

#---------------- (0) Convert to electrons ----------------------
def convert_to_electrons(data, pedestal, gain, flatten=True):
    if flatten:
        data = np.array(data).flatten()
    data_electrons = (data - pedestal) / gain  # Subtract pedestal (mean ADU of zero electron peak) and divide by gain
    return data_electrons

#---------------- (1) Calculate noise/gain ----------------------
# Function finds noise and gain from input pixel charge data
# Uses the two largest peaks in the data (zero and one electron peaks)
def calculate_noise_gain(data, n=200, fit_bounds='default'):

    data = np.array(data).flatten()
    
    # Create histogram of data up to max=100 to find peaks (zero and one electron peaks)
    # This avoids detecting spurious peaks at high charge values
    hist_max = 100
    hist_range = (np.min(data), hist_max)
    nbins = int(n * (hist_max - np.min(data)))
    counts, edges = np.histogram(data, bins=nbins, range=hist_range)
    centers = 0.5 * (edges[:-1] + edges[1:])
    
    # Find the two largest peaks
    peaks, properties = scipy_find_peaks(counts, height=0)
    
    if len(peaks) < 2:
        raise ValueError("Did not find two distinct peaks in the data")
    
    # Get indices of the two peaks with highest counts
    peak_heights = counts[peaks]
    top_two_indices = np.argsort(peak_heights)[-2:]
    top_two_peaks = peaks[top_two_indices]
    
    # Sort by position (not by height) - first peak is zero electrons, second is one electron
    top_two_peaks = np.sort(top_two_peaks)
    
    zero_peak_index = top_two_peaks[0]
    one_peak_index = top_two_peaks[1]
    
    zero_peak_charge = centers[zero_peak_index]
    one_peak_charge = centers[one_peak_index]
    
    # Define range around both peaks for fitting
    # Use a window that includes both peaks with some margin
    margin = (one_peak_charge - zero_peak_charge) * 0.5
    fit_left = zero_peak_charge - margin
    fit_right = one_peak_charge + margin
    fit_range = [fit_left, fit_right]
    
    # Create histogram for the fitting range
    fit_nbins = int(n * (fit_right - fit_left))
    fit_counts, fit_edges = np.histogram(data, bins=fit_nbins, range=fit_range)
    fit_centers = 0.5 * (fit_edges[:-1] + fit_edges[1:])
    
    # Set tight fit bounds around the detected peak positions
    if fit_bounds == 'default':
        # Tight bounds around zero peak (mu0)
        mu0_lower = zero_peak_charge - margin * 0.3
        mu0_upper = zero_peak_charge + margin * 0.3
        
        # Tight bounds around one peak (mu1)
        mu1_lower = one_peak_charge - margin * 0.3
        mu1_upper = one_peak_charge + margin * 0.3
        
        fit_bounds = (
            [0.00001, mu0_lower, 0.00001, mu1_lower, 1, 1],
            [1, mu0_upper, 1, mu1_upper, max(fit_counts), max(fit_counts)]
        )
    
    popt, pcov = curve_fit(double_gauss, fit_centers, fit_counts, maxfev=20000, bounds=fit_bounds)
    
    # Extract pedestal, noise, gain, and rest of double gaussian coefficients from curve fit
    pedestal = tuple(popt)[1]  # Pedestal is mean of zero electron peak
    noise = tuple(popt)[0]  # Noise is standard deviation of zero electron peak
    gain = tuple(popt)[3] - tuple(popt)[1]  # Gain is difference between mean of one and zero electron peaks

    return fit_counts, fit_edges, pedestal, noise, gain, popt, fit_range


#---------------- (2) Find peaks ----------------------------
# Finds all electron peaks in flattened pixel charge array per extension
# First converts data from ADU to electrons (if not already in electrons)
# Input is charge data (in ADU or electrons) from one extension
# bins is the number of bins given for initial charge histogram
# bin_factor is the multiple of the length of range used in fitting all peaks (number of bins per peak essentially)
# bin_factor is also used to define the distance parameter given to scipy_find_peaks 
# for the min distance between peaks (with buffer given by buffer) 
def find_all_peaks(data, 
                   width, 
                   buffer, 
                   pedestal,
                   noise,
                   gain,
                   bins='default',
                   flatten=True,
                   do_convert_to_electrons=True,
                   range_left='left_of_zero', 
                   range_right=2500, 
                   bin_factor=8):
    
    if flatten:
        data=np.array(data).flatten()
    
    if do_convert_to_electrons:
        data = convert_to_electrons(data, pedestal, gain, flatten=False)

    if range_left=='default':
        if do_convert_to_electrons:
            range_left = -0.5
        else:
            range_left=pedestal-3*noise-1 # If finding electron peaks in ADU, start fitting peaks at the left side of the zero electron peak

    hist_range = (range_left, range_right)

    if bins=='default':
       bins=math.floor((hist_range[1]-hist_range[0])*bin_factor)

    counts, edges = np.histogram(data, bins=bins,range=hist_range)
    centers = 0.5 * (edges[1:] + edges[:-1])
    peaks, properties = scipy_find_peaks(counts, height=0, width=width, distance=bin_factor-buffer)

    print(f'Width = {width}, buffer = {buffer}, bin factor = {bin_factor} used for peak finder')

    return counts, edges, peaks, centers, properties, hist_range


#---------------- (3) Fit nonlinearity ----------------------------
def fit_nonlinearity(peaks, centers, pedestal, gain, fit_range_right, do_convert_to_electrons=False, fit_bounds_low=-100, fit_bounds_high=100):
    # If specified, convert to electrons: subtract the pedestal (mean of zero electron peak) from all charge values, 
    # Then divide by the gain (difference between zero and 1 electron peak)
    if do_convert_to_electrons:
        peak_charge_e = np.array([(centers[p]-pedestal)/gain for p in peaks])

    # If peaks were given in electrons, do not convert to electrons. Just pick out electron values of peak locations.
    else:
        peak_charge_e = np.array([centers[p] for p in peaks])

    charge_minus_npeak = [(peak_charge_e[i] - i) for i in range(len(peaks))]
    parabola_coeff, parabola_pcov = curve_fit(parabola, peak_charge_e[:fit_range_right], charge_minus_npeak[:fit_range_right],
                           maxfev=2000, bounds=(fit_bounds_low, fit_bounds_high))
    return parabola_coeff, parabola_pcov, peak_charge_e, charge_minus_npeak

#---------------- (4) Get nonlinearity at charge value ----------------------------
# use coefficients of parabola fit to nonlinearity curve from fit_nonlinearity to estimate the nonlinearity
# at a single value or list of values
# required arguments: q (charge value that parabola equation is evaluated at, can be int, float or list), parabola_coeff (the curve_fit optimal parabola coefficients found in fit_nonlinearity)
# optional arguments: parabola_pcov, fit_range_right (used to ascertain how confident we should be in the nonlinearity value); print_values (if you want to print what the function outputs)
def get_nonlinearity_at(q, parabola_coeff, parabola_pcov=None, fit_range_right=None, print_values=True):
    a = parabola_coeff[0]
    b = parabola_coeff[1]
    c = parabola_coeff[2]

    if type(q)==float or type(q)==int:
        nonlinearity_at_q = a*q**2 + b*q + c
    
    if type(q)==list:
        nonlinearity_at_q = [(a*q_i**2 + b*q_i + c) for q_i in q]
    
    if print_values==True:
        print(f'Interpolated nonlinearity at {q} e- = {nonlinearity_at_q}')

        if fit_range_right!=None:
            print(f'Parabola fit was performed up to {fit_range_right} e-')
        if type(parabola_pcov)!=None:
            print(f'Covariance of fit = {parabola_pcov}')

    return nonlinearity_at_q

#---------------- PLOTTING FUNCTIONS ----------------------------

#---------------- Plot zero-one peaks  -----------------------
# Usage: for plotting zero-one electron peaks from each extension on same subplot or individually by extension.
# Input data is list of 2D pixel charge arrays from all 4 extensions.
# xlim can be 'default', 'none', or tuple(left, right)
# ylim can be 'none' or tuple(bottom, top)
def plot_zero_one_peaks(data_ext, 
                        zero_one_counts_ext,
                        zero_one_edges_ext,
                        pedestals, 
                        gains, 
                        double_gauss_popts, 
                        zero_one_ranges,
                        individual_figsize=(6,5), 
                        subplots_figsize=(9,7),
                        xlim='default', ylim='none',
                        fontsize=7.5,
                        yscale='linear',
                        n=200, 
                        do_convert_to_electrons=False,
                        plot_individual=False,
                        plot_together=True,
                        save_plots=False,
                        fig_path='./', file='zero_one_peaks', 
                        dpi=350):

    fig_path = Path(fig_path)
    if file != 'zero_one_peaks':
        base_name = file[:-5] + '_zero_one_peaks'
    else:
        base_name = file
    fig_name = fig_path / base_name

    if plot_individual:
        for ext, data in enumerate(data_ext):
            data = np.array(data).flatten()

            zero_one_counts=zero_one_counts_ext[ext]
            zero_one_edges=zero_one_edges_ext[ext]
            pedestal=pedestals[ext] 
            gain=gains[ext]
            double_gauss_popt=double_gauss_popts[ext]
            zero_one_range=zero_one_ranges[ext]

            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            fig.suptitle('Double-Gaussian Fit to Zero and One Electron Peaks in ADU')

            double_gauss_coeff = tuple(double_gauss_popt)+(gain,)
            data_window = data[(data > zero_one_range[0]) & (data < zero_one_range[1])]
            nbins=int(n*(zero_one_range[1] - zero_one_range[0]))

            bin_width = zero_one_edges[1] - zero_one_edges[0]
            zero_one_centers = 0.5 * (zero_one_edges[:-1] + zero_one_edges[1:])

            if yscale=='log':
                zero_one_counts = np.maximum(zero_one_counts, 1) #need in order to prevent empty bars in histogram if there are any bins that have 0 counts
                ax.set_yscale('log')
            elif yscale!='linear':
                ax.set_yscale(yscale)
            ax.bar(zero_one_edges[:-1], zero_one_counts, edgecolor='none', linewidth=0, align='edge', width=np.diff(zero_one_edges))

            ax.set_xlabel('Charge (ADU)')
            ax.set_ylabel('N')
            ax.set_title(f'EXT {ext}')
            
            ax.plot(zero_one_centers, double_gauss(zero_one_centers, *double_gauss_popt), 'r',
                label=r'$\sigma_0$ = %5.3f, $\mu_0$ = %5.3f, $\sigma_1$ = %5.3f, $\mu_1$ = %5.3f,'%double_gauss_coeff[0:4]
                +'\n'+'$N_0$ = %5.3f, $N_1$ = %5.3f, gain = %5.3f ADU/$e^{–}$'%double_gauss_coeff[4:])
            ax.legend(loc="upper right", fontsize=fontsize)
            
            if xlim=='default':
                ax.set_xlim(zero_one_range[0],zero_one_range[1])
            elif xlim!='none':
                ax.set_xlim(xlim)
            
            if ylim!='none':
                ax.set_ylim(ylim)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plots to {output_path}')
            plt.show()

        if do_convert_to_electrons:
            fig, axs = plt.subplots(2, 2, figsize=individual_figsize, constrained_layout=True)
            fig.suptitle(r'Double-Gaussian Fit to Zero and One Electron Peaks in $e^-$')
            axs = axs.flatten()
            
            for ext, data in enumerate(data_ext):
                data=np.array(data).flatten()
                ax = axs[ext]
                zero_one_range = zero_one_ranges[ext]
                pedestal = pedestals[ext]
                gain = gains[ext] 

                data_window=data[(data > zero_one_range[0]) & (data < zero_one_range[1])]

                data_window_e = convert_to_electrons(data_window, pedestal, gain)
                zero_one_range_e = convert_to_electrons(zero_one_range, pedestal, gain) 
                nbins = int(n * (zero_one_range_e[1] - zero_one_range_e[0]))

                zero_one_counts_e, zero_one_edges_e = np.histogram(data_window_e, bins=nbins, range=zero_one_range_e)
                zero_one_centers_e = 0.5 * (zero_one_edges_e[:-1] + zero_one_edges_e[1:])
                bin_width_e = zero_one_edges_e[1] - zero_one_edges_e[0]
            
                double_gauss_popt_e, double_gauss_pcov_e = curve_fit(double_gauss, zero_one_centers_e, zero_one_counts_e, maxfev=2000, 
                                                                     bounds=([0.0001, -1, 0.0001, 0.5, 0, 0], 
                                                                             [1.0,  1,  1.0,  2.0, np.inf, np.inf]))
                
                if yscale=='log':
                    zero_one_counts_e = np.maximum(zero_one_counts_e, 1) #need in order to prevent empty bars in histogram if there are any bins that have 0 counts
                    ax.set_yscale('log')
                elif yscale!='linear':
                    ax.set_yscale(yscale)

                ax.bar(zero_one_edges_e[:-1], zero_one_counts_e, align='edge', edgecolor='none', linewidth=0, width=np.diff(zero_one_edges_e))
                ax.set_xlabel(r'Charge ($e^–$)')
                ax.set_ylabel('N')
                ax.set_yscale(yscale)
                ax.set_title(f'EXT {ext}')
                ax.plot(zero_one_centers_e, double_gauss(zero_one_centers_e, *double_gauss_popt_e), 'r',
                    label=r'$\sigma_0$ = %5.3f $e^{–}$, $\mu_0$ = %5.3f $e^{–}$, $\sigma_1$ = %5.3f $e^{–}$, $\mu_1$ = %5.3f $e^{–}$'%tuple(double_gauss_popt_e)[0:4])
                ax.legend(loc="upper right", fontsize=fontsize)
                
                if xlim=='default':
                    ax.set_xlim(zero_one_range_e[0], zero_one_range_e[1])
                elif xlim!='none':
                    ax.set_xlim(xlim)

                if ylim!='none':
                    ax.set_ylim(ylim)

                if save_plots:
                    output_path = fig_name.with_stem(fig_name.stem + f'_electrons_EXT{ext}').with_suffix('.jpeg')
                    plt.savefig(str(output_path), dpi=dpi)
                    print(f'Saved plot to {output_path}')
                plt.show()

    if plot_together:

        fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True)
        fig.suptitle('Double-Gaussian Fit to Zero and One Electron Peaks in ADU')
        axs = axs.flatten()

        for ext, data in enumerate(data_ext):
            data = np.array(data).flatten()
            zero_one_counts=zero_one_counts_ext[ext]
            zero_one_edges=zero_one_edges_ext[ext]
            pedestal=pedestals[ext] 
            gain=gains[ext]
            double_gauss_popt=double_gauss_popts[ext]
            zero_one_range=zero_one_ranges[ext]

            ax = axs[ext]
            double_gauss_coeff = tuple(double_gauss_popt)+(gain,)
            data_window = data[(data > zero_one_range[0]) & (data < zero_one_range[1])]
            nbins=int(n*(zero_one_range[1]-zero_one_range[0]))

            zero_one_centers = 0.5 * (zero_one_edges[:-1] + zero_one_edges[1:])
            bin_width = zero_one_edges[1] - zero_one_edges[0]

            if yscale=='log':
                zero_one_counts = np.maximum(zero_one_counts, 1) #need in order to prevent empty bars in histogram if there are any bins that have 0 counts
                ax.set_yscale('log')
            elif yscale!='linear':
                    ax.set_yscale(yscale)

            ax.bar(zero_one_centers, zero_one_counts, align='center', edgecolor='none', linewidth=0, width=bin_width)
            
            ax.plot(zero_one_centers, double_gauss(zero_one_centers, *double_gauss_popt), 'r',
                label=r'$\sigma_0$ = %5.3f, $\mu_0$ = %5.3f, $\sigma_1$ = %5.3f, $\mu_1$ = %5.3f,'%double_gauss_coeff[0:4]
                +'\n'+'$N_0$ = %5.3f, $N_1$ = %5.3f, gain = %5.3f ADU/$e^{–}$'%double_gauss_coeff[4:])
            
            ax.set_xlabel('Charge (ADU)')
            ax.set_ylabel('N')
            ax.set_title(f'EXT {ext}')
            ax.legend(loc="upper right", fontsize=fontsize)

            if xlim=='default':
                ax.set_xlim(zero_one_range[0],zero_one_range[1])
            elif xlim!='none':
                ax.set_xlim(xlim)

            if ylim!='none':
                ax.set_ylim(ylim)
        if save_plots:
            output_path = fig_name.with_suffix('.jpeg')
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        plt.show()

        if do_convert_to_electrons:
            fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True)
            fig.suptitle(r'Double-Gaussian Fit to Zero and One Electron Peaks in $e^-$')
            axs = axs.flatten()

            for ext, data in enumerate(data_ext):
                ax = axs[ext]
                data = np.array(data).flatten()
                zero_one_range = zero_one_ranges[ext]
                pedestal = pedestals[ext]
                gain = gains[ext]

                data_window = data[(data > zero_one_range[0]) & (data < zero_one_range[1])]

                data_window_e = convert_to_electrons(data_window, pedestal, gain)
                zero_one_range_e = convert_to_electrons(zero_one_range, pedestal, gain) 
                nbins = int(n * (zero_one_range_e[1] - zero_one_range_e[0]))

                zero_one_counts_e, zero_one_edges_e = np.histogram(data_window_e, bins=nbins, range=zero_one_range_e)
                zero_one_centers_e = 0.5 * (zero_one_edges_e[:-1] + zero_one_edges_e[1:])
                bin_width_e = zero_one_edges_e[1] - zero_one_edges_e[0]
                
                double_gauss_popt_e, double_gauss_pcov_e = curve_fit(double_gauss, zero_one_centers_e, zero_one_counts_e, maxfev=2000, bounds=([0.0, -1, 0.0, 0.5, 0.0, 0.0], 
                                                                                                                                               [1.0,  1,  1.0,  2.0, 1e7, 1e7]))
                
                if yscale=='log':
                    zero_one_counts_e = np.maximum(zero_one_counts_e, 1) #need in order to prevent empty bars in histogram if there are any bins that have 0 counts
                    ax.set_yscale('log')

                elif yscale!='linear':
                    ax.set_yscale(yscale)

                ax.bar(zero_one_centers_e, zero_one_counts_e, align='center', edgecolor='none', linewidth=0, width=bin_width_e)

                ax.set_title(f'EXT {ext}')
                ax.plot(zero_one_centers_e, double_gauss(zero_one_centers_e, *double_gauss_popt_e), 'r',
                    label=r'$\sigma_0$ = %5.3f $e^{–}$, $\mu_0$ = %5.3f $e^{–}$, $\sigma_1$ = %5.3f $e^{–}$, $\mu_1$ = %5.3f $e^{–}$'%tuple(double_gauss_popt_e)[0:4])
                ax.legend(loc="upper right", fontsize=fontsize)
                ax.set_xlabel(r'Charge ($e^–$)')
                ax.set_ylabel('N')

                ax.set_title(f'EXT {ext}')
                
                if xlim=='default':
                    ax.set_xlim(zero_one_range_e[0], zero_one_range_e[1])
                elif xlim!='none':
                    ax.set_xlim(xlim)

                if ylim!='none':
                    ax.set_ylim(ylim)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + '_electrons').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            plt.show()


#---------------- Plot all electron peaks ----------------------------
# Input is list of data from each of four extensions
# ylim can be 'none' or tuple=(ylim_bottom, ylim_top)
def plot_all_peaks(counts_ext, 
                   peaks_ext, 
                   centers_ext, 
                   xlim, ylim='none', 
                   yscale='log', 
                   plot_individual=True, plot_together=False, 
                   draw_lines=True, linecolor='r', linestyle='--',
                   individual_figsize=(6,5), subplots_figsize=(9,7),
                   suptitle='Peaks in Pixel Charge Distribution',
                   save_plots=False,
                   fig_path='./', file='peak_finder', 
                   dpi=350):

    fig_path = Path(fig_path)
    if file != 'peak_finder':
        base_name = file[:-5] + '_peak_finder'
    else:
        base_name = file
    fig_name = fig_path / base_name

    if plot_individual:
        for ext, counts in enumerate(counts_ext):
            peaks=peaks_ext[ext]
            centers=centers_ext[ext]
            bin_width = centers[1] - centers[0]
            
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            fig.suptitle(suptitle+f': EXT {ext}')
            ax.bar(centers, counts, align='center', edgecolor='none', linewidth=0, width=bin_width)
            ax.set_xlabel(r'Charge ($e^-$)')
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            ax.set_xlim(xlim)
            if ylim!='none':
                ax.set_ylim(ylim)

            # draw vertical lines and labels at each peak
            if draw_lines:
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    peak_y = counts[p]

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    ax.text(peak_x,
                        peak_y,
                        f"{i}",
                        verticalalignment='bottom',
                        horizontalalignment='center',
                        color=linecolor,
                        fontsize=10)
                    
            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            plt.show()

    if plot_together:
        fig, axs = plt.subplots(2,2,figsize=subplots_figsize,constrained_layout=True)
        axs=axs.flatten()
        fig.suptitle(suptitle)

        for ext, counts in enumerate(counts_ext):
            peaks=peaks_ext[ext]
            centers=centers_ext[ext]
            bin_width = centers[1] - centers[0]
            ax = axs[ext]

            ax.bar(centers, counts, align='center', edgecolor='none', linewidth=0, width=bin_width)
            ax.set_xlabel(r'Charge ($e^-$)')
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            ax.set_xlim(xlim)
            if ylim!='none':
                ax.set_ylim(ylim)
            ax.set_title(f'EXT {ext}')

            # draw vertical lines and labels at each peak
            if draw_lines:
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    peak_y = counts[p]

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    ax.text(peak_x,
                        peak_y,
                        f"{i}",
                        verticalalignment='bottom',
                        horizontalalignment='center',
                        color=linecolor,
                        fontsize=10)
        
        if save_plots:
            output_path = fig_name.with_suffix('.jpeg')
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        plt.show()


#---------------- Plot nonlinearity ----------------------------
# xlim and ylim can be 'default', 'none', or tuple(ylim_bottom, ylim_top)
def plot_nonlinearity(peaks_ext,
                      parabola_coeffs, 
                      peak_charge_e_ext, 
                      charge_minus_npeak_ext,
                      fit_range_right_ext,
                      xlim='default', ylim='default',
                      individual_figsize=(6,5), subplots_figsize=(9,7),
                      suptitle='Pixel Charge Nonlinearity Curve',
                      line_color='r', 
                      scatter_color='b', 
                      s=2, 
                      alpha=0.5,
                      plot_individual=False, 
                      plot_together=True, 
                      save_plots=False, 
                      fig_path='./', file='nonlinearity_curve', 
                      dpi=350):

    fig_path = Path(fig_path)
    if file != 'nonlinearity_curve':
        base_name = file[:-5] + '_nonlinearity'
    else:
        base_name = file
    fig_name = fig_path / base_name

    if plot_individual:
        for ext, peaks in enumerate(peaks_ext):
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            fig.suptitle(suptitle+f': EXT {ext}')
            ax.grid()

            parabola_coeff=parabola_coeffs[ext]
            peak_charge_e=peak_charge_e_ext[ext]
            charge_minus_npeak=charge_minus_npeak_ext[ext]
            fit_range_right=fit_range_right_ext[ext]

            ax.plot(peak_charge_e, parabola(peak_charge_e, *parabola_coeff), color=line_color,
                        label=r'$%5.6f x^2 + %5.3f x + %5.3f$' %tuple(parabola_coeff))
            ax.scatter(peak_charge_e, charge_minus_npeak, c=scatter_color, s=s, alpha=alpha)
            ax.legend(loc="upper right", fontsize=8)
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim=='default':
                ax.set_ylim(min(charge_minus_npeak)-10, max(charge_minus_npeak[:fit_range_right])+5)
            elif ylim!='none':
                ax.set_ylim(ylim)

            if xlim=='default':
                ax.set_xlim(-100, fit_range_right+500)
            elif xlim!='none':
                ax.set_xlim(xlim)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            plt.show()
            
    if plot_together:
        fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True)
        axs=axs.flatten()
        fig.suptitle(suptitle)
        for ext, peak_charge_e in enumerate(peak_charge_e_ext):
            ax = axs[ext]
            ax.grid()

            parabola_coeff=parabola_coeffs[ext]
            charge_minus_npeak=charge_minus_npeak_ext[ext]
            fit_range_right=fit_range_right_ext[ext]

            ax.plot(peak_charge_e, parabola(peak_charge_e, *parabola_coeff), color=line_color,
                        label=r'$%5.6f x^2 + %5.3f x + %5.3f$' %tuple(parabola_coeff))
            ax.scatter(peak_charge_e, charge_minus_npeak, c=scatter_color, s=s, alpha=alpha)
            ax.legend(loc="upper right", fontsize=8)
            ax.set_title(f'EXT {ext}')
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim=='default':
                ax.set_ylim(min(charge_minus_npeak)-10, max(charge_minus_npeak[:fit_range_right])+5)
            elif ylim!='none':
                ax.set_ylim(ylim)
  
            if xlim=='default':
                ax.set_xlim(-100, fit_range_right+500)
            elif xlim!='none':
                ax.set_xlim(xlim)

    if save_plots:
        output_path = fig_name.with_suffix('.jpeg')
        plt.savefig(str(output_path), dpi=dpi)
        print(f'Saved plot to {output_path}')
    plt.show()


#---------------- UTILITY FUNCTIONS ----------------------------

#---------------- Get Fits ----------------------------
def get_fits(file_input):
    """
    Load FITS extensions from a file.

    Parameters
    ----------
    file_input : str or Path
        Path to the FITS file (absolute or relative to current working directory)

    Returns
    -------
    ext_charge : list
        List of data arrays from extensions 1–4
    """
    file_path = Path(file_input).resolve()
    
    # Check if file exists
    if not file_path.exists():
        raise FileNotFoundError(f"FITS file not found: {file_path}")

    # Load FITS file
    with fits.open(str(file_path)) as hdu_list:
        ext_charge = [hdu_list[i].data for i in range(1, 5)]

    return ext_charge


#---------------- Return data for each extensions in a list from pixel charge data for all extensions
def get_zero_one_peaks_ext(data_ext, fit_bounds='default'):
    zero_one_counts_ext = []
    zero_one_edges_ext = []
    pedestals = []
    gains = []
    double_gauss_popts = []
    zero_one_ranges = []
    for data in data_ext:
        data=np.array(data).flatten()

        zero_one_counts, zero_one_edges, pedestal, noise, gain, double_gauss_popt, zero_one_range = calculate_noise_gain(data, fit_bounds=fit_bounds)
        zero_one_counts_ext.append(zero_one_counts)
        zero_one_edges_ext.append(zero_one_edges)
        pedestals.append(pedestal)
        gains.append(gain)
        double_gauss_popts.append(double_gauss_popt)
        zero_one_ranges.append(zero_one_range)

    return zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges
        

def get_all_peaks_ext(data_ext, widths, buffers, pedestals, double_gauss_popts, gains, bins='default', flatten=True, do_convert_to_electrons=True, range_left='default', range_right=2000, bin_factor=8):
    counts_ext = []
    edges_ext = []
    peaks_ext = []
    centers_ext = []
    hist_ranges = []
    for ext, data in enumerate(data_ext):
        width = widths[ext]
        buffer = buffers[ext]
        pedestal = pedestals[ext]
        noise = double_gauss_popts[ext][0]
        gain = gains[ext]

        print(f'\nEXT {ext}:')
        counts, edges, peaks, centers, properties, hist_range = find_all_peaks(data, 
                                                                               width, 
                                                                               buffer, 
                                                                               pedestal,
                                                                               noise,
                                                                               gain,
                                                                               bins=bins,
                                                                               flatten=flatten,
                                                                               do_convert_to_electrons=do_convert_to_electrons,
                                                                               range_left=range_left, 
                                                                               range_right=range_right, 
                                                                               bin_factor=bin_factor)
    
        counts_ext.append(counts)
        edges_ext.append(edges)
        peaks_ext.append(peaks)
        centers_ext.append(centers)
        hist_ranges.append(hist_range)

    return counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges

def get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, gains, fit_range_right_ext, do_convert_to_electrons=False, fit_bounds_low=-100, fit_bounds_high=100, print_values=True):
    
    peak_charge_e_ext = []
    charge_minus_npeak_ext = []
    parabola_coeffs = []
    parabola_pcovs = []

    for ext, peaks in enumerate(peaks_ext):

        centers=centers_ext[ext]
        pedestal=pedestals[ext]
        gain=gains[ext]
        fit_range_right=fit_range_right_ext[ext]

        parabola_coeff, parabola_pcov, peak_charge_e, charge_minus_npeak = fit_nonlinearity(peaks,
                                                                             centers,
                                                                             pedestal, 
                                                                             gain, 
                                                                             fit_range_right, 
                                                                             do_convert_to_electrons,
                                                                             fit_bounds_low, 
                                                                             fit_bounds_high)
        
        print(f'\nEXT {ext}:')
        nonlinearity_at_500 = get_nonlinearity_at(500, parabola_coeff, parabola_pcov, fit_range_right, print_values)

        peak_charge_e_ext.append(peak_charge_e)
        charge_minus_npeak_ext.append(charge_minus_npeak)
        parabola_coeffs.append(parabola_coeff)
        parabola_pcovs.append(parabola_pcov)
    print('\n')

    return peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, parabola_pcovs, nonlinearity_at_500

def get_nonlinearity_at_ext(q, parabola_coeffs, parabola_pcovs, fit_range_right_ext, print_values=True):
    
    nonlinearity_at_q_ext = []
    for ext, parabola_coeff in enumerate(parabola_coeffs):
        parabola_pcov = parabola_pcovs[ext]
        fit_range_right = fit_range_right_ext[ext]
        nonlinearity_at_q = get_nonlinearity_at(q, parabola_coeff, parabola_pcov, fit_range_right, print_values=print_values)
        nonlinearity_at_q_ext.append(nonlinearity_at_q)

    return nonlinearity_at_q_ext

#---------------- Curves ----------------------------
def double_gauss(x, s0, m0, s1, m1, N0, N1):
    return N0 * np.exp(-(x-m0)**2/(2*s0**2)) + N1 * np.exp(-(x-m1)**2/(2*s1**2))

def parabola(x, a, b, c):
    return a*x**2 + b*x + c