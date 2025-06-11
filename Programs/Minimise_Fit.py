"""
Created on 14/09/2023
@author: Leigh Page, lp472@sussex.ac.uk

Description:
Program to find the best fit variables to fit a combined triangle wave and Voigt absorption profile to a given dataset,
then use those to calculate the pressure inside of a Rubidium vapour-cell.
The fitting is accomplished by minimising the Chi^2 between the model and the data.
Inputs:
- save_data: wether to save plot to a folder or not
- dirname_ext: directory address to read data from and/or save files to
- resolution: resolution of saved plots
- data: data to be modeled
- cell_num: the identification number of the data being modeled
- temp: temperature of the cell
- smoothing: range of moving average applied to the data to smooth it out for easier modeling
Outputs:
- variables: final best fit variables
- constants: other constants calculated for the model
- calibration: conversion values used to change data from time to frequency dependance
"""

from sys import exit
import numpy as np
from scipy import signal
from scipy.special import wofz
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd

######################################### Global Variables #########################################
period = 1.0
phaseshift = 1.0
Gamma_G = 1.0
Gamma_L = 1.0
Triangle_amp = 1.0
Voigt_amp = 1.0
nu_0 = 1.0
res = 1.0
Chi2 = 1.0
data_error = 1.0
Gaussian_nu_0 = 1.0

######################################### Functions ################################################

# Triangle wave function
def TriangleWave(x, period, phaseshift):
    return (2/np.pi)*np.arcsin(np.sin((2*np.pi*x/period) + phaseshift))

# Gaussian function array generator
def Gaussian(x, x_0, stdv):
    y_gaussian = (1/(stdv*np.sqrt(2*np.pi)))*np.exp(-0.5*(x - x_0)**2/stdv**2)
    # normalise
    y_gaussian = (y_gaussian - np.amin(y_gaussian))/(np.amax(y_gaussian) - np.amin(y_gaussian))
    return y_gaussian

# Voigt function
# Voigt line shape at x with mean nu_0, Lorentzian component HWHM Gamma_L/2and Gaussian component HWHM Gamma_G/2.
def Voigt(x, nu_0, Gamma_G, Gamma_L):
    sigma = (Gamma_G/2.0)/np.sqrt(2*np.log(2.0))
    z = ((x - nu_0) + 1j*(Gamma_L/2.0))/(sigma*np.sqrt(2.0))
    w = wofz(z)
    return (np.real(w)/(sigma*np.sqrt(2.0*np.pi)))

def convolution(arr1, arr2, x):
    #calculate average spacing
    dx = (x[-1] - x[0])/len(x) # This only works if the spacing is more or less regular
    arr3 = np.zeros(len(arr1))  
    for i in range(len(arr1)):
        for j in range(len(arr1)):
            arr3[i-int(len(arr1)/2) + 4] += arr1[j]*arr2[i-j-1]*dx #+4 keeps the centre in the write place - not sure why
    
    return arr3

def Voigt_convolution(x, nu_0, Gamma_G, Gamma_L):
    L = (Gamma_L/(2*np.pi))/((x - nu_0)**2 + (Gamma_L/2)**2)
    G = ((2*np.sqrt(np.log(2)/np.pi))/Gamma_G)*np.exp((-4*np.log(2)*(x - nu_0)**2)/(Gamma_G**2))
    V = convolution(L, G, x)
    return V

# Combined
def Combined(x, Triangle_amp, Voigt_amp, period, phaseshift, nu_0, Gamma_G, Gamma_L, res):
    #return Triangle_amp*TriangleWave(x, period, phaseshift) - Voigt_amp*Voigt_convolution(x, nu_0, Gamma_G, Gamma_L) + res
    return Triangle_amp*TriangleWave(x, period, phaseshift) - Voigt_amp*Voigt(x, nu_0, Gamma_G, Gamma_L) + res

# Chi^2 of functions to measure the quality of the fits
def Combined_Chi2(variables, data):
    # constants
    global period
    global phaseshift
    global Gamma_G
    global Gamma_L
    global data_error
    global Gaussian_nu_0
    
    # Variables
    Triangle_amp = variables[0]
    Voigt_amp = variables[1]
    nu_0 = variables[2]
    res = variables[3]
    
    # Upack data
    x = data[0]
    y_data = data[1]
    #y_errs = np.std(data[1]) + 0*data[1]
    y_errs = data_error + 0*data[1]
    
    # Get theory prediction for y from x, in other words f(x)
    y_theory = Combined(x, Triangle_amp, Voigt_amp, period, phaseshift, nu_0, Gamma_G, Gamma_L, res)
    
    # Check for values that are too high (would slow system down or even lead to infinities)
    max_theory = max(y_data)*1e9
    for i in range(len(y_theory)):
        if abs(y_theory[i]) > max_theory:
            y_theory[i] = max_theory*np.sign(y_theory[i])
    
    # Squared differences between data and theory, weighted for errors in y
    y_diff = (y_data - y_theory)**2 / (y_errs**2)
    
    # Multiply by gaussian for precise fit around particular area
    if Gaussian_nu_0 != 1.0:
        y_gaussian = Gaussian(x, Gaussian_nu_0, 0.5*np.std(x))
        y_diff = y_diff*y_gaussian
    
    # work out Chi^2
    z = np.sum(y_diff)
    #print('z = ', z)

    return z

# Reduced Chi2
def Combined_Reduced_Chi2(variables, data):
    z = Combined_Chi2(variables, data)
    return z/(len(data[0]) - (len(variables) +1)) # +1 for Gamma_L


######################################### Fitting ################################################

### Program to fit sections of measured triangle wave data to Voigt, Laurentzian and Gaussian functions
def Fit(data, temp, smoothing=1, cell_num=1, plotting=False, save_data=False, dirname_ext='', document_name='', resolution=600):
    global period
    global phaseshift
    global Gamma_G
    global Gamma_L
    global Triangle_amp
    global Voigt_amp
    global nu_0
    global res
    global Chi2
    global data_error
    global Gaussian_nu_0
    
    print('\nFitting Cell {}...'.format(str(cell_num)))
    
    ###############################################################################
    # Normalise data and save to array
    print('Normalising data...')
    for i in range(len(data)-1):
        data[i+1] = (data[i+1] - np.amin(data[i+1]))/(np.amax(data[i+1]) - np.amin(data[i+1]))
    print('* Data normalised *')
    
    ##############################################################################
    # Find period and phaseshift of triangle wave
    period, phaseshift = TriangleWave_Period_Phaseshift(data[0], data[2], smoothing)
    print('Cell ', cell_num, ' Triangle wave constants: period = ', period, ' s, phaseshift = ', phaseshift, 'rad')
    
    ################################################################################
    # Cut data to only that which we wish to fit
    data = Triangle_cut(data, period, phaseshift)
    
    ###############################################################################
    # Normalise data again and save to array
    print('Normalising data...')
    for i in range(len(data)-1):
        data[i+1] = (data[i+1] - np.amin(data[i+1]))/(np.amax(data[i+1]) - np.amin(data[i+1]))
    print('* Data normalised *')
    
    ###############################################################################
    # Flip calibration data - wires were plugged in the wrong way around
    ### !REMOVE WHEN WIRES PLACED CORRECTLY! ###
    data[3] = -1*(data[3] - 0.5) + 0.5
    
    ##############################################################################
    # Convert from seconds to Gigahertz and centering in middle
    print('Converting from s to GHz...')
    conversion, shift = Convert_t_to_f(data[0], data[2], data[3], smoothing, cell_num)
    
    # convert time data to frequency
    data[0] = data[0] - shift
    data[0] = data[0]*conversion
    
    # Convert triangle wave constants
    phaseshift = (phaseshift + 2*np.pi*shift/period) % (2*np.pi)
    period  = period*conversion
    print('Converted period =\t', period, ' GHz')
    print('converted phaseshift =\t', phaseshift, ' rad')
    print('*Conversion complete*')
    
    
    ###############################################################################
    # Data error - standard deviation of peak to peak noise
    # Requires a high pass filter to remove low frequency noise
    HP_frequency = 1 # Hz
    data_error = stdv_error(data[0], data[1], HP_frequency, period)
    print('data_error = ', data_error)
    
    
    ##############################################################################
    # Calculate Gamma_G from temperature for natural Rb
    kB = 1.381*10**(-23) #m^2.kg.s^-2.K^-1
    D1 = 794.8*10**(-9) #m
    m = 85.47*1.661*10**(-27) #kg
    T_CtoK = 273.15 #K
    Gamma_G = (2/D1)*np.sqrt(2*kB*(temp+T_CtoK)*np.log(2)/m)*10**(-9) #GHz
    print('Gamma_G = ', Gamma_G)
    # np.log() is ln()
    
    
    ##############################################################################
    # Fit modified data to a series of Voigt with various Gamma_L values
    # Initial coarse scan
    print('Coarse Voigt fitting...')
    
    # Add gaussian to fit to focus around the middle
    # where the absorption line should be if the physical experiment was done correctly
    Gaussian_nu_0 = data[0][0] + (data[0][-1] - data[0][0])/2
    print('Gaussian_nu_0, = ', Gaussian_nu_0)

    # values of Gamma_L to fit over
    n_points = 100
    Gamma_L_bnds = [1e-9, (data[0][-1] - data[0][0])/2]

    # Fitting variable limits
    Triangle_amp_bnds = [0.01, 1]
    Voigt_amp_bnds = [0.01, (np.pi/2)*(data[0][-1] - data[0][0])/4]
    nu_0_bnds = [data[0][0], data[0][-1]]
    res_bnds = [-1, 1]
    
    Gamma_L_array, Chi2_array, variables_array = fit_Gamma_L(data[0], data[1], n_points, Gamma_L_bnds, Triangle_amp_bnds, Voigt_amp_bnds, nu_0_bnds, res_bnds)
    
    # Best fit is the for the lowest value of Chi2
    #Chi2_array = Moving_Average(Chi2_array, smoothing)
    Chi2_i = np.where(Chi2_array == min(Chi2_array))[0][0]
    Gamma_L = Gamma_L_array[Chi2_i]
    Chi2 = Chi2_array[Chi2_i]
    variables = variables_array[Chi2_i]
    R_chi2 = Combined_Reduced_Chi2(variables, data)
    print('variables = ', variables)
    print('Gamma_L = ', Gamma_L)
    
    # Gamm_L error calculation
    Gamma_L_positive_std, Gamma_L_negative_std = fit_error(Gamma_L_array, Chi2_array, save_data, dirname_ext, document_name, cell_num, resolution)
    
    print('*** Coarse Voigt fitting complete ***')
    
    ##############################################################################
    # Fit modified data to a series of Voigt with various Gamma_L values
    # Second precise scan
    print('Precise Voigt fitting...')
    
    # update gaussian to focus on absorption line
    #Gaussian_nu_0 = variables[2]
    #print('Gaussian_nu_0, = ', Gaussian_nu_0)

    # values of Gamma_L to fit over
    n_points = 100
    Gamma_L_width = (Gamma_L_bnds[1] - Gamma_L_bnds[0])/5
    Gamma_L_bnds[0] = Gamma_L - Gamma_L_width
    Gamma_L_bnds[1] = Gamma_L + Gamma_L_width
    print('Gamma_L_bnds = ', Gamma_L_bnds)
    
    # Fitting variable limits
    Triangle_amp_bnds = [0.01, 1]
    Voigt_amp_bnds = [0.01, (np.pi/2)*(data[0][-1] - data[0][0])/4]
    nu_0_bnds = [data[0][0], data[0][-1]]
    res_bnds = [-1, 1]
    
    Gamma_L_array, Chi2_array, variables_array = fit_Gamma_L(data[0], data[1], n_points, Gamma_L_bnds, Triangle_amp_bnds, Voigt_amp_bnds, nu_0_bnds, res_bnds)
    
    # Best fit is the for the lowest value of Chi2
    #Chi2_array = Moving_Average(Chi2_array, smoothing)
    Chi2_i = np.where(Chi2_array == min(Chi2_array))[0][0]
    Gamma_L = Gamma_L_array[Chi2_i]
    Chi2 = Chi2_array[Chi2_i]
    variables = variables_array[Chi2_i]
    R_chi2 = Combined_Reduced_Chi2(variables, data)
    print('variables = ', variables)
    print('Gamma_L = ', Gamma_L)
    print('Combined_Reduced_Chi2 = ', R_chi2)
    
    # Gamm_L error calculation
    Gamma_L_positive_std, Gamma_L_negative_std = fit_error(Gamma_L_array, Chi2_array, save_data, dirname_ext, document_name, cell_num, resolution)
    
    print('*** Precise Voigt fitting complete ***')
    
    ##############################################################################
    # Calculate Pressure
    T_0 = 353 #K
    gamma_N2_T0 = 17.8 #GHz/amg
    T_standard = 25+T_CtoK #K
    amg = 44.615 #mol.m-3.amg
    R = 8.3145 #J.mol-1.K-1
    Torr = 101325/760 #Pa/Torr
    gamma_N2 = gamma_N2_T0*((temp+T_CtoK)/T_0)**0.3 #GHz/amg
    Pressure = amg*(Gamma_L/gamma_N2)*(temp+T_CtoK)*R/Torr #Torr
    #Calculate P at standard temp
    Pressure = T_standard*Pressure/(temp+T_CtoK) #Torr
    # Update errors
    Pressure_positive_std = Gamma_L_positive_std*(Pressure/Gamma_L)
    Pressure_negative_std = Gamma_L_negative_std*(Pressure/Gamma_L)
    print('Pressure = {} (+{} -{}) Torr'.format(Pressure, Pressure_positive_std, Pressure_negative_std))

    ##############################################################################
    # Plot best fit
    if plotting==True:
        plot(save_data, dirname_ext, document_name, resolution, cell_num, data, variables, Pressure, Pressure_positive_std, Pressure_negative_std)
    
    ##############################################################################
    # Return results
    constants = np.array([period, phaseshift, Gamma_G])
    calibration = np.array([conversion, shift])
    
    results = [cell_num, Pressure, Pressure_positive_std, Pressure_negative_std, R_chi2, Gamma_L]
    for i in range(len(variables)):
        results.append(variables[i])
    for i in range(len(constants)):
        results.append(constants[i])
    for i in range(len(calibration)):
        results.append(calibration[i])
    
    return results

# Calculate error from high frequency noise of data
def stdv_error(x, y, HP_frequency, period):
    Filtered_Data = HP_filter([x, y], HP_frequency)
    peak_data = [[], []]
    for i in range(len(x)):
        # Ignore first 1 GHz of data due to filter noise
        if np.abs(x[i]-x[0]) > 1:
            peak_data[0].append(Filtered_Data[0][i])
            peak_data[1].append(Filtered_Data[1][i])
    #Smoothed_Data = Moving_Average(Filtered_Data[1], 100)
    #peak_ids = find_peaks(Smoothed_Data,  prominence = (np.std(Smoothed_Data)))[0]
    peak_ids = find_peaks(peak_data[1],  prominence = 10*(np.std(peak_data[1])))[0]
    print('peak_ids = ', peak_ids)
    if len(peak_ids) > 0:
        #widths = peak_widths(Smoothed_Data, peak_ids, rel_height=0.5)[0]
        widths = peak_widths(peak_data[1], peak_ids, rel_height=0.5)[0]
        base_left_ids = []
        base_right_ids = []
        for i in range(len(peak_ids)):
            base_left_ids.append(peak_ids[i]-int(widths[i]+0.5))
            base_right_ids.append(peak_ids[i]+int(widths[i]+0.5))
            # limits to ids
            if base_left_ids[i] < 0:
                base_left_ids[i] = 0
            if base_right_ids[i] < 0:
                base_right_ids[i] = 0
            if base_left_ids[i] >= len(peak_data[0]):
                base_left_ids[i] = len(peak_data[0]) - 1
            if base_right_ids[i] >= len(peak_data[0]):
                base_right_ids[i] = len(peak_data[0]) - 1
        # Use data only from where there are no peaks
        no_peak_data = [[], []]
        counter = 1
        current_left = base_left_ids[0]
        current_right = base_right_ids[0]
        for i in range(len(peak_data[0])):
            # Update peak limits
            if i > current_right and len(base_left_ids) > counter:
                current_left = base_left_ids[counter]
                current_right = base_right_ids[counter]
                counter += 1
            if i < current_left or i > current_right:
                no_peak_data[0].append(peak_data[0][i])
                no_peak_data[1].append(peak_data[1][i])
        
        # peaks = [[],[]]
        # base_left = [[],[]]
        # base_right = [[],[]]
        # for i in range(len(peak_ids)):
        #     peaks[0].append(peak_data[0][peak_ids[i]])
        #     peaks[1].append(peak_data[1][peak_ids[i]])
        #     base_left[0].append(peak_data[0][base_left_ids[i]])
        #     base_left[1].append(peak_data[1][base_left_ids[i]])
        #     base_right[0].append(peak_data[0][base_right_ids[i]])
        #     base_right[1].append(peak_data[1][base_right_ids[i]])
        # plt.scatter(peaks[0], peaks[1], label='peaks', color='r')
        # plt.scatter(base_left[0], base_left[1], label='base_left', color='b')
        # plt.scatter(base_right[0], base_right[1], label='base_right', color='g')
        
        #Check if too much of data was discounted as peaks
        if len(peak_data[0]) > 2*len(no_peak_data[0]):
            noise_std = np.std(peak_data[1])
        else:
            noise_std = np.std(no_peak_data[1])
    else:
        noise_std = np.std(peak_data[1])
    # plt.plot(x, y, label='raw data', linewidth=0.1)
    # plt.plot(Filtered_Data[0], Filtered_Data[1], label='filtered data', linewidth=0.1)
    # plt.plot(peak_data[0], peak_data[1], label='no peak data', linewidth=0.1)
    # #plt.plot(no_peak_data[0], no_peak_data[1], label='no peak data', linewidth=0.1)
    # plt.title('HP filtering')
    # plt.xlabel('Frequency (GHz)')
    # plt.ylabel('Intensity')
    # #plt.ylim(-0.01, 0.01)
    # plt.legend(loc='best')
    # #plt.show()
    # plt.savefig('Batch1_Cell_2_error_plot.svg', dpi = 300)
    # plt.close()
    
    # Add error from overall drift
    #wavelength = np.array([data[0][i_start:i_end], best_fit[i_start:i_end]])
    # create 2 arrays, one wavelength long
    # wavelength_1 = []
    # wavelength_2 = []
    # for i in range(len(x)):
    #     if x[i] < x[0] + period:
    #         wavelength_1.append(y[i])
    #     elif x[i] < x[0] + 2*period:
    #         wavelength_2.append(y[i])
    # # Make them the same length
    # if len(wavelength_1) > len(wavelength_2):
    #     del wavelength_1[len(wavelength_2):-1]
    #     del wavelength_1[-1]
    # if len(wavelength_2) > len(wavelength_1):
    #     del wavelength_2[len(wavelength_1):-1]
    #     del wavelength_2[-1]
    # wavelength_1 = np.array(wavelength_1)
    # wavelength_2 = np.array(wavelength_2)
    # # calculate the standard deviation of the magnitude of the differences
    # differences = abs(wavelength_1 - wavelength_2)
    # difference_std = np.std(differences)
    
    return noise_std# + difference_std
    
# Find period and phaseshift of triangle wave data
def TriangleWave_Period_Phaseshift(x, y, smoothing):
    # Find period and phaseshift of triangle wave
    print('Finding laser scan period...')
    period_data = np.array([x, y])
    period_data[1] = (period_data[1] - ((np.amax(period_data[1]) - np.amin(period_data[1])))/2)**2 # square to get peaks and troughs
    period_data[1] = Moving_Average(Moving_Average(period_data[1], smoothing), smoothing) # Smooth to avoid false peaks
    
    # Select peaks only above 1 standard deviation
    period_data_std = np.std(period_data[1])
    period_data_peaks = find_peaks(period_data[1], prominence = (period_data_std/2))

    if len(period_data_peaks[0]) < 2:
        print('Unable to find at least 1 peak and trough in reference data')
        exit()
    period = 0
    for j in range(len(period_data_peaks[0])-1):
        period += period_data[0][period_data_peaks[0][j+1]] - period_data[0][period_data_peaks[0][j]]
    period = 2*period/(len(period_data_peaks[0])-1)
    if y[period_data_peaks[0][0]] > y[period_data_peaks[0][1]]:
        phaseshift = (0.5 - (x[period_data_peaks[0][0]]/(period/2)))*np.pi
    else:
        phaseshift = (0.5 - (x[period_data_peaks[0][1]]/(period/2)))*np.pi
    print('* Period found *')
    return period, phaseshift

# Cut data to central 1/2 wavelength of triangle wave
def Triangle_cut(data, period, phaseshift):
    print('Cutting data to only what is to be fit...')
    # Check for minimum and maximum points of triangle wave over a single wavelength in the middle of the data
    t_step = (data[0][len(data[0])-1] - data[0][0])/len(data[0])    # Average time step
    i_start = int((len(data[0]) - period/t_step)/2)
    i_end = int((len(data[0]) + period/t_step)/2)
    best_fit = 0.5*TriangleWave(data[0], period, phaseshift) + 0.5
    wavelength = np.array([data[0][i_start:i_end], best_fit[i_start:i_end]])  # Wavelength to check over
    i_min = np.where(wavelength[1] == min(wavelength[1]))[0][0]
    i_max = np.where(wavelength[1] == max(wavelength[1]))[0][0]
    if i_min < i_max:
        i_end = i_start + i_max
        i_start = i_start + i_min
    else:
        i_end = i_start + i_min
        i_start = i_start + i_max
    
    # Single section required
    data_selection = []
    for i in range(len(data)):
        data_selection.append(data[i][i_start:i_end])
    data = np.array(data_selection)
    print('*Data cut*')
    return data

# Determine time to frequency conversion factors
def Convert_t_to_f(t, y_r, y_c, smoothing, cell_num):
    # Using absorption lines a and b from calibration data of Rb87 reference cell:
    #/\                /\
    #  |              |  |
    #  |/\          /\|  |/\
    #     \        /  b     \
    #      |      |          |
    #      |/\  /\|          |/\  /
    #         \/   a            \/
    
    # Find peaks
    # Select peaks only above 1 standard deviation
    calibration_data = Moving_Average((y_r - y_c), smoothing)
    calibration_peaks = find_peaks(calibration_data, prominence = (np.std(calibration_data)))
    print('calibration_peaks = ', calibration_peaks)
    
    # Select the first and last peaks
    try:
        # Sort integers to select the two largest
        sorted_prominences = sorted(calibration_peaks[1]['prominences'])
        id_a = calibration_peaks[0][(np.where(calibration_peaks[1]['prominences'] == sorted_prominences[-1]))[0][0]]
        id_b = calibration_peaks[0][(np.where(calibration_peaks[1]['prominences'] == sorted_prominences[-2]))[0][0]]
        dt = np.abs(t[id_a] - t[id_b]) #s
    except:
        print('\n!!! No or only 1 peak found for cell {} !!!\n'.format(cell_num))
        exit()
        # dt = 0.01 # s
    print('dt = ', dt)
    
    # calibration_peaks = [[t[id_a], t[id_b]], [calibration_data[id_a], calibration_data[id_b]]]
    # plt.scatter(calibration_peaks[0], calibration_peaks[1], label='calibration_peaks', color='red')
    # plt.plot(t, calibration_data, label='Data', color='green')
    # plt.show()
        
    # conversion factors
    df = 7.647 # GHz
    conversion = df/dt
    shift = t[int((len(t)-0.5)/2)] #-0.5=-1+0.5, due to truncating
    print('df/dt =\t', conversion, ' GHz/s')
    print('shift =\t', shift, ' s')
    return conversion, shift 

# calculate minimum Chi squared over multiple Gamma_L values
def fit_Gamma_L(x, y, n_points, Gamma_L_bnds, Triangle_amp_bnds, Voigt_amp_bnds, nu_0_bnds, res_bnds):
    global Gamma_L
    
    # values of Gamma_L to fit over
    Gamma_L_array = np.arange(Gamma_L_bnds[0], Gamma_L_bnds[1], (Gamma_L_bnds[1] - Gamma_L_bnds[0])/n_points)

    # Fitting variable limits
    variable_limits = np.array([Triangle_amp_bnds, Voigt_amp_bnds, nu_0_bnds, res_bnds])
    
    # Fitting
    variable_init = np.zeros(len(variable_limits), dtype=float)
    for i in range(len(variable_init)):
        # Start with variables halfway between limits
        variable_init[i] = variable_limits[i][0] + (variable_limits[i][1] - variable_limits[i][0])/2
    variables_array = np.zeros((len(Gamma_L_array), len(variable_init)), dtype=float)
    Chi2_array = np.zeros(len(Gamma_L_array), dtype=float)
    for i in range(len(Gamma_L_array)):
        Gamma_L = Gamma_L_array[i]
        theory = optimize.minimize(Combined_Chi2,  variable_init, [x, y], method="Powell", bounds=variable_limits)
        variables_array[i] = theory["x"]
        Chi2_array[i] = Combined_Chi2(variables_array[i], [x, y])
        #print('Chi2 = ', Chi2_array[i])
    
    return Gamma_L_array, Chi2_array, variables_array

# fit error using using Chi squared
def fit_error(Gamma_L_array, Chi2_array, save_data, dirname_ext, document_name, cell_num, resolution):
    # Gamm_L error calculation
    
    #        Chi2(Gamma_L_array)
    #             ^
    #             |     |           |
    #             |     |           |
    #  min(Chi2)+1|.....|...........|
    #             |     :\         /:
    #             |     : |       | :
    #             |     : |       | :
    #             |     :  \     /  :
    #    min(Chi2)|.....:...\___/   :
    #             |     :     :     :
    #             |_____:_____:_____:___________> Gamma_L_array
    #                   :  Gamma_L  :
    #                   :           :
    #                   :           :
    #                   :<--------->:
    #                   Chi2 interval
    
    Gamma_L_negative_std = -1
    Gamma_L_positive_std = -1
    # Run through Chi2 values until first values under Chi2+1 found
    for i in range(int(len(Gamma_L_array))):
        if Chi2_array[i] < (Chi2 + 1) and Gamma_L_negative_std == -1:
            if i == 0:
                Gamma_L_negative_std = Gamma_L - Gamma_L_array[i]
            else:
                # Find  exact location of intercept assuming stright lines between Chi2 points
                Gamma_L_negative = Gamma_L_array[i] - (Gamma_L_array[i] - Gamma_L_array[i-1])*(Chi2_array[i] - (Chi2 + 1))/(Chi2_array[i] - Chi2_array[i-1])
                Gamma_L_negative_std = Gamma_L - Gamma_L_negative
        if Chi2_array[-i-1] < (Chi2 + 1) and Gamma_L_positive_std == -1:
            if i == 0:
                Gamma_L_positive_std = Gamma_L_array[-i-1] - Gamma_L
            else:
                # Find  exact location of intercept assuming stright lines between Chi2 points
                Gamma_L_positive = Gamma_L_array[-i-1] + (Gamma_L_array[-i] - Gamma_L_array[-i-1])*((Chi2 + 1) - Chi2_array[-i-1])/(Chi2_array[-i] - Chi2_array[-i-1])
                Gamma_L_positive_std = Gamma_L_positive - Gamma_L
    print('Gamma_L_positive_std = ', Gamma_L_positive_std)
    print('Gamma_L_negative_std = ', Gamma_L_negative_std)
    
    plt.scatter(Gamma_L_array, Chi2_array, marker='+')
    plt.title('Chi2 plot')
    plt.xlabel('Gamma_L')
    plt.ylabel('Chi2')
    std_neg = Gamma_L - Gamma_L_negative_std
    std_pos = Gamma_L + Gamma_L_positive_std
    plt.plot([min(Gamma_L_array), max(Gamma_L_array)], [Chi2, Chi2], linestyle=':', color='r')
    plt.plot([Gamma_L, Gamma_L], [0, np.max(Chi2_array)], linestyle=':', color='r')
    plt.plot([min(Gamma_L_array), max(Gamma_L_array)], [Chi2 + 1, Chi2 + 1], linestyle=':', color='g')
    plt.plot([std_pos, std_pos], [0, np.max(Chi2_array)], linestyle=':', color='g')
    plt.plot([std_neg, std_neg], [0, np.max(Chi2_array)], linestyle=':', color='g')
    #plt.xlim(std_neg - 3*Gamma_L_negative_std, std_pos + 3*Gamma_L_positive_std)
    #plt.xlim(4.5, 6.5)
    plt.ylim(Chi2-1, Chi2+5)
    if save_data==True:
        plt.savefig(dirname_ext+'\\'+document_name+'_Cell_{}_Chi2_plot.svg'.format(cell_num), dpi = resolution)
    else:
        plt.show()
    plt.close()
    
    return Gamma_L_positive_std, Gamma_L_negative_std
    
######################################## Extras #####################################################

def plot(save_data, dirname_ext, document_name, resolution, cell_num, data, variables, Pressure, Pressure_positive_std, Pressure_negative_std):

    
    Triangle_amp = variables[0]
    Voigt_amp = variables[1]
    nu_0 = variables[2]
    global Gamma_L
    res = variables[3]

    x = data[0]
    y_data = data[1]

    plt.plot(x, y_data, label='Data')
    plt.plot(x, Combined(x, Triangle_amp, Voigt_amp, period, phaseshift, nu_0, Gamma_G, Gamma_L, res), linestyle='--', label='Best fit')

    plt.legend(loc='best')
    plt.title('Vars: T_amp = {:.3f}, V_amp = {:.3f}, $ν_0$ = {:.3f}, $Γ_L$ = {:.3f}, res = {:.3f}'.format(Triangle_amp, Voigt_amp, nu_0, Gamma_L, res))
    plt.suptitle('P = {:.3f} (+{:.3f} -{:.3f}) Torr\tConsts: λ = {:.3f}, θ = {:.3f}, $Γ_G$ = {:.3f}'.format(Pressure, Pressure_positive_std, Pressure_negative_std, period, phaseshift, Gamma_G))
    if save_data==True:
        plt.savefig(dirname_ext+'\\'+document_name+'_Cell_{}_Voigt_fit_plot.svg'.format(cell_num), dpi = resolution)
    else:
        plt.show()
    plt.close()

# Moving average function of data. Averaging over number of data points given by average_range
def Moving_Average(data, average_range):
    # Only works averaging over an odd number of points, so add 1 to even numbers
    if (average_range % 2) == 0:
        average_range += 1

    # Loop to add up data within range and divide by range for every point
    new_data = data.copy()
    for i in range(len(data)):
        data_sum = 0
        counter = 0
        for j in range(average_range):
            # Exclude impossible indeces
            half_average = int((average_range-1)/2)
            if (i+j-half_average) >= 0 and (i+j-half_average) < len(data):
                data_sum += data[i+j-half_average]
                counter += 1
        new_data[i] = data_sum/counter
    
    return new_data

# high pass filter
def HP_filter(Data, HP_frequency):
    fN = (len(Data[0])/(2*(Data[0][-1] - Data[0][0]))) #Nyquist frequency (twice the interval between samples)
    ### Set up filters
    # Create high pass filter Butterworth filter
    N_HP = 2
    Wn_HP = HP_frequency/fN #Normalised from 0 to 1, where 1 is the Nyquist frequency
    sos_HP = signal.butter(N_HP, Wn_HP, 'highpass', output='sos')
    
    ### Filter Data
    Filtered_Data = signal.sosfilt(sos_HP, Data)
    Filtered_Data[0] = Data[0]
    
    # Plot filtered data
    # if plotting == True:
    #     title = 'Filtered Data'
    #     File_name = 'Filtered_plot'
    #     plot_filter(Filtered_Data, dirname_ext, HP_frequency, fN, sos_HP, N_HP, title, File_name)
   
    return Filtered_Data

### Plot filter
def plot_filter(Data, dirname_ext, HP_frequency, fN, sos, N, title, File_name):
    # plot HP fliter
    w, h = signal.sosfreqz(sos, worN=1024)
    plt.semilogx(fN*w/(np.pi), 20 * np.log10(np.maximum(abs(h), 1e-5)), label='N={} order HP'.format(N), color='blue')
    plt.axvline(HP_frequency, color='purple', label='{}Hz HP Cutoff'.format(HP_frequency))
    
    Fourier_Data = list(Data)
    print('Fourier_Data = ', Fourier_Data)
    Fourier_Data[0] = np.fft.fftfreq(Data[0].shape[-1], d=((Data[0][-1] - Data[0][0])/len(Data[0])))
    Fourier_Data[1] = np.fft.fft(Fourier_Data[1])
    #plt.semilogx(Fourier_Data[0],Fourier_Data[1], color='cyan')
    # Plot settings
    plt.title(title)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (dB)')
    #plt.margins(0, 0.1)
    plt.grid(which='both', axis='both')
    plt.legend(loc='best')
    #plt.savefig(dirname_ext+'\\'+File_name+'.svg', dpi = 300)
    plt.show()
    plt.close()
    return 0


############################################ Main ##################################################
# Only run if not called from another file
if __name__ == '__main__':
    
    ###### Parameters ######
    ### Overall
    cell_num = 71 # Number of cells to scan or files to read if fitting with no scanning
    save_data=True # Save the data to file
    plotting=True # Plot data
    resolution = 600 #dpi
    smoothing = 20
    local_path = 'C:\\Users\\lp472\\Desktop\\School\\Uni\\Doctorate\\Lab\\Cell_Scanner\\Scans\\'
    #local_path = '/home/leightpage/' # Path to saved data
    document_name = 'Measurements'
    #dirname_descriptor = '2023-07-13_16-18-22_Batch1'
    dirname_descriptor = '2025-01-16_11-56-35_Measurements'
    dirname_ext = local_path + dirname_descriptor

    # Read in data
    cell_data = np.array(pd.read_csv(dirname_ext+'\\'+document_name+'_Cell_'+str(cell_num)+'.csv', header=None, skiprows=11, engine='python', delimiter='\t'))
    temp = np.array([pd.read_csv(dirname_ext+'\\'+document_name+'_Cell_'+str(cell_num)+'.csv', header=None, skiprows=4, skipfooter=len(cell_data)+6, engine='python', delimiter='\t')[1][0],
                     pd.read_csv(dirname_ext+'\\'+document_name+'_Cell_'+str(cell_num)+'.csv', header=None, skiprows=5, skipfooter=len(cell_data)+5, engine='python', delimiter='\t')[1][0],
                     pd.read_csv(dirname_ext+'\\'+document_name+'_Cell_'+str(cell_num)+'.csv', header=None, skiprows=6, skipfooter=len(cell_data)+4, engine='python', delimiter='\t')[1][0]])
    # Transform data from line arrays to column arrays
    data = []
    for i in range(len(cell_data[0])):
        data.append([])
        for j in range(len(cell_data)):
            data[i].append(cell_data[j][i])
    data = np.array(data)
    temp = np.array(temp, dtype=float)
    
    # Fit data
    results = Fit(data, temp[0], smoothing, cell_num, plotting, save_data, dirname_ext, document_name, resolution)
