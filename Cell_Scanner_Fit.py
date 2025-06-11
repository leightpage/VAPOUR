"""
Created on 13/07/2022
@author: Leigh Page, lp472@sussex.ac.uk

Description:
Program to automate the measure the relative causes of absorption broadening in multiple Rb vapour-cells,
by scanning of an area with a Rasberry Pi camera and Labjack voltmeter using a COSI moving them, and then fitting the data to models.
Works in the following steps:
- Open COSI, Rasberry Pi and Labjack connections
- Home COSI
- Allow manual movement or automatic
- Callibrate corner coordinates with user feedback
- Move point to point in a loop
- Record laser transmissions with photodiodes and Labjack
- Take picures of each cell with Raspberry Pi
- Save all data in same folder
- Close connections
- Perform function fitting for each cell in turn (bounds and initial values of function constants may need to be changed according to the data)
- Save fitting results and plots
This program is modular, allowing many options to be turned off and on in the main() functuon, to suit the needs of the user

Requirement:
Use Paramiko-expect v2.9

Required in Rasberry Pi, to take image:
RPi_Camera.py

Editting from the following:
COSI_SSH_v9_micah.py
COSI-MEASURE MICAHworking.py
LabJack_write_read_loop_with_config.py

raspberry-pi-camera-guide.pdf
raspberry_pi_eduroam_2019.pdf
"""
import time
import numpy as np
import pandas as pd
import os
from Cell_Scanner import scan, save_file
from Minimise_Fit import Fit


######################################### MAIN RUN FUNCTION ############################################
def main():
    print('Start of program')

    ###### Parameters ######
    ### Overall
    save_data=True # Save the data to file
    plotting=True # Plot data
    cell_num = 2 # Number of cells to scan or files to read if fitting with no scanning
    resolution = 600 #dpi
    local_path = 'C:\\Users\\lp472\\Desktop\\School\\Uni\\Doctorate\\Lab\\Cell_Scanner\\Scans\\'
    document_name = 'Measurements'
    #dirname_descriptor = '2023-07-13_16-18-22_Batch1' # Folder to read in data from if there is no scanning
    dirname_dates = ['2025-01-23_17-05-42'] # Date and time of data measurements
    
    for a in range(len(dirname_dates)):
        dirname_descriptor = dirname_dates[a]+'_'+document_name
        ### Scan parameters
        # COSI parameters
        heating=False    #Allow heating of plate under COSI
        scanning=False          #Allow scanning over designated area
        RPi=False # Pictures with Rasberry Pi
        calibrating=False # If calibrating of the COSI coordinates needed
        # Labjack scan channels
        scan_channels = [0, 2, 1] # Channels AINx (cell data, calibration data and reference data)
        temp_channels = [3]
        ##########################################temp_channels = [3, ljm.constants.GND] # AIN3 to Ground
        # Measurement settings
        measure_time = 0.4    # seconds - Scan time
        scanRate = 10000 # Hz
        set_temp = '80' # degrees Celcius
        # Corner coordinates (x,y,z)
        #   ----C
        #   ¦   ¦
        #   B---A
        A = np.array([268.5, 147, 33.5])
        B = np.array([402.5, 146.5, 34])
        C = np.array([269, 279.5, 32.5])
        corners = [A,B,C]
        RPi_A = np.array([262, 312, 42])
        # calibrating maximum range around starting point allowed (mm)
        xy_range = 10
        z_range = 2
        # calibrating min gap between sampled points (mm)
        xy_step = 2
        z_step = 1
        probe_specs = [xy_range, xy_step, z_range, z_step]
        
        
        ### Fitting parameters
        fitting=True
        smoothing = 20 # Moving average range over data for better fitting
        
        
        
        
        ###### Operation ######
        ### Data Location
        if save_data==True and scanning==True:
            print('Create file directory')
            today_datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
            dirname_ext = local_path + today_datetime + '_' + document_name
            ### Create directory for save files
            os.makedirs(dirname_ext)
        else:
            dirname_ext = local_path + dirname_descriptor
    
        
        ### Scan
        if scanning==True:
            data, coords, temp = scan(dirname_ext, document_name, scan_channels, temp_channels, save_data, heating, calibrating, RPi, cell_num, measure_time, scanRate, set_temp, corners, RPi_A, probe_specs)
        ### Read in data
        else:
            print('Reading in data...')
            data = []
            temp = []
            for i in range(cell_num):
                cell_data = np.array(pd.read_csv(dirname_ext+'/'+document_name+'_Cell_'+str(i+1)+'.csv', header=None, skiprows=11, engine='python', delimiter='\t'))
                temp.append([pd.read_csv(dirname_ext+'/'+document_name+'_Cell_'+str(i+1)+'.csv', header=None, skiprows=4, skipfooter=len(cell_data)+6, engine='python', delimiter='\t')[1][0],
                            pd.read_csv(dirname_ext+'/'+document_name+'_Cell_'+str(i+1)+'.csv', header=None, skiprows=5, skipfooter=len(cell_data)+5, engine='python', delimiter='\t')[1][0],
                            pd.read_csv(dirname_ext+'/'+document_name+'_Cell_'+str(i+1)+'.csv', header=None, skiprows=6, skipfooter=len(cell_data)+4, engine='python', delimiter='\t')[1][0]])
                # Transform data from line arrays to column arrays
                data.append([])
                for j in range(len(cell_data[0])):
                    data[i].append([])
                    for k in range(len(cell_data)):
                        data[i][j].append(cell_data[k][j])
            data = np.array(data)
            temp = np.array(temp, dtype=float)
            print('* Data read *')
    
    
        ### Fit data to functions
        if fitting==True:
            #triangle_fit, fitting_data, fit_data, const = Fit(data, Voigt_bnds, Voigt_init, cell_num, temp)
            results = []
            for i in range(cell_num):
                # Fit data
                results.append(Fit(data[i], temp[i][0], smoothing, i+1, plotting, save_data, dirname_ext, document_name, resolution))
    
            # Save results
            # constant names
            result_names = ['cell_num', 'Pressure(Torr)', 'Pressure_positive_std(Torr)', 'Pressure_negative_std(Torr)', 'Reduced_Chi2', 'Gamma_L(GHz)', 'Triangle_Amp', 'Voigt_Amp', 'nu_0(GHz)', 'res', 'period(GHz)', 'Phaseshift(rad)', 'Gamme_G(GHz)', 'Conversion(GHz/s)', 'shift(s)']
            ### Save fitting constants
            if save_data==True:
                header = result_names[0]
                for i in range(len(result_names)-1):
                    header = header+'\t'+result_names[i+1]
                header = header+'\n'
                # Change shape of constant array for saving
                new_results = []
                for i in range(len(results[0])):
                    column = []
                    for j in range(len(results)):
                        column.append(results[j][i])
                    new_results.append(column)
                save_file(dirname_ext, document_name+'_Fitting_Constants', new_results, header)
        
    print('End of program')
    

############################################ Call run function ########################################

main()
