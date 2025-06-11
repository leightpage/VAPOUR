"""
Created on 15/09/2023
@author: Leigh Page, lp472@sussex.ac.uk

Description:
Program to automate the scanning of multiple Rb vapour-cells, by scanning of an area with a Rasberry Pi camera and Labjack voltmeter
using a COSI measure to move them.
Works in the following steps:
- Open COSI, Rasberry Pi and Labjack connections
- Home COSI
- Calibrate corner coordinates with user feedback
- Move point to point in a loop
- Record laser transmissions with photodiodes and Labjack
- Take picures of each cell with Raspberry Pi
- Save all data in same folder
- Close connections

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

import os
import time
import traceback
import paramiko
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from labjack import ljm
from paramiko_expect import SSHClientInteraction


########################################### DATA MANAGEMENT #############################################

### Saves data to .csv file
def save_file(dirname_ext, document_name, data, header):
    doc = document_name+'.csv'
    save_file = open(dirname_ext + '/' + doc, 'w')
    save_file.write(header)
    # write data
    for i in range(len(data[0])):
        for j in range(len(data)):
            save_file.write(str(data[j][i]))
            # add commas between values except end of line
            if j < (len(data)-1):
                save_file.write('\t')
        save_file.write('\n')
    save_file.close()
    
def plot(x, y, xlabel='', ylabel='', title='Data', colour='red', suptitle=False, save_data=False, dirname_ext='', document_name='Data', resolution=600):
    plt.plot(x, y, color=colour)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best')
    if suptitle!=False:
        plt.suptitle(suptitle)
    plt.title(title)
    if save_data==True:
        plt.savefig(dirname_ext+'/'+document_name, dpi = resolution)
    else:
        plt.show()
    plt.close()

### heatmap of signal quality
def heatmap(dirname_ext, save_data, probes, heat, title, subtitle, name):
    x = []
    y = []
    z = heat
    for i in range(len(probes)):
        x.append(round(probes[i][0], 3))
        y.append(round(probes[i][1], 3))
    heatmap_data = pd.DataFrame({'X (mm)': x, 'Y (mm)': y, 'V': z})
    #data_pivoted = heatmap_data.pivot('X (mm)', 'Y (mm)', 'V')
    data_pivoted = heatmap_data.pivot(index='X (mm)', columns='Y (mm)', values='V')
    ax = sns.heatmap(data_pivoted)#, vmin=0.8)
    plt.suptitle(title)
    plt.title(subtitle)
    if save_data==True:
        plt.savefig(dirname_ext + '/Corner_'+name+'_xy_Heatmap_.svg')
    else:
        plt.show()
    plt.close()


########################################## COSI FUNCTIONS ##############################################

### The COSI move function
def g1_move(g1_string, interact_COSI):
    print('Moving to '+g1_string+'...')

    interact_COSI.send(g1_string)
    interact_COSI.expect('\nMovement Complete\n', timeout=200)
    interact_COSI.current_output_clean

    print('* Moved to '+g1_string+' *')

### Creates an array of coordinates to scan over based on the number of cells to scan
def GCODE(cell_num, corners):
    #Scan paremeters (millimeters):
    #
    #         /\  z
    #         |
    #         |
    #         |
    #         |
    #       A |_ _ _ _ _ C_\  y
    #        /             /
    #       /
    #      /
    #  B |/_  x
    #
    # Corner coordinates (x,y,z)
    A = corners[0]
    B = corners[1]
    C = corners[2]
    max_cell_line = 12
    cell_gaps = max_cell_line - 1

    # Check if too many cells
    if cell_num > (max_cell_line*max_cell_line):
        print('\n\n !!! CELL NUMBER TOO BIG !!! \n\n')

    coords = A  #initialise the array
    j = 1
    k = 0
    for i in range(cell_num-1):
        # End of the line assuming a square setup - insert specific number otherwise
        if j > cell_gaps:
            k += 1
            j = 0
        # Add coordinates to array
        sign_x = np.sign(np.sign(B[0]-A[0])*(j*(B[0]-A[0]))**2 + np.sign(C[0]-A[0])*(k*(C[0]-A[0]))**2)
        sign_y = np.sign(np.sign(B[1]-A[1])*(j*(B[1]-A[1]))**2 + np.sign(C[1]-A[1])*(k*(C[1]-A[1]))**2)
        sign_z = np.sign(np.sign(B[2]-A[2])*(j*(B[2]-A[2]))**2 + np.sign(C[2]-A[2])*(k*(C[2]-A[2]))**2)
        x = sign_x*np.sqrt(((j*(B[0]-A[0]))**2 + (k*(C[0]-A[0]))**2)/cell_gaps**2) + A[0]
        y = sign_y*np.sqrt(((j*(B[1]-A[1]))**2 + (k*(C[1]-A[1]))**2)/cell_gaps**2) + A[1]
        z = sign_z*np.sqrt(((j*(B[2]-A[2]))**2 + (k*(C[2]-A[2]))**2)/cell_gaps**2) + A[2]
        coords = np.vstack((coords,np.array([x, y, z])))
        j += 1

    return coords

### Safely homes and sets the COSI arm to the starting position
def home_xyz(interact_COSI):
    print('homing ...')
    
    # Starting values to move arm safely around COSI structure while homing
    x = 200
    y = 150
    z = 10

    #Home and move z
    g1_move('g161z', interact_COSI)
    g1_move('g1z{}'.format(z), interact_COSI)

    #Home and move y
    g1_move('g161y', interact_COSI)
    g1_move('g1y{}'.format(y), interact_COSI)
    
    #Home and move x
    g1_move('g161x', interact_COSI)
    g1_move('g1x{}'.format(x), interact_COSI)
    print('*** Homing Complete ***')

### Open connection to COSI
def open_connect_COSI():
    print('Connecting to COSI ...')
    #SSH Parameters for COSI
    #Set login credentials and the server prompt
    HOSTNAME = '10.42.0.11'
    USERNAME = 'root'
    PASSWORD = 'root'
    PROMPT_COSI = 'root@beaglebone:~# '
    
    # Use SSH client to login
    try:
        # Create a new SSH client object
        client_COSI = paramiko.SSHClient()

        # Set SSH key parameters to auto accept unknown hosts
        client_COSI.load_system_host_keys()
        client_COSI.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        # Connect to the host
        client_COSI.connect(hostname=HOSTNAME, username=USERNAME, password=PASSWORD)
        print('*** COSI Connection opened ***')
    except Exception:
        traceback.print_exc()
    
    return PROMPT_COSI, client_COSI

#### COSI close connection
def close_connect_COSI(client_COSI):
    try:
        client_COSI.close()
        print('*** COSI Connection closed ***')
    except Exception:
        pass

### Scan area with COSI, saving images from Rasberry Pi and returning data from Labjack
def scan(dirname_ext, document_name, scan_channels, temp_channels, save_data, heating, calirating, RPi, cell_num, measure_time, scanRate, set_temp, corners, RPi_A, probe_specs):
    temp_channels.append(ljm.constants.GND) # Add the ground channel for temperature
    
    ### Open Connections
    print('SETTING UP CONNECTIONS')
    
    # Open connection to RPi
    if RPi==True and save_data==True:
        client_pi = open_connect_RPi()
    
    # Open connection to Labjack
    handle = open_connect_LJ()
    
    # Open connections to COSI
    PROMPT_COSI, client_COSI = open_connect_COSI()

    ### COSI moving
    print('COSI MOVING')
    
    # Create a client interaction class which will interact with the host
    with SSHClientInteraction(client_COSI, timeout=30, display=True) as interact_COSI:
        interact_COSI.expect(PROMPT_COSI)

        # Run the first command and capture the cleaned output, if you want
        # the output without cleaning, simply grab current_output instead.
        interact_COSI.send('cd /home/machinekit/BeBoPr-master')
        PROMPT_COSI='root@beaglebone:/home/machinekit/BeBoPr-master# '
        interact_COSI.expect(PROMPT_COSI, timeout=60)
        #cmd_output_uname = interact_COSI.current_output_clean
        interact_COSI.send('./mendel.elf')
        print('mendel.elf opened ...')
        interact_COSI.expect([PROMPT_COSI, 'Starting main loop...'])
        if interact_COSI.last_match == 'Starting main loop...':
            
            ### Scanning ###
            # Heat the plate under the COSI
            if heating==True:
                HeatPlate(set_temp, interact_COSI, handle, temp_channels)
            
            # Automatically home COSI
            home_xyz(interact_COSI)
            
            # Prime Labjack with initial run (first labjack measurement always has dodgy data around t=0)
            LJM_array = LJ_scan(handle, scan_channels, scanRate, 0.5)

            # Find accurate corner positions
            if calirating==True:
                # Probe the corners
                callibration=False
                while callibration==False:
                    print('Callibrating the position coordinates...')
                    probes = []
                    for i in range(len(corners)):
                        probe = prober(dirname_ext, save_data, corners[i], probe_specs, handle, scan_channels, scanRate, measure_time, interact_COSI, i)
                        print('probe = ', probe)
                        probes.append(probe)
                    print('\nOriginal corner coordinates:\nA = {}\nB = {}\nC = {}'.format(corners[0], corners[1], corners[2]))
                    print('\nProbing corner coordinates:\nA = {}\nB = {}\nC = {}'.format(probes[0], probes[1], probes[2]))
                    print('\nEnter desired corner coordinates (mm):')
                    names = ['A','B','C']
                    cartesian = ['x','y','z']
                    for i in range(len(corners)):
                        for j in range(len(corners[i])):
                            answered=False
                            while answered==False:
                                answer = float(input('{}[{}] = '.format(names[i], cartesian[j])))
                                if answer>0 and answer<600:
                                    corners[i][j] = answer
                                    answered=True
                                else:
                                    print('Xx Invalid answer xX')
                    
                    # Do you wish to repeat calirating?
                    answered=False
                    while answered==False:
                        answer = input('Do you wish to repeat calirating? (y/n):')
                        if answer=='y' or answer=='yes' or answer=='Y' or answer=='YES':
                            answered=True
                        elif answer=='n' or answer=='no' or answer=='N' or answer=='NO':
                            answered=True
                            callibration=True
                        else:
                            print('Xx Invalid answer xX')
                    


                print('*** Callibration Complete ***')
            # Generate GCODE/Gcode parameters
            coords = GCODE(cell_num, corners)

            # Scanning
            # Prepare 3D array to record all data and 2D array for coordinates & set_temperatures
            data = []
            temp = []
            print('Scanning ...')
            for i in range(len(coords)):
                # Move to gcode position
                g1_move('g1x{}g1y{}g1z{}'.format(str(coords[i][0]), str(coords[i][1]), str(coords[i][2])), interact_COSI)
                
                #Scan with Labjack
                time.sleep(5)  #Wait for vibrations from movement to subside
                LJM_array = LJ_scan(handle, scan_channels, scanRate, measure_time)
                # rearrange and invert data
                data.append([])
                # append time
                data[i].append([])
                for k in range(len(LJM_array)):
                    data[i][0].append(LJM_array[k][len(LJM_array[0])-1])
                    # append and invert data
                for j in range(len(LJM_array[0])-1):
                    data[i].append([])
                    for k in range(len(LJM_array)):
                        data[i][j+1].append(-1*LJM_array[k][j])
                data[i] = np.array(data[i])
                print('\nData = ', data[i], '\n')
                
                #Measure set_temperature from COSI and Labjack
                temp.append([[], [], []])
                temp[i][0] = float(TempPlate(interact_COSI))
                temp[i][1], temp[i][2] = LJ_temp(handle, temp_channels)
                
                #Convert lists to arrays for easier accessibility
                #temp = np.array(temp, dtype=float)
                
                print('*** Scanning Complete ***')
                
                #### save scan data
                if save_data==True:
                    settings = 'Cell_number=\t'+str(i+1)+'\nSample_Spacing(ms)=\t'+str(1/scanRate)+'\nSettle_Time(s)=\t'+str(measure_time)+'\n'
                    if heating==False:
                        set_temp = 'N/A'
                    header = settings+'Set_temperature(C)=\t{}\nBed_temperature(C)=\t{}\nThermocouple_temperature(C)=\t{}\nTemp_std(C)=\t{}\n'.format(set_temp, temp[i][0], temp[i][1], temp[i][2])
                    header = header+'Coordinate_X(mm):\t{}\nCoordinate_Y(mm):\t{}\nCoordinate_Z(mm):\t{}\n'.format(coords[i][0], coords[i][1], coords[i][2])
                    header = header+'Time(ms)'
                    for j in range(len(scan_channels)):
                        header = header+'\tAIN{}(V)'.format(j)
                    header = header+'\n'
                    save_file(dirname_ext, document_name+'_Cell_'+str(i+1), data[i], header)
                
                ### Plot data
                if plotting==True:
                    xlabel = 'Time (s)'
                    ylabel = 'Voltage (V)'
                    labels =  ['Scan', 'Reference', 'Calibration']
                    colours = ['red', 'blue', 'green']
                    suptitle = 'Cell {} Data'.format(i+1)
                    title = 'X={}mm, Y={}mm, Z={}mm'.format(coords[i][0], coords[i][1], coords[i][2])
                    
                    for j in range(len(data[i])-1):
                        plot(data[i][0], data[i][j+1], xlabel, ylabel, title, colours[j], suptitle, save_data, dirname_ext, document_name+'_Cell_{}_{}_plot.png'.format(i+1, labels[j]), resolution)
                
            
            ### Imaging ###
            if RPi==True and save_data==True:
                print('Imaging ...')
                # Convert scanning coordinates to imaging coordinates
                difference = RPi_A - corners[0]
                RPi_coords = coords + difference
                
                for i in range(len(RPi_coords)):
                    # Move to rough gcode position
                    g1_move('g1x{}g1y{}g1z{}'.format(str(RPi_coords[i][0]), str(RPi_coords[i][1]), str(RPi_coords[i][2])), interact_COSI)
                    
                    #Take picture with RPi
                    image_RPi(dirname_ext, document_name, client_pi, (i+1))
                print('*** Imaging Complete ***')
        else:
            # If problem with ./mendel.elf - close program
            print('./mendel.elf could not be started!')
            print('Ending Connection')
        
    # Close mendel.elf
    #interact_COSI.send('^C') ################################################################################################################
    #print('***mendel.elf closed ***')

    ### Close connections
    print('CLOSING CONNECTIONS')

    # Close COSI connection
    close_connect_COSI(client_COSI)
    
    # Close connection to RPi
    if RPi==True and save_data==True:
        close_connect_RPi(client_pi)
    
    # Close connection to Labjack
    close_connect_LJ(handle)

    # Return data
    return data, coords, temp


###################################### COSI TEMPERATURE FUNCTIONS ##########################################

### Heat the plate below the COSI
def HeatPlate(set_temp, interact_COSI, handle, temp_channels):
    print('Heating plate...')

    interact_COSI.send('m190 s'+set_temp)
    interact_COSI.expect("temperature for 'temp_bed' has stabilized", timeout=500)

    interact_COSI.current_output_clean

    print('*** Heating Complete ***')

### Returns the COSI plate temperature
def TempPlate(interact_COSI):
    print('Measuring plate temperature...')

    interact_COSI.send('m105')
    interact_COSI.expect('T:777.0 B:.*')    # Maybe use a delay time if this takes too long
    ConsoleOutput = interact_COSI.current_output
    a = ConsoleOutput.split(':')
    b = a[-1].rstrip('\n')
    interact_COSI.current_output_clean  #Clears current output

    print('*** temperature Measuring Complete ***')

    return str(b)


###################################### RASBERRY PI FUNCTIONS ##############################################

#### Rasberry Pi open connection
def open_connect_RPi():
    print('Connecting to Raspberry Pi ...')
    ######### SSH Terminal for RPi ##########
    HOSTNAME = '10.42.0.14'
    USERNAME = 'pi'
    PASSWORD = 'Physics1234'
    PROMPT = 'pi@10.42.0.14:~# '

    # Use SSH client to login
    try:
        # Create a new SSH client object
        client_pi = paramiko.SSHClient()

        # Set SSH key parameters to auto accept unknown hosts
        client_pi.load_system_host_keys()
        client_pi.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        # Connect to the host
        client_pi.connect(hostname=HOSTNAME, username=USERNAME, password=PASSWORD)
        print('*** RPi Connection opened ***')

    except Exception:
        traceback.print_exc()
    ######### SSH Terminal for RPi ##########
    
    return client_pi

#### Rasberry Pi close connection
def close_connect_RPi(client_pi):
    try:
        client_pi.close()
        print('*** RPi Connection closed ***')
    except Exception:
            pass

### Take picture with RPi and save to folder
def image_RPi(dirname_ext, document_name, client_pi, i):
    # Create a client interaction class which will interact with the host
    with SSHClientInteraction(client_pi, timeout=60, display=True) as interact_RPi:
        interact_RPi.send('python RPi_Camera.py')
        interact_RPi.expect('End of program', timeout=8)
        print('Picture{} taken'.format(i))

        # Open an STFP file transfer channel
        sftp = client_pi.open_sftp()
        print('SFTP open')
        pipath = '/home/pi/Desktop/image.jpg'
        localpath = dirname_ext+'/'+document_name+'_Cell_{}.jpg'.format(i)
        print('localpath')
        sftp.get(pipath,localpath)
        print('SFTP transferred')

###################################### LABJACK FUNCTIONS ###############################################

#### Labjack T7 open connection
def open_connect_LJ():
    print('Connecting to Labjack ...')
    #openS instead of oen as using string inputs instead of integers
    handle = ljm.openS("T7", "ETHERNET", "ANY") # used to have 470018689 for the last one   # Can replace any of these with "ANY" to connect to any device, connection or identifier
    info = ljm.getHandleInfo(handle)
    print("*** Labjack Connection Opened ***\nA LabJack with Device type: %i, Connection type: %i,\n"
        "Serial number: %i, IP address: %s, Port: %i,\nMax bytes per MB: %i" %
        (info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))
    
    return handle

#### Labjack T7 close connection
def close_connect_LJ(handle):
    try:
        ljm.close(handle)
        print('*** Labjack Connection closed ***')
    except Exception:
            pass

## Labjack T7 scan multiple channels
def LJ_scan(handle, scan_channels, scanRate, measure_time):

    # Channel names
    channel_names = []
    for i in range(len(scan_channels)):
        channel_names.append('AIN{}'.format(scan_channels[i]))

    # Set up scan variables
    numAddresses = len(scan_channels)
    aScanList = ljm.namesToAddresses(numAddresses, channel_names)[0]
    max_requests = 8000
    if (scanRate*measure_time) >= max_requests:
        scansPerRead = max_requests
    else:
        scansPerRead = int(scanRate*measure_time)
    #scansPerRead = int(scanRate/2)

    # Configure and start stream
    ljm.eWriteName(handle, 'STREAM_TRIGGER_INDEX', 0)       # Disable triggered stream
    ljm.eWriteName(handle, 'STREAM_CLOCK_SOURCE', 0)        # Enabling internally-clocked stream
    ljm.eWriteName(handle, "STREAM_SETTLING_US", 0)      # stream settling is 0 (default).
    ljm.eWriteName(handle, "STREAM_RESOLUTION_INDEX", 0) # and stream resolution index is 0 (default).
    ljm.eWriteName(handle, "AIN_ALL_NEGATIVE_CH", ljm.constants.GND) # All negative channels are single-ended.

    for i in range(len(scan_channels)):
        ljm.eWriteName(handle, 'AIN{}_RANGE'.format(scan_channels[i]), 10.0)              # AINx is +/-10 V
    
    # Loop incase Labjack buffer fills
    if (scanRate*measure_time) >= max_requests:
        n_loops = int((scanRate*measure_time)/max_requests + 0.5)
    else:
        n_loops = 1
    i=0
    while i < n_loops:
        scanRate = ljm.eStreamStart(handle, scansPerRead, numAddresses, aScanList, scanRate)
        ret = ljm.eStreamRead(handle)
        aData = ret[0] #The data streamed is the first output of the eStreamRead function
        DataStream = np.array([aData])
        DataArranged = np.reshape(DataStream,(scansPerRead,numAddresses))
        
        #Time between readings is assumed to be given directly by the scan rate
        Time = np.zeros((scansPerRead,1)) 
        
        for j in range(scansPerRead):
            Time[j]=(j/scanRate + i*scansPerRead/scanRate) #seconds
        ljm.eStreamStop(handle)
        LJM_array = np.append(DataArranged,Time,axis=1)

        # Save data from every loop
        if i==0:
            Data = LJM_array
        else:
            Data = np.vstack((Data, LJM_array))

        i+=1
    
    return Data

### Measure temperature in Celsius with thermocouple through Labjack
def LJ_temp(handle, temp_channels):

    ljm.eWriteName(handle, 'AIN{}_NEGATIVE_CH'.format(temp_channels[0]), temp_channels[1])           # AINx negative channel is AINx
    ljm.eWriteName(handle, 'AIN{}_RANGE'.format(temp_channels[0]), 0.01)              # AINx is +/-0.01 V
    ljm.eWriteName(handle, 'AIN{}_SETTLING_US'.format(temp_channels[0]), 0)           # AINx stream settling is 0 (default)
    ljm.eWriteName(handle, 'AIN{}_RESOLUTION_INDEX'.format(temp_channels[0]), 8)      # AINx stream resolution index is 0 (default)
    ljm.eWriteName(handle, 'AIN{}_EF_INDEX'.format(temp_channels[0]), 24)             # Type of thermocouple: T (I think)
    ljm.eWriteName(handle, 'AIN{}_EF_CONFIG_A'.format(temp_channels[0]), 1)           # Temperature unit: C


    # Read multiple temperatures and average
    numFrames = 1
    names = ['AIN{}_EF_READ_A'.format(temp_channels[0])]
    #temp = ljm.eReadNames(handle, numFrames, names)[0]
    loopAmount = 100
    intervalHandle = 1
    ljm.startInterval(intervalHandle, 1000)  # Delay between readings (in microseconds)
    
    temp = []
    for i in range(loopAmount):
            temp.append(ljm.eReadNames(handle, numFrames, names)[0])
            ljm.waitForNextInterval(intervalHandle)
    temp_std = np.std(temp)
    temp = np.average(temp)

    return temp, temp_std


######################################### CALLIBRATION FUNCTIONS ##########################################################

### Accurate calirating based on Labjack input by scanning various points within a square
def prober(dirname_ext, save_data, corner, probe_specs, handle, scan_channels, scanRate, measure_time, interact_COSI, n):
    print('Probing ...')
    ### Probing specifications ###
    names = ['A', 'B', 'C']
    # maximum range around starting point allowed (mm)
    xy_range = probe_specs[0]
    z_range = probe_specs[2]
    # min gap between sampled points (mm)
    xy_step = probe_specs[1]
    z_step = probe_specs[3]
    # Cut off off useless data (0 to 1)
    cut_off = 0.75
    xy_probe_line_number = int(2*xy_range/xy_step + 1)


    ### Find centre of cell in xy position ###
    # Sarting points
    x = corner[0] - xy_range
    y = corner[1] - xy_range
    # Initialise the array
    xy_probes = np.array([x,y])
    # Create list of coordinates to cycle through
    i = 1
    j = 1
    while i < (xy_probe_line_number**2):
        # Move forward
        x += xy_step
        # End of the line for square scanning
        if j >= xy_probe_line_number:
            x = corner[0] - xy_range
            y += xy_step
            j = 0
        # Add coordinates to array
        xy_probes = np.vstack((xy_probes,np.array([x,y])))
        i +=1
        j +=1

    # Scan each coordinate with Labjack
    xy_probe_data = []
    xy_probe_V = []
    for i in range(len(xy_probes)):
        xy_probe_data.append([])
        xy_probe_V.append([])
        # Move to loction
        g1_move('g1x{}g1y{}'.format(str(xy_probes[i][0]), str(xy_probes[i][1])), interact_COSI)
        # Measure strength of beam
        LJM_array = LJ_scan(handle, scan_channels, scanRate, measure_time)
        #Add Labjack data
        for j in range(len(LJM_array)):
            # invert data
            xy_probe_data[i].append(-1*LJM_array[j][0])
    # Convert list to array for easier accessibility
    for i in range(len(xy_probe_data)):
        xy_probe_data[i] = np.array(xy_probe_data[i], dtype=object)
    # Calculate signal strength
    for i in range(len(xy_probe_data)):
        xy_probe_V[i] = np.average(xy_probe_data[i])

    # Turn data into lines and rows to derive in 1D
    # Initialise lines of data
    x_lines = np.zeros(xy_probe_line_number)
    V_lines = np.zeros(xy_probe_line_number)
    y_rows = np.zeros(xy_probe_line_number)
    V_rows = np.zeros(xy_probe_line_number)
    for i in range(xy_probe_line_number-1):
        x_lines = np.vstack((x_lines, np.zeros(xy_probe_line_number)))
        V_lines = np.vstack((V_lines, np.zeros(xy_probe_line_number)))
        y_rows = np.vstack((y_rows, np.zeros(xy_probe_line_number)))
        V_rows = np.vstack((V_rows, np.zeros(xy_probe_line_number)))

    # Save data into lines
    i = 0
    for j in range(xy_probe_line_number):
        for k in range(xy_probe_line_number):
            x_lines[j][k] = xy_probes[i][0]
            V_lines[j][k] = xy_probe_V[i]
            y_rows[j][k] = xy_probes[k*xy_probe_line_number+j][1]
            V_rows[j][k] = xy_probe_V[k*xy_probe_line_number+j]
            i += 1

    # Remove lines and rows with low voltage (=low reflection), and replace all values with 1 or 0
    x_lines_useful = []
    V_lines_useful = []
    y_rows_useful = []
    V_rows_useful = []
    cut_off_V = min(xy_probe_V) + cut_off*(max(xy_probe_V) - min(xy_probe_V))
    for i in range(xy_probe_line_number):
        if(np.max(V_lines[i]) > cut_off_V):
            for j in range(len(V_lines[i])):
                if(V_lines[i][j] > cut_off_V):
                    V_lines[i][j] = 1
                else:
                    V_lines[i][j] = 0
            x_lines_useful.append(x_lines[i])
            V_lines_useful.append(V_lines[i])
        if(np.max(V_rows[i]) > cut_off_V):
            for j in range(len(V_rows[i])):
                if(V_rows[i][j] > cut_off_V):
                    V_rows[i][j] = 1
                else:
                    V_rows[i][j] = 0
            y_rows_useful.append(y_rows[i])
            V_rows_useful.append(V_rows[i])
    x_lines_useful = np.array(x_lines_useful)
    V_lines_useful = np.array(V_lines_useful)
    y_rows_useful = np.array(y_rows_useful)
    V_rows_useful = np.array(V_rows_useful)
    
    # Differentiate
    x_d_lines, d_lines = differentiate(x_lines_useful, V_lines_useful)
    y_d_rows, d_rows = differentiate(y_rows_useful, V_rows_useful)
    # Find mid points between minimum and maximum differentials
    x_mid = []
    for i in range(len(d_lines)):
        xi_min = (np.where(d_lines[i] == min(d_lines[i])))[0][0]
        xi_max = (np.where(d_lines[i] == max(d_lines[i])))[0][0]
        x_mid.append(x_d_lines[i][int(xi_max + (xi_min - xi_max)/2 + 0.5)])
    y_mid = []
    for i in range(len(d_rows)):
        yi_min = (np.where(d_rows[i] == min(d_rows[i])))[0][0]
        yi_max = (np.where(d_rows[i] == max(d_rows[i])))[0][0]
        y_mid.append(y_d_rows[i][int(yi_max + (yi_min - yi_max)/2 + 0.5)])
    probe = [np.average(x_mid), np.average(y_mid)]

    # Plot data on heatmap
    heatmap(dirname_ext, save_data, xy_probes, xy_probe_V, 'Mean xy Voltage Heatmap', '({},{}) centre'.format(probe[0], probe[1]), names[n])

    # Save xy probing data
    if save_data == True:
        x_probes = []
        y_probes = []
        for i in range(len(xy_probes)):
            x_probes.append(xy_probes[i][0])
            y_probes.append(xy_probes[i][1])
        data = [x_probes, y_probes, xy_probe_V]
        save_file(dirname_ext, 'Corner_'+names[n]+'_xy', data, 'X\tY\tV\n')
    

    ### Find ideal z position ###
    # Sarting points
    x = corner[0]
    y = corner[1]
    z = corner[2] - z_range
    # Initialise the array
    z_probes = [z]
    # Create list of coordinates to cycle through
    i = 1
    j = 1
    for i in range(int(2*z_range/z_step + 0.5)):
        # Move forward
        z += z_step
        z_probes.append(z) 
    z_probes = np.array(z_probes)

    #Prepare 3D array to record all probing data and 2D array for mean signal strength
    z_probe_data = []
    z_probe_V = []
    print('\nprobe_V = ', z_probe_V, '\n')
    for i in range(len(z_probes)):
        z_probe_data.append([])
        z_probe_V.append([])

    # Scan each coordinate with Labjack
    for i in range(len(z_probes)):
        # Move to loction
        g1_move('g1x{}g1y{}g1z{}'.format(str(x), str(y), str(z_probes[i])), interact_COSI)
        # Measure strength of beam
        LJM_array = LJ_scan(handle, scan_channels, scanRate, measure_time)
        #Add Labjack data
        for j in range(len(LJM_array)):
            # invert data
            z_probe_data[i].append(-1*LJM_array[j][0])
        
    # Convert list to array for easier accessibility
    for i in range(len(z_probe_data)):
        z_probe_data[i] = np.array(z_probe_data[i], dtype=object)
    print('\nz_probe_data[0] = ', z_probe_data[0], '\n')
    # Calculate signal strength
    for i in range(len(z_probe_data)):
        z_probe_V[i] = np.average(z_probe_data[i])
    print('\nz_probe_V = ', z_probe_V, '\n')
    # Find highest value along z
    zi = (np.where(z_probe_V == max(z_probe_V)))[0][0]
    probe.append(z_probes[zi])

    # Save z probing data
    if save_data == True:
        data = []
        print('len(z_probes) = ', len(z_probes))
        for i in range(len(z_probes)):
            data_line = np.array([z_probes[i], z_probe_V[i]])
            data.append(data_line)
        save_file(dirname_ext, 'Corner_'+names[n]+'Z', data, 'Z\tV\n')

    # Plot signel strength over z
    plt.plot(z_probes, z_probe_V, color='blue', label='z_probing')#, linestyle=':')
    plt.xlabel('Z (mm)')
    plt.ylabel('Voltage (V)')
    plt.legend(loc='best')
    plt.title('Max at z = {}'.format(z_probes[zi]))
    Title = 'Probing Data'
    plt.suptitle(Title)
    if save_data==True:
        plt.savefig(dirname_ext + '/Corner_'+names[n]+'_z_plot_.svg')
    else:
        plt.show()
    plt.close()
    
    print('Probing complete ...')

    # Return centre coordinates
    return probe

# Function to return differential of 2D arrays of data Y with respect to X 
def differentiate(X, Y):
    # Lines
    d_Y = []
    d_X = []
    for i in range(len(Y)):
        d_Y_0 = []
        d_X_0 = []
        for j in range(len(Y[i])-1):
            d_Y_0.append((Y[i][j+1] - Y[i][j])/(X[i][j+1] - X[i][j]))
            d_X_0.append(X[i][j])
        d_Y.append(d_Y_0)
        d_X.append(d_X_0)
    d_Y = np.array(d_Y)
    d_X = np.array(d_X)

    return d_X, d_Y



############################################ Main ##################################################
# Only run if not called from another file
if __name__ == '__main__':
    print('Start of program')

    ###### Parameters ######
    ### Overall
    cell_num = 2 # Number of cells to scan or files to read if fitting with no scanning
    save_data=True # Save the data to file
    plotting=True # Plot data
    resolution = 600 #dpi
    local_path = 'C:\\Users\\lp472\\Desktop\\School\\Uni\\Doctorate\\Lab\\Cell_Scanner\\Scans\\'
    document_name = 'Measurements'
    dirname_descriptor = '2023-07-12_16-03-22_Batch1' # Folder to read in data if there is no scanning
    
    
    ### Scan parameters
    # COSI parameters
    heating=False    #Allow heating of plate under COSI
    RPi=False # Pictures with Rasberry Pi
    calirating=True # If calibrating of the COSI coordinates needed
    # Labjack scan channels
    scan_channels = [0, 1, 2] # Channels AINx (cell data, reference data and calibration data)
    temp_channels = [3]
    ##########################################temp_channels = [3, ljm.constants.GND] # AIN3 to Ground
    # Measurement settings
    measure_time = 0.4    # seconds - Scan time
    scanRate = 10000 # Hz
    set_temp = '80' # degrees Celcius
    # Corner coordinates (x,y,z)
    A = np.array([268.5, 147, 33.5])
    B = np.array([402.5, 146.5, 34])
    C = np.array([269, 278, 32.5])
    corners = [A,B,C]
    RPi_A = np.array([262, 312, 42])
    # calibrating maximum range around starting point allowed (mm)
    xy_range = 8
    z_range = 10
    # calibrating min gap between sampled points (mm)
    xy_step = 0.5
    z_step = 0.5
    probe_specs = [xy_range, xy_step, z_range, z_step]
    
    
    
    ###### Operation ######
    ### Data Location
    if save_data==True:
        print('Create file directory')
        today_datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
        dirname_ext = local_path + today_datetime + '_' + document_name
        ### Create directory for save files
        os.makedirs(dirname_ext)
    else:
        dirname_ext = local_path + dirname_descriptor

        
    ### Scan
    data, coords, temp = scan(dirname_ext, document_name, scan_channels, temp_channels, save_data, heating, calirating, RPi, cell_num, measure_time, scanRate, set_temp, corners, RPi_A, probe_specs) 
    #Print results to terminal
    print('\nData:', data, '\n\n')
    print('\nTemp:', temp, '\n\n')

    print('End of program')
