"""
Created on 8/09/2022
@author: Leigh Page, lp472@sussex.ac.uk

Description:
Program for simple manual operation of the COSI measure, with a few automations added in

Requirement:
Use Paramiko-expect v2.9
"""

import traceback
import paramiko
from paramiko_expect import SSHClientInteraction


######################################### MAIN RUN FUNCTION ############################################
def main():
    print('Start of program')
    
    ### Parameters ###
    auto_homing=True    #Allow autmomatic homing of COSI
    heating=False    #Allow heating of plate under COSIg1
    temp = '80' # degrees Celcius

    ### Operation ###
    
    #Open connections to COSI
    PROMPT_COSI, client_COSI = open_connect_COSI()
    
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

            #Heat the plate under the COSI
            if heating==True:
                HeatPlate(temp, interact_COSI)
                b = TempPlate(interact_COSI)
                print('temp =', b)
            
            #Automatically home COSI
            if auto_homing==True:
                home_xyz(interact_COSI, PROMPT_COSI)

            # Manual Movement
            move_manually(interact_COSI, PROMPT_COSI) #Manual Gcode Commands

        else:
            # If problem with ./mendel.elf - close program
            print('./mendel.elf could not be started!')
            print('Ending Connection')
        
    #Close mendel.elf
    #interact_COSI.send('^C')
    print('***mendel.elf closed ***')

    #Close COSI connection
    close_connect_COSI(client_COSI)

    print('End of program')


###################################### COSI CONNECTION FUNCTIONS ##########################################

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


########################################## COSI MOVEMENT FUNCTIONS ##############################################

### The COSI move function
def g1_move(g1_string, interact_COSI, PROMPT_COSI):
    print('Performing '+g1_string+'...')

    interact_COSI.current_output_clean
    interact_COSI.send(g1_string)
    interact_COSI.expect('\nMovement Complete\n', timeout=60) # Wait for output
    interact_COSI.current_output_clean # Clear output

    print('* Performed '+g1_string+' *')

### Safely homes and sets the COSI arm to the starting position
def home_xyz(interact_COSI, PROMPT_COSI):
    print('homing ...')
    
    # Starting values to move arm safely around COSI structure while homing
    x = 200
    y = 150
    z = 10

    #Home and move z
    g1_move('g161z', interact_COSI, PROMPT_COSI)
    g1_move('g1z{}'.format(z), interact_COSI, PROMPT_COSI)

    #Home and move y
    g1_move('g161y', interact_COSI, PROMPT_COSI)
    g1_move('g1y{}'.format(y), interact_COSI, PROMPT_COSI)
    
    #Home and move x
    g1_move('g161x', interact_COSI, PROMPT_COSI)
    g1_move('g1x{}'.format(x), interact_COSI, PROMPT_COSI)
    print('*** Homing Complete ***')

### Allows manual movement from user input
def move_manually(interact_COSI, PROMPT_COSI):
    command = ''
    while command != 'q':
        print('Enter a manual gcode command (type q to escape):')
        command = input()
        if command != 'q':
            g1_move(command, interact_COSI, PROMPT_COSI)


###################################### COSI TEMPERATURE FUNCTIONS ##########################################

### Heat the plate below the COSI
def HeatPlate(temp, interact_COSI):
    print('Heating plate...')

    interact_COSI.send('m190 s'+temp)
    interact_COSI.expect("temperature for 'temp_bed' has stabilized", timeout=300)
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

    print('*** Temperature Measuring Complete ***')

    return str(b)


############################################ Call run function ########################################

main()
