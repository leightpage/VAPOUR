# VAPOUR
Automated spectroscopic characterisation of alkali vapour cells for magnetometry


This code was written as part of a project to automate the characterisation of alkali vapour cells used in such applications as optically pumped magnetomers (OPMs).

Automation was the key goal of this project, with two Python programs created to control everything, Cell_Scanner.py and Minimise_Fit.py, and third program to call upon both, Cell_Scanner_Fit.py. These programs were designed to be modular, with different features that can be turned on or off in the code depending on the needs. The RPi_Camera.py function was uploaded to a Raspberry Pi operating a microscope camera, for the Cell_Scanner.py program to communicate with. The Manual_move.py function was used when wanting to manually move the 3-axis robot arm of the COSI Measure during experimental setup.


The Scanning program, Cell_Scanner.py, directs the COSI Measure robot arm to move sequentially over the rubidium vapour cells while recording the outputs of the 3 photodiodes and thermocouple via the DAQ, and also capturing pictures using a single board computer (Raspberry Pi 4, Raspberry Pi Ltd, UK) operated camera.

The Fitting program, Minimise_Fit.py, fits a spectroscopic linewidth model to the data from the Scanning program to calculate the pressure of the rubidium vapour cells, using the minimise chi square method.
