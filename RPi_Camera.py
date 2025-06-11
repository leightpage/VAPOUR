"""
Program to run on Raspberry Pi, take and save a picture while powering an LED for lighting

Created on 30/05/2023
@author: Leigh Page, L.Page@sussex.ac.uk

raspberry-pi-camera-guide.pdf
https://raspberrytips.com/install-camera-raspberry-pi/
https://codingandfun.com/python-scripts-in-raspberry-pi/
https://thepihut.com/blogs/raspberry-pi-tutorials/27968772-turning-on-an-led-with-your-raspberry-pis-gpio-pins?
"""
import picamera
import RPi.GPIO as GPIO
import time
import sys

def run():
    print('Start of program')
    LED_on()
    rpi_camera_take_image()
    LED_off()
    print('End of program')

def rpi_camera_take_image():
    print('Raspberry pi camera take image')
    with picamera.PiCamera() as camera:
        camera.start_preview()
        time.sleep(1)
        camera.capture('/home/pi/Desktop/image.jpg')
        camera.stop_preview()
    
def LED_on():
    GPIO.setmode(GPIO.BCM)
    GPIO.setwarnings(False)
    GPIO.setup(18,GPIO.OUT)
    print("LED on")
    GPIO.output(18,GPIO.HIGH)

def LED_off():
    print("LED off")
    GPIO.output(18,GPIO.LOW)

run()
