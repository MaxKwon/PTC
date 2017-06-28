#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 11:29:07 2017

@author: Max K. Kwon 
"""

import numpy as np 
import matplotlib.pyplot as plt # plotting package
import matplotlib.animation as animation
from astropy.io import fits
from tkinter.filedialog import askopenfilename
from scipy import fftpack
import pylab as py

file = askopenfilename(filetypes=[('.fit','.fits', '.fts')],title='Pick bias 1')
image_data = fits.getdata(file)
header = fits.getheader(file)

print('Number of Frames', len(image_data)) # shows the number of frames in the data cube (1, 2 or many )
print('Type: ', type(image_data))
print('Shape', image_data.shape) #dimensions of array
print(repr(header))

plt.subplot(1, 6, 1*1)
plt.title("Initial Image")
plt.imshow(image_data[0], origin='lower') #showing the first frame of the cube and shifting the origin to the bottom left side of the plot 
 
          
med = np.median(image_data, axis=0) # will take the median of the entire cube with respect to the depth axis and display it as a frame (median frame)

plt.subplot(1, 6, 1*2)
plt.title("Median of All")
plt.imshow(med, origin='lower')


num_frames = len(image_data) #find the length of the first dimension of the 3D array which is the number of frames
x_dim = len(image_data[0]) #set the dimension of the x axis, returns the length of that portion of the array 
y_dim = len(image_data[0][0]) #set the dimension of the y axis, returns the length of that portion of the array             


diff_image = image_data[num_frames-1]-image_data[0] #takes the difference between the last and the first images in the 3D Array

plt.subplot(1, 6, 1*3)
plt.title("Difference")
plt.imshow(diff_image, origin='lower')

std_devs = np.zeros(shape=(y_dim, x_dim))
stack = np.zeros(shape=(num_frames))
#this for loop is really slow and somehow needs to become more efficient 
for i in range(0, y_dim):
    for j in range(0, x_dim):
       for k in range(0, num_frames):
            stack[k] = image_data[k][i][j] #makes a 1D array of the depth of the cube at a certain point on the 2D array 
            if k==(num_frames-1):
               std_devs[i][j] = np.std(stack) #after each 1D array is filled, the SD is taken (once at each pixel going through each frame)

plt.subplot(1, 6, 1*4)
plt.title("Standard Dev")
plt.imshow(std_devs, origin='lower')

num_pixels = x_dim * y_dim
sd_pixel = np.zeros(shape=(num_pixels))
pixel_index = range(num_pixels) #Ascending array from 0 to number of pixels - 1

index = 0 #keeps track of the pixel number in the for loop
#puts the pixel SD data into a 1d array in order to show it as a line graph where the x axis is the pixel number starting from the top left to tthe bottom right
for i in range(0, len(std_devs)): #goes through each row
    for j in range(0, len(std_devs[0])):#goes through each column
        sd_pixel[index] = std_devs[i][j]
        index = index + 1
        
plt.subplot(1, 6, 1*5)
plt.title("SD Per Pixel")
plt.plot(pixel_index, sd_pixel) 

print('Median SD: ', np.median(std_devs)) #the median value of the 2D array full of the SDs for each pixel 
     

def plotSDTime(data): #will plot the median standard deviation over time to see if the noise in the CCD increases with time and exposures done 
    
    mean_data = np.mean(data)          
              
    med_std_dev_frame = np.zeros(shape = (num_frames))
        
    for i in range(0, num_frames):
        med_std_dev_frame[i] = np.median((np.mean((data[i] - mean_data)**2))**(1/2.0)) # SD = (mean((value-mean_of_data)^2))^.5
                     
    frame_index = range(num_frames) #creates an array from 0 to num_frames - 1 to be used as the X axis of the plot 
                       
    print(med_std_dev_frame)
                       
    plt.subplot(1, 6, 1*6)
    plt.title("SD Over Time")
    plt.plot(frame_index, med_std_dev_frame) 
    
if (num_frames > 2): #only want SD over time if it is a data cube with more than 2 frames
    plotSDTime(image_data) 
    


