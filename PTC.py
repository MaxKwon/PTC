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
from tkinter import *
from scipy import fftpack, misc
import pylab as py

#choose which data you want (0 = no print, 1 = print)
#Params:
#info, initial, median, diff, SD pixel 2D, Sd pixel 1D, print median, SD per frame, plot LC, fourier transform
selection_array = [1,1,1,1,1,1,1,1,1,1] # 0 = no, 1 = add

#user input the number of fits files to be analyzed 
num_files = 3                
                  
files = ["" for x in range(num_files)]  #creates an array for the number of fits file

for i in range(0, num_files): #fills an array with all file names that are desired 

    files[i] = askopenfilename(filetypes=[('.fit','.fits', '.fts')], title='Pick bias')                      
    

image_data = fits.getdata(files[0]) #the image data for the first file that is input, starts with the initial file 

print("image data size: ", image_data.shape)
                         
num_frames = len(image_data) #find the length of the first dimension of the 3D array which is the number of frames
x_dim = len(image_data[0]) #set the dimension of the x axis, returns the length of that portion of the array 
y_dim = len(image_data[0][0]) #set the dimension of the y axis, returns the length of that portion of the array     
num_pixels = x_dim * y_dim #number of total pixels in the image (height * width)                  
                  
subplot_index = 1; # the position of a graph in the subplot (1 = most left), incremented everytime a plot is plotted in the subplot                                                           
subplot_rows = num_files # number of rows in the subplot
subplot_columns = 8 # number of columns in subplot 

#############################################
linear_fit_start = 0 #inclusive starts at 0
linear_fit_stop = 5 #this value is exclusive     

#Protects the code if the upper bound is longer than the array
if (linear_fit_stop > len(files)):
    linear_fit_stop = len(files)  

####some stuff to make the plots bigger####
plt.close('all')
fig = plt.figure(figsize=(12, 8.5))

###########################################                                                 
           
def printInfo(file, data):
    
    header = fits.getheader(file)
    
    print('Number of Frames', len(data)) # shows the number of frames in the data cube (1, 2 or many )
    print('Type: ', type(data))
    print('Shape', data.shape) #dimensions of array
    print(repr(header))
    
    #print("Process Complete printInfo")

def plotInitialImage(data, subplot_pos):
    
    plt.subplot(subplot_rows, subplot_columns, subplot_pos)
    plt.title("Init")
    plt.imshow(data[0], origin='lower') #showing the first frame of the cube and shifting the origin to the bottom left side of the plot 
    
    #print("Process Complete plotInitialImage")
   
def plotMedian(data, subplot_pos):     
    
    med = np.median(data, axis=0) # will take the median of the entire cube with respect to the depth axis and display it as a frame (median frame)
    
    plt.subplot(subplot_rows, subplot_columns, subplot_pos)
    plt.title("Med")
    plt.imshow(med, origin='lower')
    
    #print("Process Complete plotMedian")
    
    
def plotDiffImage(data, subplot_pos):
    
    diff_image = data[num_frames-1]-data[0] #takes the difference between the last and the first images in the 3D Array
    
    plt.subplot(subplot_rows, subplot_columns, subplot_pos)
    plt.title("Diff")
    plt.imshow(diff_image, origin='lower')
    
    #print("Process Complete plotDiffImage")


def plotSDPerPixel2D(data, subplot_pos):
    
    std_devs = np.zeros(shape=(y_dim, x_dim))

    std_devs = data.std(axis=0)

    plt.subplot(subplot_rows, subplot_columns, subplot_pos)
    plt.title("SD")
    plt.imshow(std_devs, origin='lower')
    
    #print("Process Complete plotSDPerPixel2D")
    
    return std_devs


def plotSDPerPixel1D(std_arr, subplot_pos):
    
    sd_pixel = np.zeros(shape=(num_pixels))
    pixel_index = range(num_pixels) #Ascending array from 0 to number of pixels - 1
    
    sd_pixel = std_arr.flatten()
            
    plt.subplot(subplot_rows, subplot_columns, subplot_pos)
    plt.title("SD Pix")
    plt.plot(pixel_index, sd_pixel) 
    
    #print("Process Complete plotSDPerPixel1D")


def printMedianSD(std_arr):
    print('Median SD: ', np.median(std_arr)) #the median value of the 2D array full of the SDs for each pixel
    
   #print("Process Complete printMedianSD")
     
def plotSDTime(data, subplot_pos): #will plot the median standard deviation over time to see if the noise in the CCD increases with time and exposures done 
    
    mean_data = np.mean(data)
    med_std_dev_frame = np.zeros(shape = (num_frames))
        
    for i in range(0, num_frames):
        med_std_dev_frame[i] = np.median((np.mean((data[i] - mean_data)**2))**(1/2.0)) # SD = (mean((value-mean_of_data)^2))^.5
                     
    frame_index = range(num_frames) #creates an array from 0 to num_frames - 1 to be used as the X axis of the plot 
   
    plt.subplot(subplot_rows, subplot_columns, subplot_pos)
    plt.title("SD OT")
    plt.plot(frame_index, med_std_dev_frame)   
    
   #print("Process Complete plotSDTime")
    
#Take many data cubes and plot the median counts of each cube by the median SD of the cube 
def plotPTC(files, subplot_pos):
    
    med_counts = np.zeros(shape=(num_files))
    med_std_devs = np.zeros(shape=(num_files))
    
    for i in range (0, len(files)):
        
        data = fits.getdata(files[i])
        med_counts[i] = np.median(np.median(data))
        med_std_devs[i] = np.median(data.std(axis=0))
        
    
    med_std_devs_sorted = [x for y, x in sorted(zip(med_counts, med_std_devs))] #sorts the standard deviation array in the same order as the median count array
    med_counts.sort()
    
    linear_fit_counts = np.zeros(shape=(linear_fit_stop - linear_fit_start))
    linear_fit_std_devs = np.zeros(shape=(linear_fit_stop - linear_fit_start))
    
    index = 0
    
    for i in range(linear_fit_start, linear_fit_stop):
        
        linear_fit_counts[index] = med_counts[i]
        linear_fit_std_devs[index] = med_std_devs_sorted[i]

        index = index + 1

    plt.subplot(subplot_rows, subplot_columns, subplot_pos)
    plt.title("PTC")
    plt.plot(med_counts, med_std_devs_sorted)
    plt.plot(np.unique(linear_fit_counts), np.poly1d(np.polyfit(linear_fit_counts, linear_fit_std_devs, 1))(np.unique(linear_fit_counts))) #plots linear regression model of the data
    
    slope, y_intercept = np.polyfit(linear_fit_counts, linear_fit_std_devs, 1) #returns the coeeficients of the polynomial that fits the curve, ours is of order 1 so only slope and y intercept (mx+b)
                                         
    print('Gain', slope)
    print('Y Intercept', y_intercept)
    print('counts', med_counts)
    print('std', med_std_devs_sorted)
    
   #print("Process Complete plotPTC")
    
def plotFourierTransform(data, subplot_pos):

    plt.subplot(subplot_rows, subplot_columns, subplot_pos)  # show the power spectrum
    plt.title("PS")
   
    file_F1 = fftpack.fft2(data)  # Take the fourier transform of the image.

    file_F2 = fftpack.fftshift(file_F1) # Now shift the quadrants around so that low spatial frequencies are in the center of the 2D fourier transformed image.
    
    file_psd2D = np.abs(file_F2) ** 2# Calculate a 2D power spectrum
    py.imshow(np.log10(file_psd2D), origin='lower')
    
    plt.savefig('PTC.pdf') #saves the subplot into a pdf in the folder where the program is 
    plt.savefig('PTC.jpg') #saves the subplot into a jpg in the folder where the program is 
    
   #print("Process Complete plotFourierTransform")

#checks if an array uses numbers other than 0 or 1
def notGoodNumber(arr): 
    
    for i in range(0, len(arr)):
        if ( arr[i] != 0 or arr[i] != 1):    
            return True
    return False


###########################################################

if (notGoodNumber(selection_array)): #sends out a warning (error) if a number other than 0 or 1 is used
    print("ONLY USE 1 OR 0 FOR SELECTION ARRAY")


if (selection_array[0] == 1):    
    printInfo(files[0], image_data)   

### Loop plots for all files ###

row_index = 0  #adds 8 every time to switch rows (length of one row)

for i in range (0, num_files):
    
    image_data = fits.getdata(files[i])
    
    subplot_index = 1
    
    if (selection_array[1] == 1):    
        plotInitialImage(image_data, subplot_index + row_index)    
        subplot_index = subplot_index + 1

    if (selection_array[2] == 1):
        plotMedian(image_data, subplot_index + row_index)
        subplot_index = subplot_index + 1
    
    if (selection_array[3] == 1):
        plotDiffImage(image_data, subplot_index + row_index)
        subplot_index = subplot_index + 1

    if (selection_array[4] == 1):
        sd_plot = plotSDPerPixel2D(image_data, subplot_index + row_index)
        subplot_index = subplot_index + 1

    if (selection_array[5] == 1):
        if (selection_array[4] == 1):
            plotSDPerPixel1D(sd_plot, subplot_index + row_index)
            subplot_index = subplot_index + 1
    
        else: # error given if the 2D array isnt used and the 1D standard deviation is attempted tobe used
            print("The 2D standard Deviation must be aquired in order to view the 1D Standard Deviation")

    if (num_frames > 2 and selection_array[7] == 1): #only want SD over time if it is a data cube with more than 2 frames
       plotSDTime(image_data, subplot_index + row_index)
       subplot_index = subplot_index + 1

    if (selection_array[9] == 1):
        plotFourierTransform(np.require(image_data[0], dtype=np.float32), subplot_index + row_index) #image_data[0] needs to be converted from an endian type, the encapsulating method does that 
        subplot_index = subplot_index + 1
    
    row_index = row_index + 8    
   
        
#### END FOR LOOP #####

if (selection_array[6] == 1):
    if (selection_array[4] == 1):
        printMedianSD(sd_plot)
        
    else: # error given if the 2D array isnt used and the median is attempted to be used
        print("The 2D standard Deviation must be aquired in order to view the Median Standard Deviation")

if (selection_array[8] == 1):
    plotPTC(files, 8)
    subplot_index = subplot_index + 1
    
    
fig.tight_layout()
plt.show()
