#! /usr/bin/env python
# -*- coding: utf-8 -*-
#This small script is done by Simona Lombardo
#Please write below modifications and updates

import numpy              as np
from astropy.io import fits
import numpy.ma as ma
from os import listdir
from os.path import isfile,isdir, join
import optparse
import matplotlib.pyplot             as plt
from scipy.optimize import curve_fit
import pickle
from random import randint
from matplotlib.backends.backend_pdf import PdfPages


class CMOSData:


    def __init__(self,path_file):
    
    
        self.path_file = path_file+'/'
        self.Dict_data = { f : fits.open(self.path_file+f) for f in listdir(str(self.path_file))}
        self.key_files = list(self.Dict_data.keys())[:30]
        if self.path_file.split("_")[0][-4:] == 'dark' or self.path_file.split("_")[0][-4:] == 'flat':
            self.exp_time =  float(self.path_file.split("_")[1])#+'.'+self.path_file.split("_")[2])#temporary fix but need to be in header
        else:
            print('Error!!')
                
    def get_data(self):
    
    
        self.data_mat = self.Dict_data.copy()
        self.data_mat = { f : self.Dict_data[f][0].data for f in self.key_files}
       
            
    def get_header(self):
    
    
        self.header_mat = self.Dict_data.copy()
        for k, v in self.Dict_data.items():
            self.header_mat[k] = v[0].header
            
            
    def StatsOnStack(self):
    
        self.file_numb = len(self.data_mat.values())
        self.stak_mat = []
        for v in self.data_mat.values():
            self.stak_mat.append(v)
            
        self.med_mat = np.median(self.stak_mat, axis=0)
        self.med_std = np.std(self.stak_mat, axis=0)
 
      
        
def dark_current(file_path, file_name,linelength):
 
    path_dir = str(file_path)
    
    list_path_file = listdir(str(path_dir))
    exp_time,med_dark,std_dark = [],[],[]
    for i in list_path_file:
        print(i)
        data = CMOSData(path_dir+i)
        data.get_data()
        data.StatsOnStack()
        exp_time.append(data.exp_time*linelength*50*1e-9)
        #print('exp_time',exp_time)
        #plt.imshow(data.med_mat)
        #plt.show()
        med_dark.append(np.median(data.med_mat))
        std_dark.append(np.std(data.med_mat))
        
    index=np.argsort(exp_time)
    exp_time = np.array(exp_time)[index]
    med_dark = np.array(med_dark)[index]
    std_dark = np.array(std_dark)[index]
        
    func = lambda x,m,c: x *m  + c
    popt, pcov = curve_fit(func, exp_time, med_dark)
    np.savetxt('dark_%s.txt' %file_name, np.asarray((exp_time,med_dark)))
        
    return exp_time, med_dark, std_dark, popt,pcov
    
    
    
def compute_gain(file_path, file_name,linelength):
 
    path_dir = str(file_path)
    
    list_path_file = listdir(str(path_dir))
    exp_time,med_flat,std_flat = [],[],[]
    
    for i in list_path_file:
        print(i)
        data = CMOSData(path_dir+i)
        data.get_data()
        data.StatsOnStack()
        exp_time.append(data.exp_time*linelength*50*1e-9)
        med_flat.append(np.median(data.med_mat))
        std_flat.append(np.median(data.med_std))
        
    index=np.argsort(exp_time)
    exp_time = np.array(exp_time)[index]
    med_flat = np.array(med_flat)[index]
    std_flat = np.array(std_flat)[index]
        
    func = lambda x,m,c: x *m  + c
    popt, pcov = curve_fit(func, exp_time, med_flat)
    np.savetxt('flatfield_%s.txt' %file_name, np.asarray((exp_time,med_flat)))
        
    return exp_time, med_flat, std_flat, popt,pcov
    
    
    
        
        
