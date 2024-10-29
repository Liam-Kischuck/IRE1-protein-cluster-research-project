# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 16:23:51 2022

@author: Michael
"""
import numpy as np
import matplotlib.pyplot as plt
from Read_MATLAB_file import read_matlab_data
import os

analytic_threshold = [40,70,100,140]     #The threshold particle number around which the wrapping is tested.
scales = [10,10,10,10] #Each of the Density folders increase particle number by this number.
min_radius_num = 1     #The radius number 'm' must start from this number.

lattice_site_num = [1200,1500,1800]

empty_instances_num = 0

plt.rcParams.update({'font.size': 16})
cmap = plt.get_cmap('viridis',4)
line_styles = ['solid', 'dotted', 'dashed']
labels = [r'$d_{tube} = 3.8 _ d_{IRE1}$',r'$d_{tube} = 4.8 _  d_{IRE1}$',r'$d_{tube} = 5.7 _ d_{IRE1}$']
date = 'Aug31'
for m in range(1,4):
    
    meta_folder = 'Transitions_1/Transitions_Data/Radius_' + str(m) + '/'
    Wrap_list = []
    Circular_list =[]
    standard_error_wrap = []
    standard_error_circle = []
    k_list = range(1,13)
    if m == 3:
        k_list = range(1,11)
    
    
    for k in k_list: #[1,2,3,4,5...12]
        wrap_rec = []
        circular_rec = []
        for i in range(1,5):
            instance = str(i)
            
            print(instance)
            if (k != 11) & (k != 12):
                folder = os.listdir( meta_folder + "Density_" + str(k) + "/" + date + "_subfolder_" + instance)
            elif (k == 11) :
                folder = os.listdir(meta_folder + "Density_" + '3_2' + "/" + "Jan24_" + instance + "/Jan24" + "_subfolder_" + instance)
            elif (k == 12):
                folder = os.listdir(meta_folder + "Density_" + '3_3' + "/" + "Jan24_" + instance + "/Jan24" + "_subfolder_" + instance)
            
            
            
            
            if len(folder) == 0:  #check if folder is empty
                empty_instances_num += 1
            
            
            else:
                folderfileT = meta_folder + "Density_" + str(k) + "/" +  date + "_subfolder_" + instance + "/" + date + "_time_data_" + instance + ".txt"
                folderfileWRAP = meta_folder + "Density_" + str(k) + "/" + date + "_subfolder_" + instance + "/" + date + "_wrap_data_" + instance + ".txt"
            
                if (k == 11) & (m != 3):   #Tag the new data to the end of the list
                    folderfileT = meta_folder + "Density_" + '3_2' + "/" + "Jan24_" + instance + "/Jan24" + "_subfolder_" + instance + "/" + "Jan24" + "_time_data_" + instance + ".txt"
                    folderfileWRAP = meta_folder + "Density_" + '3_2' + "/" + "Jan24_" + instance + "/Jan24" + "_subfolder_" + instance + "/" + "Jan24" + "_wrap_data_" + instance + ".txt"
                
                if (k == 12) & (m != 3):   #Tag the new data to the end of the list
                    folderfileT = meta_folder + "Density_" + '3_3'  + "/" +  "Jan24_" + instance + "/Jan24" + "_subfolder_" + instance + "/" + "Jan24" + "_time_data_" + instance + ".txt"
                    folderfileWRAP = meta_folder + "Density_" + '3_3' + "/" + "Jan24_" + instance + "/Jan24" + "_subfolder_" + instance + "/" + "Jan24" + "_wrap_data_" + instance + ".txt"
                
                
            
                y_file = open(folderfileWRAP, "r")
                x_file = open(folderfileT, "r")
            
                y_list = read_matlab_data(y_file)
                x_list = read_matlab_data(x_file)
                x_title = 'Time (s)'
                
                for j_index in range(0,6):
            
                    w_data = y_list[j_index]
                    t_data = x_list[j_index]
                    
                    counting_wrap_time = 0
                    wrapped = False
                    if w_data[0] == 1:
                        start_wrapping_at = t_data[0]
                        wrapped = True  #In case it is already wrapped.
                        
                    for n in range(len(w_data) - 1):
                        if (w_data[n] == 1) &  (wrapped == False):  #The cluster becomes a wrapped cluster. Record the thius time
                            start_wrapping_at = t_data[n] 
                            wrapped = True
                        elif (w_data[n] == 0) &  (wrapped == True):  #Record when the cluster goes back to being a circle.
                            end_wrap_at = t_data[n - 1]              #Take the point before the ciruclar cluster as the last record of a wrapped cluster.
                            wrapped =False
                            counting_wrap_time += (end_wrap_at - start_wrapping_at)
                            
                    if w_data[-1] == 1: #If you don't evaluate this, this wrapped time never gets added to the counting_wrap_time.
                        end_wrap_at = t_data[n - 1]              
                        counting_wrap_time += (end_wrap_at - start_wrapping_at)
                    
                    wrap_rec.append(counting_wrap_time/t_data[-1])
                    circular_rec.append((t_data[-1] - counting_wrap_time)/t_data[-1])
                        
                 
        Wrap_list.append(np.average(wrap_rec))
        Circular_list.append(np.average(circular_rec))
        standard_error_wrap.append(np.std(wrap_rec) / np.sqrt(np.size(wrap_rec)))
        standard_error_circle.append(np.std(circular_rec) / np.sqrt(np.size(circular_rec)))
        #/lattice_site_num[m-1]
    x_range = [(analytic_threshold[m - min_radius_num] + scales[m - min_radius_num]*h)/lattice_site_num[m-1] for h in range(-1,9) ]
    
    if m == 1:
        x_range = x_range[0:3] + [53/lattice_site_num[m-1], 56/lattice_site_num[m-1]] + x_range[3:] #For the tube with Li = 12, these were the new particle density at which we got more data.
        Wrap_list = Wrap_list[0:3] +  Wrap_list[10:12]  + Wrap_list[3:10]
        
    if m == 2:
        x_range = x_range[0:3] + [83/lattice_site_num[m-1], 86/lattice_site_num[m-1]] + x_range[3:] #For the tube with Li = 12, these were the new particle density at which we got more data.
        Wrap_list = Wrap_list[0:3] +  Wrap_list[10:12]  + Wrap_list[3:10]
        
    
    plt.plot(x_range, Wrap_list, color = cmap(m),linestyle = line_styles[m-1], label = labels[m-1], linewidth =2)
    plt.scatter(x_range, Wrap_list, color = cmap(m),s  = 30)
    # plt.scatter(x_range, Circular_list, color = 'blue', s = 20)
    # plt.plot(x_range, Wrap_list, label ='wrapped', color = 'red')
    # plt.plot(x_range, Circular_list, label = 'circle', color = 'blue')
    
    plt.errorbar(x_range, Wrap_list, yerr = standard_error_wrap, color = cmap(m), capsize= 4, elinewidth= 2.5, linestyle = '')
    #plt.errorbar(x_range, Circular_list, yerr = standard_error_circle, color = 'blue', capsize= 3)
plt.xlabel(r"Concentration $(\times 10^{-2} /d_{IRE1}^2)$ ", fontsize = 18)
plt.xticks([0.04,0.06,0.08,0.10], [4.0,6.0,8.0,10.0])


plt.ylabel("Wrapping probability", fontsize = 18)
plt.legend(fontsize = 14, frameon=False)

plt.savefig('shape_probability_Scaled.svg')
plt.show()


print('empty folder number', empty_instances_num)