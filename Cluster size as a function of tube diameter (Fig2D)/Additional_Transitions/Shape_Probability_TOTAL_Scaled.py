# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 16:23:51 2022

@author: Michael
"""
import numpy as np
import matplotlib.pyplot as plt
from Read_MATLAB_file import read_matlab_data
import os

#The 2 lists below are for the final graph
Radius = []
Transition_Sizes = []
Std_Error = []



#For the first radius  with lots of samples
analytic_threshold = [20]     #The threshold particle number around which the wrapping is tested.
scales = [2] #Each of the Density folders increase particle number by this number.
min_radius_num = 1     #The radius number 'm' must start from this number.
lattice_site_num = [1000]

date = 'May27'
for m in range(1,2):
    #This is from MassiveMatlab2/Transitions_2/Radius_4' data
    #meta_folder = 'Radius_' + str(m) + '/'
    Wrap_list = []
    Wrap_list_stdev = []
    Circular_list =[]
    for k in range(1,41):
        for i in range(1,2):
            instance = str(k)  #This line is different from the other cases, becuasse
            #The May 27 data was organized differentlys
            wrap_rec = []
            circular_rec = []
            print(instance)
            #folderfileT = meta_folder + "Density_" + str(k) + "/" +  date + "_subfolder_" + instance + "/" + date + "_time_data_" + instance + ".txt"
            #folderfileWRAP = meta_folder + "Density_" + str(k) + "/" + date + "_subfolder_" + instance + "/" + date + "_wrap_data_" + instance + ".txt"
            folderfileT = date + "_ft_data/" + date + "_subfolder_" + instance + "/" + date + "_time_data_" + instance + ".txt"
            folderfileWRAP = date + "_ft_data/" + date + "_subfolder_" + instance + "/" + date + "_wrap_data_" + instance + ".txt"

        
            y_file = open(folderfileWRAP, "r")
            x_file = open(folderfileT, "r")
        
            y_list = read_matlab_data(y_file)
            x_list = read_matlab_data(x_file)
            x_title = 'Time (s)'
            
            for j_index in range(0,25):
        
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
        Wrap_list_stdev.append(np.std(wrap_rec))
        Circular_list.append(np.average(circular_rec))
    
    x_range = [(analytic_threshold[m - min_radius_num] + scales[m - min_radius_num]*h)/lattice_site_num[m-min_radius_num] for h in range(0,40) ]
    
    
    
    #The section for finding 0.5 point
    diff_from_half = np.array(Wrap_list) - 0.5
    negatives = [ow for ow in diff_from_half if ow < 0]
    before_half_index = len(negatives) - 1
    after_half_index = before_half_index + 1
    
    #For the transition graph-----
    delta_particles = scales[m - min_radius_num]
    m_slope = (Wrap_list[after_half_index] - Wrap_list[before_half_index])/delta_particles
    n1 = analytic_threshold[m - min_radius_num] + before_half_index*scales[m - min_radius_num]
    n2 = analytic_threshold[m - min_radius_num] + after_half_index*scales[m - min_radius_num]
    
    n = [n1, n2]
    w = [Wrap_list[before_half_index], Wrap_list[after_half_index]]
    sig_w = [Wrap_list_stdev[before_half_index], Wrap_list_stdev[after_half_index]]
    
    #a, b = np.polyfit(n,w,1)
    
    a = (w[1] - w[0])/(n2 - n1)
    b = w[1] - a*n2
    
    print('slope comparison ', m_slope, a)
    
    transition_cluster_size =  (0.5 - b)/a
    
    Transition_Sizes.append(transition_cluster_size)
    Radius.append(m)
    #------------------------
    #Standard Error------------
    sig_a = (1/(n2-n1))*np.sqrt(sig_w[0]**2 + sig_w[1]**2)
    sig_b= np.sqrt(sig_w[1]**2 + (sig_a*n2)**2)
    sig_Transition_Sizes = np.sqrt((sig_b/b)**2 + (sig_a/a)**2)
    Std_Error.append(sig_Transition_Sizes)
    #--------------------------
    
    #plt.scatter(x_range, Wrap_list, color = 'red', s = 30) For the purpose of the graph, commented out.
    #plt.scatter(x_range, Wrap_list)
 #   plt.errorbar(x_range, Wrap_list, yerr = standard_error_wrap, color = 'red', capsize= 3, linewidth = 2.5, elinewidth= 2.5)    
    
#=============================================================


analytic_threshold = [40,70,100,140]     #The threshold particle number around which the wrapping is tested.
scales = [10,10,10,10] #Each of the Density folders increase particle number by this number.
min_radius_num = 1     #The radius number 'm' must start from this number.

lattice_site_num = [1200,1500,1800]

empty_instances_num = 0

plt.rcParams.update({'font.size': 16})



date = 'Aug31'
for m in range(1,4):
    
    meta_folder = 'Transitions_Data/Radius_' + str(m) + '/'
    Wrap_list = []
    Wrap_list_stdev = []
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
        Wrap_list_stdev.append(np.std(wrap_rec))
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
       
    #For the transition graph-----
    
    for pi in range(len(Wrap_list)):
        if Wrap_list[pi] < 0.5 and max(Wrap_list[:(pi + 1)]): #This is the element before the 0.5 proability.
            before_half_index = pi
        elif Wrap_list[pi] > 0.5 and min(Wrap_list[pi:]):
            after_half_index = pi
            
    diff_from_half = np.array(Wrap_list) - 0.5
    negatives = [ow for ow in diff_from_half if ow < 0]
    before_half_index = len(negatives) - 1
    after_half_index = before_half_index + 1
    
    
    delta_particles = scales[m - min_radius_num]
    m_slope = (Wrap_list[after_half_index] - Wrap_list[before_half_index])/delta_particles
    n1 = analytic_threshold[m - min_radius_num] + (before_half_index - 1)*scales[m - min_radius_num] #-1 because these data starts 1 scale size below the analytic threshold
    n2 = analytic_threshold[m - min_radius_num] + (after_half_index - 1)*scales[m - min_radius_num]
    
    n = [n1, n2]
    w = [Wrap_list[before_half_index], Wrap_list[after_half_index]]
    sig_w = [Wrap_list_stdev[before_half_index], Wrap_list_stdev[after_half_index]]
    
    #a, b = np.polyfit(n,w,1)
    
    a = (w[1] - w[0])/(n2 - n1)
    b = w[1] - a*n2
    
    print('slope comparison ', m_slope, a)
    
    transition_cluster_size =  (0.5 - b)/a
    
    Transition_Sizes.append(transition_cluster_size)
    Radius.append(m + 1)
    #------------------------
    #Standard Error------------
    sig_a = (1/(n2-n1))*np.sqrt(sig_w[0]**2 + sig_w[1]**2)
    sig_b= np.sqrt(sig_w[1]**2 + (sig_a*n2)**2)
    sig_Transition_Sizes = np.sqrt((sig_b/b)**2 + (sig_a/a)**2)
    Std_Error.append(sig_Transition_Sizes)
    #--------------------------
    
    plt.scatter(x_range, Wrap_list, color = 'red', s = 30)
   # plt.scatter(x_range, Circular_list, color = 'blue', s = 20)
    # plt.plot(x_range, Wrap_list, label ='wrapped', color = 'red')
    # plt.plot(x_range, Circular_list, label = 'circle', color = 'blue')
    
    plt.errorbar(x_range, Wrap_list, yerr = standard_error_wrap, color = 'red', capsize= 5, linewidth = 2.5, elinewidth= 2.5)
    #plt.errorbar(x_range, Circular_list, yerr = standard_error_circle, color = 'blue', capsize= 3)
    





analytic_threshold = [140]     #The threshold particle number around which the wrapping is tested.
scales = [15] #Each of the Density folders increase particle number by this number.
min_radius_num = 5     #The radius number 'm' must start from this number.
lattice_site_num = [2100]

date = 'Sep4'
for m in range(5,6):
    #This is from MassiveMatlab2/Transitions_2/Radius_4' data
    meta_folder = 'Radius_' + str(m) + '/'
    Wrap_list = []
    Wrap_list_stdev = []
    Circular_list =[]
    for k in range(1,11):
        for i in range(1,13):
            instance = str(i)
            wrap_rec = []
            circular_rec = []
            print(instance)
            folderfileT = meta_folder + "Density_" + str(k) + "/" +  date + "_subfolder_" + instance + "/" + date + "_time_data_" + instance + ".txt"
            folderfileWRAP = meta_folder + "Density_" + str(k) + "/" + date + "_subfolder_" + instance + "/" + date + "_wrap_data_" + instance + ".txt"
        
            y_file = open(folderfileWRAP, "r")
            x_file = open(folderfileT, "r")
        
            y_list = read_matlab_data(y_file)
            x_list = read_matlab_data(x_file)
            x_title = 'Time (s)'
            
            for j_index in range(0,2):
        
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
        Wrap_list_stdev.append(np.std(wrap_rec))
        Circular_list.append(np.average(circular_rec))
    
    x_range = [(analytic_threshold[m - min_radius_num] + scales[m - min_radius_num]*h)/lattice_site_num[m-min_radius_num] for h in range(-1,9) ]
    
    
    
    #The section for finding 0.5 point
    diff_from_half = np.array(Wrap_list) - 0.5
    negatives = [ow for ow in diff_from_half if ow < 0]
    before_half_index = len(negatives) - 1
    after_half_index = before_half_index + 1
    
    #For the transition graph-----
    delta_particles = scales[m - min_radius_num]
    m_slope = (Wrap_list[after_half_index] - Wrap_list[before_half_index])/delta_particles
    n1 = analytic_threshold[m - min_radius_num] + before_half_index*scales[m - min_radius_num]
    n2 = analytic_threshold[m - min_radius_num] + after_half_index*scales[m - min_radius_num]
    
    n = [n1, n2]
    w = [Wrap_list[before_half_index], Wrap_list[after_half_index]]
    sig_w = [Wrap_list_stdev[before_half_index], Wrap_list_stdev[after_half_index]]
    
    #a, b = np.polyfit(n,w,1)
    
    a = (w[1] - w[0])/(n2 - n1)
    b = w[1] - a*n2
    
    print('slope comparison ', m_slope, a)
    
    transition_cluster_size =  (0.5 - b)/a
    
    Transition_Sizes.append(transition_cluster_size)
    Radius.append(m)
    #------------------------
    #Standard Error------------
    sig_a = (1/(n2-n1))*np.sqrt(sig_w[0]**2 + sig_w[1]**2)
    sig_b= np.sqrt(sig_w[1]**2 + (sig_a*n2)**2)
    sig_Transition_Sizes = np.sqrt((sig_b/b)**2 + (sig_a/a)**2)
    Std_Error.append(sig_Transition_Sizes)
    #--------------------------
    
    #plt.scatter(x_range, Wrap_list, color = 'red', s = 30)
    #plt.scatter(x_range, Wrap_list)
    #plt.errorbar(x_range, Wrap_list, yerr = standard_error_wrap, color = 'red', capsize= 3, linewidth = 2.5, elinewidth= 2.5)    For the purpose of the graph, this  part is commented out.
    





    
analytic_threshold = [180, 250]     #The threshold particle number around which the wrapping is tested.
scales = [20,20] #Each of the Density folders increase particle number by this number.
min_radius_num = 6     #The radius number 'm' must start from this number.

lattice_site_num = [2400,2700]

empty_instances_num = 0

empty_folder_num = 0

date = 'Feb10'
for m in range(6,8):
    
    meta_folder = 'Additional_Transitions_Data/Radius_' + str(m) + '/'
    Wrap_list = []
    Wrap_list_stdev = []
    Circular_list =[]
    standard_error_wrap = []
    standard_error_circle = []
    k_list = range(1,11)
    if m == 3:
        k_list = range(1,11)
    
    
    for k in k_list: #[1,2,3,4,5...12]
        wrap_rec = []
        circular_rec = []
        for i in range(1,25):
            instance = str(i)
            
            print(instance)

            
            # if len(folder) == 0:  #check if folder is empty
            #     empty_instances_num += 1
            if len(os.listdir( meta_folder + "Density_" + str(k) + "/" +  date + "_subfolder_" + instance)) == 0:  #check if folder is empty
                empty_folder_num += 1
                #empty_files.append([m, k, i])
            
            
            else:
                folderfileT = meta_folder + "Density_" + str(k) + "/" +  date + "_subfolder_" + instance + "/" + date + "_time_data_" + instance + ".txt"
                folderfileWRAP = meta_folder + "Density_" + str(k) + "/" + date + "_subfolder_" + instance + "/" + date + "_wrap_data_" + instance + ".txt"
            

                
            
                y_file = open(folderfileWRAP, "r")
                x_file = open(folderfileT, "r")
            
                y_list = read_matlab_data(y_file)
                x_list = read_matlab_data(x_file)
                x_title = 'Time (s)'
                
                for j_index in range(0,1):
            
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
        Wrap_list_stdev.append(np.std(wrap_rec))
        Circular_list.append(np.average(circular_rec))
        standard_error_wrap.append(np.std(wrap_rec) / np.sqrt(np.size(wrap_rec)))
        standard_error_circle.append(np.std(circular_rec) / np.sqrt(np.size(circular_rec)))
        #/lattice_site_num[m-1]
    x_range = [(analytic_threshold[m - min_radius_num] + scales[m - min_radius_num]*h)/lattice_site_num[m-min_radius_num] for h in range(0,10) ]
    
    for pi in range(len(Wrap_list)):
        if Wrap_list[pi] < 0.5 and max(Wrap_list[:(pi + 1)]): #This is the element before the 0.5 proability.
            before_half_index = pi
        elif Wrap_list[pi] > 0.5 and min(Wrap_list[pi:]):
            after_half_index = pi
            
    diff_from_half = np.array(Wrap_list) - 0.5
    negatives = [ow for ow in diff_from_half if ow < 0]
    before_half_index = len(negatives) - 1
    after_half_index = before_half_index + 1
    
    #For the transition graph-----
    delta_particles = scales[m - min_radius_num]
    m_slope = (Wrap_list[after_half_index] - Wrap_list[before_half_index])/delta_particles
    n1 = analytic_threshold[m - min_radius_num] + before_half_index*scales[m - min_radius_num]
    n2 = analytic_threshold[m - min_radius_num] + after_half_index*scales[m - min_radius_num]
    
    n = [n1, n2]
    w = [Wrap_list[before_half_index], Wrap_list[after_half_index]]
    sig_w = [Wrap_list_stdev[before_half_index], Wrap_list_stdev[after_half_index]]
    
    #a, b = np.polyfit(n,w,1)
    
    a = (w[1] - w[0])/(n2 - n1)
    b = w[1] - a*n2
    
    print('slope comparison ', m_slope, a)
    
    transition_cluster_size =  (0.5 - b)/a
    
    Transition_Sizes.append(transition_cluster_size)
    Radius.append(m)
    #------------------------
    #Standard Error------------
    sig_a = (1/(n2-n1))*np.sqrt(sig_w[0]**2 + sig_w[1]**2)
    sig_b= np.sqrt(sig_w[1]**2 + (sig_a*n2)**2)
    sig_Transition_Sizes = np.sqrt((sig_b/b)**2 + (sig_a/a)**2)
    Std_Error.append(sig_Transition_Sizes)
    #--------------------------
    
    #plt.scatter(x_range, Wrap_list, color = 'red', s = 30)
   # plt.scatter(x_range, Circular_list, color = 'blue', s = 20)
    # plt.plot(x_range, Wrap_list, label ='wrapped', color = 'red')
    # plt.plot(x_range, Circular_list, label = 'circle', color = 'blue')
    
    #plt.errorbar(x_range, Wrap_list, yerr = standard_error_wrap, color = 'red', capsize= 3, linewidth = 2.5, elinewidth= 2.5)     For the purpose of the graph, commented out.
    
    
    
    
plt.xlabel("Particle density")
#plt.xticks(10*np.array(range(2,7)), 10*np.array(range(2,7))/(1000))

plt.ylabel("probability")
plt.savefig('shape_probability.svg')
plt.show()


print('empty folder number', empty_instances_num)

radius_list = np.array([10,12,15,18,21,24,27])*10/(2*np.pi)

diameter_list = radius_list*2

cmap = plt.get_cmap('viridis', 4)

no_shapeTransition = (np.pi**3) * (radius_list**2) / 100
plt.plot(radius_list, no_shapeTransition, label = r'$d_{round _ cluster} =$ tube circumference', color =cmap(0), linestyle = 'dotted', linewidth = 2)
plt.scatter(radius_list, no_shapeTransition, color =cmap(0))
#plt.scatter((6 + 3*np.array(Radius))*10/(2*np.pi), Transition_Sizes)
#plt.plot((6 + 3*np.array(Radius))*10/(2*np.pi), Transition_Sizes)

#plt.scatter(radius_list, Transition_Sizes, label = 'Simulations', color = 'cyan')
#plt.plot(radius_list, Transition_Sizes,  color = 'cyan')

plt.scatter(radius_list, Transition_Sizes, color = cmap(1))
plt.plot(radius_list, Transition_Sizes, label = 'Simulation wrapping transition', color = cmap(1), linestyle = 'dashed', linewidth = 2)
plt.errorbar(radius_list, Transition_Sizes, yerr = Std_Error, linestyle = '', color = cmap(1), capsize = 7, ecolor = cmap(1))






analytic = 4*np.pi*(radius_list**2)/100
plt.plot(radius_list, analytic, label = '$E_{round} = E_{wrap}$', color = cmap(2), linewidth = 2)
plt.scatter(radius_list, analytic,  color = cmap(2))

plt.xticks(radius_list, [int(ele)/10 for ele in np.round(diameter_list, 0)])
plt.ylabel('Cluster size (proteins)', fontsize = 18)
plt.xlabel('Tube diameter ($d_{IRE1}$)', fontsize = 18)
plt.legend(frameon = False, fontsize = 12)
plt.savefig('closedsys_wrapping_size_expanded_Scaled.svg')
plt.show()

[1.5845787080869582,
 0.20806521299004174,
 0.738773900246354,
 2.0094946835982683,
 0.800271403707987,
 0.52650065366247,
 2.8207971034596824]