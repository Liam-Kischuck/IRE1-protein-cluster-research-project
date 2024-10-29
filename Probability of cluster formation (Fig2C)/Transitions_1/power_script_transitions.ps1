$date = 'Aug31'
$prefix = 'KMC_Transitions_'
$original_file = 'KMC_Transitions_0'
$start_density_list = @(40,70,100,140)


$start_index = 1
$final_index = 4

#This is to make higher-up directories. This for loop ends at the end of the file below
for ($m = 1 ; $m -le 4 ; $m++){$folder_name =  'Radius_' + $m; mkdir $folder_name; cp KMC_Transitions_0.m $folder_name/KMC_Transitions_0.m ; 
cp matlab_slurm_0.sl $folder_name/matlab_slurm_0.sl ;  cd $folder_name; 

#This makes the secondary directories that groups together the same particle dnesities.
for ($k = 1 ; $k -le 10 ; $k++){$secondfolder_name =  'Density_' + $k; mkdir $secondfolder_name; cp KMC_Transitions_0.m $secondfolder_name/KMC_Transitions_0.m ; 
cp matlab_slurm_0.sl $secondfolder_name/matlab_slurm_0.sl ;  cd $secondfolder_name; 

#Make all the sub directories

for ($num = $start_index ; $num -le $final_index ; $num++){$target_folder = $date + '_' + $num; mkdir $target_folder; cd $target_folder; $subfolder_name = $date + '_subfolder_' + $num; mkdir $subfolder_name ; cd ..}


for ($num = $start_index ; $num -le $final_index ; $num++){$target_folder = $date + '_' + $num; cp matlab_slurm_0.sl $target_folder/matlab_slurm_$num.sl}


#Changes in matlab_slurm.sl. Change the copied sl file's execution file to the correct file name.
for ($num = $start_index ; $num -le $final_index ; $num++){$target_folder = $date + '_' + $num; cd $target_folder; (Get-Content matlab_slurm_$num.sl).replace($original_file , $prefix + $num) | Set-Content matlab_slurm_$num.sl; cd ..}

#copy the matlab file to each directory.
for ($num = $start_index ; $num -le $final_index ; $num++){$target_folder = $date + '_' + $num; $make_file = $prefix + $num; cp KMC_Transitions_0.m $target_folder/$make_file.m}

#Changes to the matlab file. Parameter changes. You can add changes like to the name of the file, and other parameters.
for ($num = $start_index ; $num -le $final_index ; $num++){$target_folder = $date + '_' + $num; cd $target_folder; $make_file = $prefix + $num + '.m'; (Get-Content $make_file).replace('instance = "0"', 'instance = "' + $num + '"') | Set-Content $make_file; cd ..}
for ($num = $start_index ; $num -le $final_index ; $num++){$target_folder = $date + '_' + $num; cd $target_folder; $make_file = $prefix + $num + '.m'; (Get-Content $make_file).replace('par_num = 30', 'par_num = ' + ($start_density_list[$m - 1] + 10*($k - 2)) ) | Set-Content $make_file; cd ..}
for ($num = $start_index ; $num -le $final_index ; $num++){$target_folder = $date + '_' + $num; cd $target_folder; $make_file = $prefix + $num + '.m'; (Get-Content $make_file).replace('Li = 9', 'Li = ' + (3*($m) + 9) ) | Set-Content $make_file; cd ..}

 Set-Location ..}

 Set-Location ..}