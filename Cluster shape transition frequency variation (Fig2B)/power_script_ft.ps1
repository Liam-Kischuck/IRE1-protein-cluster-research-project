$date = 'May27'
$prefix = 'KMC_ft_'
$original_file = 'KMC_ft_0'

#Make all the sub directories
for ($num = 0 ; $num -le 40 ; $num++){$folder_name = $date + '_' + $num; mkdir $folder_name; cd $folder_name; $subfolder_name = $date + '_subfolder_' + $num; mkdir $subfolder_name ; cd ..}

#Copy matlab_slurm_ft.sl to those sub directories.
for ($num = 0 ; $num -le 40 ; $num++){$target_folder = $date + '_' + $num; cp matlab_slurm_ft.sl $target_folder/matlab_slurm_ft$num.sl}


#Changes in matlab_slurm.sl. Change the copied sl file's execution file to the correct file name.
for ($num = 0 ; $num -le 40 ; $num++){$target_folder = $date + '_' + $num; cd $target_folder; (Get-Content matlab_slurm_ft$num.sl).replace($original_file , $prefix + $num) | Set-Content matlab_slurm_ft$num.sl; cd ..}

#copy the matlab file to each directory.
for ($num = 0 ; $num -le 40 ; $num++){$target_folder = $date + '_' + $num; $make_file = $prefix + $num; cp KMC_ft_0.m $target_folder/$make_file.m}

#Changes to the matlab file. Parameter changes. You can add changes like to the name of the file, and other parameters.
for ($num = 0 ; $num -le 40 ; $num++){$target_folder = $date + '_' + $num; cd $target_folder; $make_file = $prefix + $num + '.m'; (Get-Content $make_file).replace('instance = "0"', 'instance = "' + $num + '"') | Set-Content $make_file; cd ..}
for ($num = 0 ; $num -le 40 ; $num++){$target_folder = $date + '_' + $num; cd $target_folder; $make_file = $prefix + $num + '.m'; (Get-Content $make_file).replace('par_num = 20', 'par_num = ' + ($num + 20) ) | Set-Content $make_file; cd ..}
