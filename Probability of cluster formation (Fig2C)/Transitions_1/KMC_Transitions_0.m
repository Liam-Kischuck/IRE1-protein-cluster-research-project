%This file adds to the kMC with collective motion file
%"Collective_Fast_kMC_Endgame3.m", an added property where new particles
%are introduced to the tube with a certain proability, as well as the
%possibility for particles to exit the tube. Clusters, in this 
%implementation can partially exit from the tube. (when a part of the
%cluster crosses the boundary of the tube, those particle that belong to
%the clutser that cross are considered to have exitted, and due to this,
%the cluster shrinks. 

%This file is trying to speed up the "Collective_kMC_Tube_quicker.m" by
%taking out the for loop in count_neighbor_nf. -> it turns out it was
%taking up a lot of time. 10 times faster code!

%The LOC length can vary, and LOOK is a 1x Li*Lj list that is adjusted
%everytime LOC length is extended to include more particles. --> I found
%this to give about 1/3 times faster than the Collective_kMC_Tube_finale.m
%code.

%The unnecesary A matrix has been taken out of the computation. I still use
%it at the beginning when setting up the C matrix.


%Maybe only use this file when you have to check for just a few phase
%diagram values. It generates a lot of files of Dacccccta.

%Feb 15th 2022, A cluster exit probability was added that scales as 1/n,
%where n is the 'width' of the cluster. So if a cluster has a latitudinal
%length (width) n and is made up of N particles, then the exit probability
%is D/(4*n*sqrt(N)).

%March 7th, this version allow multiple points to be computed per phase
%point. cluster_phase_terminateNUM is the number of times the algorithm
%must return 'cluster' phase before terminating. Usually I set this to 2,
%but here I have set it to 4 (March 7th).

%Mar19: This version contains corrections to the lifetime of dimers (where
%improper cluster diffusion was breaking the clusters up too quick), and to
%the particles that would easily leave just because the systemenergy
%doesn't change.

%April 4th: This is a corrected version of the code coming from
%KMCo_Saffman_finale_2.m in the KMCo folder. The exit particle function had
%an issue in that file, so it has been corrected.

%April 13th: This version has corrected for the rare-case bug that would
%make the code go on forever. However, the E_long E_lati 3 when it should
%be class 0 issue has not been corrected yet! This is just to get some
%results for tmr, thus the name interim.

%April 19th: This fixed the E_long,lati issues where the particles
%surounded by particles showed some values ofescaping when the E_long,lati
%values should be 0. The fix partly involved cluster not getting near the
%edge, so even if you allow cluster exit, in this version clusters won't be
%able to exit unless you change the function in collective_shifter to allow
%clusters to get near the edges of the tube.

%May30th: To be able to save all data to one folder, to make it easier to
%transfer data from ComputeCanada, all data is now stored in a single
%folder '"date"_subfolder_"instance"'.

%Making this a closed system for the shape transition simulations. The
%cluster movement is also shut-off.


Li = 9;
Lj = 100;

date = "Aug31";
instance = "0";
par_num = 30;
fixed_J_value = 3;
take_record = true;

exit_cutoff = 2;   %Any cluster below this size (not including this size clusters) could exit the tube.
return_prob = 1; %The probability of a exitting particle to return backinto the tube. To be real, this should be 1.
final_iteration = 70000000;     %The maximum number of iterations allowed for the algorithm (usually 25,000,000 or less is sufficient)
terminal_comp_time =  5*3600;    %Approximately the maximum amount of time you let the algorithm run for. In second.
measure_interval = 20000;        %Specify how often you make measurements. 150,000 measurements makes a csv file of 300kB. Try and not go over this if reasonable.
                                %The number of measurements is given by the
                                %(total number of kMC_iterations)/(measure_interval).
max_run_time = 1;
                                
rng shuffle

% rho_initial = 0.00001;
% rho_end = 0.00001;

[J_hat_min, J_hat_max, J_hat_num] = deal(fixed_J_value, fixed_J_value, 6)
[rho_min, rho_max, rho_num] = deal(par_num/(Li*Lj), par_num/(Li*Lj), 1)
%14.4f
snapshots = false;

if J_hat_num > 1
    J_hat_scale = (J_hat_max - J_hat_min)/(J_hat_num - 1);
elseif (J_hat_num > 1) & (J_hat_max == J_hat_min)
    J_hat_scale = 0;     % if repeating the same J_hat value, no need for scaling.
elseif J_hat_num == 1
    J_hat_scale = 1;
elseif J_hat_min == J_hat_max
    J_hat_scale = 0;
end

%For testing----
%J_hat_scale = 1.18   %1.18^20 = 27.4 which is the range of values we are looking at.
%---------------

if rho_num > 1
    rho_scale = (rho_max - rho_min)/(rho_num - 1);
elseif rho_num == 1
    rho_scale = 1;
elseif rho_min == rho_max
    rho_scale = 0;
end

x = [];
y = [];
colour = [];

N_REC = {};      %The full record as a function of time of the number, andweight average cluster size.
TIME_REC = {};
WEIGHT_REC = {};
WRAP_REC = {};

N_rec = [];
Time_rec = [];
Weight_rec = [];

phase = ["monomer", "dimer", "oligomer", "nanocluster", "cluster", "intermediate"];
phase_num = [1,2,3,4,5,6];
Phase_to_num_dict = containers.Map(phase,phase_num);
Phase_to_num_dict;


CLUS_AVG_SIZE = {};
nanoCLUS_AVG_SIZE = {};
OLIG_AVG_SIZE = {};
DIM_AVG_SIZE = {};

CLUS_NUM = {};
nanoCLUS_NUM = {};
OLIG_NUM = {};  
DIME_NUM = {};
MONO_NUM = {};


AIS_Clus_avg_size = [];     %Recording the average size of the cluster island size ,at the end of the simulation (the AIS prefix is added).
AIS_nanoClus_avg_size = [];
AIS_Olig_avg_size = [];  
AIS_Dim_avg_size = [];

AIS_Clus_num = [];          %Recording the number of cluster for a given configuration, at the end of the simulation (the AIS prefix is added)
AIS_nanoClus_num = [];
AIS_Olig_num = [];
AIS_Dime_num = [];
AIS_Mono_num = [];

LI = [];
LJ = [];
KMC_COUNTER = [];

toc_time = [];
counter = 0;

DPE_list = [];



for rho_int = 0:(rho_num - 1)
    for J_hat_int = 0:(J_hat_num - 1)
        counter = counter + 1
        Rho = rho_min + rho_scale*rho_int
        x = [x, rho_min + rho_scale*rho_int];
        y = [y, J_hat_min + J_hat_scale*J_hat_int];
        J_hat = J_hat_min + J_hat_scale*J_hat_int
        %For testing----------------------------------
        %J_hat =  J_hat_scale^J_hat_int
        %---------------------------------------------
        
        N_particles = (rho_min + rho_scale*rho_int)*Li*Lj
        tic
        [DPE, wr, tr, nr, wrapr, AIS, island_info, found, kMC_iterations] = Phase_Calculator(J_hat, N_particles, Li, Lj, snapshots, exit_cutoff, return_prob, final_iteration, terminal_comp_time, measure_interval, max_run_time);     

        toc_time = [toc_time, toc/3600];
        
        DPE_list = [DPE_list, DPE]
        DPE = Phase_to_num_dict(DPE);
        colour = [colour, DPE];
        
        AIS_Clus_avg_size(end + 1) = AIS('avg_cluster');
        AIS_nanoClus_avg_size(end + 1) = AIS('avg_nanocluster');
        AIS_Olig_avg_size(end + 1) = AIS('avg_oligomer');
        AIS_Dim_avg_size(end + 1) = AIS('avg_dimer');

        AIS_Clus_num(end + 1) = AIS('num_cluster');
        AIS_nanoClus_num(end + 1) = AIS('num_nanocluster');
        AIS_Olig_num(end + 1) = AIS('num_oligomer');  
        AIS_Dime_num(end + 1) = AIS('num_dimer');
        AIS_Mono_num(end + 1) = AIS('num_monomer');
        
        N_rec(end + 1) = nr(end);
        Time_rec(end + 1) = tr(end);
        Weight_rec(end + 1) = wr(end);
        
        %Add a section to store the Clus_num, Olig_avg_size etc cell array
        %lists. The Phase_Claculator function must return all these arrays!
        %The Data is to be stored in a file below.
        
        N_REC = [N_REC; nr];
        WEIGHT_REC = [WEIGHT_REC; wr];
        TIME_REC = [TIME_REC; tr];
        WRAP_REC = [WRAP_REC; wrapr];
        
        %island_info has stored all the average island size and the number of islands in this order of 1,2,...,9. 
        CLUS_AVG_SIZE = [CLUS_AVG_SIZE; island_info{1}];
        nanoCLUS_AVG_SIZE = [nanoCLUS_AVG_SIZE; island_info{2}];
        OLIG_AVG_SIZE = [OLIG_AVG_SIZE; island_info{3}];
        DIM_AVG_SIZE = [DIM_AVG_SIZE; island_info{4}];

        CLUS_NUM = [CLUS_NUM; island_info{5}];
        nanoCLUS_NUM = [nanoCLUS_NUM; island_info{6}];
        OLIG_NUM = [OLIG_NUM; island_info{7}];  
        DIME_NUM = [DIME_NUM; island_info{8}];
        MONO_NUM = [MONO_NUM; island_info{9}];
        
        LI = [LI, Li];
        LJ = [LJ, Lj];
        KMC_COUNTER = [KMC_COUNTER, kMC_iterations];
        
        
%         if found == true
%             beep
%             return
%         end
    end
end
toc_time
f = figure()
scatter(x,y, 240, colour, 'filled')
colormap(parula(6));
colorbar('Ticks', linspace(1,6, 6) , 'TickLabels', ["monomer", "Dimer", "Oligomer", "Nano-Cluster", "Cluster", "Intermediate"])
xlabel('Number Density');
ylabel('J (kT)');
shg        % temporarily not show phase diagram


%-------------------writing information to file
if take_record == true
    folder = date + "_subfolder_" + instance + "/";
    A = [y; x; DPE_list; AIS_Mono_num; AIS_Dime_num; AIS_Olig_num; AIS_nanoClus_num; AIS_Clus_num; AIS_Dim_avg_size; AIS_Olig_avg_size; AIS_nanoClus_avg_size; AIS_Clus_avg_size; N_rec; Time_rec; Weight_rec; LI; LJ; toc_time; KMC_COUNTER];
    fileID = fopen(folder + date +  "_KMCf_data_" + instance + ".txt",'w');
    fprintf(fileID,'%14s %14s  %15s %14s   %14s   %14s   %14s   %14s   %14s   %14s   %14s   %14s   %14s   %14s   %14s  %14s   %14s  %14s  %14s\n','J', 'Rho', 'Dominant_Phase', 'Mono_num', 'Dime_num', 'Olig_num', 'nanoClus_num', 'Clus_num', 'Dim_avg_size', 'Olig_avg_size', 'nanoClus_avg_size', 'Clus_avg_size', 'N_record', 'Time_record', 'Weight_record', 'Li', 'Lj', 'toc_time (hrs)', 'kMC_iterations');
    fprintf(fileID,'%14.4f %14.8f %15s %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f %14.5f  %14.5f %14.5f\n', A);
    fclose(fileID);
    %-------------------

    %Store all the Clus_num, Olig_avg_sizes cell arrays into a txt file here.


    % WEIGHT_REC
    % TIME_REC
    % N_REC
    writecell(WEIGHT_REC,folder + date + "_weight_data_" + instance + ".txt")
    writecell(TIME_REC,folder + date +  "_time_data_" + instance + ".txt")
    writecell(N_REC,folder + date + "_N_data_" + instance + ".txt")
    writecell(WRAP_REC,folder + date + "_wrap_data_" + instance + ".txt")

    writecell(CLUS_AVG_SIZE,folder + date + "_clus_avg_data_" + instance + ".txt")
    writecell(nanoCLUS_AVG_SIZE,folder + date + "_nanoclus_avg_data_" + instance + ".txt")
    writecell(OLIG_AVG_SIZE,folder + date + "_olig_avg_data_" + instance + ".txt")
    writecell(DIM_AVG_SIZE,folder + date + "_dim_avg_data_" + instance + ".txt")

    writecell(CLUS_NUM,folder + date + "_clus_num_data_" + instance + ".txt")
    writecell(nanoCLUS_NUM,folder + date + "_nanoclus_num_data_" + instance + ".txt")
    writecell(OLIG_NUM,folder + date + "_olig_num_data_" + instance + ".txt")
    writecell(DIME_NUM,folder + date + "_dime_num_data_" + instance + ".txt")
    writecell(MONO_NUM,folder + date + "_mono_num_data_" + instance + ".txt")
end
%=====================================================================================================
%Below (before the next 'double line', is the whole calculation for kMC for
%a given J_hat, and density (N_particle/(Li x Lj))

function [DPE, weight_record, time_record, N_record, wrap_record, AIS, island_info, found, kMC_iterations] = Phase_Calculator(J_hat, N_particles, Li, Lj, snapshots, exit_cutoff, return_prob, final_iteration, terminal_comp_time, measure_interval, max_run_time)
    %DPE stands for Dominant Phase Estimate.
    %parameters such as , termination_length, kMCS (the number of sweeps),
    %how often the measurments are taken etc, at the moment (Nov 2021) must be 
    %adjusted wihtin the function.
    
    A = zeros(Li,Lj);
    %J_hat;      This is J/(k_B*T), so is J_hat = 2, then J = 2*k_B*T.

    FPC = {};   %stands for Fraction of Phases Collection.
    DoP_termination_record = [];   %stores the record of the last five phase determined
                                   %by the phase analysis (checking what phase
                                   %the configuration is in.
    termination_length = 12;
    cluster_phase_terminateNUM = 8;    %The number of termination records needed for a cluster phaser system to terminate. Usually this is 2.
%     if J_hat > 6
%         termination_length = 6;
%     end

    %just checking the number of particles that are added and exitted.
    added = 0;
    exitted = 0;
    c_exit = 0;
    
    %Run time parameters
    rho = N_particles/(Li*Lj);  %This is the initial density of IRE1 over the whole ER per tube.
    kMC = 0;
    kMC_0 = 1000;
    continue_monte = true;
    
    %YOU MAY WANT TO CHANGE THIS DEPENDING ON HOW LONG YOU ARE RUNNING THE
    %PROGRAM!=============================================================
    measure_interval_list = [5,5,5,10,20,30,40,40,40];  %if runnning 3million 
    %times this give 3000,000/5= 600,000 measurements for the first three J_hat values, and so on.
    %=====================================================================
    round_Jhat = round(J_hat);
    
    measure_interval = measure_interval;%measure_interval_list(round_Jhat);  %usually = 100 etc.
    initial_analysis_interval = 1000;
    max_analysis_interval = 100000;
    final_terminate = final_iteration;
    terminal_computation_time = terminal_comp_time;
    interval_scale = 1.5;
    analysis_interval = initial_analysis_interval; %this is updated as it runs.
    previous_analysis_kMC = kMC_0;
    weight_record = [];
    time_record = [];
    N_record = [];
    wrap_record = [];
    termination = false;   %This becomes 'true' when the algorithm tries to end once, but is not cluster phase. It runs a little longer then exits.
    
    Clus_avg_size = [];     %Recording the average size of the cluster island size.
    nanoClus_avg_size = [];
    Olig_avg_size = [];  
    Dim_avg_size = [];

    Clus_num = [];          %Recording the number of cluster for a given configuration
    nanoClus_num = [];
    Olig_num = [];
    Dime_num = [];
    Mono_num = [];
    
    %Initialise the lattice. 
%    NumPoints = N_particles;
%     while nnz(A) < NumPoints
%         a = mod(round( Li*rand() ) ,Li) + 1;
%         b = mod(round( Lj*rand() ), Lj) + 1;
% 
%         A(a,b) = 1;
%     end

%For circular cluster shape to start with---------
    Initial_N_particles = N_particles;
    radius = sqrt((Initial_N_particles)/pi);   %The plus 10 gives a buffer for menough particles to be included.
    num_of_particles = 0;
    t = 0;
    Initial_N_fits = false
    %Iteratively try to make a cluster of size 200. This is tricky since
    %some tube sizes are to small so simply setting a radius will not give a
    %proper clsuter of size 20. You need to increase the radius to fit 200
    %particles in a cylindrical wrapped cluster.
    while Initial_N_fits == false
        A = [];
        A = zeros(Li,Lj);
        num_of_particles = 0;
        for i = 1:Li
            for j = 1:Lj
                if (sqrt(((i - (Li/2))^2 + (j - (Lj/2))^2)) < radius) & (num_of_particles < Initial_N_particles)
                    A(i, j) = 1;
                    num_of_particles = num_of_particles + 1;
                    %NumPoints = NumPoints + 1;
                end
            end
        end
        t = t + 1;
        if num_of_particles < Initial_N_particles 
            radius = sqrt((Initial_N_particles + 0.01*t)/pi); 
        elseif num_of_particles >= Initial_N_particles
            Initial_N_fits = true;
        end  
        
    end
    t
    num_of_particles
    
    disp('confirme number...')
    nnz(A)

    if (sum(A(1,:)) > 0) & (sum(A(Li,:)) > 0)  %The case where the cluster pretty much wraps around the tube.
        A = []; %initialize lattice to make the cluster perfectly cylindrical.
        A = zeros(Li,Lj);
        cluster_width = floor(Initial_N_particles/Li); %The quotient.
        remaining_particles = mod(Initial_N_particles, Li) %The remainder.
        center = Lj/2;


        i_pl_max = Li
        j_pl_max = cluster_width
        for i_pl = 1:i_pl_max
            for j_pl = 1:j_pl_max
                A( i_pl, Lj/2 - floor((j_pl_max/2)) +j_pl)  = 1;
            end
        end
        %Add the remaining particles to the right of the cylindrical cluster.
        right_edge = Lj/2 - floor((j_pl_max/2)) + j_pl_max;
        left_edge = Lj/2 - floor((j_pl_max/2)) + 1;
        disp(remaining_particles)
        remaining_particles = round(mod(Initial_N_particles, Li)) %The remainder.
        for p = 1:remaining_particles
            if mod(p,2) == 0  %If even, add to the left.
                A(p, left_edge - 1) = 1;
            elseif mod(p,2) == 1 %If odd, add to the right.
                A(p, right_edge + 1) = 1;
            end
            disp('addedddddd')
            disp(p)
            disp(remaining_particles)
        end
    end
    disp('confirme cylindrical number...')
    nnz(A)

    %Make the concentration around the cluster right.
    %This is taken out for the closed system.
%     for i_pl = 1:Li
%         for j_pl = 1:Lj
%             if A(i_pl, j_pl) ~= 1
%                 R = rand();
%                 if R < rho
%                     A( i_pl, j_pl)  = 1;
%                 end
%             end
%         end
%     end



%A(Li/2 - floor((i_pl_max/2)) + i_pl, Lj/2 - floor((j_pl_max/2)) +j_pl)  = 1;
            

    % A(3,3) = 1;
    % A(4,3) = 1;
    % A(3,4) = 1;
    % A(4,4) = 1;
    % A(1,3) = 1;
    % A(4,1) = 1;
    % A(1,2) = 1;
    % A(3,2) = 1;
    % 
    % 
%     A = [0,1,0,0,0,1,1;
%         0,1,0,0,0,0,0;
%         0,0,0,1,0,0,0;
%         0,1,0,0,0,0,0;
%         0,0,0,0,0,0,1;
%         0,0,0,0,0,0,0]
    
%     A = [0,  0,  1,  1;
%          1,  0,  1,  1;
%          1,  1,  1,  1;
%          1,  1,  1,  1;
%          0,  1,  1,  1;];
%      
%     A = [0,  0,  0,  1,  0,  0;
%          0,  1,  0,  1,  0,  0;
%          0,  0,  0,  0,  1,  0;
%          0,  0,  0,  1,  0,  0;
%          1,  1,  0,  1,  0,  0;
%          0,  0,  0,  0,  0,  1;];
%      
%     C =  [0,  0,  0,  3,  0,  0;
%           0,  5,  0,  2,  0,  0;
%           0,  0,  0,  0,  4,  0;
%           0,  0,  0,  6,  0,  0;
%           9,  8,  0,  1,  0,  0;
%           0,  0,  0,  0,  0,  7;];
     
     
    N_particles = nnz(A);
    NumPoints = nnz(A)
           %nnz(A) gives you the number of particles.
    ONES0 = find(A==1);

    ONES = transpose(ONES0);

    site =[];

    for i = 1:length(ONES)
        j = 1 + floor((ONES(i) - 1)/size(A,1));
        site =[site; ONES(i) - (j - 1)*size(A,1),  j ];    %This 'site' code was corrected Oct 27 2021. Everything coded
    end                                                    %using this 'site' function should be updated to this corrected one.

    
    
    
    %This little section is a totally different kind of 'site' function
    %compared to usual where it lists the elements in
%     %C.==============================================================
%     ONES0 = find(C~=0)
% 
%     ONES = transpose(ONES0)
% 
%     site = zeros(length(ONES), 2);
% 
%     for i = 1:length(ONES)
%         j = 1 + floor((ONES(i) - 1)/size(C,1));
%         site(C(ONES(i)),:) = [ ONES(i) - (j - 1)*size(C,1),  j ];    %This 'site' code was corrected Oct 27 2021. Everything coded
%     end 
    %================================================================

    site_dim = size(site,1);
    
    N_particles = site_dim;                 %The random positioning doesn't always give the same number of particles as the initial N_particle #.

    C = A;                 % for each 'site' index (1 to ProteinNum), C keeps track
                           % of the position of that site index on the actual
                           % lattice. 

    for i = 1:site_dim
        C(site(i,1),site(i,2)) = i;
    end

    C;
%temporarily----------------------------------    
    if snapshots == true
        x = [];
        y = [];

        %Graph
        for i =1:Li
            for j = 1:Lj
                if A(i, j) == 1
                    x = [x, j];
                    y = [y, (Li + 1) - i];    %making the graph reflect matrix positions.
                end
            end
        end
        f1 = figure()

        scatter(x,y, 40,'blue','filled');
        xlim([0 Lj]) ;
        ylim([0 Li]);
        shg
    end
%--------------------------------------------

    %The following block is to calculate the modulos
    Cdim = size(C);
    Li = Cdim(1);
    for i = 1: Li
        ip_i(i) = i + 1;
        im_i(i) = i - 1;
    end
    ip_i(Li) = 1;
    im_i(1) = Li;

    Lj = Cdim(2);
    for j = 1: Lj
    ip_j(j) = j + 1;
    im_j(j) = j - 1;
    end
    ip_j(Lj) = 0;          %no-flux b.c.
    im_j(1) = 0;

    res = [];
    


    E_lati = zeros(Li, Lj - 1);  % Li by (Lj - 1)  No-Flux Boundary
    E_long = zeros(Li, Lj);  % Li by Lj   Periodic Boundary

    LOC_LATI = zeros(7,4*N_particles);   %The seven fold way LOC  array.      4*N_particles is the most 'relevant' bonds that there can possibly be.
    LOC_LONG = zeros(7,4*N_particles);   %The seven fold way LOC  array.
    CEP_lati = zeros(1,7);        %stands for Class End Points. It stores the index of the
    CEP_long = zeros(1,7);        %last index of some class. Very important for speeding up and not
                                  %needing to use cell arrays. Bonds that are 0, AND do
                                  %not involve a relevant spin-exchange,
                                  %have indice greater than CEP_lati or
                                  %long(7).
    LOOK_LATI = zeros(1, Li*Lj);      %The index of LOOk_LATI is the (latitudinal) bond number and the value the position in the LOC_LATI list.
    LOOK_LONG = zeros(1, Li*Lj);
    
    [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = prepare_interaction(E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, Li, Lj);



    
    


%========================================================================
%Introducing the collective moves
%Make CTM
N_particles = length(ONES);
num_column = 0;
for i = 2:(floor(N_particles/2) + 1)
    if i < floor(N_particles/2) + 1
        num_column = num_column + floor(N_particles/i);
    elseif i == floor(N_particles/2) + 1
        num_column = num_column + (N_particles - floor(N_particles/2));
    end    
end
CTM = zeros(N_particles, num_column);
max_grp_num = 0;

nbrs_info = zeros(N_particles, 4);  
for particle = 1:N_particles
    nbrs_near = neighbors(particle, site, C, im_i, ip_j, ip_i, im_j);
    nbrs_info(particle, 1:length(nbrs_near)) = nbrs_near;
end
nbrs_info;

%makes list of clusters.
CTM_site = zeros(N_particles, 2);   %Keeps track of location ofparticle i in CTM matrix.
prep = true;
for i  = 1:N_particles
    Cluster_elements = [i];
    [CTM, CTM_site, max_grp_num] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
end
prep = false;

%I have no idea why this below prepare_interaction is used..., it was
%skrewing up the energies in E_long and E_lati.
% disp("step 1:")
% E_long
% [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = prepare_interaction(E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, Li, Lj, CTM_site);
% 
% E_long
% 
% disp("End prep")

CTM;
CTM_site;
E_long;
E_lati;
C;
LOC_LONG;
LOC_LATI;


%Get rid of the need to calculate exponential every time.
exp_values = [];
for i = -3:3
    exp_values = [exp_values, min(exp(-J_hat*i), 1)];
end
exp_values;

tic
%From here start the kMC step!

%Delta_t = 0.000026;  %You get this from the diffusion coefficient of D = 0.24 um^2/s, and Delta x = 5nm.

Delta_t = 0.0001;  %You get this from the diffusion coefficient of D = 0.24 um^2/s, and Delta x = 10nm.
%%R_mem = 50;
%%D = 0.24*10^  6; what is this 10^6?
%R_ire1 = 2.8;
%A = D/log(R_mem/R_ire1);

%R_shrunk = Li*5/(2*pi) + 4;  %R_shrunk depends on the circumference of the tube. If it is a typical tube 100nm diameter, then R_shrunk = R_mem = 50.
%max_cluster_n = (R_shrunk/R_ire1)^2;   %Any cluster with a larger number of particles will give a negative D and also 
                                       %cause problems when picking the rates. So just keep the diffusion constant for cluster larrger than
                                       %this size.
% if round(A*(log(R_shrunk/R_ire1) + log(1/sqrt(max_cluster_n)) )/D, 9) <= 0
%     max_cluster_n = max_cluster_n - 1; %If the D is insanely small, make the max_cluster_n one particle smaller.
% end


run_time = 0;

found = false;
tic
%From here start the kMC step!
CEP_lati
CEP_long
J_hat
exp_values

collective_count = 0;
exitters_table = [];



%for tryen = 1:1            %To switch between a deterministic running time (setting a fixed number of iterations
while continue_monte == true    %to own using the while loop that
                                     %terminates, you just need to comment out one of these line. (Pick
                                     %either the for loop, of the while loop; it's that simple.
        kMC = kMC + 1;
       

        Q = [];
        for i = 1:7             %The seven fold method cumulative probability, times two since I separate longitudinal and latitudinal.
            if i == 1
                Q = [CEP_lati(i)*exp_values(i)/4];
            else
            Q = [Q, CEP_lati(i)*exp_values(i)/4 + Q(end)];
            end
        end
        for i = 1:7             %The seven fold method cumulative probability, longitudinal case.
            Q = [Q, CEP_long(i)*exp_values(i)/4 + Q(end)];
        end

%Stationary Cluster decay commented out--------------------------------------------
%         for i = 1:max_grp_num
%             size_of_cluster = nnz(CTM(:, i));
% %             D_mem = A*(log(R_shrunk/R_ire1) +
% %             log(1/sqrt(size_of_cluster)) )/D;
% % 
% %             if size_of_cluster > max_cluster_n
% %                 size_of_cluster = max_cluster_n; D_mem =
% %                 A*(log(R_shrunk/R_ire1) + log(1/sqrt(size_of_cluster))
% %                 )/D;
% %             end
%             %Q = [Q,  D_mem + Q(end)];    %Don't divide by 4, since this
%             %move includes all possible directions of a cluster moving.
%                 %normally 1/(sqrt(size_of_cluster))
%             Q = [Q,  1/(sqrt(size_of_cluster)) + Q(end)];
%         end    
%----------------------------------------------------------------------------------

        %particles that can exit.

        %The below is commented for the closed system==============
%         N_exitter_list = [];
%         exitters_table = [];
%         particle = 0;    %I don't think this will interfere with any other functions. The other functions that do use this variable
%                         %name properly define it within their function code.
% 
%         [N_exitter_list, exitters_table] = particles_at_edges(CTM_site, C, site, im_i, ip_j, ip_i, im_j);
%          for i = 1:4
%              Q = [Q,  N_exitter_list(i)*exp_values(i + 3)*return_prob/4 + Q(end)];      %i + 3 is because we want the 0,1,2,3 values of exp_values
% %                                                                            %which are stored in index 4,5,6,7 of the list. so for i = 1,2,3,4
% %                                                                            %i + 3 give 4,5,6,7
%          end 
        %======================================================================
        %Chances of a particle entering.
        %Q = [Q,  chances_of_addition(C, Li, Lj, N_particles, rho) + Q(end)] ;
        
        %Chances of a Cluster (partially) exitting.
        partial_exit_clusters = clusters_at_edges(CTM_site, C, CTM, exit_cutoff);
%         if length(partial_exit_clusters) > 0
%             for i = partial_exit_clusters
%                 size_of_cluster = nnz(CTM(:, i));
%                 
%                 width = cluster_width_calculator(i, CTM, site, Lj) ; %Calculate the latitudinal width of the cluster.
%                 if width ~= false        %Clusters with width longer than half the length of the tube are rejected.
%                     Q = [Q,  1/(4*width*sqrt(size_of_cluster)) + Q(end)];
%                 end
%             end  
%         end
        
        Q_end = Q(end);
         Q;
%          if sum(Q) == 0
%             C
%             E_long
%             E_lati
%             LOC_LONG
%             LOC_LATI
%         end
        
        Q = Q/Q(end);
        
        run_time = run_time - Delta_t*log(rand)/Q_end;    %testing log(rand) run_time.

        R = round(rand,9);
        while R == 0                 %It can be 0 twice in a row!
            R = round(rand,9);       %A patch-up to avaoid error with R = 0, could think more about how to fix this.
        end
        
        
        
        count = 1;
        Qclass = 0;
        while Qclass == 0
            if count == 1
                if R <= Q(count)
                    Qclass = count;
                else
                    count = count + 1;
                end

            else
                if (R > Q(count - 1)) & (R <= Q(count))
                    Qclass = count;
                else
                    count = count + 1;
                end
            end
        end
        %0.4911    0.4917    0.2954    0.3025    0.4683    1.0374    0.5070
        %J_hat
        %R
        Q;
        
       

        
        
        
        if Qclass <= 14        %A single particle move.
            %NEXT, CODE THE RANDI SECTION WHERE YOU PICK THE ACTUAL EXCHNGE THAT WILL
            %TAKE PLACE. Use randi([1,6],1,1), which gives a random integer between 1 and 6.
            %disp("single")
%             Q
%             C
%             CTM
%             LOC_LONG
%             LOOK_LONG
%             LOC_LATI
%             LOOK_LATI
            Qclass;
            if Qclass < 8
                Chosen_class = LOC_LATI(Qclass, 1:CEP_lati(Qclass));
            elseif Qclass >= 8 & Qclass < 15
                Chosen_class = LOC_LONG(Qclass - 7, 1:CEP_long(Qclass - 7));
            end
            E_long;
            E_lati;
            Chosen_class;
            numel(Chosen_class);
            bond_num = Chosen_class(randi([1,numel(Chosen_class)], 1,1));   %Choose a bond thin the chosen class.


            %For testing puposes-----
    %             bond_num = 30
    %             Qclass = 4 + 7
            %------------------------

            if Qclass < 8
                j = ceil(bond_num/size(E_lati, 1));
                i = bond_num - size(E_lati,1)*(j - 1);
                interact = [i, j; i, j + 1];
            elseif (Qclass >= 8) & (Qclass < 15)
                j = ceil(bond_num/size(E_long, 1));
                i = bond_num - size(E_long,1)*(j - 1);
                interact = [i, j; ip_i(i), j];    
            end

            Qclass;
            bond_num;
            interact;

%             temp_1 = A(interact(1,1),interact(1,2));
%             temp_2 = A(interact(2,1), interact(2,2));
%             A(interact(1,1),interact(1,2)) = temp_2;
%             A(interact(2,1), interact(2,2)) = temp_1;


            % %update the position in the site of the ONE values. +1  to the
            % %right.
            % site(random, 2) = ip_j(j); 

            %Update positions on the C matrix
            temp_1 = C(interact(1,1),interact(1,2));
            temp_2 = C(interact(2,1), interact(2,2));
            if temp_1 < 0 | temp_2 < 0
                C
                interact
                temp_1
                temp_2
                site
            end

            C(interact(1,1),interact(1,2)) = temp_2;
            C(interact(2,1), interact(2,2)) = temp_1;
            
            

            if temp_1 ~= 0 
                site(temp_1,:) = [interact(2,1), interact(2,2)];
            end
            if temp_2 ~= 0
                site(temp_2,:) = [interact(1,1), interact(1,2)];
            end    
            %fprintf('From here')
            %disp("From here to the moon and back")
            
            
            [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            
            %CTM update for single moves==================================
            if sum(0 == [temp_1, temp_2]) == 1 %If either one of the interact sites is a particle and the other is a vacancy (it can't be both particles).    
                Cluster_elements = nonzeros([temp_1, temp_2]);

                [CTM, CTM_site, max_grp_num, remaining_cluster_elements] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
                %remaining_cluster_elements = nonzeros(CTM(:,1));
                interact;
                temp_1;
                temp_2;
                remaining_cluster_elements;
                if length(remaining_cluster_elements) > 0
                    prep = true;
                    [CTM, CTM_site, max_grp_num] = CTM_updater_2(CTM, CTM_site, max_grp_num, interact, remaining_cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
                    prep = false;
                end
                    CTM_site;
                CTM;
                
%             C
%             CTM
%             LOC_LONG
%             LOOK_LONG
%             LOC_LATI
%             LOOK_LATI
            end
            

%Stationary Cluster Decay commented out --------------------------------        
%         elseif (Qclass > 14)  & (Qclass <= 14 + max_grp_num) %A collective move.
%             %disp("collective")
%             collective_count = collective_count + 1;
%             Cluster_elements = transpose(nonzeros(CTM(:, Qclass - 14)));
%-----------------------------------------------------------------------

%------------------------------------------------
%Only these two lines are commented out for the no-cluster movement version
%that still accounts for the time increment of potential cluster moves but
%doesn't actually move anything.

%             [C, site, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = Collective_Shifter(C, Cluster_elements, site, im_i, ip_j, ip_i, im_j, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, N_particles, Li, CTM_site);
%             [CTM, CTM_site, max_grp_num, remaining_cluster_elements] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
%------------------------------------------------

            
        %elseif (Qclass > 14 + max_grp_num)  &  (Qclass <= 14 + max_grp_num + 4)           %The case of particle exit from tube.

        
%commented for closed system        
%         elseif (Qclass > 14 )  &  (Qclass <= 14 + 4)     %Use this elseif statement for the case where cluster is completely stationary
%             disp("EXIT !")
% 
%             
%             exitted= exitted + 1;
%             %exit_class = Qclass - (14 + max_grp_num);
%             exit_class = Qclass - (14);  %use this if stationary cluster.
%             exitters = nonzeros(exitters_table(exit_class, :));
%           
%             
%             if length(exitters) > 1
%                 particle = randsample(exitters,1);
%             elseif length(exitters) == 1
%                 particle = exitters;
%                
%             end
%             
%             
%             [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = particle_exit(particle, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j);
%===============================================================            
        %elseif (Qclass > 14 + max_grp_num + 4) & (Qclass <=  14 + max_grp_num + 4 + 1)     %The case of particle addition to tube.                                                                         
        %Below is commented to make it a closed system, for the shape
        %transition simulations.===================================
%         elseif (Qclass > 14  + 4) & (Qclass <=  14  + 4 + 1)   %Use this elseif statement for the case where cluster is completely stationary.
%             added = added + 1;
%             location = choose_location(C, Li, Lj);
%             
%             [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = add_particle(location, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j);
          %============================================================    
%         elseif Qclass > 14 + max_grp_num + 4 + 1    
%             %disp('collective exit');
%             cluster_index = Qclass - (14 + max_grp_num + 4 + 1);
%             exit_cluster_label = partial_exit_clusters(cluster_index);
%             c_exit = c_exit + length(nonzeros(CTM(:, exit_cluster_label)));
%             [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = complete_cluster_exit(exit_cluster_label, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j, Lj, prep);
          end
        
        %This is to check if the information for each class in E_long/lati
        %is the ssame as the information in LOC (except for class 0). It
        %only checks if the number of bonds in each class are the same.
        
%         for checlass = 1:7
%             if ((CEP_lati(checlass) == 0) & (LOC_LATI(checlass, 1) ~= 0))  |  ((CEP_lati(checlass) ~= 0 & LOC_LATI(checlass, 1) == 0))
%                 C
%                 LOC_LATI
%                 CEP_lati
%                 return
%             elseif ((CEP_long(checlass) == 0) & (LOC_LONG(checlass, 1))) ~= 0  |  ((CEP_long(checlass) ~= 0) & (LOC_LONG(checlass, 1) == 0))
%                 C
%                 LOC_LONG
%                 CEP_long
%                 return
%             end
%         end
%             
%                 
%         for glass = -3:1:3
%             if glass ~= 0
%                 if length(nonzeros(LOC_LATI(glass + 4, :))) ~= CEP_lati(glass + 4)
%                     C
%                     E_lati
%                     LOC_LATI
%                     CEP_lati
%                     return
%                 end
% 
%                 if length(nonzeros(LOC_LONG(glass + 4, :))) ~= CEP_long(glass + 4)
%                     C
%                     E_long
%                     LOC_LONG
%                     CEP_long
%                     return
%                 end
%             end
%         end        
%         
%         for glass = -3:1:3
%             if glass ~= 0
%                 if length(nonzeros(LOC_LATI(glass + 4, :))) ~= sum(sum(E_lati == glass))
%                     E_lati
%                     LOC_LATI
%                     C
%                     CEP_lati
%                     return
%                 end
% 
%                 if length(nonzeros(LOC_LONG(glass + 4, :))) ~= sum(sum(E_long == glass))
%                     E_long
%                     LOC_LONG
%                     C
%                     CEP_long
%                     return
%                 end
%             end
%         end
%                 
%         
%         len = size(CTM,2);
%         for ai = 1:max_grp_num
%             if CTM(:,ai) == zeros(N_particles,1)
%                 kMC
%                 return
%             end
%         end
%         
% 
%         if nnz(C) ~= N_particles % | (length(All_cluster_elements) ~= nnz(CTM))
%             found = true;
% 
%              C
% 
%              site
% 
%              CTM(1:N_particles,1:max_grp_num)
% 
%              Qclass
%              temp_1
%              temp_2
% %             
%              Cluster_elements
% %             remaining_cluster_elements
%              max_grp_num
% 
%             
%         end

        %=============================================================
        % if temp_2 ~= 0
        %     %update the y-value of the exchanged protein.
        %     site(temp_2, 2) = j;
        % end


        if floor(kMC/measure_interval) == ceil(kMC/measure_interval)

            
            
            site_dim = size(site,1);
            N_particles;
            if N_particles ~= 0
                
                number_wrap_clusters = wrapping_status(CTM, max_grp_num, site, Li);
                wrap_record(end + 1) = number_wrap_clusters;

                DC = nbr_distribution(site_dim, site, C, im_i, ip_j, ip_i, im_j);
                weight_record(end + 1) = weight_avg_mean_cluster_size(DC, N_particles);
                time_record(end + 1)  = run_time;
                N_record(end + 1) = N_particles;
                AIS = Average_Island_Size(DC);
                Clus_avg_size(end + 1) = AIS('avg_cluster');
                nanoClus_avg_size(end + 1) = AIS('avg_nanocluster');
                Olig_avg_size(end + 1) = AIS('avg_oligomer');
                Dim_avg_size(end + 1) = AIS('avg_dimer');

                Clus_num(end + 1) = AIS('num_cluster');
                nanoClus_num(end + 1) = AIS('num_nanocluster');
                Olig_num(end + 1) = AIS('num_oligomer');  
                Dime_num(end + 1) = AIS('num_dimer');
                Mono_num(end + 1) = AIS('num_monomer');
                %fprintf('The kMC number is %i %i %i %i %i\n', kMC,  run_time, weight_record(end), Olig_num(end), N_record(end));
            elseif N_particles == 0
                disp("No particles")
                fprintf('The kMC number is %i %i \n', kMC,  run_time);
            end
        end
        
%         %--------------------
%         %check this rare issue thing!
%         if mod(kMC_test, 10000) == 0
%             kMC_test
%             fprintf('The kMC information is %i %i %i %i %i \n', kMC,  N_particles, size(site,1), previous_analysis_kMC + analysis_interval, final_terminate);
%         end
%         %--------------------
        
        
    %    if (N_particles ~= 0) & ((kMC == kMC_0) | (kMC >= previous_analysis_kMC + analysis_interval) | kMC >= final_terminate | N_particles == 0)
        if ((kMC == kMC_0) | (kMC >= previous_analysis_kMC + analysis_interval) | kMC >= final_terminate | N_particles == 0)
   
            %The "kMC >= previous_analysis_kMC + analysis_interval" part is
            %what allows it to get a particle if N_particles == 0, then
            %update to the next analysis interval.
            previous_analysis_interval = analysis_interval;
            analysis_interval = min(floor(interval_scale*analysis_interval), max_analysis_interval);
            if kMC == kMC_0
                analysis_interval = initial_analysis_interval;
            end
            previous_analysis_kMC = kMC;
            site_dim = size(site,1);
            N_particles = site_dim;
            
            fprintf('The run time is %i \n', run_time);
            fprintf('The kMCS number is %i %i \n', J_hat, kMC);
            fprintf('Number of particles is %i  \n', N_particles);
            %fprintf('The kMC number is %i %i\n', Olig_num(end), N_record(end));
            

            if N_particles > 0
                DC = nbr_distribution(site_dim, site, C, im_i, ip_j, ip_i, im_j);
                %Record everything-----------------------------------------------------
                weight_record(end + 1) = weight_avg_mean_cluster_size(DC, N_particles);
                time_record(end + 1)  = run_time;
                N_record(end + 1) = N_particles;
    
                number_wrap_clusters = wrapping_status(CTM, max_grp_num, site, Li);
                wrap_record(end + 1) = number_wrap_clusters;
                
                AIS = Average_Island_Size(DC);
                Clus_avg_size(end + 1) = AIS('avg_cluster');
                nanoClus_avg_size(end + 1) = AIS('avg_nanocluster');
                Olig_avg_size(end + 1) = AIS('avg_oligomer');
                Dim_avg_size(end + 1) = AIS('avg_dimer');
    
                Clus_num(end + 1) = AIS('num_cluster');
                nanoClus_num(end + 1) = AIS('num_nanocluster');
                Olig_num(end + 1) = AIS('num_oligomer');  
                Dime_num(end + 1) = AIS('num_dimer');
                Mono_num(end + 1) = AIS('num_monomer');
            end
            %---------------------------------------------------------------------
            
%             All_cluster_elements = [];
%             for u = 1:length(DC)
%                 if length(DC{u}) > 1
%                     All_cluster_elements = [All_cluster_elements,DC{u}];
%                 end
%             end
%             All_cluster_elements;
%             if nnz(C) ~= N_particles | (length(All_cluster_elements) ~= nnz(CTM))
%                 found = true;
%                  C
%                  site
%                  CTM(1:N_particles,1:max_grp_num)
%                  Qclass
%                  temp_1
%                  temp_2
%                  Cluster_elements
%                  max_grp_num
%             end
            
            %Phase Analysis  -------------------------
%             FP_i = Fraction_of_phases(DC);
%             FPC{end + 1} = FP_i;    %The history of the Fraction of Phases at every analysis.
%             keys(FP_i);
%             values(FP_i);
%             DoP_i = Dominant_Phase(FP_i);
%-----------------------------------------------------          
            
            if snapshots == true
                x = [];
                y = [];

                %Graph
                for i =1:Li
                    for j = 1:Lj
                        if C(i, j) ~= 0
                            x = [x, j];
                            y = [y, (Li + 1) - i];    %making the graph reflect matrix positions.
                        end
                    end
                end
                f1 = figure()

                scatter(x,y, 40,'blue','filled');
                xlim([0 Lj]) ;
                ylim([0 Li]);
                shg
            end
            

            
            
%             if DoP_i == "cluster"         %Don't scale up the interval if you got "cluster" state.
%                 analysis_interval = previous_analysis_interval;
%             end
%-------------------------------            
%             if kMC > kMC_0
%                 DoP_termination_record  =  [DoP_termination_record, DoP_i];
%             end
%             
%             if length(DoP_termination_record) > termination_length
%                 DoP_termination_record(1) = [];      %Keeping only the record of the past five DoP analyses
%             end
%             DPE = DoP_i;       %For now, this is the DPE.
%-------------------------------          
            %percentage_weight_difference < 0.05 |
%             if length(weight_record) >= 10  & (length(DoP_termination_record) > 2)
%             
%                 percentage_weight_difference = abs((sqrt(weight_record(end)) - sqrt(weight_record(end - 1)))/sqrt(weight_record(end)))
%                 if  percentage_weight_difference < 0.3 & (sum(DoP_termination_record == DoP_i) == termination_length) | (sum(DoP_termination_record(end-1:end) == "cluster") == 200) | kMC == final_terminate   %The Cluster should be "cluster") == 2 if Cluster is all you want
%                     continue_monte = false
%                 end
%             end
            %if (length(DoP_termination_record) > 2)| (N_particles == 0)
%                 if (DoP_termination_record(end) ~= "monomer") & (DoP_termination_record(end) ~= "cluster") & (termination == false) & ((kMC > final_terminate) | (toc > terminal_computation_time)) & (run_time >= max_run_time) %Only when not cluster, and it is the first time the algorithm is trying to stop. This tries to checks if the 'oligomer' is stable or not
%                     previous_analysis_kMC = kMC;
%                     analysis_interval = round(2*exp(J_hat));    %Run it for twice the average lifetime of a dimer at this interaction strength.
%                     termination = true;                         %This is so it recognizes that it has already done this final check.
                    
            if   (kMC > final_terminate) | (toc > terminal_computation_time) | (N_particles == 0) %(run_time >= max_run_time) %This is just for testing perposes
                 continue_monte = false
            end
            %end
        end
        
        C;
    end

   toc 

    x = [];
    y = [];

    %Graph
    for i =1:Li
        for j = 1:Lj
            if C(i, j) ~= 0
                x = [x, j];
                y = [y, (Li + 1) - i];    %making the graph reflect matrix positions.
            end
        end
    end
    %f = figure()
    
   
    
    kMC
    CTM(:,1:max_grp_num)

    
    
    CTM;
    CTM_site;
    site;
    collective_count;
    N_record;
    C;
    %scatter(x,y, 40,'blue','filled');
    %xlim([0 Lj]) ;
    %ylim([0 Li]);
    %shg      %   temporarily not show phase diagram
    if N_particles > 0
        DC = nbr_distribution(N_particles, site, C, im_i, ip_j, ip_i, im_j);
        AIS = Average_Island_Size(DC);
                    
        
        
        FP_i = Fraction_of_phases(DC);
        FPC{end + 1} = FP_i;    %The history of the Fraction of Phases at every analysis.
        keys(FP_i);
        values(FP_i);
        DoP_i = Dominant_Phase(FP_i)
        DoP_termination_record = [DoP_termination_record, DoP_i];
        % if length(DoP_termination_record) > termination_length
        %     DoP_termination_record(1) = [];      %Keeping only the record of the past five DoP analyses
        % end
        DPE = DoP_i;       %For now, this is the DPE.
        
    elseif N_particles == 0
        DPE = "monomer"
    end
    island_info = {Clus_avg_size, nanoClus_avg_size, Olig_avg_size, Dim_avg_size, Clus_num, nanoClus_num, Olig_num, Dime_num, Mono_num};
    
    added
    exitted
    c_exit
    J_hat
   
    kMC_iterations = kMC;                   %Total number of kMC iterations for this run. Fomr this, and the number of measurements, you can infer the measure_interval.
end






%==========================================================================================================
%Below this line is the important stuff for kMC.



%Initialise matrices that contains bond information------------------
function [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = prepare_interaction(E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, Li, Lj, CTM_site)
    if exist('CTM_site', 'var')    
        %make E_lati
        for i = 1:Li
            for j = 1:(Lj-1)
                interact = [i, j; i, j + 1];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            end
        end

        %make it for E_long
        for i = 1:Li
            for j = 1:Lj
                interact = [i, j; ip_i(i), j];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            end
        end

        %make E_lati
        for i = 1:Li
            for j = 1:(Lj-1)
                interact = [i, j; i, j + 1];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            end
        end

        %make it for E_long
        for i = 1:Li
            for j = 1:Lj
                interact = [i, j; ip_i(i), j];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            end
        end
    else
        disp("prep 1")
        %make E_lati
        for i = 1:Li
            for j = 1:(Lj-1)
                interact = [i, j; i, j + 1];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles);
            end
        end
        disp("prep 2")
        %make it for E_long
        for i = 1:Li
            for j = 1:Lj
                interact = [i, j; ip_i(i), j];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles);
            end
        end
        disp("prep 3")
        %make E_lati
        for i = 1:Li
            for j = 1:(Lj-1)
                interact = [i, j; i, j + 1];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles);
            end
        end
        disp("prep 4")
        %make it for E_long
        for i = 1:Li
            for j = 1:Lj
                interact = [i, j; ip_i(i), j];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles);
            end
        end
        disp("end prep 4")
    end
end

%--------------------------------------------------------------------


%--------------------------------------------------------------------
%below are functions for adding particles.

function location = choose_location(C, Li, Lj)
    Done = false;
    while Done == false
        position = randi(Li*2);
        if position <= Li
            if C(position, 1) == 0
                location = [position, 1];
                Done = true;
            end
        else
            if C(position - Li, Lj) == 0
                location = [position - Li, Lj];
                Done = true;
            end
        end
    end
end


function chances =  chances_of_addition(C, Li, Lj, N_particles, rho)
    fringe_particles = [transpose(nonzeros(C(:,1))), transpose(nonzeros(C(:,end)))];
    n_edgesites_unoccupied = Li*2 - length(fringe_particles);
    %The chances of an arbitrary unoccuoied site (not just at the edges) in
    %C being filled is P.
    P = rho/4;
    chances = n_edgesites_unoccupied*P;
end


function [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = add_particle(location, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j)
    %This function adds a particle labelled 'particle' and
    %appropriately updates all other lists of information.
    Lj = size(C,2);
    %Add the particle to C, A.
    %location is the 2-element [i,j] vector that give the edge location of
    %C and A that the new particle will be added to.
    new_particle = N_particles + 1;
 
    LL = location;
    
    C(LL(1), LL(2)) = new_particle;
    site(new_particle, :) = [LL(1), LL(2)];
    CTM_site(new_particle, :) = [0, 0];  
    %Update the size of CTM.
    N_particles = N_particles + 1;
    num_column = 0;
    for i = 2:(floor(N_particles/2) + 1)
        if i < floor(N_particles/2) + 1
            num_column = num_column + floor(N_particles/i);
        elseif i == floor(N_particles/2) + 1
            num_column = num_column + (N_particles - floor(N_particles/2));
        end    
    end
    previous_num_column = size(CTM, 2);
    CTM(N_particles, :) = zeros(1, previous_num_column);  %Add a bottom row.
    n_add_columns = num_column - previous_num_column;   %number of new columns to add.
    CTM(:, previous_num_column + 1: num_column) = zeros(N_particles, n_add_columns);
    
    previous_LOC_length = size(LOC_LONG,2);
    if 4*N_particles > previous_LOC_length
        LOC_LATI = [LOC_LATI, zeros(7, 4*N_particles - previous_LOC_length)];
        LOC_LONG = [LOC_LONG, zeros(7, 4*N_particles - previous_LOC_length)];
        
        LOC_length = size(LOC_LONG,2);
        
        
        for i = 1:7
            bonds_in_LOC_row = nonzeros(LOC_LATI(i,:));
            for bond_index = 1:length(bonds_in_LOC_row)
                bond = bonds_in_LOC_row(bond_index);
                LOOK_LATI(bond) = LOC_number([i, bond_index], LOC_length) ;
            end
            bonds_in_LOC_row = nonzeros(LOC_LONG(i,:));
            for bond_index = 1:length(bonds_in_LOC_row)
                bond = bonds_in_LOC_row(bond_index);
                LOOK_LONG(bond) = LOC_number([i, bond_index], LOC_length) ;
            end
        end
    end
    
    class = 0 ;   %The class is zero.
    
    
    if LL(2) == 1
%         for the side bond, up bond, and bottom bond if applicable
%         check LOOK and if they are in LOC, then you've gotta remove it.
%         Then you also need to ' shorten the entire LOC matrix by four
%         columns since the length of LOC is given by 4*N_particles.
%         Also delete the LOOK value.
%         For every deletion, decrease CEP by 1 for given class.
        
        ignore = false;
        previous_class = E_lati(LL(1), 1);
        bond_num = bond_number([LL(1), 1], size(E_lati,1));
        initial_LOOK = LOOK_LATI(bond_num);
        E_lati(LL(1), 1) = 0;
        if (initial_LOOK == 0 & class == previous_class) | (class ~= previous_class)
        [CEP_lati, LOC_LATI, LOOK_LATI] = add_update(bond_num, class, N_particles, CEP_lati, LOC_LATI, LOOK_LATI);
        end
        if initial_LOOK ~= 0 & (class ~= previous_class)
            [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
        end
        
        previous_class = E_long(LL(1), 1);
        bond_num = bond_number([LL(1), 1], size(E_long,1));
        initial_LOOK = LOOK_LONG(bond_num);
        E_long(LL(1), 1) = 0;
        if (initial_LOOK == 0 & class == previous_class) | (class ~= previous_class)
        [CEP_long, LOC_LONG, LOOK_LONG] = add_update(bond_num, class, N_particles, CEP_long, LOC_LONG, LOOK_LONG);
        end
        if initial_LOOK ~= 0 & (class ~= previous_class)
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
        
        
        previous_class = E_long(im_i(LL(1)), 1);
        bond_num = bond_number([im_i(LL(1)), 1], size(E_long,1));
        initial_LOOK = LOOK_LONG(bond_num);
        E_long(im_i(LL(1)), 1) = 0;
        if (initial_LOOK == 0 & class == previous_class) | (class ~= previous_class)
        [CEP_long, LOC_LONG, LOOK_LONG] = add_update(bond_num, class, N_particles, CEP_long, LOC_LONG, LOOK_LONG);
        end
        if initial_LOOK ~= 0 & (class ~= previous_class)
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
        
        
        Immediate_Neighbors = [im_i(LL(1)), 1;    LL(1), 2;    ip_i(LL(1)), 1];
        
    elseif LL(2) == Lj
    
        ignore = false;
        previous_class = E_lati(LL(1), Lj - 1);
        bond_num = bond_number([LL(1), Lj - 1], size(E_lati,1));
        initial_LOOK = LOOK_LATI(bond_num);
        E_lati(LL(1), Lj - 1) = 0;
        if (initial_LOOK == 0 & class == previous_class) | (class ~= previous_class)
        [CEP_lati, LOC_LATI, LOOK_LATI] = add_update(bond_num, class, N_particles, CEP_lati, LOC_LATI, LOOK_LATI);
        end
        if initial_LOOK ~= 0 & (class ~= previous_class)
            [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
        end
        
        previous_class = E_long(LL(1), Lj);
        bond_num = bond_number([LL(1), Lj], size(E_long,1));
        initial_LOOK = LOOK_LONG(bond_num);
        E_long(LL(1), Lj) = 0;
        if (initial_LOOK == 0 & class == previous_class) | (class ~= previous_class)
        [CEP_long, LOC_LONG, LOOK_LONG] = add_update(bond_num, class, N_particles, CEP_long, LOC_LONG, LOOK_LONG);
        end
        if initial_LOOK ~= 0 & (class ~= previous_class)
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
        
        
        previous_class = E_long(im_i(LL(1)), Lj);
        bond_num = bond_number([im_i(LL(1)), Lj], size(E_long,1));
        initial_LOOK = LOOK_LONG(bond_num);
        E_long(im_i(LL(1)), Lj) = 0;
        if (initial_LOOK == 0 & class == previous_class) | (class ~= previous_class)
        [CEP_long, LOC_LONG, LOOK_LONG] = add_update(bond_num, class, N_particles, CEP_long, LOC_LONG, LOOK_LONG);
        end
        if initial_LOOK ~= 0 & (class ~= previous_class)
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
        
        Immediate_Neighbors = [im_i(LL(1)), Lj;    LL(1), Lj - 1;    ip_i(LL(1)), Lj];
        
    end


    for x = 1:3
        INi = Immediate_Neighbors(x,1);
        INj = Immediate_Neighbors(x,2);
        neighbors_ofIN = neighbors_generic([INi, INj], site, C,  im_i, ip_j, ip_i, im_j);
        for no = neighbors_ofIN
            interact = [site(no, 1), site(no, 2); INi, INj];
            [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
        end
    end

    
    Cluster_elements = neighbors_generic(LL, site, C, im_i, ip_j, ip_i, im_j);
    if isempty(Cluster_elements) == false
        Cluster_elements(end + 1) = new_particle;
        %Is the 'prep = true' term needed here ?????
        prep = false;
        [CTM, CTM_site, max_grp_num, remaining_cluster_elements] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
        
    end
    

    
end


%--------------------------------------------------------------------



%--------------------------------------------------------------------
%below are the functions for partial cluster exit.
function partial_exit_clusters = clusters_at_edges(CTM_site, C, CTM, exit_cutoff)
%This function returns the cluster that are touching the edges of the
%lattice tube, and thus have the potential to move outside the tube
%boundaries and loose some particles.
    exitters = [transpose(nonzeros(C(:,1))), transpose(nonzeros(C(:,end)))];
    cluster_labels = [];
    
    if isempty(exitters) == false
        for ex = 1:length(exitters)
            CTM_location = CTM_site(exitters(ex),:);
        
            if (CTM_location(1) ~= 0) & (CTM_location(2) ~= 0)  &   (length(nonzeros(CTM(:, CTM_location(2)))) < exit_cutoff)  %This last condition was added Jan 30th to keep larger clusters in the tube.
                
                cluster_labels(end + 1) = CTM_location(2);
            end
        end
    end
    partial_exit_clusters = unique(cluster_labels);
end


function cluster_width = cluster_width_calculator(Cluster_index, CTM, site, Lj)
%calculates the latitudinal width of a Cluster. if the cluster is longer
%than half of the tube length, than the move will be rejected, and the
%cluster_width will be returned as 'false', so that a move won't happen.
    Cluster_Elements = transpose(nonzeros(CTM(:,Cluster_index)));
    width = 1;
    for elem = Cluster_Elements
        lati_position = site(elem, 2);
        distance_to_edge = min(lati_position, Lj - lati_position);
        width = max(width, distance_to_edge);
        if width == round(Lj/2) - 1
            width = false;
            break
        end
    end
    cluster_width = width;
end



function [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = complete_cluster_exit(cluster_label, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j, Lj, prep)
    %Completely gets rid of the selected cluster. 
    Cluster_Elements = sort(transpose(nonzeros(CTM(:, cluster_label))));
    
    num_el = length(Cluster_Elements);
    for i = num_el:-1:1
        subparticle = Cluster_Elements(i);
        [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = particle_exit(subparticle, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j);
    end
    
end



function [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = cluster_exit(cluster_label, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j, Lj, prep)
%This function allows part of a cluster to exit when it crosses the border
%of the lattice.
    exitters = [transpose(nonzeros(C(:,1))), transpose(nonzeros(C(:,end)))];
    exit_particles = [];
    Cluster_elements = transpose(nonzeros(CTM(:, cluster_label)));
    
    %This is to determine the direction of motion. There is a small chance
    %that the cluster spans the whole tube and that cluster elements are
    %included at both edges of the tube. Shift to the side with more edge
    %particels (a bit of an arbitrary rule). 
    fringe_particles_sites = [];   %This list can include particles from both sides.
    left_or_right = [0,0];
    for ex = length(exitters):-1:1
        CTM_location = CTM_site(exitters(ex),:);
        
        if (CTM_location(2) == cluster_label)
            particle_label = exitters(ex);
            if site(particle_label,2) == 1
                left_or_right(1) = left_or_right(1) + 1;
            elseif site(particle_label, 2) == Lj
                left_or_right(2) = left_or_right(2) + 1;
            end
            fringe_particles_sites(end + 1, :) = site(particle_label, :);
        end
    end
    
    %Choose direction
    if left_or_right(1) > left_or_right(2)
        direction = 4;
        edge = 1;
    elseif left_or_right(1) < left_or_right(2)
        direction = 2;
        edge = Lj;
    elseif left_or_right(1) == left_or_right(2)
        direction = randsample([2,4],1);
        if direction == 4
            edge = 1;
        elseif direction == 2
            edge = Lj;
        end
    end
    
    %After the single particle exits, the elements in the list
    %'Cluster_elements' will be labelled differently, so I need to keep a
    %record of the locations of the Cluster_element site.
    Cluster_elements_site = zeros(length(Cluster_elements), 2);
    for index_element = length(Cluster_elements):-1:1
        particle_label = Cluster_elements(index_element);
        if site(particle_label,2) == edge
            Cluster_elements_site(index_element, :) = [];
        else
            Cluster_elements_site(index_element, :) = site(particle_label,:);
        end
    end
    
    %Delete the edge particles and erase then from the Cluster_elements.
    for ex = 1:size(fringe_particles_sites, 1)
       
        location = fringe_particles_sites(ex, :);
        if (location(2) == edge)
            particle_label = C(location(1), location(2));
            [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = particle_exit(particle_label, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j);
            
            
        end
    end
    
    %Remake the CLuster element list.
    N_remaining = size(Cluster_elements_site,1);
    Cluster_elements = [];
    for C_index = 1:N_remaining
        Cluster_elements(end + 1) = C(Cluster_elements_site(C_index,1),Cluster_elements_site(C_index,2));  %Retrive the update particle label at that site.
    end
    
    remaining_clusters = [];
    remaining_singles = [];
    %The remaining cluster after deleting the edge particles may be
    %fragmented, so we need to find all the current labels of those
    %fragmented clusters.
    
    for shifting_element = Cluster_elements
        CTM_location = CTM_site(shifting_element,:);
        if CTM_location ~= [0,0]
            remaining_clusters(end + 1) = CTM_location(2);   %The Cluster number the particle belongs to.
        elseif CTM_location == [0,0]   %The case of a single particle.
            remaining_singles(end + 1) = shifting_element;
        end
    end
    remaining_clusters = unique(remaining_clusters);
    
    %For each remaining_cluster cluster, shift it to the its respective
    %edge direction. (left or right).
    if length(remaining_clusters) > 0
        for cluster = remaining_clusters
            
            Cluster_elements = transpose(nonzeros(CTM(:, cluster)));
            [C, site, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = Collective_Shifter(C, Cluster_elements, site, im_i, ip_j, ip_i, im_j, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, N_particles, Lj, direction);
           
            [CTM, CTM_site, max_grp_num, remaining_cluster_elements] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
        end
    end
    
    %For each remaining single particle, shift them to the specified
    %direction.
    if length(remaining_singles) > 0
        for single_particle = remaining_singles
            %a horizontal interaction, so the row is the same. edge give
            %the empty site beside it in the appropriate direction.
            interact = [site(single_particle,1), site(single_particle,2); site(single_particle,1), edge];
            %---
%             temp_1 = A(interact(1,1),interact(1,2));
%             temp_2 = A(interact(2,1), interact(2,2));
%             A(interact(1,1),interact(1,2)) = temp_2;
%             A(interact(2,1), interact(2,2)) = temp_1;


            % %update the position in the site of the ONE values. +1  to the
            % %right.
            % site(random, 2) = ip_j(j); 

            %Update positions on the C matrix
            temp_1 = C(interact(1,1),interact(1,2));
            temp_2 = C(interact(2,1), interact(2,2));
            C(interact(1,1),interact(1,2)) = temp_2;
            C(interact(2,1), interact(2,2)) = temp_1;

            if temp_1 ~= 0 
                site(temp_1,:) = [interact(2,1), interact(2,2)];
            end
            if temp_2 ~= 0
                site(temp_2,:) = [interact(1,1), interact(1,2)];
            end    
            %fprintf('From here')
            [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);


            %CTM update for single moves==================================
            if sum(0 == [temp_1, temp_2]) == 1 %If either one of the interact sites is a particle and the other is a vacancy (it can't be both particles).    
                Cluster_elements = nonzeros([temp_1, temp_2]);

                [CTM, CTM_site, max_grp_num, remaining_cluster_elements] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
                %remaining_cluster_elements = nonzeros(CTM(:,1));
                
                if length(remaining_cluster_elements) > 0
                    prep = true;
                    [CTM, CTM_site, max_grp_num] = CTM_updater_2(CTM, CTM_site, max_grp_num, interact, remaining_cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
                    prep = false;
                end
 
            end
            %----
        end
    end


end


%--------------------------------------------------------------------





%---------------------------------------------------------------------
%The next three function are for deleteing particles.
function [N_exitter_list, exitter_table] = particles_at_edges(CTM_site, C, site, im_i, ip_j, ip_i, im_j)
    exitters = [transpose(nonzeros(C(:,1))), transpose(nonzeros(C(:,end)))];
    N_exit = length(exitters);
    exitter_table = zeros(4, N_exit);     
    class_1 = [];
    class_2 = [];
    class_3 = [];
    %The first row in exitter_table is for particles exitting that don't change 
    %energy. The ones on the second row has one neighbor, the ones one the
    %third row have two neighbors.. and so on for the fourth row.
    
    if isempty(exitters) == false
        for ex = length(exitters):-1:1
            CTM_location = CTM_site(exitters(ex),:);
        
            if (CTM_location(1) ~= 0) & (CTM_location(2) ~= 0)
                particle_label = exitters(ex);
                exitters(ex) = [];
                location = site(particle_label, :);
                nbrs = neighbors_generic(location, site, C, im_i, ip_j, ip_i, im_j);  
                N_nbrs = length(nbrs);  %This is either 1, 2 or 3.
                if N_nbrs == 1
                    class_1(end + 1) = particle_label;
                elseif N_nbrs == 2
                    class_2(end + 1) = particle_label;
                elseif N_nbrs == 3
                    class_3(end + 1) = particle_label;
                end
            end
        end
    end
    exitter_table(1, 1:length(exitters)) = exitters;
    exitter_table(2, 1:length(class_1))  = class_1;
    exitter_table(3, 1:length(class_2))  = class_2;
    exitter_table(4, 1:length(class_3))  = class_3;
    N_exitter_list = [length(exitters), length(class_1), length(class_2), length(class_3)];
end

function nbrs = neighbors_generic(xx, site, C, im_i, ip_j, ip_i, im_j)
    %For site xx, which is a vector like [i, j], determine the site index number of the neighbors. 
    %The position xx does not need to be occupied, but this function find
    %the occupied neighbors of the location xx, occupied or not.
    %site is a Nx2 array where N is the number of proteins (ones in the matrix 
    %used for the spin exchange. C is the LxL matrix (L is the lattice
    %dimensions) where the positions of 1's are replaced by their index in
    %the array site (the column index of site) 
    %look up table for mod
    iv = xx(1);
    jv = xx(2);
    x_site = [iv, jv];
    nbr_site = [iv, ip_j(jv);  ip_i(iv) , jv;  iv, im_j(jv);  im_i(iv), jv];
    
    %Remove those that are not included due to nf b.c. in the j-direction.
    for ii = 0:3
        if (nbr_site(4 - ii, 2) == 0) % |  (C(beside_site(4 - ii, 1), beside_site(4 - ii, 2)) == 0)
            nbr_site(4 - ii,:) = [];     %remove row that is not included due to nfbc.
        end
    end
    
    nbrs = [];
    for i = 1:length(nbr_site)
        nbr_index = C(nbr_site(i, 1), nbr_site(i, 2));
        if nbr_index ~= 0
            nbrs = [nbrs, nbr_index];
            nbrs;
        end
    end
end

function [N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num] = particle_exit(particle, N_particles, site, C, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, CTM, CTM_site, max_grp_num, ip_i, ip_j, im_i, im_j)
    %This function removes the particle labelled 'particle' and
    %appropriately updates all other lists of information.
    Lj = size(C,2);
    %Remove the particle from C, A and update all labells in C greater than
    %'particle'
    prep = false;
    if particle > size(site,1)  %This is to find the error.
        particle
        site
        
        nonzeros(C)
        CTM
        
    end
    site_of_p = site(particle, :);
    C(site_of_p(1), site_of_p(2)) = 0;
    
    
    if (CTM_site(particle, 1) ~= 0) & (CTM_site(particle, 2) ~= 0)
        Cluster_elements = [0, particle];
        [CTM, CTM_site, max_grp_num, remaining_clusters] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep);
        if length(remaining_clusters) > 0
        prep = true;
        interact = [0,0; site_of_p(1), site_of_p(2)];
        [CTM, CTM_site, max_grp_num] = CTM_updater_2(CTM, CTM_site, max_grp_num, interact, remaining_clusters, site, C, im_i, ip_j, ip_i, im_j, prep);
        prep = false;
        end
    end
    
    N = N_particles;
    N_change_label = N - particle;
    for n = particle + 1:N
        ni = site(n,1);
        nj = site(n,2);
        C(ni, nj) = C(ni, nj) - 1;
        ctm_i = CTM_site(n,1);
        ctm_j = CTM_site(n,2);
        if (ctm_i ~= 0) & (ctm_j ~= 0)
            CTM(ctm_i, ctm_j) = CTM(ctm_i, ctm_j) - 1;
        end
    end
    %Delete the the row in site. 
    %CTM_site delete the row for the particle.
    site(particle,:) = [];
    CTM_site(particle,:) = [];
    N_particles = N_particles - 1;
    num_column = 0;
    for i = 2:(floor(N_particles/2) + 1)
        if i < floor(N_particles/2) + 1
            num_column = num_column + floor(N_particles/i);
        elseif i == floor(N_particles/2) + 1
            num_column = num_column + (N_particles - floor(N_particles/2));
        end    
    end
    %The template of CTM matrix is: CTM = zeros(N_particles, num_column);
    %To adjust the size of CTM, the following is done.
    CTM(N_particles + 1, :) = [];    %Remove the lowest row.
    previous_num_column = size(CTM, 2);
    %Delete right-most columns of CTM.
    CTM(:, num_column + 1: previous_num_column) = [];
    
    CEP_long;
    LOC_LONG;
    
    class = 0;

    if site_of_p(2) == 1
%         for the side bond, up bond, and bottom bond if applicable
%         check LOOK and if they are in LOC, then you've gotta remove it.
%         Then you also need to ' shorten the entire LOC matrix by four
%         columns since the length of LOC is given by 4*N_particles.
%         Also delete the LOOK value.
%         For every deletion, decrease CEP by 1 for given class.
        ignore = true;
        previous_class = E_lati(site_of_p(1), 1);
        bond_num = bond_number([site_of_p(1), 1], size(E_lati,1));
        E_lati(site_of_p(1), 1) = 0;
        initial_LOOK = LOOK_LATI(bond_num);
        if initial_LOOK ~= 0 
            [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
        end
            
        previous_class = E_long(site_of_p(1), 1);
        bond_num = bond_number([site_of_p(1), 1], size(E_long,1));
        E_long(site_of_p(1), 1) = 0;
        initial_LOOK = LOOK_LONG(bond_num);
        if initial_LOOK ~= 0 
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
            
        previous_class = E_long(im_i(site_of_p(1)), 1);
        bond_num = bond_number([im_i(site_of_p(1)), 1], size(E_long,1));
        E_long(im_i(site_of_p(1)), 1) = 0;
        initial_LOOK = LOOK_LONG(bond_num);
        if initial_LOOK ~= 0 
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
            
        Immediate_Neighbors = [im_i(site_of_p(1)), 1;    site_of_p(1), 2;    ip_i(site_of_p(1)), 1];
        
        for x = 1:3
            INi = Immediate_Neighbors(x,1);
            INj = Immediate_Neighbors(x,2);
            neighbors_ofIN = neighbors_generic([INi, INj], site, C,  im_i, ip_j, ip_i, im_j);
            for no = neighbors_ofIN
                interact = [site(no, 1), site(no, 2); INi, INj];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            end
        end
        
    elseif site_of_p(2) == Lj
        ignore = true;
        previous_class = E_lati(site_of_p(1), Lj - 1);
        bond_num = bond_number([site_of_p(1), Lj - 1], size(E_lati,1));
        E_lati(site_of_p(1), Lj - 1) = 0;
        initial_LOOK = LOOK_LATI(bond_num);
        if initial_LOOK ~= 0 
            [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
        end
            
        previous_class = E_long(site_of_p(1), Lj);
        bond_num = bond_number([site_of_p(1), Lj], size(E_long,1));
        E_long(site_of_p(1), Lj) = 0;
        initial_LOOK = LOOK_LONG(bond_num);
        if initial_LOOK ~= 0 
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
            
        previous_class = E_long(im_i(site_of_p(1)), Lj);
        bond_num = bond_number([im_i(site_of_p(1)), Lj], size(E_long,1));
        E_long(im_i(site_of_p(1)), Lj) = 0;
        initial_LOOK = LOOK_LONG(bond_num);
        if initial_LOOK ~= 0 
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
            
        Immediate_Neighbors = [im_i(site_of_p(1)), Lj;    site_of_p(1), Lj - 1;    ip_i(site_of_p(1)), Lj];
        
        for x = 1:3
            INi = Immediate_Neighbors(x,1);
            INj = Immediate_Neighbors(x,2);
            neighbors_ofIN = neighbors_generic([INi, INj], site, C,  im_i, ip_j, ip_i, im_j);
            for no = neighbors_ofIN
                interact = [site(no, 1), site(no, 2); INi, INj];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            end
        end
    

        
    else   %Cases where the particle is not at the edge, but I'm deleting it if it belongs to a cluster that is exitting.
        
        ignore = true;
        previous_class = E_lati(site_of_p(1), im_j(site_of_p(2)));
        bond_num = bond_number([site_of_p(1), im_j(site_of_p(2))], size(E_lati,1));
        E_lati(site_of_p(1), im_j(site_of_p(2))) = 0;
        initial_LOOK = LOOK_LATI(bond_num);
        if initial_LOOK ~= 0 
            [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
        end
        
        ignore = true;
        previous_class = E_lati(site_of_p(1), site_of_p(2));
        bond_num = bond_number([site_of_p(1), site_of_p(2)], size(E_lati,1));
        E_lati(site_of_p(1), site_of_p(2)) = 0;
        initial_LOOK = LOOK_LATI(bond_num);
        if initial_LOOK ~= 0 
            [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
        end
            
        previous_class = E_long(site_of_p(1), site_of_p(2));
        bond_num = bond_number([site_of_p(1), site_of_p(2)], size(E_long,1));
        E_long(site_of_p(1), site_of_p(2)) = 0;
        initial_LOOK = LOOK_LONG(bond_num);
        if initial_LOOK ~= 0 
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
            
        previous_class = E_long(im_i(site_of_p(1)), site_of_p(2));
        bond_num = bond_number([im_i(site_of_p(1)), site_of_p(2)], size(E_long,1));
        E_long(im_i(site_of_p(1)), site_of_p(2)) = 0;
        initial_LOOK = LOOK_LONG(bond_num);
        if initial_LOOK ~= 0 
            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
        end
        
        
        Immediate_Neighbors = [im_i(site_of_p(1)), site_of_p(2);    site_of_p(1), im_j(site_of_p(2));    ip_i(site_of_p(1)), site_of_p(2); site_of_p(1), ip_j(site_of_p(2))];
        for x = 1:4
            INi = Immediate_Neighbors(x,1);
            INj = Immediate_Neighbors(x,2);
            neighbors_ofIN = neighbors_generic([INi, INj], site, C,  im_i, ip_j, ip_i, im_j);
            for no = neighbors_ofIN
                interact = [site(no, 1), site(no, 2); INi, INj];
                [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            end
        end
         
        
    end
end

%Functions for deleting particles end here.
%-----------------------------------------------------------









function [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site, cluster_move)
    
    Li = size(C, 1);
    Lj = size(C, 2);
    
    skip_long_updates = false;
    
    
    %Make a 'cheep way' to calculate the surrounding nerighborhood of the
    %particles to see if it is isolated -> would be able to ignore all the
    %other updates.
    int_maxj = max(interact(1,2), interact(2,2));
    int_minj = min(interact(1,2), interact(2,2));
    int_maxi = max(interact(1,1), interact(2,1));
    int_mini = min(interact(1,1), interact(2,1));
    if ~(  (((int_maxj + 2) > Lj) |  ((int_minj - 2) < 1)) |  (((int_maxi + 2) > Li) |  ((int_mini - 2) < 1))  )  %If not near the boundaries at the tube exits (also for the top and botton, temporary condition).
                                                                                                          %Hope to allow cases where I can calculate the neighborhood when it circles back over the top and bottom
                                                                                                   %but that might make this part more expensive to execute for every iteration.
        SUM = sum(sum(C((int_mini - 2):(int_maxi + 2), (int_minj - 2): (int_maxj + 2))));
        if int_maxi == int_mini   %laterral interaction
            particle_label = C(int_mini, int_maxj) + C(int_mini, int_minj);
            if SUM == particle_label
                skip_long_updates = true;
            end
        elseif int_maxj == int_minj   %longitudinal interaction
            particle_label = C(int_mini, int_minj) + C(int_maxi, int_minj);
            if SUM == particle_label
                skip_long_updates = true;
            end
        end
    end
                
        
    
    
    %The following gives the six occupied neighbors of some pair of sites that interacted.
    %Let interact(1,:) be the (i, j)_A corrdinates for particle A and interact(2,:) be
    %the (i, j)_B coordinates of site (partice) B.
    if (C(interact(1,1), interact(1,2)) == 0 & C(interact(2,1), interact(2,2)) == 0) 
                  %In this case the bonds don't change
        
    
    elseif ((C(interact(1,1), interact(1,2)) ~= 0) & (C(interact(2,1), interact(2,2)) ~= 0))
                  %In this case the bonds don't change
                  %This is the reason why we need to make the clusters
                  %shift and each particle must be updated in order!!!!
%          if  interact(1,1) == interact(2,1)  %This implise latitudinal interaction
%              lati_bond = min(interact(1,2), interact(2,2));
%              E_i = interact(1,1);
%              E_j = min(interact(1,2), interact(2,2));
%              E_lati(E_i, E_j) = 0;
%          elseif  interact(1,2) == interact(2,2)  %This implise longitudinal interaction
%             long_bond = min(interact(1,1), interact(2,1));
%             E_j = interact(1,2);
%             E_i = min(interact(1,1), interact(2,1));
%             E_long(E_i, E_j) = 0;
%          end
                  
                  
%     %Avoid update all 23 bonds if the particle is isolated.            
%     elseif skip_long_updates == true
%         if int_maxi == int_mini   %laterral interaction
%             if C(int_mini, int_maxj) ~= 0  %A move from left to right. 
%                 
%                 %Delete the latitudianl bond to the left
%                 E_i = int_mini;
%                 E_j = int_minj - 1;
%                 bond_num = bond_number([E_i, E_j], Li);   
%                 initial_LOOK = LOOK_LATI(bond_num);
%                 [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
%                 %Delete the two longitudinal bonds
%                 delete_bond_indices
%                 
%                 %Add update to two long and one lati bond.
%                 
%             %Repeat this for the case in the other direction of movement
%             
%         %Repeat this for the case of longitudinal interaction
        
        
        
    else
        close_six = [];
        if interact(1,1) == interact(2,1)               %This implise latitudinal interaction on the lattice.
            lati_bond = min(interact(1,2), interact(2,2));%whichever of the two with smaller j-th index determins which 
                                                          %hop (bond) it is assiciated to in the
                                                          %E_lati matrix.
            src = [interact(1,1), lati_bond]  ;          %src short for 'source'
            
            %Avoid update all 23 bonds if the particle is isolated. 
            if skip_long_updates == false
                close_six = [src(1), src(2); src(1), ip_j(src(2));
                             im_i(src(1)), src(2); im_i(src(1)), ip_j(src(2));          %start from above source, go clockwise.
                             src(1), ip_j(ip_j(src(2))); ip_i(src(1)), ip_j(src(2));
                             ip_i(src(1)), src(2); src(1), im_j(src(2))];
             elseif skip_long_updates == true
                 close_six = [src(1), src(2); src(1), ip_j(src(2))];
             end
        end
        if interact(1,2) == interact(2,2)
            long_bond = min(interact(1,1), interact(2,1));
            if (long_bond == 1) & ((interact(1,1) == Li) | (interact(2,1) == Li))
                        long_bond = Li;     %At the periodic boundary, the long_bond = min() argument doesn't work...
            end
            src = [long_bond, interact(1,2)];
            if skip_long_updates == false
                close_six = [src(1), src(2); ip_i(src(1)), src(2);
                             im_i(src(1)), src(2); src(1), ip_j(src(2));
                             ip_i(src(1)), ip_j(src(2)); ip_i(ip_i(src(1))), src(2);
                             ip_i(src(1)), im_j(src(2)); src(1), im_j(src(2))];
             elseif skip_long_updates == true   
                 close_six = [src(1), src(2); ip_i(src(1)), src(2)];
             end
        end
        
        
        if skip_long_updates == false
            for ii = 0:7
                if close_six(8 - ii, 2) == 0
                    close_six(8 - ii,:) = [];     %remove row that is not included due to nfbc.
                end
            end
        end

        num_close_six = size(close_six,1);     %The number of sites remaining in close_six.


        %close_six
        
        %site_dim(1)
        for n = 1:size(close_six,1)     %This for loop checks all bonds that are affected by the spin-exchange of the particles of the interacting bond.
            beside_site = [im_i(close_six(n,1)), close_six(n,2);         %The four neighbors of site n in the close_six list.
                           close_six(n,1), ip_j(close_six(n,2));
                           ip_i(close_six(n,1)), close_six(n,2); 
                           close_six(n,1), im_j(close_six(n,2))];    
                       
            occupancy_status = [C(close_six(n,1), close_six(n,2))];
            %Take out parts that are not included due to the no flux b.c.
            %(nfbc)
            for ii = 0:3
                if (beside_site(4 - ii, 2) == 0) % |  (C(beside_site(4 - ii, 1), beside_site(4 - ii, 2)) == 0)
                    beside_site(4 - ii,:) = [];     %remove row that is not included due to nfbc.
                
               
                elseif (n > 2) & (isequal(beside_site(4 - ii, :), interact(1,:)) |  isequal(beside_site(4 - ii, :), interact(2,:)) )
                    beside_site(4 - ii,:) = [];       %when n>2, the bonds of the interact sites are already update, so we can remove 
                                                      %the interact site from the 'beside_site' of this particular 'center_site'.
                else
                    occupancy_status(end + 1) = C(beside_site(4 - ii, 1), beside_site(4 - ii, 2));
                end
            end

            
            
            
            num_neighbors = size(beside_site, 1);
            %num_neighbors
            if sum(occupancy_status) ~= 0
            for m =  1:num_neighbors   %For each of the six neighbors of the two sites that exchanged particle occupancy, updating the bonds that were affected. 
                ignore = false;
                center_site = close_six(n,:);
                particle_label = nonzeros([C(center_site(1), center_site(2)),  C(beside_site(m,1), beside_site(m,2))]);  %The particle's label number. The center site is determined by the smaller (i,j) value, and maynot be the site that is occupied.
                center_neighbors = size(count_neighbors_nf(center_site, beside_site(m,:), C, Li, Lj, im_i, ip_j, ip_i, im_j), 1) ;   %Just the info of the number of neighbors
                beside_site_neighbors = size(count_neighbors_nf(beside_site(m,:), center_site, C, Li, Lj, im_i, ip_j, ip_i, im_j), 1);
                if (C(beside_site(m,1), beside_site(m,2)) ~= 0) & (C(center_site(1), center_site(2)) == 0)
                    class = -(center_neighbors - beside_site_neighbors);     %The negative sign is because H = -J*sigma_i*sigma_j.
                
                elseif (C(beside_site(m,1), beside_site(m,2)) == 0) & (C(center_site(1), center_site(2)) ~= 0)
                    class = -( beside_site_neighbors - center_neighbors);     %The negative sign is because H = -J*sigma_i*sigma_j.
                
                
                %'class' is basically (delta)H of the potential spin-exchange.
                elseif (C(beside_site(m,1), beside_site(m,2)) ~= 0) & (C(center_site(1), center_site(2)) ~= 0)
                    class = 0 ;   
                    ignore = true;
                
                
                elseif (C(beside_site(m,1), beside_site(m,2)) == 0) & (C(center_site(1), center_site(2)) == 0)
                    class = 0;
                    ignore = true;
                    
                end

%The following is I'm not sure why it works, but it was making particles that were surrounded by 4 neighbors have E_long,lati values of 2 or 3, messing
%up the whole energy landscape.
%                 if class == 0   %This is so that particles that belong to a cluster will not exchange easily, even if the energy of the system doesn't increase.
%                     if ~exist('cluster_move', 'var') & exist('CTM_site', 'var')
%                         
%                         if (CTM_site(particle_label, 1) ~= 0) & (CTM_site(particle_label, 2) ~= 0)   %Meaning it is in a cluster.
%                             class = center_neighbors     
%                         end
%                     elseif exist('cluster_move', 'var')  %If it is a cluster move, just base 'class' on your immediate neighbors.
% %                         neighbor_set = neighbors(particle_label, site, C, im_i, ip_j, ip_i, im_j);
% %                         neighbor
% %                         if C(center_site(1), center_site(2)) == particle_label;
% %                             xx =  [beside_site(m,1), beside_site(m,2)];       %This would be the unoccupied site beside the particle.
% %                         elseif  C(beside_site(m,1), beside_site(m,2)) == particle_label
% %                             xx = [center_site(1), center_site(2)];            %This would be the unoccupied site beside the particle.
% %                         end
% %                         exchange_neighbor_set = neighbors_generic(xx, site, C, im_i, ip_j, ip_i, im_j);
%                         
%                         
%                         %if   %Meaning even if you move, you still belong to the same cluster.
%                         if (CTM_site(particle_label, 1) ~= 0) & (CTM_site(particle_label, 2) ~= 0)   %Meaning it is in a cluster.
%                             class = center_neighbors;     
%                         end
% %                             %what about ignore?
% %                             
% %                         end
%                         
%                     end
%                 end
                
                
                if (center_site(1) == beside_site(m, 1))   %This implise latitudinal interaction on the lattice. 
                    %Below two blocks are the updating of a single bond, where the bond is between the 'center_site' and the 'beside_site'.
                    E_i = center_site(1);
                    E_j = min(center_site(2), beside_site(m, 2));
                    previous_class = E_lati(E_i, E_j);
                    
                    
                    E_lati(E_i, E_j) = class;                  % This should be the class of the interaction
                    bond_num = bond_number([E_i, E_j], Li);    %Li is the same as the row dimension of the E_lati matrix
                    initial_LOOK = LOOK_LATI(bond_num);
                    
                    %updating LOC_LATI
                    if (ignore == false)
                       if class == previous_class
                           if class == 0 & previous_class == 0 & initial_LOOK == 0
                               [CEP_lati, LOC_LATI, LOOK_LATI] = add_update(bond_num, class, N_particles, CEP_lati, LOC_LATI, LOOK_LATI);
                           elseif initial_LOOK ~= 0
                           end
                       elseif class ~= previous_class
                           if initial_LOOK == 0 
                               [CEP_lati, LOC_LATI, LOOK_LATI] = add_update(bond_num, class, N_particles, CEP_lati, LOC_LATI, LOOK_LATI);
                           elseif initial_LOOK ~= 0
                               [CEP_lati, LOC_LATI, LOOK_LATI] = add_update(bond_num, class, N_particles, CEP_lati, LOC_LATI, LOOK_LATI);
                               [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
                           end
                       end
                        
                        
                    
                    
                    elseif (ignore == true)         
                        if initial_LOOK == 0
                            
                        elseif initial_LOOK ~= 0
                            [CEP_lati, LOC_LATI, LOOK_LATI] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_lati, LOC_LATI, LOOK_LATI);
                        end
                        
                    end
                

                elseif center_site(2) == beside_site(m, 2)     %This implise longitudinal interaction on the lattice.
                    E_j = center_site(2);
                    E_i = min(center_site(1), beside_site(m, 1));
                    if (E_i == 1) & ((center_site(1) == Li) | (beside_site(m, 1) == Li))
                        E_i = Li;     %At the periodic boundary, the E_i = min() argument doesn't work...
                    end
                    previous_class = E_long(E_i, E_j);
                    E_long(E_i, E_j) = class;
                    bond_num = bond_number([E_i, E_j], Li);   
                    initial_LOOK = LOOK_LONG(bond_num);
                    
                    %update LOC_LONG
                    if (ignore == false)
                       if class == previous_class
                           if class == 0 & previous_class == 0 & initial_LOOK == 0
                               [CEP_long, LOC_LONG, LOOK_LONG] = add_update(bond_num, class, N_particles, CEP_long, LOC_LONG, LOOK_LONG);
                           elseif initial_LOOK ~= 0
                           end
                       elseif class ~= previous_class
                           if initial_LOOK == 0 
                               [CEP_long, LOC_LONG, LOOK_LONG] = add_update(bond_num, class, N_particles, CEP_long, LOC_LONG, LOOK_LONG);
                           elseif initial_LOOK ~= 0
                               [CEP_long, LOC_LONG, LOOK_LONG] = add_update(bond_num, class, N_particles, CEP_long, LOC_LONG, LOOK_LONG);
                               [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
                           end
                       end
                        
                        
                    
                    
                    elseif (ignore == true) 
                        if initial_LOOK == 0
                            
                        elseif initial_LOOK ~= 0
                            [CEP_long, LOC_LONG, LOOK_LONG] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP_long, LOC_LONG, LOOK_LONG);
                        end
                        
                        
                    end

                  
                end
               
                
            end
            end
           
        end
        
    end
end

    




   

function bdn = bond_number(bond, E_row_dim)
%bond is the site of the bond (interaction) as define for E_lati and E_long
%E_row_dim is the row dimensions of the bond (interaction) matrix.
    bdn = bond(1) + (bond(2) - 1)*E_row_dim;
end


function Ln = LOC_number(LOC_site, LOC_length)
    %This gives the number associated to some site in the LOC matrix.
    %The numbering reads as, for a (7,10) matrix for example, 1 -> 10 for the 
    %elements in the first row, then 11 -> 20 for those in the second row, 
    %and so on until 61 -> 70. So LOC_length in this case would be 10.
    %LOC_site is a tuple, like [3,4].
    Ln = LOC_site(2) + (LOC_site(1) - 1)*LOC_length;
end


function j_remain = j_axis_LOCation(bond_index_in_LOC, LOC_length)
    %This gives, for a given class row in LOC_LATI or LONG, the j-axis index 
    %of the element 'bond_num'.
    j_remain = rem(bond_index_in_LOC, LOC_length);
    if j_remain == 0
        j_remain = LOC_length;    %If the j axis index is == LOC_length, the rem give 
    end                           %0 instead of the value of LOC_length.
end



%Super important two update functions
function [CEP, LOC, LOOK] = add_update(bond_num, class, N_particles, CEP, LOC, LOOK)
    %This fundtion appropriately adds a bond number to a new class row in
    %the LOC array. It properly updates The values in CEP and the LOCation
    %value in the array LOOK.
    LOC_length = size(LOC,2);
    last_in_class_plusone = CEP(class + 4) + 1;
    LOC(class + 4, last_in_class_plusone) = bond_num;
    LOOK(bond_num) = LOC_number([class + 4, last_in_class_plusone], LOC_length) ;
    CEP(class + 4) = last_in_class_plusone;
end

function [CEP, LOC, LOOK] = delete_update(initial_LOOK, previous_class, N_particles, bond_num, ignore, CEP, LOC, LOOK)
    %Deletes a bond information from CEP LOC and LOOK.
    LOC_length = size(LOC,2);
    j_index_prev = j_axis_LOCation(initial_LOOK, LOC_length);
    
    class_end_point_bond = LOC(previous_class + 4, CEP(previous_class + 4));
    if class_end_point_bond ~= bond_num
        LOC(previous_class + 4, j_index_prev) = class_end_point_bond;
        %Must update the bond at 'LOC(CEP(previous_class + 4))' LOCation in the LOOK array.
        LOOK(class_end_point_bond) = LOC_number([previous_class + 4, j_index_prev], LOC_length);
    end
   
    
    
    LOC(previous_class + 4, CEP(previous_class + 4)) = 0;
    CEP(previous_class + 4) = CEP(previous_class + 4) - 1;
    if (ignore == true)          %Only execute if it has newly become a 'ignore' case.
        LOOK(bond_num) = 0;
    end
end



%=======================================================================================================================================
%Below is all the updater stuff for collective moves.

%Constructing CTM_updater_2() December 9th, 2021.======================
function [CTM, CTM_site, max_grp_num] = CTM_updater_2(CTM, CTM_site, max_grp_num, interact, remaining_cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep)
    %This function is just to update the three CTM variables and nothing
    %else. interact is used instead of Cluster_elements, since this
    %function is only useful for single particle moves. 
   % disp("This is #2");
   particle_exitted = false;
    N_particles = size(site,1);
    if (interact(1,1) ~= 0) & (interact(1,2) ~= 0)
        quasi_particle_1 = C(interact(1,1), interact(1,2)); 
        quasi_particle_2 = C(interact(2,1), interact(2,2));
    elseif (interact(1,1) == 0) & (interact(1,2) == 0) %The case where a partcle is leaving is when interact = [0,0; x,y] where x,y is the site the particle was before exitting.
        particle_exitted = true;
        
    end
    remaining_cluster_elements;
    if isempty(remaining_cluster_elements)
        %do nothing
    elseif (particle_exitted == true) | ((quasi_particle_1 ~= 0) & (quasi_particle_2 == 0)) | ((quasi_particle_1 == 0) & (quasi_particle_2 ~= 0)) 
        %When exchanged between particle and a vacancy, an update of the previous cluster may be needed.
        %First you gotta delete everything related to
        %remaining_cluster_elements.
        CTM;
        CTM_site;
        leader_particle = remaining_cluster_elements(1);
        group_number = CTM_site(leader_particle,2)  ;  %The group number is the column number the group is stored in CTM.
        
        for i = 1:length(remaining_cluster_elements ) %delete CTM_site history for the remaining guys.
%             CTM_site(remaining_cluster_elements(i),:) = [0,0];
%         end
%             %-----------------------------------
           % fprintf('lalala %d', remaining_cluster_elements(i));
            particle = remaining_cluster_elements(i);
      %      disp(CTM);
            begin_row = CTM_site(particle,1);
            former_group_below = nonzeros(CTM(begin_row + 1:N_particles,group_number));
            CTM(begin_row:begin_row + length(former_group_below),group_number) = zeros(length(former_group_below) + 1,1);
            CTM(begin_row:begin_row + length(former_group_below) -1, group_number) = [former_group_below];
            CTM_site(particle,:) = zeros(1,2);
            %The whole column must update
            current_group = nonzeros(CTM(:,group_number)) ;
            CTM_site = CTM_site_update(group_number, CTM_site, current_group);
          %  fprintf('lalala Ending %d', remaining_cluster_elements(i));
            CTM;
            CTM_site;
        end
%             %-------------------------------

%         %----------------------------------
        
        if CTM(:,group_number) == zeros(N_particles,1)
            if group_number ~= max_grp_num
                
                CTM(:,group_number) = CTM(:, max_grp_num);
                CTM_site = CTM_site_update(group_number, CTM_site, nonzeros(CTM(:, max_grp_num)));
                CTM(:,max_grp_num) = zeros(N_particles, 1);
                
            elseif group_number == max_grp_num
                CTM(:,group_number) = zeros(N_particles, 1);
            end
            max_grp_num = max_grp_num - 1;
        end
        %----------------------------------
        

        CTM;
        CTM_site;
        %each remaining cluster element must be updated independently.
        length(remaining_cluster_elements);
        if length(remaining_cluster_elements) > 1
            for eli = 1:length(remaining_cluster_elements)
                c_element = [remaining_cluster_elements(eli)];
                [CTM, CTM_site, max_grp_num] = CTM_updater_1(CTM, CTM_site, max_grp_num, c_element, site, C, im_i, ip_j, ip_i, im_j, prep);
            end
        end
   
        
    end        
        
        
end
%======================================================================



function [CTM, CTM_site, max_grp_num, remaining_cluster] = CTM_updater_1(CTM, CTM_site, max_grp_num, Cluster_elements, site, C, im_i, ip_j, ip_i, im_j, prep)
    %This function updates everything, except for the case where the single
    %particle detachment leaves the cluster divided. That will be
    %separately done in CTM_updater_2 ... ! (2 means it is secondary, CTM_updater_1
    %does most of the heavy lifting in terms of updating CTM and CTM_site for
    %single and collective movement)
    remaining_cluster = [];
    Np = length(Cluster_elements);
    N_particles = size(CTM,1);
    nbrs_info = zeros(Np, 4);
    if Cluster_elements(1) ~= 0     %This means the cluster_element list is a usual case.
    for particle_index = 1:Np
        particle = Cluster_elements(particle_index);
        nbrs_near = neighbors(particle, site, C, im_i, ip_j, ip_i, im_j);
        nbrs_info(particle_index, 1:length(nbrs_near)) = nbrs_near;
    end
    elseif Cluster_elements(1) == 0    %if Cluster_elements = [0, particle], then 0 indicates an exit move.
        Np = 1;
        nbrs_info = [];
    end
        
    %update CTM cluster elements
    for i = 1:Np
        if Cluster_elements(1) ~= 0
            small_nbr_grp = [Cluster_elements(i), transpose(nonzeros(nbrs_info(i,:)))];
        elseif Cluster_elements(1) == 0
            small_nbr_grp = [Cluster_elements(2)];    %For the case of a particle exit, the 2ndindex is the particle label that exitted.
            Cluster_elements = [Cluster_elements(2)]; %take away the 0, since it is not needed any more.
        end
        belong_to_grps = [];
        if length(small_nbr_grp) == 1   %This case would be a single particle detachment from cluster, with no reattachment to another.
            particle = Cluster_elements(1);
            part_of = CTM_site(particle,2) ; %The j index represents the group this particle used belongs to. 0 if the particle is a monomer.
            if part_of > 0
                begin_row = CTM_site(particle,1);
                former_group_below = nonzeros(CTM(begin_row + 1:N_particles,part_of));
                CTM(begin_row:begin_row + length(former_group_below),part_of) = zeros(length(former_group_below) + 1,1);
                CTM(begin_row:begin_row + length(former_group_below) -1, part_of) = [former_group_below];
                CTM_site = CTM_site_update(part_of, CTM_site, former_group_below);
                CTM_site(particle,:) = zeros(1,2);
                remaining_cluster = transpose(nonzeros(CTM(:,part_of)));
            end
            
           
            
            
            
            
        elseif length(small_nbr_grp) > 1
            for member = small_nbr_grp
                if member < 0
                    member
                    CTM_site
                    CTM
                    site
                    small_nbr_grp
                    Cluster_elements
                end
                part_of = CTM_site(member,2);  %The j index represents the group this particle already belongs to. 0 if the particle is a monomer.
                if part_of > 0
                    belong_to_grps = [belong_to_grps, part_of];
                end
            end
            
            if Np == 1       %This is a deleting the history is the case is a single particle detachment and reattach to another cluster.
                particle = Cluster_elements(1);
                previous_cluster = CTM_site(particle,2);
                any_previous = sum(previous_cluster == belong_to_grps(2:end));    %Check if any new neighbors belong to the previous cluster.
                special = false;
                if any_previous == length(belong_to_grps(2:end)) & ~(isempty(belong_to_grps(2:end)))
                    special = true;
                end
                %The following if statement can be made more efficient if
                %you new the previous nbrs of particle. YOu would need to
                %use any_previous
                if (previous_cluster ~= 0)  %length(previous_cluster) is to check if
                     if prep == false                                                                  %any_previous is < the length, then some are missing
                       % fprintf('hello %d', particle)   ;                                                                %potentially leading to a fragmented cluster -> need to 
                        small_nbr_grp  ;                                                                %update it!
                        for indy2 = length(small_nbr_grp):-1:1
                            nbr_p = small_nbr_grp(indy2);
                            if CTM_site(nbr_p,2) == previous_cluster
                                small_nbr_grp(indy2) = [];
                            end
                        end
                        small_nbr_grp;
                     end
                                                                                       
                    begin_row = CTM_site(particle,1);
                    former_group_below = nonzeros(CTM(begin_row + 1:N_particles,previous_cluster));
                    CTM(begin_row:begin_row + length(former_group_below),previous_cluster) = zeros(length(former_group_below) + 1,1);
                    CTM(begin_row:begin_row + length(former_group_below) -1, previous_cluster) = [former_group_below];
                    CTM_site = CTM_site_update(previous_cluster, CTM_site, former_group_below);
                    
                    CTM_site(particle,:) = zeros(1,2);
                    belong_to_grps = belong_to_grps(2:end) ;     %Only keep the groups of the neighbors, b/c the single particle no longer
                                                                % belongs to the previous cluster.  
                    %Move all previous_cluster values in belong_to_grps.
                     remaining_cluster = transpose(nonzeros(CTM(:, previous_cluster)));
%                     small_nbr_grp
                    if prep == false
                        belong_to_grps(belong_to_grps == previous_cluster) = [];
                    end
%                     if length(remaining_cluster) > 1 & length(unique(belong_to_grps)) > 1 
%                         for indy = length(belong_to_grps):-1:1
%                             if belong_to_grps(indy) == previous_cluster
%                                 belong_to_grps(indy) = [];
%                             end
%                         end
%                     for indy2 = length(small_nbr_grp):-1:1
%                         nbr_p = small_nbr_grp(indy2)
%                         if CTM_site(nbr_p,2) == previous_cluster
%                             small_nbr_grp(indy2) = [];
%                         end
%                     end
%                     end
%                     belong_to_grps
                    
                     small_nbr_grp;
%                     remaining_cluster = transpose(nonzeros(CTM(:, previous_cluster)))


                    
                end
            end
            
            
            belong_to_grps = unique(belong_to_grps )  ;  %unique also sorts the list in ascending order.
            if isempty(belong_to_grps) == true & (isempty(small_nbr_grp) == false)
                max_grp_num = max_grp_num + 1;
                column_number = max_grp_num;
                small_nbr_grp = unique([[Cluster_elements(1)], small_nbr_grp]);
                CTM(1:length(small_nbr_grp), max_grp_num) = small_nbr_grp;
                CTM_site = CTM_site_update(max_grp_num, CTM_site, unique(small_nbr_grp));
            elseif (length(belong_to_grps) == 1) & (isempty(belong_to_grps) == false) & (isempty(small_nbr_grp) == false)
                    
                    small_nbr_grp = unique([[Cluster_elements(1)], small_nbr_grp]);
                    column_number = belong_to_grps(1);
                    group = unique([transpose(nonzeros(CTM(:,column_number))), small_nbr_grp]);

                    CTM(1:length(group), column_number) = group;
                    CTM_site = CTM_site_update(column_number, CTM_site, group);
                    
            elseif length(belong_to_grps) > 1
                %MERGE!!
               
                small_nbr_grp = unique([[Cluster_elements(1)], small_nbr_grp]);
                principle_column = belong_to_grps(1);
                group = [transpose(nonzeros(CTM(:,principle_column))), small_nbr_grp];

                merge_list = belong_to_grps(2:end);
                for merge_column_index = 1:length(merge_list)
                    merge_column = merge_list(merge_column_index);
                    group = [group, transpose(nonzeros(CTM(:,merge_column)))];
                    
                end
                group = unique(group);
                CTM(1:length(group), principle_column) = group;
                CTM_site = CTM_site_update(principle_column, CTM_site, group);
                
                for merge_column_index = length(merge_list):-1:1
                    merge_column = merge_list(merge_column_index);
                    group = [group, transpose(nonzeros(CTM(:,merge_column)))];
                    if merge_column < max_grp_num
                        max_group = unique([transpose(nonzeros(CTM(:,max_grp_num)))]);
                        CTM(:,merge_column) = CTM(:, max_grp_num);
                        CTM(:, max_grp_num) = zeros(N_particles, 1);
                        CTM_site = CTM_site_update(merge_column, CTM_site, max_group);
                        max_grp_num = max_grp_num - 1;
                    elseif merge_column == max_grp_num
                        CTM(:, max_grp_num) = zeros(N_particles, 1);
                        max_grp_num = max_grp_num - 1;
                    end
                end
                
            end

        end
    end
end


function CTM_site = CTM_site_update(column, CTM_site, grp_elements)
    %the grp_elements assumes the first one will be at the top of the list,
    %and that after the list end, if the length is less than N_particles,
    %then the rest of the column in CTM is just zeros. 
    row = 1;
    for indy = 1:length(grp_elements)
        grp_elem = grp_elements(indy);
        CTM_site(grp_elem,:) = [row, column];
        row = row + 1;
    end
end

%=========================================================================
%Finished functions.

function [C, site, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = Collective_Shifter( C, Cluster_elements, site, im_i, ip_j, ip_i, im_j, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, N_particles, Li,  CTM_site, direction)
    %This function shifts the cluster in a random direction (usual setting is to have a random direction). If the direction
    %is specified, it will shift the cluster in that direction.
    cluster_move = true;
    random_direction = randi(4);
  
    if ~exist('direction', 'var')
        direction = random_direction;
    end
    
     %Relevant boundary condition direction
     if direction == 1
         rbc = im_i;
     elseif direction == 2
         rbc = ip_j;
     elseif direction == 3
         rbc = ip_i;
     elseif direction == 4
         rbc = im_j;
     end

     direction_sites = [] ;  %To store the index in the direction of movement of the elements in the cluster.
                             %The index of the direction_sites is the same as
                             %the index of Cluster_elements and can be used to
                             %retrieve the Cluster_elements number.
     for p = Cluster_elements
        if (direction == 2) | (direction == 4)   %A move to the right or left.
            direction_sites(end + 1) = site(p, 2);
        elseif (direction == 1) | (direction == 3)   %A move to the up or down.
            direction_sites(end + 1) = site(p, 1);
        end
     end
     
     %To make it in the 
     
     Cluster_elements;
     direction_sites;
     LI = size(C, 1);
     %==========================================================
            
    C;
     Cluster_elements;
     direction_sites;
     %This whole following two blocks are to to know in what order to move the cluster
     %elements. The result is a list of cluster elements in the order that they should be moved. 
     check_order = [1];
     if (direction == 1) | (direction == 4)
        for m = 2:length(direction_sites)
            included = false;
            n = 1;
            while included == false
                n;
                if (n == 1) & direction_sites(m) < direction_sites(check_order(1))
                    check_order = [m, check_order];
                    included = true;
                elseif (n < length(check_order)) & (direction_sites(check_order(n)) <= direction_sites(m)) & (direction_sites(m) < direction_sites(check_order(n + 1))) 
                    %check_order stores the index of the cluster elements.
                    %Their correponding values in direction_sites
                    %determines where they are placed in check_order.
                    check_order = [check_order(1:n), m, check_order(n+1:end)];             %The index of the cluster element is stored in it's position in check_order.
                    included = true;
                elseif n == length(check_order)
                    check_order = [check_order, m];    %This makes it work!
                    included = true;
                end
                n = n + 1;
            end
        end
     elseif (direction == 2) | (direction == 3)
         for m = 2:length(direction_sites)
            included = false;
            n = 1;
            while included == false
                if (n == 1) & direction_sites(m) >= direction_sites(check_order(1))
                    check_order = [m, check_order];
                    included = true;

                elseif (n < length(check_order)) & (direction_sites(check_order(n)) > direction_sites(m)) & (direction_sites(m) >= direction_sites(check_order(n + 1)))
                    check_order = [check_order(1:n), m, check_order(n+1:end)];
                    included = true;
                elseif n == length(check_order)
                    check_order = [check_order, m];
                    included = true;
                end
                n = n + 1;
            end
         end
     end
     
     %If direction = 3, and this is a move across periodic boundary, then
     %you need to make the right tip of the cluster moves first. This is
     %because the cluster must be moved (and updated) from the front of the
     %cluster (in the direction of motion).
     check_order;
     length_of_check_order  = length(check_order);
     for num = 1:length_of_check_order
            check_direction(num) = direction_sites(check_order(num));
     end
     
     check_direction = unique(check_direction); %We only need to know which rows are occupied. This orders the list in increasing order.
     length_of_check_dir = length(check_direction);
     %The following reorders the list so that it starts from the 'real' tip of the
     %cluster
     stick_out_row = NaN;
     if (direction == 1) | (direction == 3)
         for i  = 1:length_of_check_dir
             if i < length_of_check_dir
                 if abs((check_direction(i + 1) - check_direction(i))) > 1
                     if direction == 3
                        stick_out_row = check_direction(i);   %if direction == 3.
                     elseif direction == 1
                         stick_out_row = check_direction(i + 1);
                     break
                     end
                 end
             end
         end
     end
     
     stick_out_particles = [];
     for ele = Cluster_elements
         if site(ele, 1) == stick_out_row
            stick_out_particles(end+1)  = ele;
         end
     end
     
     
     
%      direction_sites
%      if ((direction == 1) | (direction == 3)) & (length(unique(direction_sites)) == Li)
%         check_order = sort(direction_sites)
%      end
     
     
     
   %  if (direction == 2) | (direction == 4)
         check_order;
         for num = 1:length(check_order)
            check_order(num) = Cluster_elements(check_order(num));
         end
   %  end
   
   %Flip check_order if it crosses the periodic boundary. Move the first
   %appearance of a particle in stick_out_particles to the front of the
   %list.
%    check_order;
%    disp('here')
%    stick_out_particles
   for i = 1:length_of_check_order
       if any(stick_out_particles == check_order(i) )  %if the particle belongs to the stickout row...
           check_order = [check_order(i:end), check_order(1:i-1)];
           break
       end
       
   end
   
%     check_order
%     Cluster_elements
%     rbc
    %Move everything!  This version keeps the cluster stiff. It can still merge
    %with other clusters
    Lj = ip_j(end - 1);   %The second to last element of ip_j should be Lj.
    if ( (direction == 2) | (direction == 4) ) & ((rbc(site(check_order(1),2)) == 0) | (rbc(site(check_order(1),2)) == 1) | (rbc(site(check_order(1),2)) == Lj))
        
        if ((rbc(site(check_order(1),2)) == 0) | (rbc(site(check_order(1),2)) == 1) | (rbc(site(check_order(1),2)) == Lj))
            %Do nothing. The (== 0) case is becuse it is at the edges, but the (== 1) and (== Lj) is because it messes the E_long,lati
            %matrices. But this makes it impossible to simulate cluster exit, so if considering cluster exit, you'll need to change this!
            
        end
    else 
            
       % for element = Cluster_elements
        for element = check_order
            interact = [site(element,1), site(element,2)];
            %A(interact(1,1),interact(1,2)) = 0;
            C(interact(1,1),interact(1,2)) = 0;
%             
        end

            % %update the position in the site of the ONE values. +1  to the
            % %right.
            % site(random, 2) = ip_j(j); 

        interact_list_1 = [];
        interact_list_2 = [];
        old_interact = [];
        interact_counter = 1;
        %for element = Cluster_elements
        
        for element = check_order  
            
            interact = [site(element,1), site(element,2)];
            old_interact(interact_counter,:) = interact;
            if (direction == 2) | (direction == 4)   %A move to the right or left.
                interact = [interact; site(element,1), rbc(site(element,2))];
            elseif (direction == 1) | (direction == 3)   %A move up or down.
                interact = [interact; rbc(site(element,1)), site(element,2)];
            end
            
            interact_list_1(interact_counter,:) = interact(1,:);
            interact_list_2(interact_counter,:) = interact(2,:);
            interact_counter = interact_counter + 1;
            %A(interact(2,1), interact(2,2)) = 1;
            C(interact(2,1), interact(2,2)) = element;
            
            
            
            
            
            site(element,:) = interact(2,:);
%             if temp_2 ~= 0
%                 site(temp_2,:) = interact(1,:);
%             end 

        end
        for i_interact = 1:(interact_counter - 1) 
            
            interact = [interact_list_1(i_interact,:); interact_list_2(i_interact,:)];
            
            
            
            [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
            %I will need to add some way to update 'site' in the
            %interaction_updater for single particle moves!
            

        end
        
        %This only updates with the bond that was in the direction that the particle was moved. This can lead to the
        %interaction_updater not updating the bonds at the front-end of the cluster, leaving some previous interactions in E_long,lati present when
        %there should be none. Thus, you must also update from a bond in the direction of movement from the current location. The reason the front-end of the cluster
        %isn't updated is because the bonds above use as a pivot the previous positions of the particles. When the fron-end (which can
        %have many different interactions with its surroundings, moves forward, the interaction of the new location with its previous
        %position comes out as two particles being present in both locations, resulting in the interaction_updater ignoring the
        %necessary updates. Thus, if we re-update using the interact pivot as the current location, it will see that in front of the current
        %location there is an empty space, thus the local bonds will be updated.
        
        for element = check_order 
            interact = [site(element,1), site(element,2)];
            old_interact(interact_counter,:) = interact;
            if (direction == 2) | (direction == 4)   %A move to the right or left.
                interact = [interact; site(element,1), rbc(site(element,2))];
            elseif (direction == 1) | (direction == 3)   %A move up or down.
                interact = [interact; rbc(site(element,1)), site(element,2)];
            end
            
            
            [E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI] = interaction_updater(interact, E_long, E_lati, LOC_LONG, LOC_LATI, CEP_long, CEP_lati, LOOK_LONG, LOOK_LATI, C, ip_i, ip_j, im_i, im_j, N_particles, CTM_site);
        end
        
  
        
    end
end



    




   






%========================================================================================================================================
%The following (above the next 'double line') is mainly to do with the
%phase analysis. Just counting the number of oligomers, monomers etc.


function number_wrap_clusters = wrapping_status(CTM, max_grp_num, site, Li)
    %This function determine how many cluster in the tube are wrapping the
    %tube, meaning that they include particles from every i-index of the
    %tube.
    number = 0;
    for cluster_number = 1:max_grp_num
        element_rows = [];
        cluster_elementss = nonzeros(CTM(:,cluster_number));
        particle_num = length(cluster_elementss);
        for element_index = 1:particle_num
            element_rows = [element_rows, site(cluster_elementss(element_index),1)];   %site(cluster_elementss(element_index),1) is the i-index (longitudinal index) of the element.
        end
        if length(unique(element_rows)) == Li;
            number = number + 1;
        end
    end
    number_wrap_clusters = number;
end


function AIS = Average_Island_Size(DC)
    %DC is a cell array of the form provided by nbe_distribution. Each cell
    %of the cell array DC has the constituing protein indices (where the indices are
    %determined by the site array) per cluster or monomer. ParticleNum is
    %just the total number of particles (proteins) on the lattice.
    %This function gives a map container (dictionary in python) where the
    %key is the island phase and the values are the average size of the
    %islands in that phase (Oligomer or CLuster). The number of the phases
    %are also included.
    Data_type = {'avg_dimer', 'avg_oligomer', 'avg_nanocluster', 'avg_cluster', 'num_monomer', 'num_dimer', 'num_oligomer', 'num_nanocluster', 'num_cluster'};
    ParticleNum  = 0;
    num_monomer = 0;
    num_dimer = 0;
    num_oligomer = 0;
    num_nanocluster = 0;
    num_cluster = 0;
    
    Dim_sizes = [];
    Olig_sizes = [];
    nanoClus_sizes = [];
    Clus_sizes = [];
    for i = 1:length(DC)
        intensity = length(DC{i});
        ParticleNum = ParticleNum + intensity;
        if intensity == 1
            num_monomer = num_monomer + 1;
        end
        if (intensity == 2) | (intensity == 3)
            Dim_sizes = [Dim_sizes, intensity];
        end
        if (intensity <= 10) & (intensity >= 4)
            Olig_sizes = [Olig_sizes, intensity];
        end
        if (intensity > 10) & (intensity <= 20)
            nanoClus_sizes = [nanoClus_sizes, intensity];
        end
        if (intensity > 20) 
            Clus_sizes = [Clus_sizes, intensity];
        end
    end
    num_dimer = length(Dim_sizes);
    num_oligomer = length(Olig_sizes);
    num_nanocluster = length(nanoClus_sizes);
    num_cluster = length(Clus_sizes);
    
    
    if length(Clus_sizes) == 0   %This is required since is Clis_size = [], the mean becomes NaN, instead of 0.
        Clus_sizes = 0;
    end
    if length(nanoClus_sizes) == 0
        nanoClus_sizes = 0;
    end
    if length(Olig_sizes) == 0
        Olig_sizes = 0;
    end
    if length(Dim_sizes) == 0
        Dim_sizes = 0;
    end
    
    AIS = containers.Map(Data_type, [mean(Dim_sizes), mean(Olig_sizes), mean(nanoClus_sizes), mean(Clus_sizes), num_monomer, num_dimer, num_oligomer, num_nanocluster,  num_cluster]);

end







function nacs = number_avg_mean_cluster_size(DC, N_particles)
    cluster_num = zeros(1, N_particles);
    for clus = 1:length(DC)
        size_of = length(DC{clus});
        cluster_num(size_of) = cluster_num(size_of) + 1;
    end
    nacs = N_particles/sum(cluster_num);
end


function wamcs = weight_avg_mean_cluster_size(DC, N_particles)
    cluster_num = zeros(1, N_particles);
    for clus = 1:length(DC)
        size_of = length(DC{clus});
        cluster_num(size_of) = cluster_num(size_of) + 1;
    end
    
    wamcs = 0;
    for clus_size  = 1:length(cluster_num)
        wamcs = wamcs + (clus_size^2)*cluster_num(clus_size);
    end
    wamcs = wamcs / N_particles;
end






function DoP = Dominant_Phase(FP)
    %Returns the dominant phases of the configuration (phase with maximum
    % fraction in FP). FP is Map.container of 4 cells with each cell
    % representing the phase, and the fraction of particles in that phase.
    
    Phases_list = cellstr(keys(FP));
    value_list = cell2mat(values(FP));
    max_fraction = max(value_list);
    max_phase_index = find(value_list == max_fraction);
    if length(max_phase_index) == 1
        DoP = string(Phases_list(max_phase_index));
    end
    if length(max_phase_index) > 1
        DoP = "intermediate";
    end

end





function FP = Fraction_of_phases(DC)
    %DC is a cell array of the form provided by nbr_distribution. Each cell
    %of the cell array DC has the constituing protein indices (where the indices are
    %determined by the site array) per cluster or monomer. ParticleNum is
    %just the total number of particles (proteins) on the lattice.
    %This function gives a map container (dictionary in python) where the
    %key is the island phase and the values are the percentage of particles
    %(proteins) in that phase
    phase = {'monomer', 'dimer', 'oligomer', 'nanocluster', 'cluster'};
    ParticleNum  = 0;
    phase_intensity = zeros(1,5);
    for i = 1:length(DC)
        intensity = length(DC{i});
        ParticleNum = ParticleNum + intensity;
        if intensity == 1
            phase_intensity(1) = phase_intensity(1) + intensity;
        end
        if (intensity == 2) |  (intensity == 3)
            phase_intensity(2) = phase_intensity(2) + intensity;
        end
        if (intensity <= 10) & (intensity >= 4)
            phase_intensity(3) = phase_intensity(3) + intensity;
        end
        if (intensity > 10)  &  (intensity <= 20)
            phase_intensity(4) = phase_intensity(4) + intensity;
        end
        if (intensity > 20)
            phase_intensity(5) = phase_intensity(5) + intensity;
        end
    end
    FP = containers.Map(phase, phase_intensity/ParticleNum);

end


function P = nbr_distribution(N_particles, site, C, im_i, ip_j, ip_i, im_j )
    %This function depends on the boundary conditions, but only through
    %it's use of the neighbors_nf() function, which here is set to be
    %no-flux b.c. It reterns all the groups of monomers, dimers,...
    %clusters, where the particles are labeled. {[1,2,3,7]},{[4]},{[5,6]}
    %is a cell array of a tetramer, monomer, and a dimer respectively.
%     site = zeros(N_particles, 2);
%     for id = 1:size(C,1)
%         for jd = 1:size(C,2)
%             if C(id, jd) ~= 0
%                 site(C(id, jd), 1) = id;
%                 site(C(id, jd), 2) = jd;
%             end
%         end
%     end
        

    P = {};    % the histogram
    nbr_list = {};
    site_dim = size(site, 1);
    for x = 1:site_dim
        
        x_nbrs = neighbors_nf(x, site, C, im_i, ip_j, ip_i, im_j);
        nbr_list{end + 1} =  [x, x_nbrs];
    end
    
    nbr_list;
    P{1} = nbr_list{1};
    for s = 2:site_dim

        included = [];
        for node = 1:length(nbr_list{s})
%             for cell_index = 1: length(histogram)
%                 any(histogram{cell_index} == nbr_list(s, node))
            for cell_index = 1:length(P)
                if any(P{cell_index} == nbr_list{s}(node))
                    included = [included, cell_index];
                end
            end
        end
        included = unique(included);
        if length(included) ~= 0
        P{included(1)} = cat(2, P{included(1)}, nbr_list{s});
        for n  = 2:length(included)
            P{included(1)} = cat(2, P{included(1)}, P{included(n)});
            P{included(n)} = [];
        end
        P{included(1)} = unique(P{included(1)});
        index = cellfun(@isempty, P) == 0;
        P = P(index);
        end

        if length(included) == 0
            P{end + 1} = nbr_list{s};
        end         
    end
end

%This neighbors is used for collective moves.
function nbrs = neighbors(x, site, C, im_i, ip_j, ip_i, im_j)
    %For site x, determine the site index number of the neighbors. 
    %site is a Nx2 array where N is the number of proteins (ones in the matrix 
    %used for the spin exchange. C is the LxL matrix (L is the lattice
    %dimensions) where the positions of 1's are replaced by their index in
    %the array site (the column index of site) 
    %look up table for mod
    
    iv = site(x,1);
    jv = site(x,2);
    x_site = [iv, jv];
    nbr_site = [iv, ip_j(jv);  ip_i(iv) , jv;  iv, im_j(jv);  im_i(iv), jv];
    
    %Remove those that are not included due to nf b.c. in the j-direction.
    for ii = 0:3
        if (nbr_site(4 - ii, 2) == 0) % |  (C(beside_site(4 - ii, 1), beside_site(4 - ii, 2)) == 0)
            nbr_site(4 - ii,:) = [];     %remove row that is not included due to nfbc.
        end
    end
    
    nbrs = [];
    for i = 1:length(nbr_site)
        nbr_index = C(nbr_site(i, 1), nbr_site(i, 2));
        if nbr_index ~= 0
            nbrs = [nbrs, nbr_index];
        end
    end
end


function nbrs = neighbors_nf(x, site, C, im_i, ip_j, ip_i, im_j)
    %For site x, determine the site index number of the neighbors. 
    %site is a Nx2 array where N is the number of proteins (ones in the matrix 
    %used for the spin exchange. C is the LxL matrix (L is the lattice
    %dimensions) where the positions of 1's are replaced by their index in
    %the array site (the column index of site) 
    %look up table for mod
    
%     Li = size(C,1);
%     for i = 1: Li
%         ip_i(i) = i + 1;
%         im_i(i) = i - 1;
%     end
%     ip_i(Li) = 1;
%     im_i(1) = Li;
%     
%     Lj = size(C,2);
%     for j = 1: Lj
%     ip_j(j) = j + 1;
%     im_j(j) = j - 1;
%     end
%     ip_j(Lj) = 0;
%     im_j(1) = 0;
    
    
    iv = site(x,1);
    jv = site(x,2);
    x_site = [iv, jv];
    
    
    
    nbr_site = [iv, ip_j(jv);  ip_i(iv) , jv;  iv, im_j(jv);  im_i(iv), jv];
    
    %Remove those that are not included due to nf b.c. in the j-direction.
    for ii = 0:3
        if (nbr_site(4 - ii, 2) == 0) % |  (C(beside_site(4 - ii, 1), beside_site(4 - ii, 2)) == 0)
            nbr_site(4 - ii,:) = [];     %remove row that is not included due to nfbc.
        end
    end
    
    nbrs = [];
    for i = 1:length(nbr_site)
        nbr_index = C(nbr_site(i, 1), nbr_site(i, 2));
        if nbr_index ~= 0
            nbrs = [nbrs, nbr_index];
        end
    end
end








function cnn = count_neighbors_nf(center_site, beside_site, C, Li, Lj, im_i, ip_j, ip_i, im_j)
    %for any site, it counts the number of occupied neighbors, excluding
    %site 'beside_site'. Beside_site is just a list e.g. [6,10] or [4,4].
    %Assumes no flux b.c.
    %site is a Nx2 array where N is the number of proteins (ones in the matrix 
    %used for the spin exchange. C is the LxL matrix (L is the lattice
    %dimensions) where the positions of 1's are replaced by their index in
    %the array site (the column index of site) 
    %look up table for mod
    
%     Li = size(C,1);
%     for i = 1: Li
%         ip_i(i) = i + 1;
%         im_i(i) = i - 1;
%     end
%     ip_i(Li) = 1;
%     im_i(1) = Li;
%     
%     Lj = size(C,2);
%     for j = 1: Lj
%     ip_j(j) = j + 1;
%     im_j(j) = j - 1;
%     end
%     ip_j(Lj) = 0;          %no-flux b.c.
%     im_j(1) = 0;     
    
    cnn = [];
    iv = center_site(1);
    jv = center_site(2);
    
    nbr_site_prep = [iv, ip_j(jv);  ip_i(iv) , jv;  iv, im_j(jv);  im_i(iv), jv];
    nbr_site = [];
    for i = 1:4
        if ~((nbr_site_prep(i, 1) == beside_site(1)) & (nbr_site_prep(i, 2) == beside_site(2)))
            nbr_site = [nbr_site;  nbr_site_prep(i, 1), nbr_site_prep(i, 2)];
        end
    end
    
    for i = 1:length(nbr_site)
        if (nbr_site(i, 2) ~= 0)    %if == 0, then it is not allowed by the nf b.c.
            nbr_index = C(nbr_site(i, 1), nbr_site(i, 2));
            if nbr_index ~= 0
                cnn = [cnn; nbr_site(i, 1), nbr_site(i, 2)];
            end
        end
    end
   
end
