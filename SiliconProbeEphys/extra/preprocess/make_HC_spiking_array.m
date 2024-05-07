function make_HC_spiking_array(frb4, year)
%% Combine 'Spikes_5LoomsHC_L1/L2' to make 'all_hist' array and 'all_av_hist' array. 
% 'Spikes_5Looms'.. made from 'Core_Loom_Processing.m'

% Need to add the number of frames added before the loom started in 'core_loom..' - this is 'frb4' - normally 180 or 60; 
% Also need to add year.

% Created by Burnett Apirl 2021 - adapted from 'add_REAL_DEPTH.m'

%% Contains the spiking (Hz) - over intervals of 0.0167s (1 frame). 

% if year == 2020
% Loom start = 180
% Loom end = 395 
% Each loom = 43 - (43*5 = 215)
% Array end = 575 - (180+215+180)

% elseif year == 2021
% Loom start = 180
% Loom end = 410
% Each loom = 45 - (45*5 = 225)
% Array end = 585 - (180+225+180)

%% 'all_hist' is an array that contains the spiking information, averaged over 1 frame (0.0167s) with each repetitions of the loom as a new row.
% 'all_av_hist' is the same but each row is the average response over the 5 looms for each animal.

% Each row is a cluster/cell from one recording.
% The first columns are the spiking information with each column containing
% the total spikes per 1 frame.
% Later columns contain information about:
% - the animal number
% - the animal geno
% - the real depth of that cell/cluster.
% - the loom number.

%%
files = dir('*5LoomsHC*');
num_files = numel(files);

het_animals = ["7614", "7476", "7790", "7269","1970", "1385", "1394", "1386"];

% frb4 = 60;
% year = 2021;

if year == 2020
    loomfr = 43;
    szarray = loomfr*5+frb4*2;
elseif year == 2021
    loomfr = 46;
    szarray = loomfr*5+frb4*2;
end

% This will contain every cell's response to each loom bout individually.
all_hist = [];

%This will contain the average response of each cell over the 5 bouts.
all_av_hist = [];

% This will contain the average response over the acclim.
all_acc_hist = [];



if year == 2020
    
    for i = 1:num_files
        fname = files(i).name;
        load(fname);
        
        ani = str2double(fname(10:13));
        probe_depth = str2double(fname(21:24));
        date =  str2double(fname(1:6));
        
        if contains(string(ani), het_animals)
            geno = 0;
        else
            geno = 1;
        end
        
        % L1
        n_cl = numel(all_spikes_L1(:,1));
        
        
        %% Add details in separate columns to each LOOM array individually.
        sz_acc = [numel(all_spikes_L1(1,:)), numel(all_spikes_L2(1,:)), numel(all_spikes_L3(1,:)), numel(all_spikes_L4(1,:)), numel(all_spikes_L5(1,:))];
        max_val = find(sz_acc == max(sz_acc));
        if numel(max_val>2)
            max_val = max_val(1);
        end
        
        allsp_Acc = [];
        
        % Loom1
        for di = 1:n_cl
            real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe.
            all_spikes_L1(di, szarray+1) = date;
            all_spikes_L1(di, szarray+2) = ani;
            all_spikes_L1(di, szarray+3) = 1;
            all_spikes_L1(di, szarray+4) = real_depth;
            all_spikes_L1(di, szarray+5) = geno;
            
            all_spikes_L2(di, szarray+1) = date;
            all_spikes_L2(di, szarray+2) = ani;
            all_spikes_L2(di, szarray+3) = 2;
            all_spikes_L2(di, szarray+4) = real_depth;
            all_spikes_L2(di, szarray+5) = geno;
            
            all_spikes_L3(di, szarray+1) = date;
            all_spikes_L3(di, szarray+2) = ani;
            all_spikes_L3(di, szarray+3) = 3;
            all_spikes_L3(di, szarray+4) = real_depth;
            all_spikes_L3(di, szarray+5) = geno;
            
            all_spikes_L4(di, szarray+1) = date;
            all_spikes_L4(di, szarray+2) = ani;
            all_spikes_L4(di, szarray+3) = 4;
            all_spikes_L4(di, szarray+4) = real_depth;
            all_spikes_L4(di, szarray+5) = geno;
            
            all_spikes_L5(di, szarray+1) = date;
            all_spikes_L5(di, szarray+2) = ani;
            all_spikes_L5(di, szarray+3) = 5;
            all_spikes_L5(di, szarray+4) = real_depth;
            all_spikes_L5(di, szarray+5) = geno;
            
            allsp_AV(di, szarray+1) = date;
            allsp_AV(di, szarray+2) = ani;
            allsp_AV(di, szarray+3) = 6;
            allsp_AV(di, szarray+4) = real_depth;
            allsp_AV(di, szarray+5) = geno;
            
            % Acclim -
            if max_val == 1
                acc_av = mean(all_spikes_L1(di, 1:frb4)*60);
                acc_std = std(all_spikes_L1(di, 1:frb4)*60);
                acc_var = var(all_spikes_L1(di, 1:frb4)*60);
                fr_acc = numel(all_spikes_L1(di, 1:frb4));
            elseif max_val == 2
                acc_av = mean(all_spikes_L2(di, 1:frb4)*60);
                acc_std = std(all_spikes_L2(di, 1:frb4)*60);
                acc_var = var(all_spikes_L2(di, 1:frb4)*60);
                fr_acc = numel(all_spikes_L2(di, 1:frb4));
            elseif max_val == 3
                acc_av = mean(all_spikes_L3(di, 1:frb4)*60);
                acc_std = std(all_spikes_L3(di, 1:frb4)*60);
                acc_var = var(all_spikes_L3(di, 1:frb4)*60);
                fr_acc = numel(all_spikes_L3(di, 1:frb4));
            elseif max_val == 4
                acc_av = mean(all_spikes_L4(di, 1:frb4)*60);
                acc_std = std(all_spikes_L4(di, 1:frb4)*60);
                acc_var = var(all_spikes_L4(di, 1:frb4)*60);
                fr_acc = numel(all_spikes_L4(di, 1:frb4));
            elseif max_val == 5
                acc_av = mean(all_spikes_L5(di, 1:frb4)*60);
                acc_std = std(all_spikes_L5(di, 1:frb4)*60);
                acc_var = var(all_spikes_L5(di, 1:frb4)*60);
                fr_acc = numel(all_spikes_L5(di, 1:frb4));
            end
            
            allsp_Acc(di, 1) = acc_av;
            allsp_Acc(di, 2) = acc_std;
            allsp_Acc(di, 3) = acc_var;
            allsp_Acc(di, 4) = fr_acc;
            allsp_Acc(di, 5) = date;
            allsp_Acc(di, 6) = ani;
            allsp_Acc(di, 7) = real_depth;
            allsp_Acc(di, 8) = geno;
            
        end
        
        %Vertically concatenate all the cells for each
        av_hist = vertcat(all_spikes_L1, all_spikes_L2, all_spikes_L3, all_spikes_L4, all_spikes_L5);
        all_hist = vertcat(all_hist, av_hist);
        
        %Vertically concatenate all the cells for average
        all_av_hist = vertcat(all_av_hist, allsp_AV);
        
        %Vertically concatenate all the cells for acclim
        all_acc_hist = vertcat(all_acc_hist, allsp_Acc);
    end
    
    % Muliply all_hist_WT by 60 to get frequency of spikes in Hz. - currently spikes/frame.
    all_hist(:,1:szarray) = all_hist(:,1:szarray)*60;
    all_av_hist(:,1:szarray) = all_av_hist(:,1:szarray)*60;
    
elseif year == 2021
    
    
    for i = 1:num_files
        fname = files(i).name;
        load(fname);
        
        ani = str2double(fname(10:13));
        probe_depth = str2double(fname(21:24));
        date =  str2double(fname(1:6));
        
        if contains(string(ani), het_animals)
            geno = 0;
        else
            geno = 1;
        end
        
        % L1
        n_cl = numel(allsp_L1(:,1));
        
        
        %% Add details in separate columns to each LOOM array individually.
        sz_acc = [numel(all_spikes_L1(1,:)), numel(all_spikes_L2(1,:)), numel(all_spikes_L3(1,:)), numel(all_spikes_L4(1,:)), numel(all_spikes_L5(1,:))];
        max_val = find(sz_acc == max(sz_acc));
        if numel(max_val>2)
            max_val = max_val(1);
        end
        
        allsp_Acc = [];
        
        % Loom1
        for di = 1:n_cl
            real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe.
            allsp_L1(di, szarray+1) = date;
            allsp_L1(di, szarray+2) = ani;
            allsp_L1(di, szarray+3) = 1;
            allsp_L1(di, szarray+4) = real_depth;
            allsp_L1(di, szarray+5) = geno;
            
            allsp_L2(di, szarray+1) = date;
            allsp_L2(di, szarray+2) = ani;
            allsp_L2(di, szarray+3) = 2;
            allsp_L2(di, szarray+4) = real_depth;
            allsp_L2(di, szarray+5) = geno;
            
            allsp_L3(di, szarray+1) = date;
            allsp_L3(di, szarray+2) = ani;
            allsp_L3(di, szarray+3) = 3;
            allsp_L3(di, szarray+4) = real_depth;
            allsp_L3(di, szarray+5) = geno;
            
            allsp_L4(di, szarray+1) = date;
            allsp_L4(di, szarray+2) = ani;
            allsp_L4(di, szarray+3) = 4;
            allsp_L4(di, szarray+4) = real_depth;
            allsp_L4(di, szarray+5) = geno;
            
            allsp_L5(di, szarray+1) = date;
            allsp_L5(di, szarray+2) = ani;
            allsp_L5(di, szarray+3) = 5;
            allsp_L5(di, szarray+4) = real_depth;
            allsp_L5(di, szarray+5) = geno;
            
            allsp_AV(di, szarray+1) = date;
            allsp_AV(di, szarray+2) = ani;
            allsp_AV(di, szarray+3) = 6;
            allsp_AV(di, szarray+4) = real_depth;
            allsp_AV(di, szarray+5) = geno;
            
            % Acclim -
            if max_val == 1
                acc_av = mean(all_spikes_L1(di, szarray:end)*60);
                acc_std = std(all_spikes_L1(di, szarray:end)*60);
                acc_var = var(all_spikes_L1(di, szarray:end)*60);
                fr_acc = numel(all_spikes_L1(di, szarray:end));
            elseif max_val == 2
                acc_av = mean(all_spikes_L2(di, szarray:end)*60);
                acc_std = std(all_spikes_L2(di, szarray:end)*60);
                acc_var = var(all_spikes_L2(di, szarray:end)*60);
                fr_acc = numel(all_spikes_L2(di, szarray:end));
            elseif max_val == 3
                acc_av = mean(all_spikes_L3(di, szarray:end)*60);
                acc_std = std(all_spikes_L3(di, szarray:end)*60);
                acc_var = var(all_spikes_L3(di, szarray:end)*60);
                fr_acc = numel(all_spikes_L3(di, szarray:end));
            elseif max_val == 4
                acc_av = mean(all_spikes_L4(di, szarray:end)*60);
                acc_std = std(all_spikes_L4(di, szarray:end)*60);
                acc_var = var(all_spikes_L4(di, szarray:end)*60);
                fr_acc = numel(all_spikes_L4(di, szarray:end));
            elseif max_val == 5
                acc_av = mean(all_spikes_L5(di, szarray:end)*60);
                acc_std = std(all_spikes_L5(di, szarray:end)*60);
                acc_var = var(all_spikes_L5(di, szarray:end)*60);
                fr_acc = numel(all_spikes_L5(di, szarray:end));
            end
            
            allsp_Acc(di, 1) = acc_av;
            allsp_Acc(di, 2) = acc_std;
            allsp_Acc(di, 3) = acc_var;
            allsp_Acc(di, 4) = fr_acc;
            allsp_Acc(di, 5) = date;
            allsp_Acc(di, 6) = ani;
            allsp_Acc(di, 7) = real_depth;
            allsp_Acc(di, 8) = geno;
            
        end
    
    %Vertically concatenate all the cells for each
    av_hist = vertcat(allsp_L1, allsp_L2, allsp_L3, allsp_L4, allsp_L5);
    all_hist = vertcat(all_hist, av_hist);
    
    %Vertically concatenate all the cells for average
    all_av_hist = vertcat(all_av_hist, allsp_AV);
    
    %Vertically concatenate all the cells for acclim
    all_acc_hist = vertcat(all_acc_hist, allsp_Acc);
    end 
    
end


% Muliply all_hist_WT by 60 to get frequency of spikes in Hz. - currently spikes/frame.
all_hist(:,1:szarray) = all_hist(:,1:szarray)*60;
all_av_hist(:,1:szarray) = all_av_hist(:,1:szarray)*60;


save('210706_Spiking_Looms_Ptchd1_HC.mat', 'all_hist', 'all_av_hist', 'all_acc_hist');

clear
 
end 
 
  
    
    
    
    