%% Core Loom Processing Script
% Modified by Burnett 01/03/21

% This script is for the processing of the loom data. 

% Stimulus starts with 15 looms in a row. 
% Then 'normal' 5*5 loom stimuli HIGH contrast
% Then 'normal' 5*5 loom stimuli LOW contrast
% Then Different speeds
% Then random locations
% Then flashes
% Then optomotor. 

%% Path where data is: 
folder_to_analyze = '/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Test/Test2_1650_Loom';
exp_name = folder_to_analyze(end-14:end);

%collect all needed files automatically, they should be inside folder_to_analyze
spike_cluster_file= fullfile(folder_to_analyze,'spike_clusters.npy');
spike_time_file= fullfile(folder_to_analyze,'spike_times.npy');
cluster_info_file= fullfile(folder_to_analyze,'cluster_group.tsv');
cluster_info_gen_file= fullfile(folder_to_analyze,'cluster_info.tsv');

%check if all files are there and return dat and log files
datfile= check_files_for_analysis(folder_to_analyze, spike_cluster_file, spike_time_file, cluster_info_file);
if isempty(datfile) 
    return;
end

main_dir = '/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe'; 

%for this session add folder with the functions to read npy files
addpath(fullfile(main_dir,'npy-matlab'))
addpath('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Scripts/PreProcess/Looms');

%% LOGS 

%get all stimuli log files and check repetitiions, there should be
%sum(nrep_per_stim) large gaps in red frames
%get the list of files in all folders and subfolders
%allsubfolders=['**',filesep,'*.*'];
filelist = dir(folder_to_analyze);%get list of files and folders in the given folder
filelist = filelist(~[filelist.isdir]);  %remove folders from list
%find all log files
log_files=[];
for li =1:length(filelist)
    if  contains(filelist(li).name,'.log','IgnoreCase', true)  && contains(filelist(li).name,'stimuli','IgnoreCase', true)  
        log_files=[log_files; filelist(li)];
    end
end

%Find the number of LOOM log files. 
loom_logs = []; 
for li =1:length(filelist)
    if  contains(filelist(li).name,'Loom','IgnoreCase', true)  && contains(filelist(li).name,'stimuli','IgnoreCase', true)  
        loom_logs=[loom_logs; filelist(li)];
    end
end

%sort log files by date
[~,idx] = sort([log_files.datenum]);

% Log 2 = multiloom 
% Logs 7,10, 13, 16, 19 = HC
% Logs 24, 27, 30, 33, 36 = LC
% Logs 41, 44, 47, 50, 53 = DS
% Log 58 = diffP

%% RED FRAMES

%find large gaps in the red signal that are the gaps between stimuli
%presentations or repetitions of the chirp
%get the syncing signal
[red_signal, sample_rate] = read_sync_signal_from_Intan_RHD2000_file(datfile);
% figure, plot(red_signal);

%each red frame is sampled multiple times
sec_per_frame=1/60; %0.0167;
frame_duration = int32(sec_per_frame*sample_rate);

redsum = movsum(red_signal, frame_duration); 
rf=redsum>10;
red_signal = rf & red_signal; %this will remove spurious red signal in between red frames

%remove red frames registered within <2 frames 
no_doubles = floor(1.9*frame_duration);
for i=1:length(red_signal)-no_doubles
    if red_signal(i)==1
        red_signal(i+1:i+no_doubles)=0;
    end
end
%compute the time between neighboring red frames
rftiming=find(red_signal);
rftiming2 = rftiming; 
drf=rftiming(2:end)-rftiming(1:end-1);
%find large gaps between red frames, they are gaps between stimuli

ri = find(drf>10000); 
n_gaps=length(ri);


% % % % % % % % % % % % % % % % % % DEBUG: 

  figure, plot(red_signal);
  for i = 8
    hold on, plot([rftiming(ri(i)) rftiming(ri(i))], [0, 1],'g');
    hold on, plot([rftiming(ri(i+1)) rftiming(ri(i+1))], [0, 1],'r');
  end 


%% PROCESS STIMULI
% For each will need: red_signal , istart, iend, logfile.  


%% If need to read logs:

%     filename=loom_logs(idx(i_log)).name;   
%     fullfilename_log=fullfile(folder_to_analyze, filename);
%     [nrep, ifi, tloom, tbtw, Speed, MaxR, Pos, loomCol, bgCol] = read_loom_log(fullfilename_log);
    
 %   get the folder and filename from the full path
%     [folder_to_analyze, log_name, ext]  = fileparts(fullfilename_log);


%% Process MULTILOOM
i_log = 2;
istart =  ri(1)+1;
iend = ri(2);
% Must read log file of RF gap to find out how many frames were presented
% in the gap. 
nfr_gap = 181; 

% % % % % % % % % % % % % % % % % % % 
process_multiloom(folder_to_analyze, rftiming, drf, istart, iend, nfr_gap); 


%% Process HC

%     i_log = 7; 
    istart = ri(3)+1;
    iend = ri(8);
    ri2 = ri(3:8); 
    
 process_HC_LC(folder_to_analyze, rftiming, drf, istart, iend, ri2, "HC");    
      
%      figure, plot(red_signal); %     
%      hold on, plot([rftiming(ri(3)) rftiming(ri(3))], [0, 1],'r');
%      hold on, plot([rftiming(ri(4)) rftiming(ri(4))], [0, 1],'r');
%      hold on, plot([rftiming(ri(5)) rftiming(ri(5))], [0, 1],'r');
%      hold on, plot([rftiming(ri(6)) rftiming(ri(6))], [0, 1],'r');
%      hold on, plot([rftiming(ri(7)) rftiming(ri(7))], [0, 1],'r');
%      hold on, plot([rftiming(ri(8)) rftiming(ri(8))], [0, 1],'r');
%      hold on, plot([rftiming(ri(9)) rftiming(ri(9))], [0, 1],'r');
    
%% Process LC

%     i_log = 24; 
    istart = ri(10)+1;
    iend = ri(15);
    ri2 = ri(10:15); 

 process_HC_LC(folder_to_analyze, rftiming, drf, istart, iend, ri2, "LC");    
      
%      figure, plot(red_signal); %     
%      hold on, plot([rftiming(ri(10)) rftiming(ri(10))], [0, 1],'r');
%      hold on, plot([rftiming(ri(11)) rftiming(ri(11))], [0, 1],'r');
%      hold on, plot([rftiming(ri(12)) rftiming(ri(12))], [0, 1],'r');
%      hold on, plot([rftiming(ri(13)) rftiming(ri(13))], [0, 1],'r');
%      hold on, plot([rftiming(ri(14)) rftiming(ri(14))], [0, 1],'r');
%      hold on, plot([rftiming(ri(15)) rftiming(ri(15))], [0, 1],'r');

%% Process DS

%     i_log = 41;
    istart = ri(17)+1;
    iend = ri(22);
    ri2 = ri(17:22);
    
 process_HC_LC(folder_to_analyze, rftiming, drf, istart, iend, ri2, "DS");    
     

%% Process Diff Pos

% i_log = 58;
istart =  ri(25)+1;
iend = ri(26);
ri2 = ri(25:26); 
% 
%  figure, plot(red_signal); %     
%      hold on, plot([rftiming(ri(25)) rftiming(ri(25))], [0, 1],'r');
%      hold on, plot([rftiming(ri(26)) rftiming(ri(26))], [0, 1],'r');

process_diff_positions(folder_to_analyze, rftiming, drf, istart, iend, ri2); 


%% Process Flashes

% i_log = 77;
istart =  ri(27);
iend = ri(31);
ri2 = ri(27:31); 

process_flashes(folder_to_analyze, rftiming, drf, istart, iend, ri2); 


%% Process OPTO

% i_log = 88;
istart =  ri(31);
iend = ri(35);
ri2 = ri(31:35); 

%      figure, plot(red_signal); %     
%      hold on, plot([rftiming(ri(31)) rftiming(ri(31))], [0, 1],'r');
%      hold on, plot([rftiming(ri(32)) rftiming(ri(32))], [0, 1],'r');
%      hold on, plot([rftiming(ri(33)) rftiming(ri(33))], [0, 1],'r');
%      hold on, plot([rftiming(ri(34)) rftiming(ri(34))], [0, 1],'r');
%      hold on, plot([rftiming(ri(35)) rftiming(ri(35))], [0, 1],'r');

 process_opto(folder_to_analyze, rftiming, drf, istart, iend, ri2);


%%


