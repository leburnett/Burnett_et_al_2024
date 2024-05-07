%the script loads the recording which was done during presentations of one
%or more stimuli. The start and end time points of each stimulus will be
%found and used in the analysis of the stimulus. It it assumed that stimuli are separated by relatively large pause
% 
% %Olga Symonova
% %22.09.2020

% folder_to_analyze = 'C:\DATA\EPhys\silicone_probe\Sorting280720_1479_Dreadds(-)_1\CNO23min_LISFFC';
folder_to_analyze = 'C:\Data_analysis\Silicon_Probe_Analysis\Example_Data\200925_GN7475_Test1_1250';

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

%for this session add folder with the functions to read npy files
addpath(fullfile(pwd,'npy-matlab'));
%for this session add folder with the functions to analyse specific stimuli
addpath(fullfile(pwd,'ff'));
addpath(fullfile(pwd,'checkers'));

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
%sort log files by date
[~,idx] = sort([log_files.datenum]);

%find large gaps in the red signal that are the gaps between stimuli
%presentations or repetitions of the chirp
%get the syncing signal
[red_signal, sample_rate] = read_sync_signal_from_Intan_RHD2000_file(datfile);
figure, plot(red_signal);

%each red frame is sampled multiple times
sec_per_frame=0.0167;
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
drf=rftiming(2:end)-rftiming(1:end-1);
%find large gaps between red frames, they are gaps between stimuli
ri = find(drf>100000); %Larger than 5 s. 
n_gaps=length(ri);

%save('temp.mat');
istart=1;
iend=ri(1);
igap=0;
nrep_gap=0;

%iterate thru stimuli sorted by time
for i= 3 %1:length(idx)
    filename=log_files(idx(i)).name;
    % if chirp stimulus, get the number of repetitions
    if  contains(filename,'full_field_curve','IgnoreCase', true)  
        fullfilename_log=fullfile(folder_to_analyze, filename);
        [nrep, ifi, tbefore, tafter, bgIntensity, channels, nframes] = read_full_field_curve_log(fullfilename_log);
        nrep_gap=nrep-1;
    end

    %get the first red frame of the stimulus
    if igap==0 
        istart=rftiming(1);
    else
        istart=rftiming(ri(igap)+1);
    end
    %move the current gap and get the last red frame of the stimulus
    igap=igap+nrep_gap+1;
    nrep_gap=0;
    if igap>n_gaps
        iend=rftiming(end);        
    else
        iend=rftiming(ri(igap));        
    end
    
   %analyse chirp signal
    if  contains(filename,'full_field_curve','IgnoreCase', true)  
        fullfilename_log=fullfile(folder_to_analyze, filename);   
%         analyze_chirp_stimuli(red_signal, istart, iend, datfile, fullfilename_log);
    %analyse checker stimuli
    elseif contains(filename,'checker','IgnoreCase', true)  
        fullfilename_log=fullfile(folder_to_analyze, filename);        
        analyze_checkers(red_signal, istart, iend, fullfilename_log);
    %analyse linear_intensity_steps stimulus
    elseif contains(filename,'linear_intensity_steps','IgnoreCase', true)  
        fullfilename_log=fullfile(folder_to_analyze, filename);              
        %this function does not exist yet %analyze_linear_intensity_steps(red_signal, istart, iend, datfile, fullfilename_log);
    %analyse moving_bar stimulus
    elseif contains(filename,'moving_bar_var_intensities','IgnoreCase', true)  
        fullfilename_log=fullfile(folder_to_analyze, filename);          
        %this function does not exist yet %analyze_moving_bar_var_intensities(red_signal, istart, iend, datfile, fullfilename_log);
    elseif contains(filename,'Loom','IgnoreCase', true)  
        fullfilename_log=fullfile(folder_to_analyze, filename);
        analyze_loom_stimuli(red_signal, istart, iend, fullfilename_log)
        %this function does not exist yet %analyze_moving_bar_var_intensities(red_signal, istart, iend, datfile, fullfilename_log);
    else
        disp('Did not find matching analysis script');
    end     
    
end



%% check that all needed files exist in the given folder, issue a warning msg othewise
function datfile = check_files_for_analysis(folder_to_analyze, spike_cluster_file, spike_time_file, cluster_info_file)
    %check if the files exist in the folder
    datfile='';
    logfile='';
    
    if ~exist(spike_cluster_file,'file')
        disp('File spike_clusters.npy is not in the folder. Exitting...');
        return;
    end
    
    if ~exist(spike_time_file,'file')
        disp('File spike_times.npy is not in the folder. Exitting...');
        return;
    end
    
     if ~exist(cluster_info_file,'file')
        disp('File cluster_info.tsv is not in the folder. Exitting...');
        return;
     end
    
     %find dat and log files
     d = dir(folder_to_analyze);
     nfiles = length(d);
     
     for fi=3:nfiles
         filename=d(fi).name;
         if endsWith(filename,'.dat')
             datfile=fullfile(folder_to_analyze,filename);
         end
%          if endsWith(filename,'.log')
%              logfile=fullfile(folder_to_analyze,filename);
%          end
     end
     if ~exist(datfile,'file')
        disp('Could not find *.dat file in the folder. Exitting...');
        return;
     end
%      if ~exist(logfile,'var')
%         disp('Could not find *.log file in the folder. Exitting...');
%         return;
%      end
end
%%

