function process_flashes(folder_to_analyze, rftiming, drf, istart, iend, ri2)

% Four different speeds
% 1.5s, 1s, 0.5s, 0.25s

% Frames presented (905, 1209, 1215, 1231)

exp_name = folder_to_analyze(end-14:end);

%collect all needed files automatically, they should be inside folder_to_analyze
spike_cluster_file= fullfile(folder_to_analyze,'spike_clusters.npy');
spike_time_file= fullfile(folder_to_analyze,'spike_times.npy');
cluster_info_file= fullfile(folder_to_analyze,'cluster_group.tsv');
cluster_info_gen_file= fullfile(folder_to_analyze,'cluster_info.tsv');

main_dir = '/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe'; 

%for this session add folder with the functions to read npy files
addpath(fullfile(main_dir,'npy-matlab'))
addpath('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Scripts/PreProcess/Looms');

% intra red  frame interval
ifired=median(drf);
% interval in sec*sampling_rate between 2 frames
ifi=ifired/5;

%get all good and mua clusters
    cluster_id = readNPY(spike_cluster_file);
    
    spike_t = readNPY(spike_time_file);
    goodclusters = get_clusters_by_property(cluster_info_file,'good');
    muaclusters = get_clusters_by_property(cluster_info_file,'mua');
    %join mua and good clusters
    goodclusters=[goodclusters,muaclusters];
    %get depth info of the clusters
    ids_depth = get_clusters_depth(cluster_info_gen_file);
    %number of all clusters to analyse
    n_gcl=length(goodclusters);
    

%% For each different speed of flash create a separate array. 

 % Beginning and end of each bout
    bout_start = ri2(1:end-1)+1;
    bout_end = ri2(2:end);

% 4 different speeds
    
    %% Speed = 1.5s
    
    i = 1; 
    
    fl_start = bout_start(i);
    fl_end = bout_end(i); 
    
    % Number of rfs. 
    n_rf = fl_end-fl_start; 
    
    %number of all frames
    allframesnum= n_rf*5;
    
    %all spike counts, per cluster, per rep, binned by one frame length
    all_spikes=zeros(n_gcl,allframesnum-1);
    
    rftiming_i=rftiming(fl_start:fl_end);
    
    %timing of all frames, using linear interpolation of time of red frames
    frame_timing_i=interp1(1:5:allframesnum,rftiming_i(1:end-1),1:allframesnum,'linear','extrap');

    %get spikes for all clusters
    for ci=1:size(ids_depth,1)
        %id of the cluster
        idc=ids_depth(ci,1);
        %select spikes of the current cluster
        spikes=find(cluster_id==idc);
        %get the time of spikes
        ts=spike_t(spikes);
        %bin by frame
        [N,edges] = histcounts(ts,frame_timing_i+istart-1);
        all_spikes(ci,:)=N;
    end
  
    all_spikes15 = all_spikes; 
        
%     figure; imagesc(all_spikes15); caxis([0 5])
%     hold on 
%     val =15; 
%     for i = 1:10
%         plot([val val], [-100 200], 'w')
%         val=val+90;
%     end 
%     
    
    %%  Speed = 1s
    
     i = 2; 
    
    fl_start = bout_start(i);
    fl_end = bout_end(i); 
    
    % Number of rfs. 
    n_rf = fl_end-fl_start; 
    
    %number of all frames
    allframesnum= n_rf*5;
    
    %all spike counts, per cluster, per rep, binned by one frame length
    all_spikes=zeros(n_gcl,allframesnum-1);
    
    rftiming_i=rftiming(fl_start:fl_end);
    
    %timing of all frames, using linear interpolation of time of red frames
    frame_timing_i=interp1(1:5:allframesnum,rftiming_i(1:end-1),1:allframesnum,'linear','extrap');

    %get spikes for all clusters
    for ci=1:size(ids_depth,1)
        %id of the cluster
        idc=ids_depth(ci,1);
        %select spikes of the current cluster
        spikes=find(cluster_id==idc);
        %get the time of spikes
        ts=spike_t(spikes);
        %bin by frame
        [N,edges] = histcounts(ts,frame_timing_i+istart-1);
        all_spikes(ci,:)=N;
    end
  
    all_spikes1 = all_spikes; 
    
    figure; imagesc(all_spikes1); caxis([0 5])
    hold on 
    val =0; 
    for i = 1:20
        plot([val val], [-100 200], 'w')
        val=val+60;
    end 
    
    %% Speed = 0.5s
    
     i = 3; 
    
    fl_start = bout_start(i);
    fl_end = bout_end(i); 
    
    % Number of rfs. 
    n_rf = fl_end-fl_start; 
    
    %number of all frames
    allframesnum= n_rf*5;
    
    %all spike counts, per cluster, per rep, binned by one frame length
    all_spikes=zeros(n_gcl,allframesnum-1);
    
    rftiming_i=rftiming(fl_start:fl_end);
    
    %timing of all frames, using linear interpolation of time of red frames
    frame_timing_i=interp1(1:5:allframesnum,rftiming_i(1:end-1),1:allframesnum,'linear','extrap');

    %get spikes for all clusters
    for ci=1:size(ids_depth,1)
        %id of the cluster
        idc=ids_depth(ci,1);
        %select spikes of the current cluster
        spikes=find(cluster_id==idc);
        %get the time of spikes
        ts=spike_t(spikes);
        %bin by frame
        [N,edges] = histcounts(ts,frame_timing_i+istart-1);
        all_spikes(ci,:)=N;
    end
  
    all_spikes05 = all_spikes; 
    
    figure; imagesc(all_spikes05); caxis([0 5])
    hold on 
    val =15; 
    for i = 1:40
        plot([val val], [-100 200], 'w')
        val=val+30;
    end 
    
    %% Speed = 0.25s
    
    i = 4; 
    
    fl_start = bout_start(i);
    fl_end = bout_end(i); 
    
    % Number of rfs. 
    n_rf = fl_end-fl_start; 
    
    %number of all frames
    allframesnum= n_rf*5;
    
    %all spike counts, per cluster, per rep, binned by one frame length
    all_spikes=zeros(n_gcl,allframesnum-1);
    
    rftiming_i=rftiming(fl_start:fl_end);
    
    %timing of all frames, using linear interpolation of time of red frames
    frame_timing_i=interp1(1:5:allframesnum,rftiming_i(1:end-1),1:allframesnum,'linear','extrap');

    %get spikes for all clusters
    for ci=1:size(ids_depth,1)
        %id of the cluster
        idc=ids_depth(ci,1);
        %select spikes of the current cluster
        spikes=find(cluster_id==idc);
        %get the time of spikes
        ts=spike_t(spikes);
        %bin by frame
        [N,edges] = histcounts(ts,frame_timing_i+istart-1);
        all_spikes(ci,:)=N;
    end
  
    all_spikes025 = all_spikes(:, 1:1250); 
       
    figure; imagesc(all_spikes025); caxis([0 5])
    hold on 
    val =7.8; 
    for i = 1:80
        plot([val val], [-100 200], 'w')
        q = i/2; 
        if isinteger(q)
           val= val+15.3;
        else
            val = val+15; 
        end 
    end 
    % still seems to drift. 
    
    %% 
    
   % save('Spikes_Flashes.mat', 'all_spikes15', 'all_spikes1', 'all_spikes05', 'all_spikes025')
end 



