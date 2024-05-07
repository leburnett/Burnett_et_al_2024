function process_opto(folder_to_analyze, rftiming, drf, istart, iend, ri2)

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
    
    %% 
    nrep = 4; 
   
    % count how many frames within each repetition
    bout_start = ri2(1:end-1)+1;
    bout_end = ri2(2:end);
     
    %number of red frames per rep
    nrf_rep = floor(min(rftiming(bout_end)-rftiming(bout_start))/ifired);

    %number of red frames per rep
    RFnum=nrf_rep; %-1

    %number of all frames
    allframesnum=RFnum*5;

    %all spike counts, per cluster, per rep, binned by one frame length
    all_spikes=zeros(n_gcl,nrep,allframesnum-1);

        %iterate thru all repetitiions
    for rid=1:nrep
        %same length for all repetitions
        bout_endi=bout_start(rid)+RFnum;
        rftiming_i=rftiming(bout_start(rid):bout_endi);
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
            all_spikes(ci,rid,:)=N;        
        end
    end

        L1_hist = squeeze(all_spikes(:,1,1:594));
        R1_hist = squeeze(all_spikes(:,2,1:594));
        L2_hist = squeeze(all_spikes(:,3,1:594));
        R2_hist = squeeze(all_spikes(:,4,1:594));

        
        figure; imagesc(L1_hist); caxis([0 5])
        figure; imagesc(R1_hist); caxis([0 5])
        figure; imagesc(L2_hist); caxis([0 5])
        figure; imagesc(R2_hist); caxis([0 5])

% save('Spikes_Opto.mat', 'all_spikes', 'L1_hist', 'L2_hist', 'R1_hist', 'R2_hist'); 


end 