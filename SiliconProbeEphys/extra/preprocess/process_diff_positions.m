function process_diff_positions(folder_to_analyze, rftiming, drf, istart, iend, ri2)

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

% DEBUGGING. 
%  figure, plot(red_signal_local); %
%      hold on, plot([rftiming(istart) rftiming(istart)], [0, 1],'g'); %   should be the first RF of the stimulus. 
%      hold on, plot([rftiming(iend) rftiming(iend)], [0, 1],'r'); %    should be the last RF of the stimulus. 

ifired=median(drf);

%number of red frames 
nrf = floor(min(rftiming(iend)-rftiming(istart))/ifired);
RFnum = nrf-1; %also iend-istart;
allframesnum = (RFnum*5);

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

%% Finding the timing of the frames during the stimulus using the timing of the RFs. 

rftiming_i=rftiming(istart:iend);
%timing of all frames, using linear interpolation of time of red frames
frame_timing_i=interp1(1:5:allframesnum,rftiming_i(1:end),1:allframesnum,'linear','extrap');

%% Create array of spiking over stimulus + RF gaps either side. 
 all_spikes=zeros(n_gcl,allframesnum-1);

%get spikes for all clusters
for ci=1:size(ids_depth,1)
    %id of the cluster
    idc=ids_depth(ci,1);
    %select spikes of the current cluster
    spikes=find(cluster_id==idc);
    %get the time of spikes
    ts=spike_t(spikes);
    %bin by frame
    [N,edges] = histcounts(ts,frame_timing_i+istart-1);  %istart-1
    all_spikes(ci,:)=N;
end

all_spikes_diffpos = all_spikes;

% Plot figure. 
figure; imagesc(all_spikes_diffpos); caxis([0 5])
hold on
curr_pos = 0; 
for i = 1:150
    plot([curr_pos curr_pos], [-100 200], 'w') %end of RF gap
    curr_pos = curr_pos+45; 
end 


%%
% all_spikes_diffpos %81*6889

data_array = zeros(n_gcl, 150); 

ind_to_choose = 1:46:6897;
for j = 1:150
    row = ind_to_choose(j);
    XY_array(j,1) = vbl_array(row,7);
    XY_array(j,2) = vbl_array(row,8);
end
   
FullScreen = zeros(10,15);
NormActAll = zeros(n_gcl, 10, 15);

for i = 1:n_gcl
    dt = all_spikes_diffpos(i, :);
    dt(end:6900) = 0;
    data_array(:,i) = reshape(sum(reshape(dt,[],150)),[],1);
    
    
    Activity = zeros(10,15);
    NumF = zeros(10,15);
    
    for j = 1:150
        Activity(XY_array(j, 1), XY_array(j,2)) = Activity(XY_array(j, 1), XY_array(j,2)) + data_array(j,i);
        NumF(XY_array(j, 1), XY_array(j,2)) = NumF(XY_array(j, 1), XY_array(j,2)) +1;
    end
    
    NormAct = zeros(10,15);
    NormAct = Activity./NumF;
    
    NormActAll(i, :, :) = NormAct; 
%     figure; imagesc(NormAct); colorbar
    FullScreen = FullScreen+NormAct; 
end

FullScreen = FullScreen/n_gcl; 
figure; imagesc(FullScreen); colorbar
  
save('Spikes_DiffPos.mat', 'NormActAll', 'FullScreen', 'all_spikes_diffpos');


end 
