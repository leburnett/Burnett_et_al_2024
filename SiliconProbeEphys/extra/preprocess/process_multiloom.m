function process_multiloom(folder_to_analyze, rftiming, drf, istart, iend, nfr_gap)

% i_log = 2;
% istart =  ri(1)+1;
% iend = ri(2);

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
RFnum = nrf; %also iend-istart;
allframesnum = RFnum*5;

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
frame_timing_i=interp1(1:5:allframesnum,rftiming_i(1:end-1),1:allframesnum,'linear','extrap');

%% Finding the timing of the frames during the RF gaps before and after the stimulus. 

% RF gap PRE
timingpregap = rftiming(istart-1);
timingpostgap =rftiming(istart);

%frame timing for the 3s pre gap.
frame_timing_i_pre = interp1([1,nfr_gap], [timingpregap, timingpostgap], 1:1:nfr_gap);
frame_timing_i_pre = frame_timing_i_pre(1:nfr_gap-1); % 181 is istart therefore = frame_timing_i(1);
%frame timing for the 3s post gap.

timingpregap = rftiming(iend);
timingpostgap =rftiming(iend+1);

%frame timing for the 3s pre gap.
frame_timing_i_post = interp1([1,nfr_gap], [timingpregap, timingpostgap], 1:1:nfr_gap);

%% Create array of spiking over stimulus + RF gaps either side. 

allframes_withgaps = allframesnum+(nfr_gap-1)+nfr_gap;
all_spikes=zeros(n_gcl,allframes_withgaps-1);
frame_timing_all = horzcat(frame_timing_i_pre, frame_timing_i, frame_timing_i_post);

%get spikes for all clusters
for ci=1:size(ids_depth,1)
    %id of the cluster
    idc=ids_depth(ci,1);
    %select spikes of the current cluster
    spikes=find(cluster_id==idc);
    %get the time of spikes
    ts=spike_t(spikes);
    %bin by frame
    %             [N,edges] = histcounts(ts,frame_timing_i+istart-1);  %istart-1
    [N,edges] = histcounts(ts,frame_timing_all+(istart-1));
    all_spikes(ci,:)=N;
end
%      end

all_spikes_multiloom = all_spikes;

% Plot figure. 
figure; imagesc(all_spikes_multiloom); caxis([0 5])
hold on
plot([180 180], [-100 200], 'w') %end of RF gap

plot([225 225], [-100 200], 'w')
plot([270 270], [-100 200], 'w')
plot([315 315], [-100 200], 'w')
plot([360 360], [-100 200], 'w')
plot([405 405], [-100 200], 'w')
plot([450 450], [-100 200], 'w')
plot([495 495], [-100 200], 'w')
plot([540 540], [-100 200], 'w')
plot([585 585], [-100 200], 'w')
plot([630 630], [-100 200], 'w')
plot([675 675], [-100 200], 'w')
plot([720 720], [-100 200], 'w')
plot([765 765], [-100 200], 'w')
plot([810 810], [-100 200], 'w')
plot([855 855], [-100 200], 'w')

end 
