function analyze_checkers(red_signal, istart, iend, logfile)
%this file contains a set of routines to analyse the data recorded with a
%silicone probe during the presentation of the moving checkers stimuli
%
%Date: 15.09.2020
% O. Symonova
figure, plot(red_signal);

istart = 3100000;
iend = 3170000;

%get the folder and filename from the full path
[folder_to_analyze, log_name,ext]  = fileparts(logfile);

%collect all needed files automatically
spike_cluster_file= fullfile(folder_to_analyze,'spike_clusters.npy');
spike_time_file= fullfile(folder_to_analyze,'spike_times.npy');
cluster_info_file= fullfile(folder_to_analyze,'cluster_group.tsv');
cluster_info_gen_file= fullfile(folder_to_analyze,'cluster_info.tsv');
res_folder= fullfile(folder_to_analyze,'results');

%make the folder with results
if exist(res_folder,'dir')==0
    mkdir(res_folder);
end    

%local copy of red frames relevant to the stimulus
red_signal_i=red_signal(istart:iend);

%for debugging plot red frames and start and end of the stimulus
% debug=1;
% if debug
%     figure, plot(red_signal);
%     hold on, plot([istart istart], [0, 1],'g');
%     hold on, plot([iend iend], [0, 1],'r');
% end

%clean red frames and get time stamp of each frame
[frame_timing, red_signal_i] = get_clean_frame_timing(red_signal_i, log_name);
dt=frame_timing(2:end)-frame_timing(1:end-1);
ifi=median(dt);
rftiming=frame_timing(1:5:end);

%reconstruct stimulus
[frame_array,npatt,checker_size]=reconstruct_checker_stimuli(logfile);
[imheight, imwidth, nfr]=size(frame_array);

%frame_history_window
frame_history=30;

%cluster label and time of each spike
cluster_id = readNPY(spike_cluster_file);
spike_t = readNPY(spike_time_file);

%find ids of good clusters
goodclusters = get_clusters_by_property(cluster_info_file,'good');
muaclusters = get_clusters_by_property(cluster_info_file,'mua');
goodclusters=[goodclusters,muaclusters];

n_gcl=length(goodclusters);
tstart=frame_timing(frame_history+1)+istart-1;
tend=frame_timing(end)+istart-1;

%first save all RF to visualize them in the same scale
RF=zeros(n_gcl,imheight, imwidth,frame_history);

%parallel pool
poolobj =parpool(4);

%for each cluster find correlation with the stimulus
parfor ci=1:n_gcl
    idc=goodclusters(ci);
    %select spikes per cluster
    spikes=find(cluster_id==idc);
    
    %get the time of spikes
    ts=spike_t(spikes);  
    
    %take only spikes after the start of the stimulus+frame_history
    ts=ts(ts>tstart);
    ts=ts(ts<tend);
    ts=ts-istart+1;
    
    numspikes = length(ts);     
    
    %get the RF and dynamic filter DF
    %get the frame number corresponding to each spike
    %compute how many red frames happened before each spike
    nredfr_before_spike=arrayfun(@(i) sum(red_signal_i(1:ts(i))),1:numspikes);
    %time of the last red frame 
    trf=rftiming(nredfr_before_spike)';
    %how many frames after the last red frame, round down
    %rarely there is a delay between red frames, make sure to not take more
    %than 5 frames in between
    nfr_after_lr=min(floor((double(ts)-trf)./ifi),4);  
    %frame id at the time when spike happened
    fri=((nredfr_before_spike-1)*5)'+1+nfr_after_lr;
    %frame id at which frame history stats
    fri_start_frame_history=fri-frame_history+1;

    %sum up frames corresponding to each spike    
    RFi=zeros(imheight, imwidth,frame_history);
    for is=1:numspikes
        inds=fri_start_frame_history(is):fri(is);
        npatinds=floor(inds/npatt)+1;
        RFi=RFi+double(frame_array(:,:,npatinds));
    end
    RF(ci,:,:,:)=RFi./numspikes;      
end

%flip array to match mouse's view
parfor ci=1:n_gcl
    RFtemp=squeeze(RF(ci,:,:,:));
    %PTB flips images up-down, different coordinate frame
    RFtemp=flipud(RFtemp);
    %mouse towards the projector, hence left-right are flipped
    RFtemp=fliplr(RFtemp);
    RF(ci,:,:,:)=RFtemp;    
end

% save('RF_all_temp.mat', 'goodclusters','RF', '-v7.3');

%get depth of the clusters
ids_depth = get_clusters_depth(cluster_info_gen_file);

%save variables for the RF analysis
%the indices in RF coincide with indices in 'goodclusters'
% cluster_id can be recovered from goodclusters (i.e. goodclusters)
% ids_depth is the pairwise array of cluster_id and depth
matRF = fullfile(res_folder,'RF_all.mat');
save(matRF, 'goodclusters','RF', 'ids_depth', '-v7.3');

%make crop size be twice the size of the checker
disp(['Checker size =',num2str(checker_size)]);
crop_size_half=checker_size*4;
full_crop_size=crop_size_half*2+1;
%array to visualize the summary of RF (cropped and only 1 frame)
RFsummary=zeros(n_gcl,full_crop_size, full_crop_size);
parfor ci=1:n_gcl
    RFi=squeeze(RF(ci,:,:,:));
    vari=var(RFi,0,3);
    [mval,cx]=max(max(vari,[],1));
    [mval,cy]=max(max(vari,[],2));    
    %disp([ci,mval]);
    meani=mean(RFi(:));
    %find the largest deviation at the RF in time
    [mval,ct]=max(abs(squeeze(RFi(cy,cx,:))-meani));
    RFi_tocrop=zeros(crop_size_half*2+imheight,crop_size_half*2+imwidth);
    RFi_tocrop(crop_size_half+1:crop_size_half+imheight,crop_size_half+1:crop_size_half+imwidth)=RFi(:,:,ct);
    RFsummary(ci,:,:)=RFi_tocrop(cy:cy+2*crop_size_half,cx:cx+2*crop_size_half);
end

%sort clusters by depth
ids_depth_sorted=sortrows(ids_depth,2);

nfigcol=11;
nfigrow=ceil(n_gcl/nfigcol);
figname = [log_name(1:end-4),'_RFsummary.fig'];
figure('name',figname,'units','normalized','outerposition',[0 0 1 1]);
spi=1;
%make the summary figure, where RF are sorted by depth
for i=1:length(ids_depth_sorted)    
    cluster_index=ids_depth_sorted(i,1);    
    cluster_RF_index=find(goodclusters==cluster_index);
    %if the cluster id is not 'good' cluster
    if isempty(cluster_RF_index)
        continue;
    end
    %depth of the cluster
    depthi=ids_depth_sorted(i,2);
    %create subplot
    subplot(nfigrow,nfigcol,spi);  
    imagesc(squeeze(RFsummary(cluster_RF_index,:,:)));
    axis off;  
    colormap(redblue)
    titlestr = ['c=',num2str(cluster_index), ' d=', num2str(depthi)];%, '\mum'];
    title(titlestr);  
    spi=spi+1;
end
%save figure
figname=fullfile(res_folder,figname);
savefig(figname);

%save videos of each RF
parfor ci=1:n_gcl
    idc=goodclusters(ci);
    iddi = find(ids_depth(:,1)==idc);
    depthi = ids_depth(iddi,2);        
    avi_file_RF=['cluster',num2str(idc),'_depth', num2str(depthi),'_checker_RF.avi'];
    avi_file_RF=fullfile(res_folder,avi_file_RF);    
    make_video_RF(squeeze(RF(ci,:,:,:)), avi_file_RF,[],[],0);  
end

delete(poolobj);

end
%%

%% read the file and extract cluster id and depth info
function ids_depth = get_clusters_depth(cluster_info_gen_file)
    f=fopen(cluster_info_gen_file);
    %first line is the header
    t=fgetl(f);
    ids_depth=[];
    while ~feof(f)
        t=fgetl(f);
        strline=strsplit(t);
            id=str2num(strline{1});
            depth=str2num(strline{7});
            ids_depth=[ids_depth; [id,depth]];        
    end
    fclose(f);
end

%% this function plots the average intensity at a specific location of the
%stimulus during several frames before a spike
function plot_dynamic_filters(DFlist,goodclusters,RF)
    clusters2plot = [DFlist(1:end).cluster_id];
    n_gcl=length(goodclusters);
    for ci=1:n_gcl
        idc=goodclusters(ci);
        ind = find(clusters2plot==idc);
        if isempty(ind)
            continue;
        end
        start_row = DFlist(ind).start_row;
        end_row = DFlist(ind).end_row;
        start_col = DFlist(ind).start_col;
        end_col = DFlist(ind).end_col;
        
        filterfile=['cluster',num2str(idc),'_checker_DF.fig'];        
        filterfile_png=['cluster',num2str(idc),'_checker_DF.png'];
        traces=squeeze(RF(ci,start_row:end_row,start_col:end_col,:));
        nrow = size(traces,1);
        ncol = size(traces,2);
        nframes = size(traces,3);
        traces=reshape(traces,[nrow*ncol,nframes]);
        
        legend_labels=cell(nrow*ncol,1);
   
        i=1;
        for ri=start_row:end_row
            for ci=start_col:end_col               
                legend_labels(i)={[num2str(ri),'_',num2str(ci)]};
                i=i+1;
            end
        end
            
        fig=figure;        
        plot(traces')
        legend(legend_labels,'Location','eastoutside','Interpreter', 'none');
           
        savefig(filterfile);
        saveas(fig,filterfile_png,'png');  
        close(fig);
    end
end

%% removes red frames registered multiple times and inserts missed frames
function [frame_timing, red_signal_local] = get_clean_frame_timing(red_signal_local, filename)
    
    rftiming=find(red_signal_local);

    %compute the time between neighboring red frames
    drf=rftiming(2:end)-rftiming(1:end-1);
    %intra red frame interval
    ifired=median(drf);
    
    %find red frames with 120ms in between (instead of 83)
    big_gap_duration = 1.5*ifired;
    missed_rf=find(drf>big_gap_duration);     
    %interval in sec*sampling_rate between  frames
    ifi=ifired/5;
    %check if there is a very large gap: several consecutive frames have been missed, results might not be reliable 
    large_missed_rf=find(drf>ifired*2.2); %if more consecutive frames missing issue warning
    if ~isempty(large_missed_rf)
        warning([filename, ': more than one consecutive frame is missing. The results might be wrong.']);    
    end
    if isempty(missed_rf)
        nummissed=0;
    else
        nummissed=length(missed_rf);
    end
    disp([num2str(nummissed),' frames were missed']);
    %insert the missed redfraemes
    while ~isempty(missed_rf)
        % the time where the red frame has to be
        mrf=rftiming(missed_rf)+ifired;
        red_signal_local(mrf)=1;
        rftiming=find(red_signal_local);
        drf=rftiming(2:end)-rftiming(1:end-1);
        missed_rf=find(drf>big_gap_duration); 
    end

    %number of red frames (we do not consider the last red frame and frames after it)
    RFnum=length(rftiming)-1;
    %number of all frames
    allframesnum=RFnum*5;
    %timing of all frames, computed using linear interpolation of time of red
    %frames
    frame_timing=interp1(1:5:allframesnum,rftiming(1:end-1),1:allframesnum,'linear','extrap');
end

%% read the file that encodes the cluster id and associated property
%  returns the indices of the clusters with the given property 'propstr'
function ids = get_clusters_by_property(filename, propstr)
    f=fopen(filename);
    %first line is the header
    t=fgetl(f);
    ids=[];
    while ~feof(f)
        t=fgetl(f);
        strline=strsplit(t);
        if strcmp(strline{2},propstr)
            id=str2num(strline{1});
            ids=[ids, id];
        end
    end
    fclose(f);
end


