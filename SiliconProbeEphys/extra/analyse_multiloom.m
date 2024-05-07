%% Analyse firing of cells to looming stimulus - multiloom
% Based on 'make_firing_array_loom_stimulus.m' but for newly processed data 2021. 

% Made 31/03/21 by Burnett. 

% Averaged across cells - WT / HET 
% Last columns contain, depth of cluster, animal # and geno. 

% To be used after 'combine_and_add_depth_Multiloom.m' - 


%%  1-  Load '210323_Spiking_Looms_HC_L1.mat' - 'all_hist' and 'all_av_hist'

W = 3515; 
l = 46;

% Find basic facts about data. 
n_cells = numel(all_hist(:,1));

% For the average across the 5 bouts
allWT = find(all_hist(:,W+4)==1); 
allHET = find(all_hist(:,W+4)==0); 

%%
nWT = numel(allWT);
nHET = numel(allHET);

depthsWT = all_hist(allWT, W+3); 
depthsHET = all_hist(allHET, W+3); 

avDWT = mean(depthsWT);
avDHET = mean(depthsHET);

%% FINDING RESPONSIVE CELLS. 
% Should find the cells that respond to the loom - % of cells found that are responsive to the loom in WT/HET and plot the depths of this subpopulation cells. 

% Look at the response of each cell from 0 -60 and from 61 - 100

resp_IND = zeros(n_cells,4); 

for i = 1:n_cells
    
    resp_b4 = mean(all_hist(i, 1:60));
    
%      if resp_b4 > 0.25
        
        resp_l1 = mean(all_hist(i, 61:60+l));
        delta_resp = resp_l1-resp_b4;
        
        resp_IND(i,1) = resp_b4;
        resp_IND(i,2) = resp_l1;
        resp_IND(i,3) = delta_resp;
        
        if resp_b4 == 0 && resp_l1 == 0
            resp_IND(i,4) = NaN;
            all_hist(i, W+5) = NaN;
            all_hist(i, W+6) = NaN;
        else
            val = (resp_l1-resp_b4)/(resp_l1+resp_b4);
            if val ~=Inf
                resp_IND(i,4) = val;
            else
                resp_IND(i,4) = NaN;
            end
            all_hist(i, W+5) = resp_IND(i,4);
            
            % Responsive or not responsive. 
            if all_hist(i, W+5) > 0
                all_hist(i, W+6) = 1;
            else
                all_hist(i, W+6) = 0;
            end

        end
        
%     else 
%         all_hist(i, 815) = NaN;
%         all_hist(i, 816) = NaN; 
%     end
    
end



% Plot the distribution of resp_index. 
%           figure; histogram(resp_IND(:,4), 200)


all_av_hist = all_hist; 

%% SORTING CELLS BASED ON RESPONSE TO LOOM STIMULUS

all_WT_RESP = find(all_av_hist(:,W+6)==1 & all_av_hist(:,W+4)==1); 
allWT_RESP_ZERO = find(all_av_hist(:,W+6)==0 & all_av_hist(:,W+4)==1); 
% % % % % %

all_HET_RESP = find(all_av_hist(:,W+6)==1 & all_av_hist(:,W+4)==0); 
allHET_RESP_ZERO = find(all_av_hist(:,W+6)==0 & all_av_hist(:,W+4)==0); 
% % % % % % %

n_up_WT = numel(all_WT_RESP); 
n_zero_WT = numel(allWT_RESP_ZERO);

n_up_HET = numel(all_HET_RESP); 
n_zero_HET = numel(allHET_RESP_ZERO);

num_vals = zeros(2,4);
num_vals(1,:) = [n_up_WT, n_zero_WT, n_up_HET, n_zero_HET]; 
num_vals(2, 1:2) = num_vals(1, 1:2)/nWT; 
num_vals(2, 3:4) = num_vals(1, 3:4)/nHET; 

% Plot the firing of the 'responsive' and 'non-responsive' cells. 

figure; imagesc(all_av_hist(allHET, 1:W)); title('HET'); caxis([0 100])
figure; imagesc(all_av_hist(allWT, 1:W)); title('WT'); caxis([0 100])

figure; imagesc(all_av_hist(all_HET_RESP, 1:W)); title('HET - RESP'); caxis([0 100])
figure; imagesc(all_av_hist(all_WT_RESP, 1:W)); title('WT - RESP'); caxis([0 100])

figure; imagesc(all_av_hist(allHET_RESP_ZERO, 1:W)); title('HET - Non RESP'); caxis([0 20])
figure; imagesc(all_av_hist(allWT_RESP_ZERO, 1:W)); title('WT - Non RESP'); caxis([0 20])

colormap(gray)
lf = 46+180;

hold on 
plot([60 60], [0 1500], 'w')
plot([60+lf 60+lf], [0 1500], 'w')
plot([60+lf*2 60+lf*2], [0 1500], 'w')
plot([60+lf*3 60+lf*3], [0 1500], 'w')
plot([60+lf*4 60+lf*4], [0 1500], 'w')
plot([60+lf*5 60+lf*5], [0 1500], 'w')
plot([60+lf*6 60+lf*6], [0 1500], 'w')
plot([60+lf*7 60+lf*7], [0 1500], 'w')
plot([60+lf*8 60+lf*8], [0 1500], 'w')
plot([60+lf*9 60+lf*9], [0 1500], 'w')
plot([60+lf*10 60+lf*10], [0 1500], 'w')
plot([60+lf*11 60+lf*11], [0 1500], 'w')
plot([60+lf*12 60+lf*12], [0 1500], 'w')
plot([60+lf*13 60+lf*13], [0 1500], 'w')
plot([60+lf*14 60+lf*14], [0 1500], 'w')
plot([60+lf*15 60+lf*15], [0 1500], 'w')
xticklabels('')
yticklabels('')
ax = gca;
ax.FontSize = 20; 

%% Plot per animal

% For the average across the 5 bouts
ani = 1971; 
allANI = find(all_av_hist(:,352)==ani); 
figure; imagesc(all_av_hist(allANI, 1:350)); title(string(ani)); caxis([0 100])

%% Sort by depth

av_hist_sorted = sortrows(all_av_hist, W+3); 
allWT_sort = find(av_hist_sorted(:,W+4)==1 & av_hist_sorted(:,W+6)==1); 
allHET_sort = find(av_hist_sorted(:,W+4)==0 & av_hist_sorted(:,W+6)==1); 

figure; imagesc(av_hist_sorted(allHET_sort, 1:W)); title('HET'); caxis([0 150])
% colormap(gray)
lf = 46+180;
hold on 
plot([60 60], [0 1500], 'w')
for i = 1:15
plot([60+lf*i 60+lf*i], [0 1500], 'w')
end 
xticklabels('')
yticklabels('')
ax = gca;
ax.FontSize = 20; 

figure; imagesc(av_hist_sorted(allWT_sort, 1:W)); title('WT'); caxis([0 150])
% colormap(gray)
lf = 46+180;
hold on 
plot([60 60], [0 1500], 'w')
for i = 1:15
plot([60+lf*i 60+lf*i], [0 1500], 'w')
end 
xticklabels('')
yticklabels('')
ax = gca;
ax.FontSize = 20; 


%% DEPTHS OF CELLS - STATS - HISTOGRAMS

% ALL CELLS
depthsWT = all_av_hist(allWT, W+3)*-1; 
depthsHET = all_av_hist(allHET, W+3)*-1;

figure; histogram(depthsWT, [-2000:100:-200], 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal'); 
hold on;  histogram(depthsHET, [-2000:100:-200], 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal'); 
box off
ylabel('Depth - um')
xlabel('# Cells')
ax = gca;
ax.FontSize = 20;

figure; histogram(depthsWT, [-2000:100:-200], 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal', 'Normalization', 'pdf'); 
hold on;  histogram(depthsHET, [-2000:100:-200], 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal', 'Normalization', 'pdf'); 
box off
ylabel('Depth - um')
xlabel('PDF')
ax = gca;
ax.FontSize = 20;

% Stats - depths of ALL cells
[p, h] = kstest2(depthsWT, depthsHET)
mean(depthsWT)
mean(depthsHET)




% RESPONSIVE CELLS

depth_resp_WT = all_av_hist(all_WT_RESP, W+3)*-1; 
depth_resp_HET = all_av_hist(all_HET_RESP, W+3)*-1;

% histogram PLOTS

figure; histogram(depth_resp_WT, [-2000:100:-200], 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal'); 
hold on;  histogram(depth_resp_HET, [-2000:100:-200], 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal'); 
box off
ylabel('Depth - um')
xlabel('# Cells')
ax = gca;
ax.FontSize = 20;

figure; histogram(depth_resp_WT, [-2000:100:-200], 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal', 'Normalization', 'pdf'); 
hold on;  histogram(depth_resp_HET, [-2000:100:-200], 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Orientation', 'Horizontal', 'Normalization', 'pdf'); 
box off
ylabel('Depth - um')
xlabel('PDF')
ax = gca;
ax.FontSize = 20;

% Stats - depths of RESPONSIVE cells
[p, h] = kstest2(depth_resp_WT, depth_resp_HET)
mean(depth_resp_WT)
mean(depth_resp_HET)











%% MAKE FIRING ARRAY FOR ALL CELLS. - CAN USE THE INDICES OF 'RESP CELLS' LATER. 

% Each row is a cluster. 
% Each column contains info about that cluster's response to the LOOM STIMULUS. 
% The Last columns == depth, animal and geno. Use these later for stats. 

% THIS ARRAY IS NOT NORMALISED BY THE FIRING RATE OF THE CLUSTER. 'RAW' SPIKING. 

% renaming things for ease. 
all_cells_hist = all_hist; 
all_hist = all_av_hist;

%
% lf = 43;
% figure;
% imagesc(all_av_hist(all_HET_RESP, 1:350))
% hold on 
% plot([60 60], [0 2500], 'w')
% plot([60+lf 60+lf], [0 2500], 'w')
% plot([60+lf*2 60+lf*2], [0 2500], 'w')
% plot([60+lf*3 60+lf*3], [0 2500], 'w')
% plot([60+lf*4 60+lf*4], [0 2500], 'w')
% plot([60+lf*5 60+lf*5], [0 2500], 'w')
% caxis([0 300]); title('HET-RESP')

% Make array - each row is a cluster. 
firing_array  = zeros(n_cells, 55);

st = 60; 
lf = 46+180; 
br = 180; 
    
for i =  1: n_cells
    
    % 1 - Average firing rate for the 1s before the loom and the 1s after loom
    firing_array(i,1) = mean([all_hist(i,1:st)]);
    
    
    
    % 2 - Number of spikes during loom stimulus
    all_spikes1 = mean(all_hist(i,st+1:st+lf-br)); % 45 for first loom, 46 for the rest. 
    firing_array(i,2) = all_spikes1; 
    
    all_spikes2 =  mean(all_hist(i,st+lf+1:st+lf*2-br));
    firing_array(i,3) = all_spikes2;
    
    all_spikes3 =  mean(all_hist(i,st+lf*2+1:st+lf*3-br));
    firing_array(i,4) =all_spikes3;
    
    all_spikes4 = mean(all_hist(i,st+lf*3+1:st+lf*4-br));
    firing_array(i,5) = all_spikes4;
    
    all_spikes5 = mean(all_hist(i,st+lf*4+1:st+lf*5-br));
    firing_array(i,6) =all_spikes5;
    
    all_spikes6 =  mean(all_hist(i,st+lf*5+1:st+lf*6-br));
    firing_array(i,7) = all_spikes6;
    
    all_spikes7 =  mean(all_hist(i,st+lf*6+1:st+lf*7-br));
    firing_array(i,8) =all_spikes7;
    
    all_spikes8 = mean(all_hist(i,st+lf*7+1:st+lf*8-br));
    firing_array(i,9) = all_spikes8;
    
    all_spikes9 = mean(all_hist(i,st+lf*8+1:st+lf*9-br));
    firing_array(i,10) =all_spikes9;
    
    all_spikes10 =  mean(all_hist(i,st+lf*9+1:st+lf*10-br));
    firing_array(i,11) = all_spikes10;
    
    all_spikes11 =  mean(all_hist(i,st+lf*10+1:st+lf*11-br));
    firing_array(i,12) =all_spikes11;
    
    all_spikes12 = mean(all_hist(i,st+lf*11+1:st+lf*12-br));
    firing_array(i,13) = all_spikes12;
    
    all_spikes13 = mean(all_hist(i,st+lf*12+1:st+lf*13-br));
    firing_array(i,14) =all_spikes13;
    
    all_spikes14 =  mean(all_hist(i,st+lf*13+1:st+lf*14-br));
    firing_array(i,15) = all_spikes14;
    
    all_spikes15 =  mean(all_hist(i,st+lf*14+1:st+lf*15-br));
    firing_array(i,16) =all_spikes15;

    
   
    %Average spikes during loom.
    firing_array(i,17) = mean([all_spikes1, all_spikes2, all_spikes3, all_spikes4, all_spikes5, all_spikes6, all_spikes7, all_spikes8, all_spikes9, all_spikes10, all_spikes11, all_spikes12, all_spikes13, all_spikes14, all_spikes15]);
    
    
    % 3 - Time to Peak per loom
    rows = find(all_hist(i,st+1:st+lf)==max(all_hist(i,st+1:st+lf-br)));
    rowpeak = rows(1);  % For loom one - add 9 frames
    t2peak1 = rowpeak/60; %Find time in seconds.
    firing_array(i,18) = t2peak1;
    
    rows = find(all_hist(i,st+lf+1:st+lf*2)==max(all_hist(i,st+lf+1:st+lf*2-br)));
    rowpeak = rows(1);
    t2peak2 = rowpeak/60;
    firing_array(i,19) = t2peak2 ;
    
    rows = find(all_hist(i,st+lf*2+1:st+lf*3)==max(all_hist(i,st+lf*2+1:st+lf*3-br)));
    rowpeak = rows(1);
    t2peak3 = rowpeak/60;
    firing_array(i,20) =t2peak3;
    
    rows = find(all_hist(i,st+lf*3+1:st+lf*4)==max(all_hist(i,st+lf*3+1:st+lf*4-br)));
    rowpeak = rows(1);
    t2peak4 = rowpeak/60;
    firing_array(i,21) = t2peak4;
    
    rows = find(all_hist(i,st+lf*4+1:st+lf*5)==max(all_hist(i,st+lf*4+1:st+lf*5-br)));
    rowpeak = rows(1);
    t2peak5 = rowpeak/60;
    firing_array(i,22) =t2peak5;
   
    
    rows = find(all_hist(i,st+lf*5+1:st+lf*6)==max(all_hist(i,st+lf*5+1:st+lf*6-br)));
    rowpeak = rows(1);  % For loom one - add 9 frames
    t2peak6 = rowpeak/60; %Find time in seconds.
    firing_array(i,23) = t2peak6;
    
    rows = find(all_hist(i,st+lf*6+1:st+lf*7)==max(all_hist(i,st+lf*6+1:st+lf*7-br)));
    rowpeak = rows(1);
    t2peak7 = rowpeak/60;
    firing_array(i,24) = t2peak7 ;
    
    rows = find(all_hist(i,st+lf*7+1:st+lf*8)==max(all_hist(i,st+lf*7+1:st+lf*8-br)));
    rowpeak = rows(1);
    t2peak8 = rowpeak/60;
    firing_array(i,25) =t2peak8;
    
    rows = find(all_hist(i,st+lf*8+1:st+lf*9)==max(all_hist(i,st+lf*8+1:st+lf*9-br)));
    rowpeak = rows(1);
    t2peak9 = rowpeak/60;
    firing_array(i,26) = t2peak9;
    
    rows = find(all_hist(i,st+lf*9+1:st+lf*10)==max(all_hist(i,st+lf*9+1:st+lf*10-br)));
    rowpeak = rows(1);
    t2peak10 = rowpeak/60;
    firing_array(i,27) =t2peak10;
    
    
    rows = find(all_hist(i,st+lf*10+1:st+lf*11)==max(all_hist(i,st+lf*10+1:st+lf*11-br)));
    rowpeak = rows(1);  % For loom one - add 9 frames
    t2peak11 = rowpeak/60; %Find time in seconds.
    firing_array(i,28) = t2peak11;
    
    rows = find(all_hist(i,st+lf*11+1:st+lf*12)==max(all_hist(i,st+lf*11+1:st+lf*12-br)));
    rowpeak = rows(1);
    t2peak12 = rowpeak/60;
    firing_array(i,29) = t2peak12 ;
    
    rows = find(all_hist(i,st+lf*12+1:st+lf*13)==max(all_hist(i,st+lf*12+1:st+lf*13-br)));
    rowpeak = rows(1);
    t2peak13 = rowpeak/60;
    firing_array(i,30) =t2peak13;
    
    rows = find(all_hist(i,st+lf*13+1:st+lf*14)==max(all_hist(i,st+lf*13+1:st+lf*14-br)));
    rowpeak = rows(1);
    t2peak14 = rowpeak/60;
    firing_array(i,31) = t2peak14;
    
    rows = find(all_hist(i,st+lf*14+1:st+lf*15)==max(all_hist(i,st+lf*14+1:st+lf*15-br)));
    rowpeak = rows(1);
    t2peak15 = rowpeak/60;
    firing_array(i,32) =t2peak15;
    
    % Average T2PEak
    firing_array(i,33) =mean([t2peak1, t2peak2, t2peak3, t2peak4, t2peak5, t2peak6, t2peak7, t2peak8, t2peak9, t2peak10, t2peak11, t2peak12, t2peak13, t2peak14, t2peak15]);
    
    
    
    
    
    % 5 - Max Number of spikes during loom stimulus (15-20) 
    
    % 2 - Number of spikes during loom stimulus
    m1 = max(all_hist(i,st+1:st+lf-br)); % 45 for first loom, 46 for the rest. 
    firing_array(i,34) = m1; 
    
    m2 =  max(all_hist(i,st+lf+1:st+lf*2-br));
    firing_array(i,35) = m2;
    
    m3 =  max(all_hist(i,st+lf*2+1:st+lf*3-br));
    firing_array(i,36) =m3;
    
    m4 = max(all_hist(i,st+lf*3+1:st+lf*4-br));
    firing_array(i,37) = m4;
    
    m5 = max(all_hist(i,st+lf*4+1:st+lf*5-br));
    firing_array(i,38) =m5;
    
    m6 =  max(all_hist(i,st+lf*5+1:st+lf*6-br));
    firing_array(i,39) = m6;
    
    m7 =  max(all_hist(i,st+lf*6+1:st+lf*7-br));
    firing_array(i,40) =m7;
    
    m8 = max(all_hist(i,st+lf*7+1:st+lf*8-br));
    firing_array(i,41) = m8;
    
    m9 = max(all_hist(i,st+lf*8+1:st+lf*9-br));
    firing_array(i,42) =m9;
    
    m10 =  max(all_hist(i,st+lf*9+1:st+lf*10-br));
    firing_array(i,43) = m10;
    
    m11 =  max(all_hist(i,st+lf*10+1:st+lf*11-br));
    firing_array(i,44) =m11;
    
    m12 = max(all_hist(i,st+lf*11+1:st+lf*12-br));
    firing_array(i,45) = m12;
    
    m13 = max(all_hist(i,st+lf*12+1:st+lf*13-br));
    firing_array(i,46) =m13;
    
    m14 =  max(all_hist(i,st+lf*13+1:st+lf*14-br));
    firing_array(i,47) = m14;
    
    m15 =  max(all_hist(i,st+lf*14+1:st+lf*15-br));
    firing_array(i,48) =m15;

   
    %Average max spikes during loom.
    firing_array(i,49) = mean([m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15]);
    
    % Add depth, animal and geno.
    firing_array(i,50:55) = all_hist(i, W+1:W+6);
    
end

save('210826_Firing_Array_Ptchd1_Multiloom_N8.mat', 'firing_array');



    %% Stats - - - - - one loom every 45 frames. 
    
    % 1 - Time to Peak for each loom. 
    % 2 - Time to peak average per cell.
    % 3 - Average firing rate per cell.  - look at distribution of cell firing rates. WT HET. 
    % 4 - Change in firing rate during loom. (Spikes 309-342)/ Spikes(258:300) 
  
  allWT_RESP_sSC = find(all_av_hist(all_WT_RESP,W+4)<1400);
  allWT_RESP_dSC = find(all_av_hist(all_WT_RESP,W+4)>=1400);
  
  allHET_RESP_sSC = find(all_av_hist(all_HET_RESP,W+4)<1400);
  allHET_RESP_dSC = find(all_av_hist(all_HET_RESP,W+4)>=1400);
    
  
  wt_resp = find(firing_array(:,53)==1 & firing_array(:,55)==1); 
  het_resp = find(firing_array(:,53)==0 & firing_array(:,55)==1); 
  
  % all cells 
%   wt_resp = find(firing_array(:,53)==1); 
%   het_resp = find(firing_array(:,53)==0); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  Make Stats array.
  
  stats_firing_array = zeros(8, 49);
  
  allnums = [1:1:49]; 
  
  % Can change 'allWT' to 'all_WT_RESP' or 'allWT_RESP_sSC' or 'allWT_RESP_dSC' - same for HET - to get just
  % responsive cells. 
  
  for j = allnums
  WTvals = firing_array(wt_resp, j); 
  HETvals = firing_array(het_resp, j); 
  
  nWT = numel(WTvals);
  nHET = numel(HETvals);
  
  stats_firing_array(1, j) = nanmean(WTvals);
  stats_firing_array(2, j) = nanmean(HETvals);
  
  [p1,h1] = ranksum(WTvals, HETvals);
  
  stats_firing_array(3, j) =p1 ;
  stats_firing_array(4, j) =h1; 
      
  [h2, p2] = kstest2(WTvals, HETvals);
    
  stats_firing_array(5, j) =p2; 
  stats_firing_array(6, j) =h2; 
     
   stats_firing_array(7, j) = nanstd(WTvals)/sqrt(nWT);
  stats_firing_array(8, j) =nanstd(HETvals)/sqrt(nHET);   
      
  end 
%       
  stats_firing_table = array2table(stats_firing_array, 'VariableNames', {'BaselineFiring', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12', 'N13', 'N14', 'N15', 'AvNumSpikes', 'T2P1', 'T2P2', 'T2P3', 'T2P4', 'T2P5', 'T2P6', 'T2P7', 'T2P8', 'T2P9', 'T2P10', 'T2P11', 'T2P12', 'T2P13', 'T2P14', 'T2P15', 'AvT2P', 'MaxL1', 'MaxL2', 'MaxL3', 'MaxL4', 'MaxL5', 'MaxL6', 'MaxL7', 'MaxL8', 'MaxL9', 'MaxL10', 'MaxL11', 'MaxL12', 'MaxL13', 'MaxL14', 'MaxL15', 'MaxAv'}, 'RowNames', {'MeanWT', 'MeanHET', 'RSp', 'RSh', 'KSp', 'KSh', 'SEMWT', 'SEMHET'});
 
  save('210722_StatsFiringTable_Ptchd1_Multiloom_ALLCELLS.mat', 'stats_firing_table');

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   %% 5 depths! 
    
   for jj = 1:4 % 7 different depths ranges.
       
       
       if jj ==1
           Dstart = 0;
           Dend = 750;
       elseif jj ==2
           Dstart = 750;
           Dend = 1250;
       elseif jj ==3
           Dstart = 1250;
           Dend = 1750;
       elseif jj ==4
           Dstart = 1750;
           Dend = 2500;
       end  
       
       %        WTrows = find(firing_array(:,26) == 1 & firing_array(:, 25) > Dstart & firing_array(:,25)<= Dend);
       %        HETrows = find(firing_array(:,26) == 0 & firing_array(:, 25) > Dstart & firing_array(:,25)<= Dend);
       
       % ALL CELLS
%        wt_rows = find(firing_array(:,52)> Dstart & firing_array(:,52) <= Dend  & firing_array(:,53)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells.
%        het_rows = find(firing_array(:,52)> Dstart & firing_array(:,52) <= Dend & firing_array(:,53)==0);
%        
       % RESP CELLS
       wt_rows = find(firing_array(:,52)> Dstart & firing_array(:,52) <= Dend & firing_array(:,55)==1 & firing_array(:,53)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 
       het_rows = find(firing_array(:,52)> Dstart & firing_array(:,52) <= Dend & firing_array(:,55)==1 & firing_array(:,53)==0);
%        
%        WTrows = find(firing_array(:,24) == 1);
%        HETrows = find(firing_array(:,24) == 0);
       
       stats_firing_array = zeros(8, 49);
       
       for j = 1:49
           WTvals = firing_array(wt_rows, j); 
           HETvals = firing_array(het_rows, j);
           
           nWT = numel(WTvals);
           nHET = numel(HETvals);
           
           stats_firing_array(1, j) = nanmean(WTvals);
           stats_firing_array(2, j) =nanmean(HETvals);
           
           [p1,h1] = ranksum(WTvals, HETvals);
           stats_firing_array(3, j) =p1 ;
           stats_firing_array(4, j) =h1;
           
           [h2, p2] = kstest2(WTvals, HETvals);
           stats_firing_array(5, j) =p2;
           stats_firing_array(6, j) =h2;
           
           stats_firing_array(7, j) = std(WTvals)/sqrt(nWT);
           stats_firing_array(8, j) =std(HETvals)/sqrt(nHET);
           
       end
       
  stats_firing_table = array2table(stats_firing_array, 'VariableNames', {'BaselineFiring', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12', 'N13', 'N14', 'N15', 'AvNumSpikes', 'T2P1', 'T2P2', 'T2P3', 'T2P4', 'T2P5', 'T2P6', 'T2P7', 'T2P8', 'T2P9', 'T2P10', 'T2P11', 'T2P12', 'T2P13', 'T2P14', 'T2P15', 'AvT2P', 'MaxL1', 'MaxL2', 'MaxL3', 'MaxL4', 'MaxL5', 'MaxL6', 'MaxL7', 'MaxL8', 'MaxL9', 'MaxL10', 'MaxL11', 'MaxL12', 'MaxL13', 'MaxL14', 'MaxL15', 'MaxAv'}, 'RowNames', {'MeanWT', 'MeanHET', 'RSp', 'RSh', 'KSp', 'KSh', 'SEMWT', 'SEMHET'});


       if jj ==1
           D1 = stats_firing_table;
         %  save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/DataArrays_Flipped/StatsByDepth/210120_Stats_D1_RespCELLS.mat'), 'D1');
       elseif jj ==2
           D2 = stats_firing_table;
          % save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/DataArrays_Flipped/StatsByDepth/210120_Stats_D2_RespCELLS.mat'), 'D2');
       elseif jj == 3
           D3 = stats_firing_table;
           %save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/DataArrays_Flipped/StatsByDepth/210120_Stats_D3_RespCELLS.mat'), 'D3');
       elseif jj ==4
           D4 = stats_firing_table;
           %save(strcat('//Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/DataArrays_Flipped/StatsByDepth/210120_Stats_D4_RespCELLS.mat'), 'D4');
       elseif jj ==5
           D5 = stats_firing_table;
           %save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/DataArrays_Flipped/StatsByDepth/210120_Stats_D5_RespCELLS.mat'), 'D5');
      end
       %  save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/210104_DepthStatsArrays/210104_Stats_D', string(jj), '.mat', strcat('D', string(jj))));
       
   end
   
   save('210826_Stats_RespCells_Ptchd1_Multiloom.mat', 'stats_firing_table', 'D1', 'D2', 'D3', 'D4');
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   %% 2 depths
   
    for jj = 1:2 % 7 different depths ranges.
       
       % Define depth range.
       if jj ==1
           Dstart = 400;
           Dend = 1300;
       elseif jj ==2
           Dstart = 1300;
           Dend = 3000;
       end 

       wt_rows = find(firing_array(:,25)> Dstart & firing_array(:,25) <= Dend & firing_array(:,28)==1 & firing_array(:,26)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 
       het_rows = find(firing_array(:,25)> Dstart & firing_array(:,25) <= Dend & firing_array(:,28)==1 & firing_array(:,26)==0);
       
       stats_firing_array = zeros(8, 22);
       
       for j = 1:21
           WTvals = firing_array(wt_rows, j); 
           HETvals = firing_array(het_rows, j);
           
           nWT = numel(WTvals);
           nHET = numel(HETvals);
           
           stats_firing_array(1, j) = nanmean(WTvals);
           stats_firing_array(2, j) =nanmean(HETvals);
           
           [p1,h1] = ranksum(WTvals, HETvals);
           stats_firing_array(3, j) =p1 ;
           stats_firing_array(4, j) =h1;
           
           [h2, p2] = kstest2(WTvals, HETvals);
           stats_firing_array(5, j) =p2;
           stats_firing_array(6, j) =h2;
           
           stats_firing_array(7, j) = std(WTvals)/sqrt(nWT);
           stats_firing_array(8, j) =std(HETvals)/sqrt(nHET);
           
       end
       
%        stats_firing_table = array2table(stats_firing_array, 'VariableNames', {'BaselineFiring', 'NumSpikes1', 'NumSpikes2', 'NumSpikes3', 'NumSpikes4', 'NumSpikes5', 'AvNumSpikes', 'T2P1', 'T2P2', 'T2P3', 'T2P4', 'T2P5', 'AvT2P', 'DeltaFiring'}, 'RowNames', {'MeanWT', 'MeanHET', 'RSp', 'RSh', 'KSp', 'KSh', 'SEMWT', 'SEMHET'});
         stats_firing_table = array2table(stats_firing_array, 'VariableNames', {'BaselineFiring', 'NumSpikes1', 'NumSpikes2', 'NumSpikes3', 'NumSpikes4', 'NumSpikes5', 'AvNumSpikes', 'T2P1', 'T2P2', 'T2P3', 'T2P4', 'T2P5', 'AvT2P', 'DeltaFiring', 'DeltaFiringL1', 'MaxL1', 'MaxL2', 'MaxL3', 'MaxL4', 'MaxL5', 'MaxAv', 'Depth'}, 'RowNames', {'MeanWT', 'MeanHET', 'RSp', 'RSh', 'KSp', 'KSh', 'SEMWT', 'SEMHET'});


       if jj ==1
           sSC_stats = stats_firing_table;
         %  save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/DataArrays_Flipped/StatsByDepth/210120_Stats_D1_RespCELLS.mat'), 'D1');
       elseif jj ==2
           dSC_stats = stats_firing_table;
          % save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Loom/DataArrays_Flipped/StatsByDepth/210120_Stats_D2_RespCELLS.mat'), 'D2');
       end
       
   end
   
      save('210323_Stats_Arrays_sSC-dSC_Resp_Cells.mat', 'sSC_stats', 'dSC_stats'); 

   
   
   
   
   
   
   
   
   
   
   
   
   
 %% MAKE PLOTS
 
 % For each DEPTH make a plot containing the following subplots:
 
 % 1 - Average Firing (Mean + SEM) - across all looms 
 % 2 - Baseline Firing
 % 3 - Mean all looms 
 % 4 - Boxplots of Max Spiking / Average spiking / decay value +
 % 4b - dot plots of the values above by depth. 
 
 % 5 - Errorbar plots comparing ACROSS depths. 
 
 
 %% Mean + SEM - Basline firing at each depth. 
 
%  col = 'r';
  col = [255/255 114/255 32/255]; % orange ptchd1
 
 all_HIST_LOOM = all_hist; 
 
 
st = 60; 
lf = 46+180; 
br = 180; 
 
 % IF ALL DEPTHS = 
%  WTrows = find(all_HIST_LOOM(:,603) == 1);
% HETrows = find(all_HIST_LOOM(:,603) == 0);
 
 for jj = 1:4
 
       
       if jj ==1
           Dstart = 0;
           Dend = 750;
       elseif jj ==2
           Dstart = 750;
           Dend = 1250;
       elseif jj ==3
           Dstart = 1250;
           Dend = 1700;
       elseif jj ==4
           Dstart = 1700;
           Dend = 2200;
       elseif jj == 5
           Dstart = 0;
           Dend = 2000;
       end
%        
%         if jj ==1
%            Dstart = 700;
%            Dend = 1300;
%        elseif jj ==2
%            Dstart = 1300;
%            Dend = 3000;
%         end 
%        
%        WTrows = find(all_HIST_LOOM(:,603) == 1 & all_HIST_LOOM(:, 601) > 0 & all_HIST_LOOM(:,601)<= 1000);
%        HETrows = find(all_HIST_LOOM(:,603) == 0 & all_HIST_LOOM(:, 601) > 0 & all_HIST_LOOM(:,601)<= 1000);
%
%        WTrows = find(all_HIST_LOOM(all_WT_RESP, 354) > Dstart & all_HIST_LOOM(all_WT_RESP,354)<= Dend);
%        HETrows = find(all_HIST_LOOM(all_HET_RESP, 354) > Dstart & all_HIST_LOOM(all_HET_RESP,354)<= Dend);

%        WTrows = find(all_HIST_LOOM(:,813)> Dstart & all_HIST_LOOM(:,813) <= Dend & all_HIST_LOOM(:, 814)==1 & all_HIST_LOOM(:, 816)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells.
%        HETrows = find(all_HIST_LOOM(:,813)> Dstart & all_HIST_LOOM(:,813) <= Dend & all_HIST_LOOM(:, 814)==0 & all_HIST_LOOM(:, 816)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells.

        WTrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==1 & all_HIST_LOOM(:, W+6)==1 ); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells.
        HETrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==0 & all_HIST_LOOM(:, W+6)==1 ); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells.

        nWT = numel(WTrows);
       nHET = numel(HETrows);
       
%        WTrows = all_WT_RESP;
%        HETrows = all_HET_RESP; 
       
%        WTrows = allWT;
%        HETrows = allHET; 
% 
%         WTrows = allWT_RESP_dSC;
%         HETrows = allHET_RESP_dSC; 
%         
%         
%        WTrows =  allWT_RESP_ZERO;
%        HETrows = allHET_RESP_ZERO; 
  
       %% BASELINE PLOT 
%        
%        WT_base = [all_HIST_LOOM(WTrows, 1:60), all_HIST_LOOM(WTrows, 291:350)]; 
%        HET_base = [all_HIST_LOOM(HETrows, 1:60), all_HIST_LOOM(HETrows, 291:350)]; 
%         
%        mWT = mean(WT_base);
%        mHET = mean(HET_base);
%        
%        x = (1:1:120);
%        
%        semWT = std(WT_base)/sqrt(nWT);
%        y1 = mWT+semWT;
%        y2 = mWT-semWT;
%        
%        semHET = std(WT_base)/sqrt(nHET);
%        y3 = mHET+semHET;
%        y4 = mHET-semHET;
%        
%        figure
%        plot(x, y1, 'w')
%        hold on
%        plot(x, y2, 'w')
%        patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
%        plot(mWT', 'k', 'LineWidth', 1)
%        
%        plot(x, y3, 'w')
%        hold on
%        plot(x, y4, 'w')
%        patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')
%        plot(mHET', col, 'LineWidth', 1)
%        xticks(1:60:300)
%        xticklabels({'0', '1', '2', '3', '4', '5'})
%        
% %        title(strcat('Depth - ', string(jj)))
%        box off
%        ax = gca;
%        ax.FontSize = 20; 
       
       
       %% AVERAGE OVER LOOMS 
       
       lWT = all_HIST_LOOM(WTrows, 1:W); 
       lHET = all_HIST_LOOM(HETrows, 1:W); 
       
       mlWT = (mean(lWT)); 
       mlHET = (mean(lHET)); 
       
       semWT = std(lWT)/sqrt(nWT);
       y1 = mlWT+semWT;
       y2 = mlWT-semWT;
       
       semHET = std(lHET)/sqrt(nHET);
       y3 = mlHET+semHET;
       y4 = mlHET-semHET;
       
       x = (1:1:W);
       
       lf = 46+180; 
        
       figure
       plot(x, y1, 'w')
       hold on
       plot(x, y2, 'w')
       patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mlWT', 'k', 'LineWidth', 1)
%      
       plot(x, y3, 'w')
       hold on
       plot(x, y4, 'w')
       patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mlHET', 'Color', col, 'LineWidth', 1)
       
       plot([60 60], [0 60], 'k:', 'LineWidth', 1.2)
       for i = 1:14
           plot([60+lf*i 60+lf*i], [0 60], 'k')
       end
       
%        vals = 60:60:300;
%        xticks(vals)
%        xticks([19, 49, 79, 109, 139, 169, 199, 129, 159, 189, 209, 239, 269, 299])
%        xticklabels({'0','1', '2',  '3',  '4'})
        xticks('')
%        xlabel('Time - s')
       ylabel('Activity - Hz')
       set(gca, 'FontSize', 20)
      axis([0 W 0 60])
       box off
       
%        
       
        figure
        subplot(2,1,1)
       plot(x, y1, 'w')
       hold on
       plot(x, y2, 'w')
       patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mlWT', 'k', 'LineWidth', 1)
       ylim([0 60])
       box off
       xlim([0 3500])
       subplot(2,1,2)
       plot(x, y3, 'w')
       hold on
       plot(x, y4, 'w')
       patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mlHET', 'Color', col, 'LineWidth', 1)
       ylim([0 60])
       xlim([0 3500])
       box off
       
       
       %% Average of loom response - trace for 'one' loom.  
% 
%          
%     % 2 - Number of spikes during loom stimulus
      WTvals = vertcat(all_hist(WTrows,st+1:st+lf-br), all_hist(WTrows,st+lf+1:st+lf*2-br), all_hist(WTrows,st+lf*2+1:st+lf*3-br), all_hist(WTrows,st+lf*3+1:st+lf*4-br),all_hist(WTrows,st+lf*4+1:st+lf*5-br),all_hist(WTrows,st+lf*5+1:st+lf*6-br), all_hist(WTrows,st+lf*6+1:st+lf*7-br), all_hist(WTrows,st+lf*7+1:st+lf*8-br), all_hist(WTrows,st+lf*8+1:st+lf*9-br), all_hist(WTrows,st+lf*9+1:st+lf*10-br), all_hist(WTrows,st+lf*10+1:st+lf*11-br), all_hist(WTrows,st+lf*11+1:st+lf*12-br), all_hist(WTrows,st+lf*12+1:st+lf*13-br), all_hist(WTrows,st+lf*13+1:st+lf*14-br), all_hist(WTrows,st+lf*14+1:st+lf*15-br));
      HETvals = vertcat(all_hist(HETrows,st+1:st+lf-br), all_hist(HETrows,st+lf+1:st+lf*2-br), all_hist(HETrows,st+lf*2+1:st+lf*3-br), all_hist(HETrows,st+lf*3+1:st+lf*4-br),all_hist(HETrows,st+lf*4+1:st+lf*5-br),all_hist(HETrows,st+lf*5+1:st+lf*6-br), all_hist(HETrows,st+lf*6+1:st+lf*7-br), all_hist(HETrows,st+lf*7+1:st+lf*8-br), all_hist(HETrows,st+lf*8+1:st+lf*9-br), all_hist(HETrows,st+lf*9+1:st+lf*10-br), all_hist(HETrows,st+lf*10+1:st+lf*11-br), all_hist(HETrows,st+lf*11+1:st+lf*12-br), all_hist(HETrows,st+lf*12+1:st+lf*13-br), all_hist(HETrows,st+lf*13+1:st+lf*14-br), all_hist(HETrows,st+lf*14+1:st+lf*15-br));

      % Combine means. 
      mw = mean(WTvals);
      mh = mean(HETvals);
      
       semWT = std(WTvals)/sqrt(nWT);
       y1 = mw+semWT;
       y2 = mw-semWT;
       
       semHET = std(HETvals)/sqrt(nHET);
       y3 = mh+semHET;
       y4 = mh-semHET;
       
       x3 = (1:1:46);
       
       subplot(4,1,jj)
       plot(x3, y1, 'w')
       hold on
       plot(x3, y2, 'w')
       patch([x3 fliplr(x3)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mw', 'k', 'LineWidth', 1)
     
       plot(x3, y3, 'w')
       hold on
       plot(x3, y4, 'w')
       patch([x3 fliplr(x3)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mh', 'Color', col, 'LineWidth', 1)
       
       vals = 1:6:30;
       xticks(vals)
%        xticks([19, 49, 79, 109, 139, 169, 199, 129, 159, 189, 209, 239, 269, 299])
       xticklabels({'0','0.1', '0.2',  '0.3',  '0.4', '0.5' , '0.6'})
       xlabel('Time - s')
       ylabel('Activity - Hz')
%        ylabel('zscore')
       set(gca, 'FontSize', 20)
       axis([0 30 0 80])
       box off
%        
   
 end 

 
 
 
 
 
 
%  figure
%   for jj = 1:4
%    
%        if jj ==1
%            Dstart = 0;
%            Dend = 700;
%        elseif jj ==2
%            Dstart = 700;
%            Dend = 1250;
%        elseif jj ==3
%            Dstart = 1250;
%            Dend = 1700;
%        elseif jj ==4
%            Dstart = 1700;
%            Dend = 2200;
%        elseif jj == 5
%            Dstart = 0;
%            Dend = 2500;
%        end
% 
%         WTrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. & all_HIST_LOOM(:, W+6)==1 
%         HETrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==0); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. & all_HIST_LOOM(:, W+6)==1 
% 
%         nWT = numel(WTrows);
%        nHET = numel(HETrows);
%        
%       WTvals = vertcat(all_hist(WTrows,st+1:st+lf-br), all_hist(WTrows,st+lf+1:st+lf*2-br), all_hist(WTrows,st+lf*2+1:st+lf*3-br), all_hist(WTrows,st+lf*3+1:st+lf*4-br),all_hist(WTrows,st+lf*4+1:st+lf*5-br),all_hist(WTrows,st+lf*5+1:st+lf*6-br), all_hist(WTrows,st+lf*6+1:st+lf*7-br), all_hist(WTrows,st+lf*7+1:st+lf*8-br), all_hist(WTrows,st+lf*8+1:st+lf*9-br), all_hist(WTrows,st+lf*9+1:st+lf*10-br), all_hist(WTrows,st+lf*10+1:st+lf*11-br), all_hist(WTrows,st+lf*11+1:st+lf*12-br), all_hist(WTrows,st+lf*12+1:st+lf*13-br), all_hist(WTrows,st+lf*13+1:st+lf*14-br), all_hist(WTrows,st+lf*14+1:st+lf*15-br));
%       HETvals = vertcat(all_hist(HETrows,st+1:st+lf-br), all_hist(HETrows,st+lf+1:st+lf*2-br), all_hist(HETrows,st+lf*2+1:st+lf*3-br), all_hist(HETrows,st+lf*3+1:st+lf*4-br),all_hist(HETrows,st+lf*4+1:st+lf*5-br),all_hist(HETrows,st+lf*5+1:st+lf*6-br), all_hist(HETrows,st+lf*6+1:st+lf*7-br), all_hist(HETrows,st+lf*7+1:st+lf*8-br), all_hist(HETrows,st+lf*8+1:st+lf*9-br), all_hist(HETrows,st+lf*9+1:st+lf*10-br), all_hist(HETrows,st+lf*10+1:st+lf*11-br), all_hist(HETrows,st+lf*11+1:st+lf*12-br), all_hist(HETrows,st+lf*12+1:st+lf*13-br), all_hist(HETrows,st+lf*13+1:st+lf*14-br), all_hist(HETrows,st+lf*14+1:st+lf*15-br));
% 
%       % Combine means. 
%       mw = mean(WTvals);
%       mh = mean(HETvals);
%       
%        semWT = std(WTvals)/sqrt(nWT);
%        y1 = mw+semWT;
%        y2 = mw-semWT;
%        
%        semHET = std(HETvals)/sqrt(nHET);
%        y3 = mh+semHET;
%        y4 = mh-semHET;
%        
%        x3 = (1:1:46);
%        
%        subplot(4,1,jj)
%        
%        plot(x3, y1, 'w')
%        hold on
%        plot(x3, y2, 'w')
%        patch([x3 fliplr(x3)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
%        plot(mw', 'k', 'LineWidth', 1)
%      
%        plot(x3, y3, 'w')
%        hold on
%        plot(x3, y4, 'w')
%        patch([x3 fliplr(x3)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')
%        plot(mh', 'Color', col, 'LineWidth', 1)
%        
%        xticks([])
% %        vals = 0:6:42;
% %        xticks(vals)
% %        xticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7'});
% % %        xticks([19, 49, 79, 109, 139, 169, 199, 129, 159, 189, 209, 239, 269, 299])
% %        xticklabels({'0','0.1', '0.2',  '0.3',  '0.4', '0.5' , '0.6'})
%        xlabel('Time - s')
%        ylabel('Activity - Hz')
% %        ylabel('zscore')
%        set(gca, 'FontSize', 20)
%        axis([0 30 0 30])
%        box off
%        
%   end 
%  
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 %% Subplot with Mean + SEM at different depth intervals. 
 
  figure
 for jj = 1:4
     
       if jj ==1
           Dstart = 0;
           Dend = 750;
       elseif jj ==2
           Dstart = 750;
           Dend = 1250;
       elseif jj ==3
           Dstart = 1250;
           Dend = 1750;
       elseif jj ==4
           Dstart = 1750;
           Dend = 2200;
       elseif jj == 5
           Dstart = 0; 
           Dend = 2500; 
       end      
       
          
       WTrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==1); %  & all_HIST_LOOM(:, W+6)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells.
       HETrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==0); % & all_HIST_LOOM(:, W+6)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 
       
        nWT = numel(WTrows)
       nHET = numel(HETrows)
       
       %% AVERAGE OVER LOOMS 
       
       lWT = all_HIST_LOOM(WTrows, 1:W); 
       lHET = all_HIST_LOOM(HETrows, 1:W); 
       
       mlWT = (mean(lWT)); 
       mlHET = (mean(lHET)); 
       
       semWT = std(lWT)/sqrt(nWT);
       y1 = mlWT+semWT;
       y2 = mlWT-semWT;
       
       semHET = std(lHET)/sqrt(nHET);
       y3 = mlHET+semHET;
       y4 = mlHET-semHET;
       
       x = (1:1:W);
       
       lf = 46+180; 
       
       subplot(4,1,jj)
       plot(x, y1, 'w')
       hold on
       plot(x, y2, 'w')
       patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mlWT', 'k', 'LineWidth', 1)
     
       plot(x, y3, 'w')
       hold on
       plot(x, y4, 'w')
       patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')
       plot(mlHET', 'Color', col, 'LineWidth', 1)
%        
       plot([60 60], [0 60], 'k:', 'LineWidth', 1.2)
       for i = 1:15
           plot([60+lf*i 60+lf*i], [0 60], 'k')
       end
       xticks([])
         xlabel('Time - s')
%        ylabel('Activity - Hz')
       set(gca, 'FontSize', 15)
       axis([0 W 0 50])
       box off
       
 end 
 
 
 %% Sig test: Depths of ranges 
  
 for jj = 1:4
     
       if jj ==1
           Dstart = 0;
           Dend = 700;
       elseif jj ==2
           Dstart = 700;
           Dend = 1250;
       elseif jj ==3
           Dstart = 1250;
           Dend = 1750;
       elseif jj ==4
           Dstart = 1750;
           Dend = 2200;
       end      
       
       WTrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==1 & all_HIST_LOOM(:, W+6)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 
       HETrows = find(all_HIST_LOOM(:,W+3)> Dstart & all_HIST_LOOM(:,W+3) <= Dend & all_HIST_LOOM(:, W+4)==0 & all_HIST_LOOM(:, W+6)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 

       
       %% AVERAGE OVER LOOMS 
       
       lWT = all_HIST_LOOM(WTrows, W+3); 
       lHET = all_HIST_LOOM(HETrows, W+3); 
       
        nanmean(lWT)
        nanmean(lHET)
        
        [p,h] = ranksum(lWT, lHET)
        
        [h2, p2] = kstest2(lWT, lHET)
       
       
 end 
 
 
 
 
 
 depth = W+3
 
       
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
       %%  Plot of firing to looms over depth
      
       val = 3521;
     offset = 0; 
     lf = 47;
        
 for jj = 2:4
 
      
     % Define depth range.
       if jj ==1
           Dstart = 0;
           Dend = 700;
       elseif jj ==2
           Dstart = 700;
           Dend = 1250;
       elseif jj ==3
           Dstart = 1250;
           Dend = 1700;
       elseif jj ==4
           Dstart = 1700;
           Dend = 2200;
%        elseif jj ==5
%            Dstart = 1750;
%            Dend = 2000;
       end
   
          
%        WTrows = find(all_HIST_LOOM(:,813)> Dstart & all_HIST_LOOM(:,813) <= Dend & all_HIST_LOOM(:, 814)==1 & all_HIST_LOOM(:, 816)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 
%        HETrows = find(all_HIST_LOOM(:,813)> Dstart & all_HIST_LOOM(:,813) <= Dend & all_HIST_LOOM(:, 814)==0 & all_HIST_LOOM(:, 816)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 
%        
%        % ALL CELLS 
        WTrows = find(all_HIST_LOOM(:,813)> Dstart & all_HIST_LOOM(:,813) <= Dend  & all_HIST_LOOM(:, 816)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 
        HETrows = find(all_HIST_LOOM(:,813)> Dstart & all_HIST_LOOM(:,813) <= Dend  & all_HIST_LOOM(:, 816)==1); %%%%%%%%%%%%%%%%%%%%% CHange here for all/ resp cells. 


    
      % AVERAGE OVER LOOMS - excluding change in lumin. 
       
       lWT = (all_HIST_LOOM(WTrows, 1:val)); 
       lHET = (all_HIST_LOOM(HETrows, 1:val)); 
        
       mlWT = (mean(lWT)); 
       mlHET = (mean(lHET)); 
       
       x2 = (1:1:val);
       
%        figure
       subplot(4,2,[1,3,5,7])
       plot(mlWT' - offset, 'k', 'LineWidth', 1)
       hold on
       plot([60 60], [100 -100], 'k:', 'LineWidth', 1.2)
       for i = 1:15
           plot([60+lf*i 60+lf*i], [100 -100], 'k')
       end
       axis([0 800 -100 30])
       box off

        ax = gca;
        ax.FontSize = 20; 
        yticklabels({''})
       xticklabels({''})
       
       subplot(4,2,[2,4,6,8])
       plot(mlHET'- offset, 'Color', col, 'LineWidth', 1)
       hold on
       plot([60 60], [100 -100], 'k:', 'LineWidth', 1.2)
       for i = 1:15
           plot([60+lf*i 60+lf*i], [100 -100], 'k')
       end
        axis([0 800 -100 30])
       box off
        ax = gca;
        ax.FontSize = 20; 
   
%        vals = -9:60:242;
%        xticks(vals)
%        xticklabels({'0','1', '2',  '3',  '4'})
%        xlabel('Time - s')
      yticklabels({''})
       xticklabels({''})
       
       offset = offset + 50; 
 
 end 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  
 %% Make PLOTS - DOTS VERTICALLY 
 
 figure
Y = [4,3,2,1]; 
X = shape_info_table.Decay_L1_W(1:4); 
% X = shape_info(1:7, 1); 
plot(X, Y, 'k.-', 'MarkerSize', 25)
yticklabels({''})
 hold on 
 X2 = shape_info_table.Decay_L1_H(1:4); 
% X2 = shape_info(1:7, 2); 
plot(X2, Y, 'ro-', 'MarkerSize', 15)
box off  
ax = gca;
ax.FontSize = 20; 

% xlabel('Time - s')
% xlabel('')
xlabel('Activity - Hz')
axis([0 12 0 6]) 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 3 - Errorbar - all depths! 

%   col = [255/255 114/255 32/255]; % orange ptchd1
    col = [0.95 (114/255)-0.05 (32/255)-0.05]; 
  
 % AVERAGE FIRNG 
X = [1:1:15]; 
Y = stats_firing_table{1,2:16}; 
Y2 = stats_firing_table{2,2:16}; 
errWT = stats_firing_table{7,2:16}; 
errHET = stats_firing_table{8,2:16}; 


Y = stats_firing_table{1,18:32}; 
Y2 = stats_firing_table{2,18:32}; 
errWT = stats_firing_table{7,18:32}; 
errHET = stats_firing_table{8,18:32}; 


Y = stats_firing_table{1,34:48}; 
Y2 = stats_firing_table{2,34:48}; 
errWT = stats_firing_table{7,34:48}; 
errHET = stats_firing_table{8,34:48}; 


% figure
errorbar(X,Y,errWT,'k.-', 'MarkerSize', 28, 'CapSize',3)
hold on 
 errorbar(X,Y2,errHET, 'o-', 'Color', col, 'MarkerSize', 8, 'CapSize',3, 'MarkerFaceColor', 'w')

axis([0 16 0 140])
xticks(X)
xticklabels(X)
ylabel('Max Activity - Hz')
% ylabel('Time - s')
xlabel('Loom')
set(gca, 'FontSize', 20)
box off
hold off

















% 4 - Errorbar - specific depths! 
 
X = [1:1:15]; 
Y = D2{1,2:16}; 
Y2 = D2{2,2:16}; 
errWT = D2{7,2:16}; 
errHET = D2{8,2:16}; 


Y = D2{1,18:32}; 
Y2 = D2{2,18:32}; 
errWT = D2{7,18:32}; 
errHET = D2{8,18:32}; 


Y = D2{1,34:48}; 
Y2 = D2{2,34:48}; 
errWT = D2{7,34:48}; 
errHET = D2{8,34:48}; 


figure
errorbar(X,Y,errWT,'k.-', 'MarkerSize', 28, 'CapSize',3)
hold on 
errorbar(X,Y2,errHET, 'ro-', 'MarkerSize', 8, 'CapSize',3)


 axis([0.5 5.5 0 150])
xticks(X)
xticklabels(X)
ylabel('Activity - Hz')
xlabel('Loom')
set(gca, 'FontSize', 20)
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% T2P
 
X = [1,2,3,4,5]; 
Y = stats_firing_table{1,2:6}; 
Y2 = stats_firing_table{2,2:6}; 
errWT = stats_firing_table{7,2:6}; 
errHET = stats_firing_table{8,2:6}; 

figure
errorbar(X,Y,errWT,'k.-', 'MarkerSize', 28, 'CapSize',3)
hold on 
errorbar(X,Y2,errHET, 'ro-', 'MarkerSize', 8, 'CapSize',3)
axis([0.5 5.5 0 12])
xticks(X)
xticklabels(X)
ylabel('Activity - Hz')
xlabel('Loom')
set(gca, 'FontSize', 20)
box off


 % MAX 
 
X = [1,2,3,4,5]; 
Y = stats_firing_table{1,8:12}; 
Y2 = stats_firing_table{2,8:12}; 
errWT = stats_firing_table{7,8:12}; 
errHET = stats_firing_table{8,8:12}; 

figure
errorbar(X,Y,errWT,'k.-', 'MarkerSize', 28, 'CapSize',3)
hold on 
errorbar(X,Y2,errHET, 'ro-', 'MarkerSize', 8, 'CapSize',3)
axis([0.5 5.5 0 0.3])
xticks(X)
xticklabels(X)
ylabel('Time - s')
xlabel('Loom')
set(gca, 'FontSize', 20)
box off


