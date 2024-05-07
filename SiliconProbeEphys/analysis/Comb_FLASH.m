% Combine FLASH analysis arrays:
% Cretaed by Burnett - 24/02/22

% LOAD table with information about the depth of the sSC in different animals. 
% load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Setd5_Rec_Depth_Info.mat', 'animal_depth_info_table')
% load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Setd5_Rec_Depth2.mat', 'animal_depth_info_table')
load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/220912_AniDepth_Ptchd1.mat', 'animal_depth_info_table');

% ANALYSIS RESULTS FILES STORED HERE:
% "/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/SPIKESS/Flashes/FLASH_OFO"
% /Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/FLASH/SPIKES/RESULTS_FLASH

%% For all the analysis files want to:

% het_animals = [7269, 7476, 7614];
het_animals = [1385, 1386, 1394, 2709, 4369];
% het_animals = [2833, 3557, 4124]; 

%spiking data 
% dir1 = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/SPIKESS/Flashes/FLASH_OFO';
% dir1 = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/FlashOFO/NEW_withzetaP'; 
dir1 = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/FLASH/SPIKES/RESULTS2';
% dir1 = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Cul3/FLASH/SPIKES/RESULTS';
files1 = dir(fullfile(dir1, '*FLASH_1-0*'));
nfiles =numel(files1);

%% 
Flash_Table = table(); 

for i = 1:nfiles
    
    %FIRST - load the new OLGA file with the cell ID and p value 
    fname = fullfile(files1(i).folder, files1(i).name);
    load(fname)
    expname = files1(i).name; 
    
    animal_num = expname(8:11);
    date = str2num(expname(1:6));
    test = str2num(expname(13));
    probe = str2num(expname(15:18));
    GENOO = [];
    DEPTHH = []; 
    TOTAL_SP = [];
    
    n_gcl = numel(OFO(:,1)); %number of units in this recording
   
    row = find(animal_depth_info_table.Ani == str2num(animal_num)); %  & animal_depth_info_table.Date == date
    if ~isempty(row)
    
        sSC_depth = animal_depth_info_table.Depth(row);

    for jj = 1:n_gcl
        d = (probe - (800 - OFO(jj, 125)));
        d2 = sSC_depth - d;
        DEPTHH(jj,1) = d2; % was d2 % % % % % % % % % % % CHANGED from 127
        
        if ismember(str2num(animal_num), het_animals)
            GENOO(jj, 1) = 2;
        else
            GENOO(jj, 1) = 1; 
        end 
        
        TOTAL_SP(jj,1) = sum(OFO(jj, 1:120));
    end 
   
    % Make TABLE
    Date = OFO(:,121);
    Ani = OFO(:,122);
    Exp = OFO(:,123);
    Geno = GENOO;
    ProbeDepth = OFO(:,124);
    Depth = DEPTHH;
    TotalSpikes = TOTAL_SP;
    OFO = OFO(:, 1:120); % Average spiking over 1 frame (0.0167s bins) - average over 20 REPS ( 2 x 10) 
    OFO_NORM = OFO_NORM(:, 1:120); % Normalised (zscored) version of OFO above. 
    SPT = spikes_per_time;
    SPAV = spikes_rep_av;
    SPSTIM= spikes_per_stim;
    
    tbl = table(Date, Ani, Exp, Geno, ProbeDepth, TotalSpikes, Depth, OFO, OFO_NORM, SPT, SPAV, SPSTIM); 
    
    Flash_Table = vertcat(Flash_Table, tbl);

    end 
    
end 

save('220913_FlashTable_Cul3_120-OFO.mat', 'Flash_Table'); %, 'animal_depth_info_table')

% In reality - all we want if OFO and OFO_NORM. 
% save('220620_Setd5_FLASH_table_pval_CELLID_RAWnotNORM.mat', 'Flash_Table', 'animal_depth_info_table');


%% Add variables.  
% OFO = 120 frames
% 30 OFF - 60 ON - 30 OFF 

for kk = 1:height(Flash_Table)
    
    data = Flash_Table.OFO(kk, :);
    ONvals = data(31:90);
    OFFvals = [data(1:30), data(91:120)]; 
    
    ONimmed = data(31:45);
    OFFimmed = data(91:105);
    
    maxON = max(ONvals);
    maxOFF = max(OFFvals);
    
    avON = nanmean(ONvals);
    avOFF = nanmean(OFFvals);
    
    avONimmed = nanmean(ONimmed);
    avOFFimmed = nanmean(OFFimmed);
    
    diffsum_all = (avOFF-avON)/(avON+avOFF);
    diffsum_immed = (avOFFimmed-avONimmed)/(avONimmed+avOFFimmed);
    
    Flash_Table.MaxON(kk) = maxON;
    Flash_Table.MaxOFF(kk) = maxOFF;
    
    Flash_Table.AvON(kk) = avON;
    Flash_Table.AvOFF(kk) = avOFF;
    
    Flash_Table.AvONImmed(kk) = avONimmed;
    Flash_Table.AvOFFImmed(kk) = avOFFimmed;
    
    Flash_Table.DiffSum(kk) = diffsum_all;
    Flash_Table.DifSumImmed(kk) = diffsum_immed;
    
end 


% for kk = 1:height(Flash_Table)
%     if Flash_Table.Ani(kk) == 2027
%         Flash_Table.Geno(kk) = 2;
%     end 
% end 

%%

dataWT = Flash_Table.AvOFFImmed(allWT)*60;
dataHET = Flash_Table.AvOFFImmed(allHET)*60;

nWT = numel(dataWT);
nHET = numel(dataHET);

nanmean(dataWT)
nanmean(dataHET)
[p,h] = ranksum(dataWT, dataHET)
[h, p] = kstest2(dataWT, dataHET)


% 

xWT = ones(nWT, 1);
xHET = ones(nHET, 1)*2;

d = vertcat(dataWT, dataHET);
gp = vertcat(xWT, xHET);

b = boxplot(d, gp, 'Colors', [0 0 0; 0 0 0], 'Symbol','w.');
set(b, 'linew', 1.25);
ylim([-2 50])
box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
f = gcf;
f.Position = [1053  427  172  257]; 



figure; histogram(dWT,  'Normalization', 'pdf'); hold on 
histogram(dHET,  'Normalization', 'pdf');


%%
col = [255/255 114/255 32/255]; 

%% GENERAL OVERVIEW - WT/HET

allWT = find(Flash_Table.Geno == 1 & Flash_Table.Depth<0 & Flash_Table.Depth>-350  & Flash_Table.TotalSpikes >1 & abs(Flash_Table.DifSumImmed)>0.2); %& Flash_Table.Ani~=1387 
allHET = find(Flash_Table.Geno == 2 & Flash_Table.Depth<0 & Flash_Table.Depth>-350 & Flash_Table.Ani~=1387 & Flash_Table.TotalSpikes >1 & abs(Flash_Table.DifSumImmed)>0.2);

%% QUICK GLANCE: 

data = Flash_Table.OFO;
data = Flash_Table.OFO_NORM;
data = Flash_Table.avSpPT*60;
figure; imagesc(data(allWT, :)); caxis([0 350])%caxis([-0.1 2])
figure; imagesc(data(allHET, :)); caxis([0 350]) %caxis([-0.1 2])









%% KMEANS - SORT WT and HET together - KEY PLOTS - HEATMAP AND MEAN+SEM

close all

data = [(Flash_Table.OFO)*60, Flash_Table.Geno];

% Only for cells between 0 and -400. 
allWT = find(Flash_Table.Geno == 1 & Flash_Table.Depth<0 & Flash_Table.Depth>-300  & Flash_Table.TotalSpikes >1 & Flash_Table.Ani~=1387 & abs(Flash_Table.DifSumImmed)>0.2);
allHET = find(Flash_Table.Geno == 2 & Flash_Table.Depth<0 & Flash_Table.Depth>-300 & Flash_Table.Ani~=1387 & Flash_Table.TotalSpikes >1 & abs(Flash_Table.DifSumImmed)>0.2);

all_depths = find(Flash_Table.Depth<0 & Flash_Table.Depth>-350 & Flash_Table.TotalSpikes>1 & abs(Flash_Table.DifSumImmed)>0.1); % &  & Flash_Table.Depth>-600 &  Flash_Table.Depth>-1000 % 
% all_depths = find(Flash_Table.Depth>0 | Flash_Table.P_VALUE>=0.05);

data = data(all_depths, :);
% data(all_depths, :) = []; 

%% % % % % Evaluate the optimal number of clusters:
eva = evalclusters(data, 'kmeans','CalinskiHarabasz', 'KList', [1:30]); % 'silhouette',
figure; plot(eva)

%% RUN KMEANS

% Set the elbow value as the number of clusters. 
n_k = 8;

data(:, 122) = kmeans(data(:, 1:120), n_k);
data = sortrows(data, 122);

allWT = find(data(:, 121)==1);
allHET = find(data(:, 121)==2);

nWT = numel(allWT);
nHET = numel(allHET);

% save('/Users/lauraburnett/Documents/Burnett_etal/DATA/SETD5/EPHYS/SiliconProbe/Flashes/220620_FLASH_OFO_Hz_Kmeans10.mat', 'data');

%%

v = linspace(1,0, 1000);
cmap = [v', v', v'];

%% % % % % % % % PLOT HEATMAPS
close all

spn = 15; 

figure
ax = subplot(1,spn,1:spn-1);
imagesc(data(allWT,1:120));
hold on
plot([30 30], [0 nWT], 'r:', 'LineWidth', 2)
plot([90 90], [0 nWT], 'r:', 'LineWidth', 2)
box off
ax.XTick = []; 
colormap(ax(1), cmap)
% colorbar('SouthOutside')
% caxis(ax, [-0.05 0.35])
% caxis(ax, [-0.2 0.8])
caxis(ax, [0 80])
ax.YTick = [];

ax2 = subplot(1,spn,spn);
imagesc(data(allWT, 122))
colormap(ax2, gray)
% axis off
% sgtitle('WT')
f = gcf;
f.Position = [4     1   356   804];
box on
ax2.LineWidth = 0.5;
ax2.Color = 'k';
ax2.TickLength = [0 0];
ax2.XTick = [];
ax2.YTick = []; 
f = gcf;
f.Position = [100 600 300 nWT*4];

%
figure
ax = subplot(1,spn,1:spn-1);
imagesc(data(allHET,1:120)); 
hold on
plot([30 30], [0 nHET], 'r:', 'LineWidth', 2)
plot([90 90], [0 nHET], 'r:', 'LineWidth', 2)
box off
ax.XTick = []; 
colormap(ax(1), cmap)
% colorbar('SouthOutside')
% caxis(ax, [-0.05 0.35])
% caxis(ax, [-0.2 0.8])
caxis(ax, [0 80])
ax.YTick = [];

ax2 = subplot(1,spn,spn);
imagesc(data(allHET, 122))
colormap(ax2, gray)
f = gcf;
f.Position = [699 1 356   804];
box on
ax2.LineWidth = 0.5;
ax2.Color = 'k';
ax2.TickLength = [0 0];
ax2.XTick = [];
ax2.YTick = []; 
f2 = gcf;
f2.Position = [402   171   300   nWT*4]; %[100 600 300 nHET*2];    


%% Separate mean/sem plot with WT/ HET COMBINED
% close all

col = [1 0 0.8];

x = 1:1:120;

figure
for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 122)== j & data(:,121) == 1);
    % Separate these rows into a different array
    data2 = data(all_type, 1:120);
    
    % Variables for Mean/SEM
    n_trials = numel(all_type);
    if n_trials==1
        mWT = data2;
    else
        mWT = nanmean(data2);
    end
    semWT = nanstd(data2)/sqrt(n_trials);
    y1 = mWT+semWT;
    y2 = mWT-semWT;
    
    % Make the plot
    subplot(n_k,1,j); hold on
    
    % Find Min/Max
%     maxx = max(nanmean(data2))+0.3;
%     minn = min(nanmean(data2))-0.15;
%     
%     maxx = 42 ;%1.1; % 0.55
%     minn = 0 ; % -0.25; %-0.08

     all_cl= find(data(:, 122)==j); 
     max1 =max(nanmean(data(all_cl, 1:120)));
     maxx = max1 + max1*0.75;
     
     if maxx > 15
         minn = -3;
     else 
         minn = 0;
     end 

    % Rectangle of light stim
    rectangle('Position', [0 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
%     rectangle('Position', [30 minn 60 diff([minn maxx])], 'FaceColor', [1 1 0 0.1], 'EdgeColor', 'none')
    rectangle('Position', [90 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
    
    % Plot SEM
    plot(x, y1, 'w')
    plot(x, y2, 'w')
    patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    % Plot mean
    plot(mWT, 'k', 'LineWidth', 1.2);
    box off
    ylim([minn maxx])
    ax = gca;
    ax.TickDir = 'out';
    ax.XAxis.Visible = 'off';
%     if j~= n_k
%     ax.YAxis.Visible = 'off';
%     end 
    ax.LineWidth = 1;
    
end

for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 122)== j & data(:,121) == 2);
    % Separate these rows into a different array
    data2 = data(all_type, 1:120);
    
    % Variables for Mean/SEM
    n_trials = numel(all_type);
    if n_trials == 1
        mHET = data2;
    else
        mHET = nanmean(data2);
    end
    semHET = nanstd(data2)/sqrt(n_trials);
    y1 = mHET+semHET;
    y2 = mHET-semHET;
    
    % Make the plot
    subplot(n_k,1,j); hold on
    
    % Find Min/Max
%     maxx = max(nanmean(data2))+0.3;
%     minn = min(nanmean(data2))-0.15;

%      maxx = 40;% 1.1; % 0.55
%     minn = 0; %-0.25; %-0.08

    % Plot SEM
    plot(x, y1, 'w')
    plot(x, y2, 'w')
    patch([x fliplr(x)], [y1 fliplr(y2)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    % Plot mean
    plot(mHET, 'Color', col, 'LineWidth', 1.2);
    box off
%     ylim([minn maxx])
    ax = gca;
    ax.TickDir = 'out';
    ax.XAxis.Visible = 'off';
%     if j~= n_k
%     ax.YAxis.Visible = 'off';
%     end 
    ax.LineWidth = 1;
    
end
f4 = gcf;
% f4.Position = [1056  1 256 804];
f4.Position = [ 633    23   263   782]; % [237    41   272   764]; 

% save('/Users/lauraburnett/Documents/Burnett_etal/SVG_Figs/Setd5/EPHYS/FLASHES/K7/220225_DATA_Array_Flashes_Setd5_K7.mat', 'data', 'n_k')



%% Variance /Spread in WT/HET between clusters

allWT = find(data(:, 121)==1);
allHET = find(data(:, 121)==2);

nWT = numel(allWT);
nHET = numel(allHET);

clust_WT = data(allWT, 122);
clust_HET = data(allHET, 122);

[p,h] = kstest2(clust_WT,clust_HET)
figure; histogram(clust_WT, 'Normalization', 'pdf');
hold on;
histogram(clust_HET, 'Normalization', 'pdf')

%% Percentage of Flash Responsive Cells


all_ani = unique(Flash_Table.Ani);

WTVALS = [];
HETVALS = []; 

% het_animals = [7269, 7476, 7614];
het_animals = [1385, 1386, 1394, 2709, 4369];
% het_animals = [2833, 2027, 3557, 4124]; 

       
for j = 1:numel(all_ani)
    ani= all_ani(j);
 
    ani_all =find(Flash_Table.Ani == ani & Flash_Table.Depth<=0 & Flash_Table.Depth>-350 & Flash_Table.TotalSpikes>1); % Flash_Table.Depth<-500 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500
    n_all = numel(ani_all);
    ani_rows = find(Flash_Table.Ani == ani & Flash_Table.Depth<=0 & Flash_Table.Depth>-350 & Flash_Table.TotalSpikes>1 & Flash_Table.P_VALUE>0.05); % Flash_Table.Depth<-500 & Flash_Table.Depth>-1000
    n_resp = numel(ani_rows);
    
    resp_ratio = n_resp/n_all;
    
    if ismember(ani, het_animals)
        genoo = 2;
    else
        genoo = 1;
    end
    
    if genoo ==1
        WTVALS = [WTVALS, resp_ratio];
    elseif genoo ==2
        HETVALS = [HETVALS, resp_ratio];
    end
end

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)
[p,h] = ttest2(WTVALS, HETVALS)

%% BAR CHART + ANIMAL POINTS - % OF FLASH RESPONSIVE CELLS
% Plus errorbar

close
semwt = nanstd(WTVALS)/sqrt(numel(WTVALS)); 
semhet = nanstd(HETVALS)/sqrt(numel(HETVALS)); 

figure
bar(1, nanmean(WTVALS), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
bar(2, nanmean(HETVALS),'FaceColor', col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% WT
for i = 1:numel(WTVALS)
    xval =  1;
    yval = WTVALS(i);
   
    marker = 'k.';
    scatter(xval, yval, 1000,  marker, 'jitter', 'on');
    hold on 
end 

% HET
for i = 1:numel(HETVALS)
    xval =  2;
    yval = HETVALS(i);
   
%     marker = 'r.';
    scatter(xval, yval, 1000, col, '.', 'jitter', 'on');
    hold on 
end 

% axis([0 3 -1000 0])
axis([0 3 0 1])
box off
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.02, 0.02];
ax.LineWidth= 1.5;

errorbar(1, nanmean(WTVALS), semwt, 'k', 'LineWidth', 1.5);
errorbar(2, nanmean(HETVALS), semhet, 'Color', col, 'LineWidth', 1.5);

f = gcf;
f.Position = [680   750   230   348]; 









%% % % % % % % Plots of the Mean/SEM of each cluster:

x = 1:1:90;

% WT % % % % %
figure
for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 92)== j & data(:,91) == 1);
    % Separate these rows into a different array
    data2 = data(all_type, 1:90);
    
    % Variables for Mean/SEM
    n_trials = numel(all_type);
    if n_trials == 1
        mWT = data2;
    else
        mWT = nanmean(data2);
    end
    semWT = nanstd(data2)/sqrt(n_trials);
    y1 = mWT+semWT;
    y2 = mWT-semWT;
    
    % Make the plot
    subplot(n_k,1,j); hold on
    
    % Find Min/Max
%     maxx = max(nanmean(data2))+0.075;
%     minn = min(nanmean(data2))-0.075;
% 
%     maxx = 0.55; %1.8;
%     minn = -0.075; % -0.8;
%     
    % Rectangle of light stim
%     rectangle('Position', [0 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
% %     rectangle('Position', [30 minn 60 diff([minn maxx])], 'FaceColor', [1 1 0 0.1], 'EdgeColor', 'none')
%     rectangle('Position', [90 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
%     
    % Plot SEM
    plot(x, y1, 'w')
    plot(x, y2, 'w')
    patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    % Plot mean
    plot(mWT, 'k', 'LineWidth', 1.2);
    box off
%     ylim([minn maxx])
    ax = gca;
    ax.TickDir = 'out';
    ax.XAxis.Visible = 'off';
    ax.LineWidth = 1;
    
end
f2 = gcf;
f2.Position = [365  1   256   804];


% HET % % % % %
figure
for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 92)== j & data(:,91) == 2);
    % Separate these rows into a different array
    data2 = data(all_type, 1:90);
    
    % Variables for Mean/SEM
    n_trials = numel(all_type);
    if n_trials == 1
        mHET = data2;
    else 
        mHET = nanmean(data2);
    end 
    semHET = nanstd(data2)/sqrt(n_trials);
    y1 = mHET+semHET;
    y2 = mHET-semHET;
    
    % Make the plot
    subplot(n_k,1,j); hold on
    
    % Find Min/Max
    %     maxx = max(nanmean(data2))+0.075;
    %     minn = min(nanmean(data2))-0.075;
    
%     maxx = 0.55; %1.8;
%     minn = -0.075; % -0.8;
%     
%     % Rectangle of light stim
%     rectangle('Position', [0 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
%     %     rectangle('Position', [30 minn 60 diff([minn maxx])], 'FaceColor', [1 1 0 0.1], 'EdgeColor', 'none')
%     rectangle('Position', [90 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
%     
    % Plot SEM
    plot(x, y1, 'w')
    plot(x, y2, 'w')
    patch([x fliplr(x)], [y1 fliplr(y2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    % Plot mean
    plot(mHET, 'r', 'LineWidth', 1.2);
    box off
%     ylim([minn maxx])
    ax = gca;
    ax.TickDir = 'out';
    ax.XAxis.Visible = 'off';
    ax.LineWidth = 1;
    
end
f3 = gcf;
f3.Position = [1056  1 256 804];

%% PIE CHART of 'types' of responses for WT/HET

vals = NaN(2, n_k);

for jj = 1:n_k
    
    allWT = find(data(:, 122)== jj & data(:,121) == 1);
    allHET = find(data(:, 122)== jj & data(:,121) == 2);
    
    vals(1,jj) = numel(allWT);
    vals(2,jj) = numel(allHET);

end 

figure
% lb = {'ON-OFF', 'ONsus', 'NS', 'OFFsus', 'ON', 'OFFsupp', 'OFF'};
% lb = {'None', 'On', 'Off', 'On-Off'};
% expl = [1 1 1 1 1 1 1 1 1 1];
expl = ones(n_k,1);

X = vals(1,:);
pie(X, expl)
% title(ax1,'Setd5^+^/^+');
colormap(cmap)

figure
Y = vals(2,:);
pie(Y, expl)
% title(ax2,'Setd5^-^/^-');
colormap(cmap)

cmap = [0 0 0; 0.2  0.2 0.2; 0.3  0.3 0.3; 0.5 0.5 0.5; 0.6 0.6 0.6; 0.7 0.7 0.7; 0.8 0.8 0.8; 1,1, 1];
cmapH = [0 0 0;0.2  0 0; 0.3  0 0; 0.5  0 0; 0.7 0 0; 0.8 0 0; 1,0, 0];

cmap = [0 0 0; 0.25  0 0.1; 0.35 0 0.15; 0.45 0 0.25; 0.7 0 0.45; 0.8 0 0.55; 0.9 0 0.7; 1,0, 0.8];

%% Taking into consideration OLGA's new p-value from zeta-score! 

allWT = find(Flash_Table.Geno == 1);
allHET = find(Flash_Table.Geno == 2);

pval = 0.05;

allWT_p = find(Flash_Table.Geno == 1 & Flash_Table.P_VALUE <pval & Flash_Table.Depth <0);
allHET_p = find(Flash_Table.Geno == 2 & Flash_Table.P_VALUE <pval & Flash_Table.Depth <0);

allWT_np = find(Flash_Table.Geno == 1 & Flash_Table.P_VALUE >=pval & Flash_Table.Depth <0);
allHET_np = find(Flash_Table.Geno == 2 & Flash_Table.P_VALUE >=pval & Flash_Table.Depth <0);

data = [Flash_Table.avSpPT*60];

figure; imagesc(data(allWT, :)); caxis([0 100])
figure; imagesc(data(allHET, :)); caxis([0 100])

figure; imagesc(data(allWT_p, :)); caxis([0 25])
figure; imagesc(data(allHET_p, :)); caxis([0 25])

figure; imagesc(data(allWT_np, :)); caxis([0 25])
figure; imagesc(data(allHET_np, :)); caxis([0 25])

% %
%% Plot spiking with cells sorted by depth

data = [(Flash_Table.avSpPT)*60, Flash_Table.Geno, Flash_Table.Depth];

all_depths = find(Flash_Table.Depth<0 & Flash_Table.Depth>-1000 & Flash_Table.P_VALUE<0.05); % & & Flash_Table.TotalSpikes<3500  & Flash_Table.Depth>-600 &  Flash_Table.Depth>-1000 

data = data(all_depths, :);

data = sortrows(data, 1202);

allWT = find(data(:, 1201)==1);
allHET = find(data(:, 1201)==2);

nWT = numel(allWT);
nHET = numel(allHET);

figure; imagesc(data(allWT, 1:1200)); caxis([0 100])
figure; imagesc(data(allHET, 1:1200));caxis([0 100])


%% Plot - log of pvalue versus depth. 

figure
for j = 1:height(Flash_Table)
    
    xval = log(Flash_Table.P_VALUE(j));
    yval = Flash_Table.Depth(j);
    
    genooo = Flash_Table.Geno(j);
    
    if genooo == 1
        col = 'k';
    else
        col = 'r';
    end 
    
    plot(xval, yval, 'o', 'Color', col)
    hold on 
end 


%% Percentage of all WT/HET cells that are 'responsive'


allWTp = find(data(:, 121)==1);
allHETp = find(data(:, 121)==2);

nWTp = numel(allWTp);
nHETp = numel(allHETp);

allWT = find(Flash_Table.Geno == 1 & Flash_Table.Depth < 0);
allHET = find(Flash_Table.Geno == 2 & Flash_Table.Depth < 0);

nWT = numel(allWT);
nHET = numel(allHET);

percentresp_WT = (nWTp/nWT)*100; % 35.95%
percentresp_HET = (nHETp/nHET)*100; %33.85%

% % Plot figure - PIE of % responsive cells
% figure; pie([percentresp_WT, 100-percentresp_WT], [1,1])
% colormap(gray)
% 
% % Plot figure - PIE of % responsive cells
% figure; pie([percentresp_HET, 100-percentresp_HET], [1,1])
% colormap(gray)

% STACKED BAR CHART - BETTER. 
b = bar([percentresp_WT, 100-percentresp_WT; percentresp_HET, 100-percentresp_HET], 'stacked'); 
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [1 1 1];
b(1).LineWidth = 1.5;
b(2).LineWidth = 1.5;
box off
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
yticks([0, 25,50,75,100])
f= gcf;
f.Position = [440   457   192   341];





































































% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%% without NS

% lb = {'ON-OFF', 'ONsus', 'OFFsus', 'ON', 'OFFsupp', 'OFF'};
expl = ones(n_k,1);

X = vals(1,:);
ax1 = subplot(1,2,1);
pie(ax1,X, expl)
% title(ax1,'Setd5^+^/^+');
colormap(gray)

Y = vals(2,:);
ax2 = subplot(1,2,2);
pie(ax2,Y, expl)
% title(ax2,'Setd5^-^/^-');
colormap(gray)

%% 

allWT = find(data(:, 121)==1);
allHET = find(data(:, 121)==2);

distWT = data(allWT, 122);
distHET = data(allHET,122);

% [p,h] = kstest2(distWT, distHET)

% figure
% histogram(distWT, [1:1:7], 'FaceColor', 'k')
% hold on 
% histogram(distHET, [1:1:7], 'FaceColor', 'r')




%% PLOT TRACES OF FLASHES ON / OFF - 10 REPS

n_k = 8;

data = [Flash_Table.avSpPT, Flash_Table.Geno];
data(:, 1202) = kmeans(data(:, 1:1200), n_k);
data = sortrows(data, 1202);

figure; imagesc(data(:, 1202))
figure; imagesc(data(:, 1:1200))

all_type = find(data(:, 1201)==2);
figure; plot(nanmean(data(all_type,1:1200)))
title('7')


data = data(all_type, :);

tt = 120; 
md = data(i, 1:1200);

% for i = 1:100 

    close
%     d1 = Flash_Table{i,10};
%     d2 = Flash_Table{i,11};
%     d3 = Flash_Table{i,12};
%     md = Flash_Table{i,13};
    
    if Flash_Table.Geno(i) ==1
        col = 'k';
    else
        col = 'r';
    end
    
    %% PLOT INDIVUDAL TRACES
    
    md = data(i, 1:1200);
    if range(md)>0
        close
    tt = 120;
%     md = data(i, 1:1200);
    col = 'r';
    
    maxx = max(movmean(md,3))+1;
    minn = min(md)-1;
    rng = range([maxx minn]);
    
    figure
    rectangle('Position', [0 minn 30 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*2 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*3 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*4 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*5 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*6 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*7 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*8 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*9 minn 30 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    hold on
    %     plot(d1, 'Color', [0.9 0.9 0.9])
    %     plot(d2, 'Color', [0.9 0.9 0.9])
    %     plot(d3, 'Color', [0.9 0.9 0.9])
    n = 3;
    
    
    plot(movmean(md,3), 'Color', col, 'LineWidth', 1.2)
    ylim([minn maxx])
    ax = gca;
    ax.TickDir = 'out';
    ax.XAxis.Visible = 'off';
    ax.LineWidth = 1;
    ax.FontSize = 16;
    %     title(string(i))
    
    f = gcf;
    f.Position = [100         519        1221         165];
    i = i+1;
    else 
        i = i+1;
    end 

%%
   
title('')



















%% MEAN of allWT/HET
figure
plot(nanmean(data(allWT, 1:90)), 'k'); box off
hold on
plot(nanmean(data(allHET, 1:90)), 'r');


%% Find "Flash selective" cells - VR - 

% 0-30 = 0ff
% 31-90 = ON
% 91 - 120 = OFF

for i = 1:height(Flash_Table)
    
off_av = nanmean(data(i, 1:30));
on_fast = nanmean(data(i, 31:45));
on_all = nanmean(data(i, 31:90));
off_after= nanmean(data(i, 91:105));
% rng(i,1) = range(data(i, 1:120));

if on_fast>off_av && off_after>on_all
    Flash_Table.Type(i) = 3; %on/off
elseif on_fast<off_av && off_after>on_all
    Flash_Table.Type(i) = 2; %off
elseif on_fast>off_av &&  off_after<on_all
    Flash_Table.Type(i) = 1; %on 
% elseif on_fast < nanmean([off_av,off_after]) 
%     Flash_Table.Type(i) = 4; 
else 
    Flash_Table.Type(i) = 0; % None
end 

end 

%%

all_typeW = find(Flash_Table.Type == 0 & Flash_Table.Geno ==1 & Flash_Table.TotalSpikes >10);
all_typeH = find(Flash_Table.Type == 0 & Flash_Table.Geno ==2 & Flash_Table.TotalSpikes >10);

data = Flash_Table.OFO_NORM;
figure
subplot(1,2,1)
imagesc(data(all_typeW, :))
subplot(1,2,2)
imagesc(data(all_typeH, :))

%%
figure
for k = 0:3
    
    type = k;
    
    allWT = find(Flash_Table.Geno == 1 & Flash_Table.TotalSpikes >10 & Flash_Table.Type == k & Flash_Table.Depth > -400 & Flash_Table.Depth <0);
    allHET = find(Flash_Table.Geno == 2 & Flash_Table.TotalSpikes >10 & Flash_Table.Type == k & Flash_Table.Depth > -400 & Flash_Table.Depth <0);
    
    w = nanmean(data(allWT, 1:120));
    h = nanmean(data(allHET, 1:120));
    
    nWT = numel(allWT);
    nHET = numel(allHET);
    
%     subplot(4,1,k+1)
%     plot(movmean(w,3), 'k', 'LineWidth', 1.3)
%     hold on
%     plot(movmean(h,3), 'r', 'LineWidth', 1.3)
%     box off
%     xticks('')
%     hold on
%     plot([30 30], [-0.5 1], 'k:', 'LineWidth', 1.2)
%     plot([90 90], [-0.5 1], 'k:', 'LineWidth', 1.2)
    
    
        subplot(2,4,k+1)
        imagesc(data(allWT, 1:120))
        box off
        xticks('')
        subplot(2,4,(k+1)+4)
        imagesc(data(allHET, 1:120))
        box off
        xticks('')
    
    vals(1,k+1) = numel(allWT);
    vals(2,k+1) = numel(allHET);
    
end


%% PIE 

lb = {'None', 'On', 'Off', 'On-Off'};

X = vals(1,:);
ax1 = subplot(1,2,1);
pie(ax1,X, lb)
title(ax1,'WT');

Y = vals(2,:);
ax2 = subplot(1,2,2);
pie(ax2,Y, lb)
title(ax2,'HET');



%% Indivudal Trials

close
data = Flash_Table.avSpPT; 
md = data(i, 1:1200);

    tt = 120;
    col = 'r';
    
    maxx = max(movmean(md,3))+1;
    minn = min(md)-1;
    rng = range([maxx minn]);
    
    figure
    rectangle('Position', [0 minn 30 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    hold on
    rectangle('Position', [90 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*2 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*3 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*4 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*5 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*6 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*7 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*8 minn 60 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    rectangle('Position', [90+tt*9 minn 30 rng], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'none')
    plot(md, 'Color', col)
i=i+1;
f = gcf; 
f.Position = [ 1441         544        1996         254]; 
    
% DRIFT
% 7270 = 9, 26, 28, 33
% 7269 - 204, 684
% 1002, 1018

% NO DRIFT
% i = 1010? 




%% 

data = Flash_Table.SPAV(allWT, :);

data = Flash_Table.SPT(allWT, :);
figure; plot(nanmean(data))

figure; imagesc(data)


%% Cul3


% Plotting depth verus diffsum 

n = height(Flash_Table);

figure
for i = 1:n
    
%     if Flash_Table.Depth(i) < 0
        ani=  Flash_Table.Ani(i);
        
        if abs(Flash_Table.DifSumImmed(i))>0.05
            
            % Set colour of dot for individual mouse:
            if ani == 2869 
                col = [1 0 0];
                xval = 1;
            elseif ani == 2027 
                col = [0 0 0];
                xval = 2;
            elseif ani == 3558
                col = [0 0 0];
                xval = 3;
            elseif ani == 3557
                col = [1 0 0];
                xval = 4;
            elseif ani == 4123
                col = [0 0 0];
                xval = 5;
            elseif ani == 4124 
                col = [1 0 0];
                xval = 6;
            end
            
        elseif abs(Flash_Table.DifSumImmed(i))<=0.05
            
            if ani == 2869 
                col = [0.7 0.7 0.7];
                xval = 1;
            elseif ani == 2027 
                col = [0.7 0.7 0.7];
                xval = 2;
            elseif ani == 3558
                col = [0.7 0.7 0.7];
                xval = 3;
            elseif ani == 3557
                col = [0.7 0.7 0.7];
                xval = 4;
            elseif ani == 4123
                col = [0.7 0.7 0.7];
                xval = 5;
            elseif ani == 4124 
                col = [0.7 0.7 0.7];
                xval = 6;
            end
            
        end
        
        rndnum = (2 * rand - 1)/4;
        xval = xval + rndnum;
        yval = Flash_Table.Depth(i);
        plot(xval, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 10)
        hold on
%     end
    
end

box off 
ax= gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
ax.FontSize = 18;
xlim([0 7])
ylim([-1500 100])
f = gcf;
f.Position = [70   250   461   806];
% title('Loom')



% animal_depth_table = array2table(data, 'VariableNames', {'Ani', 'Geno', 'Depth'});















