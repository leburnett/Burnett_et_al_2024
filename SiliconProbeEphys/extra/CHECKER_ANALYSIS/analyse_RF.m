
%% Analyse temporal and spatial RFs. 
% Combines all '...GRF.mat' files to make large arrays with all info from
% all animals then analyses and plots graphs. 

% Created by Burnett 19/04/21
% Used after RFtoRFsum.m 

%% 1 - COMBINE ARRAYS TO MAKE GIANT ARRAYS WITH INFO FROM ALL ANIMALS: 

files = dir('*GRF.mat'); 
files2 = dir('*ids_info.mat');
num_files = numel(files);

all_TRF = [];
all_GTRF = [];
all_GSRF_r = []; 
all_GSRF_full = []; 
all_GRF_info = [];
all_info = []; 

for i = 1:num_files
    filename = files(i).name;
    filename2 = files2(i).name;
    load(filename)
    load(filename2)
   
    GTRF = tempRFs(good_RF,:);

    all_TRF = vertcat(all_TRF, tempRFs); 
    all_GTRF = vertcat(all_GTRF, GTRF);
    all_GSRF_r = vertcat(all_GSRF_r, GRF_r);
    all_GSRF_full = vertcat(all_GSRF_full, GRF_full);
    all_GRF_info = vertcat(all_GRF_info, GRF_info);
    all_info = vertcat(all_info, ids_info); 
end

% Add WT/HET column to the 'info' arrays. 

ng = numel(all_GRF_info(:, 1));
n = numel(all_info(:,1)); 

for j = 1:ng
    if all_GRF_info(j,9)==7270 || all_GRF_info(j,9)==7788 || all_GRF_info(j,9)==7475 ||all_GRF_info(j,9)==7616 || all_GRF_info(j,9)==2832 || all_GRF_info(j,9)==2830 || all_GRF_info(j,9)==1971
        all_GRF_info(j,11) = 1; %WT
    else 
        all_GRF_info(j,11) = 0; %HET
    end 
end

for j = 1:n
    if all_info(j,9)==7270 || all_info(j,9)==7788 || all_info(j,9)==7475 ||all_info(j,9)==7616 || all_info(j,9)==2832 || all_info(j,9)==2830 || all_info(j,9)==1971
        all_info(j,11) = 1; %WT
    else 
        all_info(j,11) = 0; %HET
    end 
end

save('210520_RF_DATA_Cul3_N3.mat', 'all_TRF', 'all_GTRF', 'all_GSRF_r', 'all_GSRF_full', 'all_GRF_info', 'all_info');





%% 2 - REALIGN TEMP RF - since some RFs seem to be amiss - need to 'realign' temp RFs. 

% In theory, the 'sensitivity' should be back to baseline at 30 since this
% is when the spike happened. 

%% Make z-scored normalised version to see differences better. 
all_GTRFz = zeros(ng, 45);

for i = 1:ng
    all_GTRFz(i, 1:45) = zscore(all_GTRF(i, 1:45)); 
end 

figure; imagesc(all_GTRFz); 
hold on; plot([15 15], [0 n], 'k:')
hold on; plot([20 20], [0 n], 'k:')
hold on; plot([25 25], [0 n], 'k:')
hold on; plot([30 30], [0 n], 'k:')
hold on; plot([35 35], [0 n], 'k:')
hold on; plot([40 40], [0 n], 'k:')


%% Make new array - only 35 frames in length not 45. 

trf = zeros(ng, 35); 

% VALUES TO SHIFT: 

% Cul3 - 2021
vals_5 = [];
vals_10 = [1:11, 13:17, 30, 32:61, 80:90, 92:96, 102, 130:134, 136:140, 157:159, 161:173];
vals_15 = [18,19, 23:28, 91, 97:101, 103:119, 121:123, 129, 125:127, 141:152, 154:156];
vals_rm = [12, 20:22, 29, 31, 62:64, 120, 124, 128, 135, 153, 160, 174:179];

% Setd5 - 2020
vals_5 = [36:67, 74:84, 86:127, 157:159, 161:163, 166, 168:171, 163:186, 203:211, 229:237, 241:249, 251, 253:262];
vals_10 = [213:215, 217:225, 227, 228];
vals_15 = [35];
vals_rm = [149, 160, 164, 165, 167, 172, 187, 194, 212, 216, 226, 238:240, 250, 252];

% Setd5 - 2021

vals_5 = [5,6,8:17, 24:27, 29:39, 41:47, 49];
vals_10 = [];
vals_15 = [];
vals_rm = [18:23, 28, 40, 48, 50, 51];


for i = 1:ng
    if ismember(i, vals_5)
       trf(i, 1:35) = all_GTRFz(i, 6:40); 
    elseif  ismember(i, vals_10)
       trf(i, 1:35) = all_GTRFz(i, 11:45); 
    elseif  ismember(i, vals_15)
        xtr = zeros(1,5);
       trf(i, 1:35) = [all_GTRFz(i, 16:45), xtr]; 
    else
      trf(i, 1:35) = all_GTRFz(i, 1:35);   
    end 
end 

% Cul3 - 20 not 35

ng = numel(all_GTRFz(:,1));
trf = zeros(ng, 20);

for i = 1:ng
    if ismember(i, vals_5)
       trf(i, 1:20) = all_GTRFz(i, 6:25); 
    elseif  ismember(i, vals_10)
       trf(i, 1:20) = all_GTRFz(i, 11:30); 
    elseif  ismember(i, vals_15)
       trf(i, 1:20) = all_GTRFz(i, 16:35); 
    else
      trf(i, 1:20) = all_GTRFz(i, 1:20);   
    end 
end 

extr = zeros(ng, 10);
extr2 = zeros(ng, 5);
trf = [extr, trf, extr2];

trf2 = trf; 
trf(vals_rm, :) = []; 


% Visualise again - CHECK ALIGNMENT. 
figure; imagesc(trf)
hold on; plot([15 15], [0 n], 'k:')
hold on; plot([20 20], [0 n], 'k:')
hold on; plot([25 25], [0 n], 'k:')
hold on; plot([30 30], [0 n], 'k:')
hold on; plot([35 35], [0 n], 'k:')
hold on; plot([40 40], [0 n], 'k:')

% Remove rows from other arrays:
all_GRF_info(vals_rm, :) = []; 
all_GSRF_full(vals_rm, :, :) = []; 
all_GSRF_r(vals_rm, :) = []; 
all_GTRFz(vals_rm, :) = []; 

save('210520_DATA_RF_aligned_RM_N3_Cul3.mat', 'trf', 'trf2', 'all_GRF_info', 'all_GSRF_full', 'all_GSRF_r', 'all_GTRFz', 'vals_5', 'vals_10', 'vals_rm');

%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% 3 - TEMPORAL - Find out whether unipolar or bipolar. 

v2 = 1.5;
n = numel(all_GRF_info(:,1));

for i = 1:n
    av_2025 = mean(trf(i, 20:25));
    av_2530 = mean(trf(i, 25:30));
    maxv = max(trf(i, 20:30));
    minv = min(trf(i, 20:30));
    
    all_GRF_info(i, 12) = av_2025;
    all_GRF_info(i, 13) = av_2530;
    if av_2530 < -1
        all_GRF_info(i, 14) = 0; 
    else
        all_GRF_info(i, 14) = 1; 
    end 
    
    all_GRF_info(i, 15) = maxv;
    all_GRF_info(i, 16) = minv;
    
    % COl 17  = bimodal?
    if maxv > v2 && minv < -v2
        all_GRF_info(i, 17) = 2; %bimodal
    elseif maxv >v2 && minv >= -v2
        all_GRF_info(i, 17) = 1; % unimodal - ON
    elseif maxv <=v2 && minv < -v2
        all_GRF_info(i, 17) = 0; % unimodal - OFF
    end 
    
    % Col 18 = GROUP - bimodal and On/off. 
    if all_GRF_info(i, 14) == 1 && all_GRF_info(i, 17)==2 % Bimodal on 
        all_GRF_info(i, 18) = 1; 
    elseif all_GRF_info(i, 14) == 0 && all_GRF_info(i, 17)==2 % Bimodal off 
        all_GRF_info(i, 18) = 2;
    elseif all_GRF_info(i, 17)==1 % uni on 
        all_GRF_info(i, 18) = 3;
    elseif all_GRF_info(i, 17)==0 % uni off
        all_GRF_info(i, 18) = 4;
    end    
        
end 

% Add the columns for group/ depth. 

sr = all_GRF_info(:,18); % GROUP 
sr2 = all_GRF_info(:,6); % DEPTH
trf3 = [trf, sr, sr2];


 %% Plots of WT/HET sorted by 'group' - multi/bipolar - on/off

% WT 
allWT = find(all_GRF_info(:, 11)==1); 
trf_WT = trf3(allWT, :); 
trf_WT = sortrows(trf_WT, 37);

figure; 
imagesc(trf_WT(:, 1:35));
colormap(redblue); hold on;
plot([30 30], [1 n], 'k:', 'LineWidth', 1.5)
title('WT')


% HET 
allHET = find(all_GRF_info(:, 11)==0); 
trf_HET = trf3(allHET, :);
trf_HET = sortrows(trf_HET, 37);

figure; 
imagesc(trf_HET(:, 1:35));
colormap(redblue); hold on;
plot([30 30], [1 n], 'k:', 'LineWidth', 1.5)
title('HET')

% Adding different y label axes. 
% ylab = num2str(all_GTRF(allWT, 51)); 
% ylab2 = num2str(all_GTRF(allHET, 51)); 
% yticks(1:1:nWT)
% yticklabels({ylab})

%
grf_info_table = array2table(all_GRF_info, 'VariableNames', {'ID', 'Amp','Channel', 'D1', 'Spikes', 'Depth', 'nna', 'Date', 'Animal', 'Trial', 'Geno', 'Av2025', 'Av2530', 'OnOff', 'MaxVal', 'MinVal', 'Biphasic', 'Group'});


%% Analysis by group.

all1 = find(grf_info_table.Group ==1); % biphasic on last
all2 = find(grf_info_table.Group ==2); % biphasic off last
all3 = find(grf_info_table.Group ==3); % uniphasic on 
all4 = find(grf_info_table.Group ==4); % uniphasic off 

allwt = find(grf_info_table.Geno ==1);
allhet = find(grf_info_table.Geno ==0);

% WT
var = grf_info_table.Depth(allwt)*-1; 
group = grf_info_table.Group(allwt); 

boxplot(var, group, 'Color', 'k')
xticks([1,2,3,4])
% xticklabels({'WT-ON', 'WT-OFF', 'HET-ON', 'HET-OFF'})
xticklabels({'BI-ON', 'BI-OFF', 'UNI-ON', 'UNI-OFF'})
xtickangle(45)
ylabel('Depth from surface (um)')
title('WT')
ax=gca;
ax.FontSize = 12; 
ylim([-2000 -200])
box off

% HET
var = grf_info_table.Depth(allhet)*-1; 
group = grf_info_table.Group(allhet); 

boxplot(var, group, 'Color', 'k')
xticks([1,2,3])
% xticklabels({'WT-ON', 'WT-OFF', 'HET-ON', 'HET-OFF'})
xticklabels({'BI-ON', 'BI-OFF', 'UNI-OFF'})
xtickangle(45)
ylabel('Depth from surface (um)')
title('HET')
ax=gca;
ax.FontSize = 12; 
ylim([-2000 -200])
box off

%% Depths of visually responsive cells - histogram

depthsWT = grf_info_table.Depth(allwt)*-1; 
depthsHET = grf_info_table.Depth(allhet)*-1; 

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

[p, h] = kstest2(depthsWT, depthsHET)
mean(depthsWT)
mean(depthsHET)

        

%% Heatmaps of TempRF - sorted by spatial RF grouping. 
wtbion = find(grf_info_table.Group ==1 & grf_info_table.Geno == 1); 
wtbioff = find(grf_info_table.Group ==2 & grf_info_table.Geno == 1); 
wtunion = find(grf_info_table.Group ==3  & grf_info_table.Geno == 1); 
wtunioff = find(grf_info_table.Group ==4  & grf_info_table.Geno == 1); 

hetbion = find(grf_info_table.Group ==1 & grf_info_table.Geno == 0); 
hetbioff = find(grf_info_table.Group ==2 & grf_info_table.Geno == 0); 
hetunion = find(grf_info_table.Group ==3  & grf_info_table.Geno == 0); 
hetunioff = find(grf_info_table.Group ==4  & grf_info_table.Geno == 0); 

TRF_wtbion = trf(wtbion, 1:35); 
TRF_wtbioff = trf(wtbioff, 1:35); 
TRF_wtunion = trf(wtunion, 1:35); 
TRF_wtunioff = trf(wtunioff, 1:35); 

TRF_hetbion = trf(hetbion, 1:35); 
TRF_hetbioff = trf(hetbioff, 1:35); 
TRF_hetunion = trf(hetunion, 1:35); 
TRF_hetunioff = trf(hetunioff, 1:35); 

figure
subplot(4,2,1)
imagesc(TRF_wtbioff)
colorbar
title('WT - bi - OFF')
colormap(redblue)
% m2 = mean(mean(TRF_wton));
% caxis([m2-(m2/2), m2+(m2/2)]) 
caxis([-5 5])

subplot(4,2,2)
imagesc(TRF_hetbioff)
colorbar
title('HET - bi - OFF')
colormap(redblue)
caxis([-5 5])

subplot(4,2,3)
imagesc(TRF_wtunion)
colorbar
title('WT -uni - ON')
colormap(redblue)
caxis([-5 5])

subplot(4,2,4)
imagesc(TRF_hetunion)
colorbar
title('HET -uni - ON')
colormap(redblue)
caxis([-5 5])

subplot(4,2,5)
imagesc(TRF_wtunioff)
colorbar
title('WT - uni - OFF')
colormap(redblue)
caxis([-5 5])

subplot(4,2,6)
imagesc(TRF_hetunioff)
colorbar
title('HET - uni - OFF')
colormap(redblue)
caxis([-5 5])

subplot(4,2,7)
imagesc(TRF_wtbion)
colorbar
title('WT - bi - ON')
colormap(redblue)
caxis([-5 5])

subplot(4,2,8)
imagesc(TRF_hetbion)
colorbar
title('HET - bi - ON')
colormap(redblue)
caxis([-5 5])


%%  MEAN + SEM 

col = 'm';

TRF_wtbioff = trf(wtbioff, 1:35); 
TRF_wtunion = trf(wtunion, 1:35); 
TRF_wtunioff = trf(wtunioff, 1:35); 

TRF_hetbioff = trf(hetbioff, 1:35); 
TRF_hetunion = trf(hetunion, 1:35); 
TRF_hetunioff = trf(hetunioff, 1:35);     

nwtbioff = numel(TRF_wtbioff(:,1));
nwtunion = numel(TRF_wtunion(:,1));
nwtunioff = numel(TRF_wtunioff(:,1));

nhetbioff = numel(TRF_hetbioff(:,1));
nhetunion = numel(TRF_hetunion(:,1));
nhetunioff = numel(TRF_hetunioff(:,1));

mean_wtbioff = mean(TRF_wtbioff);
mean_wtunion = mean(TRF_wtunion);
mean_wtunioff = mean(TRF_wtunioff);

mean_hetbioff = mean(TRF_hetbioff);
mean_hetunion = mean(TRF_hetunion);
mean_hetunioff = mean(TRF_hetunioff);

x = (1:1:35);

% bioff
semWT2off = std(TRF_wtbioff)/sqrt(nwtbioff); 
y1 = mean_wtbioff+semWT2off;
y2 = mean_wtbioff-semWT2off;

%union
semWT1on = std(TRF_wtunion)/sqrt(nwtunion); 
y3 = mean_wtunion+semWT1on;
y4 = mean_wtunion-semWT1on;

%unioff
semWT1off = std(TRF_wtunioff)/sqrt(nwtunioff); 
y5 = mean_wtunioff+semWT1off;
y6 = mean_wtunioff-semWT1off;

 % HET 

% bioff
semHET2off = std(TRF_hetbioff)/sqrt(nhetbioff); 
y1b = mean_hetbioff+semHET2off;
y2b = mean_hetbioff-semHET2off;

%union
semHET1on = std(TRF_hetunion)/sqrt(nhetunion); 
y3b = mean_hetunion+semHET1on;
y4b = mean_hetunion-semHET1on;

%unioff
semHET1off = std(TRF_hetunioff)/sqrt(nhetunioff); 
y5b = mean_hetunioff+semHET1off;
y6b = mean_hetunioff-semHET1off;

%%

% COL1 - BIPHASIC 
% COL2 - UNI ON
% COL3 - UNI OFF

% ROW 1 - WT
% ROW2 - HET
% ROW3 - BOTH 

figure

%WT BI
subplot(3,3,1)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtbioff', 'k', 'LineWidth', 1.3)
title('WT - BI')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

% WT - uni On 
subplot(3,3,2)
plot(x, y3, 'w')
hold on
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunion', 'k', 'LineWidth', 1.3)
title('WT - UNI - ON')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

% WT - uni off
subplot(3,3,3)
plot(x, y5, 'w')
hold on 
plot(x, y6, 'w')
patch([x fliplr(x)], [y5 fliplr(y6)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunioff', 'k', 'LineWidth', 1.3)
title('WT - UNI - OFF')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

%%

subplot(3,3,4)
plot(x, y1b, 'w')
hold on 
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetbioff', col, 'LineWidth', 1.3)
title('HET - BI')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

subplot(3,3,5)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunion', col, 'LineWidth', 1.3)
title('HET - UNI - ON')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

subplot(3,3,6)
plot(x, y5b, 'w')
hold on 
plot(x, y6b, 'w')
patch([x fliplr(x)], [y5b fliplr(y6b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunioff', col, 'LineWidth', 1.3)
title('HET - UNI - OFF')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

%%


subplot(3,3,7)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtbioff', 'k', 'LineWidth', 1.3)
plot(x, y1b, 'w')
hold on 
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetbioff', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

subplot(3,3,8)
plot(x, y3, 'w')
hold on
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunion', 'k', 'LineWidth', 1.3)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunion', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

subplot(3,3,9)
plot(x, y5, 'w')
hold on
plot(x, y6, 'w')
patch([x fliplr(x)], [y5 fliplr(y6)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunioff', 'k', 'LineWidth', 1.3)
plot(x, y5b, 'w')
hold on 
plot(x, y6b, 'w')
patch([x fliplr(x)], [y5b fliplr(y6b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunioff', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')


%% JUST OVERLAPPING WT/HET


figure
subplot(1,3,1)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtbioff', 'k', 'LineWidth', 1.3)
plot(x, y1b, 'w')
hold on 
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetbioff', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')
title('BI - OFF')
box off
xticks([0, 15, 30])
xticklabels({'-500', '-250', '0'})
xlabel('Time before spike - ms')
ylabel('Sensitivity (z-score)')

subplot(1,3,3)
plot(x, y3, 'w')
hold on
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunion', 'k', 'LineWidth', 1.3)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunion', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')
title('UNI - ON')
box off
xticks([0, 15, 30])
xticklabels({'-500', '-250', '0'})
xlabel('Time before spike - ms')
ylabel('Sensitivity (z-score)')

subplot(1,3,2)
plot(x, y5, 'w')
hold on
plot(x, y6, 'w')
patch([x fliplr(x)], [y5 fliplr(y6)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunioff', 'k', 'LineWidth', 1.3)
plot(x, y5b, 'w')
hold on 
plot(x, y6b, 'w')
patch([x fliplr(x)], [y5b fliplr(y6b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunioff', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')
title('UNI - OFF')
box off
xticks([0, 15, 30])
xticklabels({'-500', '-250', '0'})
xlabel('Time before spike - ms')
ylabel('Sensitivity (z-score)')



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% 4 - SPATIAL ANALYSIS 

% First look at 2D section through point of highest variance. 
% Use all_GSRF_r

num = numel(all_GSRF_r(:,1));
GRF_r2 = zeros(num, 161); 
     
    for j = 1:num
        GRF_r2(j, :) = zscore(all_GSRF_r(j,:));
        
        average_val = mean(GRF_r2(j,:));
        [max_val, maxi] = max(GRF_r2(j,:));
        [min_val, mini] = min(GRF_r2(j,:));
        
        grf_info_table.AverageValCentre(j) = average_val;
        grf_info_table.MaxVal(j) = max_val;
        grf_info_table.MaxVali(j) = maxi;
        grf_info_table.MinVal(j) = min_val;
        grf_info_table.MinVali(j) = mini;
        
        average60to100 = mean(GRF_r2(j, 70:90));
%         averageall = mean(GRF_r2(j, :));
        if average60to100 >= 0
            onoff = 1;
            grf_info_table.OnOffSRF(j) = 1;
        elseif average60to100 < 0
            onoff = 0;
            grf_info_table.OnOffSRF(j) = 0;
        end
        
        grf_info_table.Average60to100(j) = average60to100;
        
        % Finding WHP
        if onoff == 1
            maxon = max(GRF_r2(j,70:90));
            data = GRF_r2(j,40:120);
            data = data-(maxon/2);
            vals = sign(data);
            vals2 = diff(vals);
            vall = find(vals2 ~=0);
            if numel(vall) == 2
                difval = diff(vall);
            elseif numel(vall)~=2
                difval = NaN;
            end
            
        elseif onoff == 0
            minoff = min(GRF_r2(j,70:90));
            data = GRF_r2(j,40:120);
            data = data+(abs(minoff)/2);
            vals = sign(data);
            vals2 = diff(vals);
            vall = find(vals2 ~=0);
            if numel(vall) == 2
                difval = diff(vall);
            elseif numel(vall)~=2
                difval = NaN;
            end
        end
        
        grf_info_table.WHP(j) = difval;
   
    end


    
allWT_on = find(grf_info_table.Geno==1 & grf_info_table.OnOffSRF==1);  % 19
allHET_on = find(grf_info_table.Geno==0 & grf_info_table.OnOffSRF==1); % 9
allWT_off = find(grf_info_table.Geno==1 & grf_info_table.OnOffSRF==0); % 121
allHET_off = find(grf_info_table.Geno==0 & grf_info_table.OnOffSRF==0); % 97



%% Heatmap plot split by genotype and On/OFF. 

% PLOT 
figure
subplot(2,2,1)
imagesc(GRF_r2(allWT_on,:))
colorbar
title('WT - ON')
colormap(redblue)
caxis([-4 4])

subplot(2,2,2)
imagesc(GRF_r2(allWT_off,:))
colorbar
title('WT - OFF')
colormap(redblue)
caxis([-4 4])

subplot(2,2,3)
imagesc(GRF_r2(allHET_on,:))
colorbar
title('HET - ON')
colormap(redblue)
caxis([-4 4])

subplot(2,2,4)
imagesc(GRF_r2(allHET_off,:))
colorbar
title('HET - OFF')
colormap(redblue)
caxis([-4 4])

%% 

WTon = GRF_r2(allWT_on,1:161); 
WToff = GRF_r2(allWT_off,1:161);
HETon = GRF_r2(allHET_on,1:161);
HEToff = GRF_r2(allHET_off,1:161); 

nWTon = numel(allWT_on); 
nWToff = numel(allWT_off); 
nHETon = numel(allHET_on); 
nHEToff = numel(allHET_off); 

mean_WTon = mean(WTon); 
mean_WToff = mean(WToff); 
mean_HETon = mean(HETon); 
mean_HEToff = mean(HEToff); 

x = (1:1:161);

semWTon = std(WTon)/sqrt(nWTon); 
y1 = mean_WTon+semWTon;
y2 = mean_WTon-semWTon;

semWToff = std(WToff)/sqrt(nWToff); 
y1b = mean_WToff+semWToff;
y2b = mean_WToff-semWToff;
     
semHETon = std(HETon)/sqrt(nHETon); 
y3 = mean_HETon+semHETon;
y4 = mean_HETon-semHETon;

semHEToff = std(HEToff)/sqrt(nHEToff); 
y3b = mean_HEToff+semHEToff;
y4b = mean_HEToff-semHEToff;


figure
subplot(3,2,1)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WTon', 'k', 'LineWidth', 1.3)
title('WT - ON')
axis([0 161 -2 3])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')

subplot(3,2,2)
plot(x, y1b, 'w')
hold on
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WToff', 'k', 'LineWidth', 1.3)
title('WT - OFF')
axis([0 161 -3 2])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')

subplot(3,2,3)
plot(x, y3, 'w')
hold on 
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HETon', col, 'LineWidth', 1.3)
title('HET - ON')
axis([0 161 -2 3])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')

subplot(3,2,4)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HEToff', col, 'LineWidth', 1.3)
title('HET - OFF')
axis([0 161 -3 2])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')

subplot(3,2,5)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WTon', 'k', 'LineWidth', 1.3)
plot(x, y3, 'w')
hold on 
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HETon', col, 'LineWidth', 1.3)
title('ON')
axis([0 161 -2 3])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')

subplot(3,2,6)
plot(x, y1b, 'w')
hold on
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WToff', 'k', 'LineWidth', 1.3)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HEToff', col, 'LineWidth', 1.3)
title('OFF')
axis([0 161 -3 2])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')

%% Boxplot  - whp for both  + Indiv cells plotted on top. 

allWT_on = find(grf_info_table.Geno==1 & grf_info_table.OnOffSRF==1);  % 19
allHET_on = find(grf_info_table.Geno==0 & grf_info_table.OnOffSRF==1); % 9
allWT_off = find(grf_info_table.Geno==1 & grf_info_table.OnOffSRF==0); % 121
allHET_off = find(grf_info_table.Geno==0 & grf_info_table.OnOffSRF==0); % 97

for k = 1:num
  if grf_info_table.Geno(k)==1 & grf_info_table.OnOffSRF(k)==1
       grf_info_table.SRFGROUP(k) =1;
  elseif grf_info_table.Geno(k)==1 & grf_info_table.OnOffSRF(k)==0
      grf_info_table.SRFGROUP(k) =2;
  elseif grf_info_table.Geno(k)==0 & grf_info_table.OnOffSRF(k)==1
     grf_info_table.SRFGROUP(k) =3;
  elseif grf_info_table.Geno(k)==0 & grf_info_table.OnOffSRF(k)==0
      grf_info_table.SRFGROUP(k) =4;
  end 
end 



var = grf_info_table.WHP;
group = grf_info_table.SRFGROUP;

figure
boxplot(var, group, 'Color', 'k')
hold on 
xticks([1,2,3,4])
xticklabels({'WT-ON', 'WT-OFF', 'HET-ON', 'HET-OFF'})
xtickangle(45)
ylabel('Width (pixels)')
title('Width at Half Maximum')
ax = gca;
ax.FontSize = 12; 
box off


for i = 1:num
    jit = rand(1)/3-0.17; 
    
    if grf_info_table.SRFGROUP(i) == 1
        xval = 1; 
        yval = grf_info_table.WHP(i); 
        plot(xval+jit, yval, 'ko', 'MarkerSize', 7)
%         hold on 
    elseif grf_info_table.SRFGROUP(i)  == 2
        xval = 2; 
        yval =  grf_info_table.WHP(i); 
        plot(xval+jit, yval, 'ko', 'MarkerSize', 7)
%         hold on 
    elseif grf_info_table.SRFGROUP(i)  == 3
        xval = 3; 
        yval =  grf_info_table.WHP(i); 
        plot(xval+jit, yval, 'mo', 'MarkerSize', 7)
%         hold on 
    elseif grf_info_table.SRFGROUP(i) == 4
        xval = 4; 
        yval =  grf_info_table.WHP(i); 
        plot(xval+jit, yval, 'mo', 'MarkerSize', 7)
%         hold on 
    end 
end 
axis([0 5 0 55])
ax = gca;
ax.FontSize = 16; 


whp_WTOFF = grf_info_table.WHP(allWT_off);
whp_HETOFF = grf_info_table.WHP(allHET_off);

[p,h] = ranksum(whp_WTOFF, whp_HETOFF) % p = 0.0012; 
[h, p] = kstest2(whp_WTOFF, whp_HETOFF) % p = 0.0151



%% Plot values - per cell - per depth. 

figure
for j = 1:n
    
    if grf_info_table.SRFGROUP(j) == 1 || grf_info_table.SRFGROUP(j) == 3
        subplot(1,2,1)
        xval = grf_info_table.WHP(j); % whp
        yval = grf_info_table.Depth(j); %depth
        
        if grf_info_table.SRFGROUP(j) == 1 
            marker = 'ko'; 
        elseif grf_info_table.SRFGROUP(j) == 3
            marker = 'mo';
        end 
        plot(xval, -yval, marker, 'MarkerSize', 8)
        hold on 
%           axis([0 50000 -1600 -500])
        title('ON')
        ylabel('Depth from surface - um')
        xlabel('Width - pixels')
        ylim([-2000 -300])
        xlim([0 50])
        box off
        ax = gca; 
        ax.FontSize = 14;
        
    elseif grf_info_table.SRFGROUP(j) == 2 || grf_info_table.SRFGROUP(j) == 4
        subplot(1,2,2)
        xval = grf_info_table.WHP(j); % whp
        yval = grf_info_table.Depth(j); %depth
        
        if grf_info_table.SRFGROUP(j) == 2
            marker = 'ko'; 
        elseif grf_info_table.SRFGROUP(j) == 4
            marker = 'mo';
        end 
        plot(xval, -yval, marker, 'MarkerSize', 8)
        hold on 
%          axis([0 50000 -1600 -500])
        title('OFF')
        ylabel('Depth from surface - um')
        xlabel('Width - pixels')
        ylim([-2000 -300])
        box off
         xlim([0 50])
         ax2 = gca; 
        ax2.FontSize = 14;
    end 
end
sgtitle('Width at half-peak by Depth')
ax3 = gca;
ax3.FontSize = 14; 
linkaxes([ax, ax2], 'xy')






%% FULL 2D Spatial RF Analysis 

for i = 1:50
data = squeeze(all_GSRF_full(i, :, :)); 
subplot(10,5,i); imagesc(data); colormap(redblue); axis off; axis square;
end 


data = squeeze(all_GSRF_full(240, :, :)); 
figure; imagesc(data); colormap(redblue); axis off; axis square;
















