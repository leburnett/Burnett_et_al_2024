%% Script for making full arrays of good temporal and spatial receptive fields. 
% Created by Burnett 14/12/20
% Used after RFtoRFsum.m 

% 0 should be at the tip bu Tomas flipped so that 0 is at the TOP of the
% probe and 800 is at the tip. eg 100 

% Extract row from TempRFs which are in 'good_RF' add the columns from
% GRF_info. 


%% MAKE TABLE WITH INFORMATION ABOUT ALL CELLS RECORDED. 

% Open all files '*GRF.mat'
files = dir('*info.mat'); 
num_files = numel(files);

% Make array with ALL The SPATIAL ARRAYS: 

all_RF_info = [];

for i = 1:num_files
    filename = files(i).name;
    load(filename)
    all_RF_info = vertcat(all_RF_info, ids_info);  
end

   all_RF_info(:,4) = [];
   all_RF_info(:,7) = [];
       
% Add WT/HET
n_all = numel(all_RF_info(:,1));

for j = 1:n_all
    if all_RF_info(j,7)==2869 || all_RF_info(j,7)==3558 || all_RF_info(j,7)==4123 ||all_RF_info(j,7)==1389 ||all_RF_info(j,7)==2710 ||all_RF_info(j,7)==4366 ||all_RF_info(j,7)==7270 ||all_RF_info(j,7)==7616 ||all_RF_info(j,7)==7788
        all_RF_info(j,9) = 1; %WT
    else 
        all_RF_info(j,9) = 0; %HET
    end 
end 

info_table = array2table(all_RF_info, 'VariableNames', {'ID', 'Amp', 'Channel', 'Spikes', 'Depth', 'Good', 'Animal', 'Trial', 'Geno'});

save('220107_Setd5_ALL_RF_info.mat', 'all_RF_info', 'info_table')

%% ALL CELLS RECORDED

channel = all_RF_info(:,3);
num_spikes = all_RF_info(:,4);
depth = all_RF_info(:,5);
good = all_RF_info(:,6);
animal = all_RF_info(:,7);
recording = all_RF_info(:,8);
geno = all_RF_info(:,9);

all_rf_table = table(channel, num_spikes, depth, good, animal, recording, geno);
save('220107_All_Cells_Recorded_Details_Setd5.mat', 'all_rf_table');

%% STATS ALL CELLS

% Depth of ALL cells

var = all_rf_table.depth*-1; 
group = all_rf_table.animal;

[p,t,stats] = anova1(var, group)
[c,m,h] = multcompare(stats, 'CType', 'bonferroni')

% %-Visually Responsive Cells - per animal

var = all_rf_table.good; 
group = all_rf_table.animal;

[p,t,stats] = anova1(var, group)
[c,m,h] = multcompare(stats, 'CType', 'bonferroni')

% %-Visually Responsive Cells - per geno

% m(:,1) = means
wtvals = [m(2,1), m(5,1), m(6,1)];
hetvals = [m(1,1), m(3,1), m(4,1)];
[p, h] = ttest(wtvals, hetvals)

% animal averages
aa(1,1) = mean(wtvals);
aa(1,2) = 1-WTVR;

aa(2,1) = mean(hetvals);
aa(2,2) = 1-HETVR;

SUP = find(all_rf_table.depth < 1250);

% all cells
allVRW = numel(find(all_rf_table.geno(SUP) == 1 & all_rf_table.good(SUP) == 1));
allW = numel(find(all_rf_table.geno(SUP) == 1));
allVRH = numel(find(all_rf_table.geno(SUP) == 0 & all_rf_table.good(SUP) == 1));
allH = numel(find(all_rf_table.geno(SUP) == 0));

aa(1,1) = allVRW/allW;
aa(2,1) = 1-(allVRW/allW);

aa(1,2) = allVRH/allH;
aa(2,2) = 1-(allVRH/allH);

% BAR PLOT

b = bar(aa', 'Stacked', 'LineWidth', 1.2);
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
box off
ylim([0 1.1])
ax=gca;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.FontSize = 16; 
xticks([1,2])
xticklabels({'WT', 'HET'})





%% TEMPORAL RECEPTIVE FIELDS.

% Open all files '*GRF.mat'
files = dir('*GRF.mat'); 
num_files = numel(files);

%% TEMPORAL ARRAY 
% Generates array 'GTRF' for each recording - saves separately. 
% Also generates vertcat version of these arrays 'all_GTRF' which contains
% all GTRF for ALL recordings. 

% THESE ARE NOT TEMPORALLY ALIGNED> 

all_GTRF = [];

for i = 1:num_files
    filename = files(i).name;
    load(filename)
    n_cells = numel(tempRFs(:,1)); 
    
    % DEPTH IS ALREADY FINE FROM RFtoRFsum.m 
    
%     %Flipping depth back to -800. 
%     for jj = 1:n_cells
%     depth_probe = str2num(filename(end-19:end-16)); 
% %      depth_probe = str2num(filename(end-11:end-8)); 
%     d1 = tempRFs(jj, 51); 
%     % TOMAS FLIPPED - 0 is at TOP, 800 is at TIP. 
%     diff1 = 800 - d1; 
%     d2 = depth_probe -diff1; 
% %     diff1 = depth_probe - d1; 
% %     diff2 = 800-diff1; 
% %     d2 = depth_probe -diff2; 
%     tempRFs(jj,51) = d2; 
%     end 
    
    %Make GTRF. 
    GTRF = tempRFs(good_RF,:);
%     figure; plot(mean(GTRF(:, 1:45)))
%     figure; plot(mean(GRF_r(:, 1:161)))
    GTRF(:,53) = GRF_info(:,2); % Amplitude
    GTRF(:,54) = GRF_info(:,3); % Channel
    GTRF(:,55) = GRF_info(:,5); % Num Spikes
    GTRF(:,56) = GRF_info(:,10); % Trial
    
%     x = input(num2str(i));
%     close all
    
    new_filename = filename(1:end-7);
    filename_to_save = strcat(new_filename,'GTRF_LONGF.mat');
    save(filename_to_save, 'GTRF')
    
    all_GTRF = vertcat(all_GTRF, GTRF);
end

n = numel(all_GTRF(:, 1));

% Add WT/HET
% for j = 1:n
%     if all_GTRF(j,50)==7270 || all_GTRF(j,50)==7788 || all_GTRF(j,50)==7475 ||all_GTRF(j,50)==7616 || all_GTRF(j,50)==2832 || all_GTRF(j,50)==2830 || all_GTRF(j,50)==1971
%         all_GTRF(j,57) = 1; %WT
%     else 
%         all_GTRF(j,57) = 0; %HET
%     end 
% end

for j = 1:n
    if all_GTRF(j,50)==1389 || all_GTRF(j,50)==2710 || all_GTRF(j,50)==4366 ||all_GTRF(j,50)==2869 ||all_GTRF(j,50)==3558 || all_GTRF(j,50)==4123 || all_GTRF(j,50)==7270 || all_GTRF(j,50)==7616 || all_GTRF(j,50)==7788    
        all_GTRF(j,57) = 1; %WT
    else 
        all_GTRF(j,57) = 0; %HET
    end 
end

save('220107_ALL_GTRF_LONG_SETD5_N6.mat', 'all_GTRF'); %

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% REALIGN TEMP RF - since some RFs seem to be amiss - need to 'realign' temp RFs. 
% In theory, the 'sensitivity' should be back to baseline at 30 since this
% is when the spike happened. 

% Find the total number of 'good' RFs found. 
n = numel(all_GTRF(:,1));

%% Make z-scored normalised version to see differences better. 
all_GTRFz = zeros(n, 45);

for i = 1:n
    all_GTRFz(i, 1:45) = zscore(all_GTRF(i, 1:45)); 
end 

figure; imagesc(all_GTRFz(:, 1:45)); 
hold on; plot([15 15], [0 n], 'k:')
hold on; plot([20 20], [0 n], 'k:')
hold on; plot([25 25], [0 n], 'k:')
hold on; plot([30 30], [0 n], 'k')
hold on; plot([35 35], [0 n], 'k:')
hold on; plot([40 40], [0 n], 'k:')


%% Make new array - only 35 frames in length not 45. 

trf = zeros(n, 35); 

%Setd5 - 05/01/22
vals_5 = [];
vals__5 = [37:63, 65:67, 74:127, 157:163, 166, 169:176, 178:190, 192:201,206:215, 232:240, 244:254, 256:260, 266:274, 276, 278:287];
vals__10 = [216:218, 220:229, 231]; 
vals_10 = [];
vals__15 = []; 
vals_15=[];
vals__20 = [35,36];

vals_rm = [149, 164, 165, 167, 168, 177, 191, 202:205,219, 230, 241:243, 255, 261:265, 275, 277]; 


%Cul3 - 2021
% vals_5 = [18:21, 23:28, 63:65, 71:93, 95:103, 115:130];
% vals_10 = [22, 94, 134];

% vals__5 = [1:33, 45:60, 64, 78:97, 175:185, 316:336, 351:363];
% vals_5 = [98:105, 111:121, 160:174, 262:273, 282:299, 345:349];
% vals_10 = []; 
% vals_15 = []; 
% 
% vals_rm = [155:159, 375, 376, 337:344, 350]; 

% Setd5 - 2020
% vals_5 = [36:67, 74:127, 142, 171:200, 217:264, 281:312, 315:336];
% vals_10 = [35, 265:280];

% Ptchd1

% vals__5 = [1:10,12:13, 124, 160, 175:183, 201:210, 310:330];
% vals_5 = [11, 14:29,40:41,32:38 53:87, 184:196, 211:214];
% vals_10 = [88:90, 94, 100:102,112, 114, 118:120, 154, 157:159,161:164, 252:272]; 
% vals_15 = [167, 168, 365:373]; 
% 
% vals_rm = [30:31,39,44, 49,91:93, 95:99, 103, 105, 106, 152,153, 155,156, 166, 169:174, 215:219, 301, 306, 360:364, 364]; 

% vals_5 = [76:86,88:129,144,173:175, 177:179, 182, 184:187, 189:202, 219:227, 229:241, 243:266,283:291, 295:305, 307:314, 317:325, 327:338];
% vals_10 = [267:280, 282];

% for i = 1:n
%     if ismember(i, vals_5)
%        trf(i, 1:35) = all_GTRFz(i, 6:40); 
%     elseif  ismember(i, vals_10)
%        trf(i, 1:35) = all_GTRFz(i, 11:45); 
%     else
%       trf(i, 1:35) = all_GTRFz(i, 1:35);   
%     end 
% end 

for i = 1:n
    if ismember(i, vals_5)
       trf(i, 6:35) = all_GTRFz(i, 1:30); 
       trf(i, 1:5) = zeros(1,5); 
    elseif  ismember(i, vals_10)
       trf(i, 11:35) = all_GTRFz(i, 1:25); 
       trf(i, 1:10) = zeros(1,10); 
    elseif ismember(i, vals_15)
       trf(i, 16:35) = all_GTRFz(i, 1:20); 
       trf(i, 1:15) = zeros(1,15); 
    elseif ismember(i, vals__5)
        trf(i, 1:35) = all_GTRFz(i, 6:40); 
    elseif ismember(i, vals__10)
        trf(i, 1:35) = all_GTRFz(i, 11:45);
    elseif ismember(i, vals__15)
        trf(i, 1:30) = all_GTRFz(i, 16:45); 
    elseif ismember(i, vals__20)
        trf(i, 1:25) = all_GTRFz(i, 21:45); 
    else
      trf(i, 1:35) = all_GTRFz(i, 1:35);   
    end 
end 

% for i = 191:193
%     trf(i, 6:35) = trf(i, 1:30); 
%     trf(i, 1:5) = zeros(1,5);
% end 
% 
% trf(161, 1:30) = trf(161, 6:35); 
% trf(161, 31:35) = zeros(1,5);

trf(vals_rm, :) = [];
% trf([62, 135], :) = [];

% Visualise again
figure; imagesc(trf)
hold on; plot([15 15], [0 n], 'k:')
hold on; plot([20 20], [0 n], 'k:')
hold on; plot([25 25], [0 n], 'k:')
hold on; plot([30 30], [0 n], 'k')
hold on; plot([35 35], [0 n], 'k:')
hold on; plot([40 40], [0 n], 'k:')

%%

% Rm vals
%Cul3 - 2021
% vals_rm = [20, 62, 64, 109, 134, 148:153]; 
% Setd5 - 2020
% % vals_rm = [147:155, 313, 314]; 
% 
 trf2 = trf;
% % vals_rm = [149:157, 165, 176, 180, 181, 183, 188, 203, 210, 228, 242, 281, 294, 306, 315, 316, 326]; 
% trf2(vals_rm, :) = []; 
% 
% figure; imagesc(trf2)

%Ptchd1
trf_data = all_GTRF(:, 46:end);
trf_data(vals_rm, :) = []; 
% trf_data([62, 135], :) = []; 


% save('210518_Setd5_aligned_tempRF.mat', 'trf', 'trf2', 'trf_data', 'all_GTRF', 'vals_5', 'vals_10', 'vals_rm', 'trf_table');

% sorttrf = kmeans(trf2(:, 25:30), 30); 
% trf_sort = sortrows(trf2, sorttrf);
% figure; imagesc(trf_sort);

%% Find out whether unipolar or bipolar. - ADD GROUP TO TRF_TABLE. 
v2 = 1.5;

n = numel(trf_data(:,1));
for i = 1:n
    av_2025 = mean(trf2(i, 20:25));
    av_2530 = mean(trf2(i, 25:30));
    maxv = max(trf2(i, 20:30));
    minv = min(trf2(i, 20:30));
    
    trf_data(i, 13) = av_2025;
    trf_data(i, 14) = av_2530;
    if av_2530 < -1
        trf_data(i, 15) = 0; 
    else
        trf_data(i, 15) = 1; 
    end 
    
    trf_data(i, 16) = maxv;
    trf_data(i, 17) = minv;
    
    if maxv > v2 && minv < -v2
        trf_data(i, 18) = 2; %bimodal
    elseif maxv >v2 && minv >= -v2
        trf_data(i, 18) = 1; % unimodal - ON
    elseif maxv <=v2 && minv < -v2
        trf_data(i, 18) = 0; % unimodal - OFF
    end 
    
    if trf_data(i, 15) == 1 && trf_data(i, 18)==2
        trf_data(i, 19) = 1; % ON BIPOLAR
    elseif trf_data(i, 15) == 0 && trf_data(i, 18)==2
        trf_data(i, 19) = 2; % OFF BIPOLAR
    elseif trf_data(i, 18)==1
        trf_data(i, 19) = 3; % UNI ON
    elseif trf_data(i, 18)==0
        trf_data(i, 19) = 4; % UNI OFF
    end    
        
end 

sr = trf_data(:,19); 
sr2 = trf_data(:,6); 
trf3 = [trf2, sr, sr2];

% trf3 = sortrows(trf3, 36);
% figure; 
% imagesc(trf3(:, 1:35));
% colormap(redblue); hold on;
% plot([30 30], [1 n], 'k:', 'LineWidth', 1.5)


%% ADD TIME TO SPIKE FROM PEAK GIVEN GROUP
n = numel(trf_data(:,1));

for p = 1:n
    if trf_table.Group(p)==1 || trf_table.Group(p)==3  % ON 
    v = trf_table.MaxVal(p);
    t = find(trf(p,:)==v);
    v2 = 30-t(1); % spike happens at frame 30. How many frames from MIN/MAX. 
    v3 = v2*(1000/60);
    trf_table.Time2Spike(p) = v3;
    elseif trf_table.Group(p)==2 || trf_table.Group(p)==4  % OFF 
    v = trf_table.MinVal(p);
    t = find(trf(p,:)==v);
    v2 = 30-t(1); % spike happens at frame 30. How many frames from MIN/MAX. 
    v3 = v2*(1000/60);
    trf_table.Time2Spike(p) = v3;
    end
end 

%% Boxplot  - Depth per geno. 

n = numel(trf_data(:,1));

gp = find(trf_table.Depth <1250);

var = trf_table.Spikes(gp); 
group2 = trf_table.Geno(gp);

figure
hold on 
for j = 1:n
    if trf_table.Depth(j)  < 1250
        if trf_table.Geno(j) == 1
            marker = 'o';
            col = [0.5 0.5 0.5];
            sz = 90;
            xval = 2;
        elseif trf_table.Geno(j) == 0
            marker = 'o';
            sz = 90;
            %         col = [1 0 1];
            %         col = [255/255 114/255 32/255];
            col = 'r';
            xval = 1;
        end
    yval = trf_table.Spikes(j); %depth
    scatter(xval, yval,sz, col, marker, 'jitter', 'on', 'jitterAmount', 0.2)
    hold on 
    end 
end
boxplot(var, group2, 'Color', 'k', 'OutlierSize', 0.01)
set(findobj(gca,'type','line'),'linew',1.5)
xticks([1,2,3,4])
xticklabels({'HET', 'WT'})
xtickangle(45)
% ylabel('Time to Spike (ms)')
ylabel('Spikes')
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.2;
box off 
ylim([0 150])


[p,t,stats] = anova1(var, group2)
[c,m,h] = multcompare(stats, 'CType', 'bonferroni')

%% Stats - # Spikes - responsive cells

allWT = find(trf_table.Geno ==1 & trf_table.Depth < 1250);
allHET = find(trf_table.Geno ==0 & trf_table.Depth < 1250);

sWT = trf_table.Spikes(allWT);
sHET = trf_table.Spikes(allHET);

[p, h] = ranksum(sWT, sHET)
mean(sWT)
mean(sHET)

%% Depth of VR cells

allWT = find(trf_table.Geno ==1);
allHET = find(trf_table.Geno ==0);

sWT = trf_table.Depth(allWT);
sHET = trf_table.Depth(allHET);

[p, h] = ranksum(sWT, sHET)
mean(sWT)
mean(sHET)
numel(sWT)
numel(sHET)

%% Plots of WT/HET sorted by 'group' - multi/bipolar - on/off

% WT 
allWT = find(trf_data(:, 12)==1); 
trf_WT = trf3(allWT, :); 
trf_WT = sortrows(trf_WT, 37);

figure; 
imagesc(trf_WT(:, 1:35));
colormap(redblue); hold on;
plot([30 30], [1 n], 'k:', 'LineWidth', 1.5)
title('WT')


% HET 
allHET = find(trf_data(:, 12)==0); 
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


%% TABLE OF DATA 
trf_table = array2table(trf_data, 'VariableNames', {'Diff', 'cx', 'cy', 'Date', 'Animal','Depth', 'ID', 'Amp','Channel', 'Spikes', 'Trial', 'Geno', 'Av2025', 'Av2530', 'OnOff', 'MaxVal', 'MinVal', 'Biphasic', 'Group'});

save('220105_Adjusted_trf_Setd5.mat', 'trf', 'trf_data', 'trf2', 'trf3', 'vals_5', 'vals_10', 'vals_15', 'vals__5', 'vals_rm', 'trf_table');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Analysis by group.

all1 = find(trf_table.Group ==1); % biphasic on last
all2 = find(trf_table.Group ==2); % biphasic off last
all3 = find(trf_table.Group ==3); % uniphasic on 
all4 = find(trf_table.Group ==4); % uniphasic off 

allwt = find(trf_table.Geno ==1);
allhet = find(trf_table.Geno ==0);

var = trf_table.Depth(allhet)*-1; 
group = trf_table.OnOff(allhet); 

figure
boxplot(var, group)
xticks([1,2,3,4])
% xticklabels({'WT-ON', 'WT-OFF', 'HET-ON', 'HET-OFF'})
% xticklabels({'BI-ON', 'BI-OFF', 'UNI-ON', 'UNI-OFF'})
xticklabels({'ON', 'OFF'})
xtickangle(45)
ylabel('Depth from surface (um)')
title('WT')
ax=gca;
ax.FontSize = 12; 
ylim([-2000 -200])
box off


[p,t,stats] = anova1(var, group)
[c,m,h] = multcompare(stats)
        


%% Boxplot  - Depth per geno. 

n = numel(trf_data(:,1));

var = trf_table.Depth*-1; 
% group1 = trf_table.OnOff; 
group2 = trf_table.Geno;

figure
hold on 
for j = 1:n
    if trf_table.Geno(j) == 1 
        marker = 'o';
        col = [0.4 0.4 0.4];
        sz = 90;
        xval = 2;
    elseif trf_table.Geno(j) == 0
        marker = 'o';
        sz = 90; 
%         col = [1 0 1];
%         col = [255/255 114/255 32/255];
        col = 'r';
        xval = 1;
    end   
    yval = trf_table.Depth(j); %depth
    scatter(xval, -yval,sz, col, marker, 'jitter', 'on', 'jitterAmount', 0.2)
    hold on 
end
boxplot(var, group2, 'Color', 'k', 'OutlierSize', 0.01)
set(findobj(gca,'type','line'),'linew',1.5)
xticks([1,2,3,4])
xticklabels({'HET', 'WT'})
xtickangle(45)
ylabel('Depth (um)')
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.75;
box off 
ylim([-1800 -400])

[p,t,stats] = anova1(var, group2)
[c,m,h] = multcompare(stats, 'CType', 'bonferroni')


%% Boxplot  - Depth per Animal. 

n = numel(trf_data(:,1));

var = trf_table.Depth*-1; 
% group1 = trf_table.OnOff; 
group2 = trf_table.Animal;

col_het = 'r';

figure
hold on
for j = 1:n
    if trf_table.Animal(j) == 7269%1394 %2833
        marker = 'o';
%         col = [1 0 1];
%         col = [255/255 114/255 32/255]; 
            col = col_het;
        sz = 75;
        xval = 1;
    elseif trf_table.Animal(j) == 7270%1389%2869
        marker = 'o';
        sz = 75; 
        col = [0.4 0.4 0.4];
        xval = 2;
    elseif trf_table.Animal(j) == 7476%2709 %3557
        marker = 'o';
%         col = [1 0 1];
%         col = [255/255 114/255 32/255];
            col = col_het;
        sz = 75;
        xval = 3;
    elseif trf_table.Animal(j) == 7788%2710 %3558
        marker = 'o';
        sz = 75; 
        col = [0.4 0.4 0.4];
        xval = 6;
    elseif trf_table.Animal(j) == 7616%4366 %4123
        marker = 'o';
        col = [0.4 0.4 0.4];
        sz = 75;
        xval = 5;
    elseif trf_table.Animal(j) == 7614%4369%4124
        marker = 'o';
        sz = 75; 
%         col = [255/255 114/255 32/255]; 
%         col = [1 0 1];
          col = col_het;
        xval = 4;
    end   
    yval = trf_table.Depth(j); %depth
    scatter(xval, -yval,sz, col, marker, 'jitter', 'on', 'jitterAmount', 0.2)
    hold on 
end
boxplot(var, group2, 'Color', 'k')
set(findobj(gca,'type','line'),'linew',1.5)
hold on 
% xticks([1,2,3,4])
% xticklabels({'HET', 'WT'})
xtickangle(45)
ylabel('Depth (um)')
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.5;
box off 
ylim([-1800 -400])

[p,t,stats] = anova1(var, group2)
[c,m,h] = multcompare(stats)






% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Plot depth of ALL cells - plot individual cells on top and make visually responsive cells coloured. 

% WT - 7270, 7788, 7616
% HET - 7269, 7476, 7614

n = numel(all_RF_info(:,1));

var = all_rf_table.depth*-1; 
group = all_rf_table.animal;

figure
hold on
for j = 1:n
    if all_rf_table.good(j) == 1
        
        if all_rf_table.animal(j) == 7269 %1394 %2833 % HET
%             col = [1 0 1];
%             col = [255/255 114/255 32/255]; 
              col = 'r';
            xval = 1;
        elseif all_rf_table.animal(j) ==  7270%1389 %2869
            col = [0.4 0.4 0.4];
            xval = 2;
        elseif all_rf_table.animal(j) ==  7476%2709 %3557 % HET
%             col = [1 0 1];
%             col = [255/255 114/255 32/255]; 
                col = 'r'; 
            xval = 3;
        elseif all_rf_table.animal(j) == 7788 %2710 %3558
            col = [0.4 0.4 0.4];
            xval = 6;
        elseif all_rf_table.animal(j) == 7616 %4366 %4123
            col = [0.4 0.4 0.4];
            xval = 5;
        elseif all_rf_table.animal(j) == 7614 %4369 %4124 %HET
%             col = [1 0 1];
%             col = [255/255 114/255 32/255]; 
                col = 'r'; 
            xval = 4;
        end
        
    elseif all_rf_table.good(j) == 0
        
        if all_rf_table.animal(j) == 7269 %1394 %2833
            col = [0.9 0.9 0.9];
            xval = 1;
        elseif all_rf_table.animal(j) == 7270 %1389 %2869
            col = [0.9 0.9 0.9];
            xval = 2;
        elseif all_rf_table.animal(j) == 7476 %2709 %3557
            col = [0.9 0.9 0.9];
            xval = 3;
        elseif all_rf_table.animal(j) == 7788 %2710 %3558
            col = [0.9 0.9 0.9];
            xval = 6;
        elseif all_rf_table.animal(j) == 7616 %4366 %4123
            col = [0.9 0.9 0.9];
            xval = 5;
        elseif all_rf_table.animal(j) == 7614 %4369 %4124
            col = [0.9 0.9 0.9];
            xval = 4;
        end
    end
    marker = 'o';
    sz = 75;
    yval = all_rf_table.depth(j); %depth
    scatter(xval, -yval,sz, col, marker, 'jitter', 'on', 'jitterAmount', 0.2)
    hold on 
end
boxplot(var, group, 'Color', 'k')
set(findobj(gca,'type','line'),'linew',1)
hold on 

xtickangle(45)
ylabel('Depth (um)')
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.5;
box off 
ylim([-2500 0])
xlim([0 7])
%% STATS

% Depth of ALL cells

var = all_rf_table.depth*-1; 
group = all_rf_table.animal;

[p,t,stats] = anova1(var, group)
[c,m,h] = multcompare(stats, 'CType', 'bonferroni')


% Depth of VR cells

var = gsrf_info(:,5)*-1; 
group = gsrf_info(:,7);

[p,t,stats] = anova1(var, group)
[c,m,h] = multcompare(stats, 'CType', 'bonferroni')



%% ADD the 'Group' column from gsrf_info - about WT/HET/ON/OFF from spatial RF - to gtrf_table. 


% From other script - removing cells with 'RFs' outside of a normal
% position in the VF. 

% % Remove rows with RF outside of normal VF place. 
% %rowstorm = find(T2(:, 34)<200 | T2(:, 34)>500 | T2(:, 35)<50 | T2(:, 35)>300); 
%  rowstorm = find(T2(:, 34)<250 | T2(:, 34)>450 | T2(:, 35)<100 | T2(:, 35)>250); 
% T2B(rowstorm, :) = []; 

%% Heatmaps of TempRF - sorted by spatial RF grouping. 

wtbioff = find(trf_table.Group ==2 & trf_table.Geno == 1); 
wtunion = find(trf_table.Group ==3  & trf_table.Geno == 1); 
wtunioff = find(trf_table.Group ==4  & trf_table.Geno == 1); 

hetbioff = find(trf_table.Group ==2 & trf_table.Geno == 0); 
hetunion = find(trf_table.Group ==3  & trf_table.Geno == 0); 
hetunioff = find(trf_table.Group ==4  & trf_table.Geno == 0); 

TRF_wtbioff = trf2(wtbioff, 1:35); 
TRF_wtunion = trf2(wtunion, 1:35); 
TRF_wtunioff = trf2(wtunioff, 1:35); 

TRF_hetbioff = trf2(hetbioff, 1:35); 
TRF_hetunion = trf2(hetunion, 1:35); 
TRF_hetunioff = trf2(hetunioff, 1:35); 

figure
subplot(3,2,1)
imagesc(TRF_wtbioff)
colorbar
title('WT - bi - OFF')
colormap(redblue)
% m2 = mean(mean(TRF_wton));
% caxis([m2-(m2/2), m2+(m2/2)]) 
caxis([-5 5])

subplot(3,2,2)
imagesc(TRF_hetbioff)
colorbar
title('HET - bi - OFF')
colormap(redblue)
caxis([-5 5])

subplot(3,2,3)
imagesc(TRF_wtunion)
colorbar
title('WT -uni - ON')
colormap(redblue)
caxis([-5 5])

subplot(3,2,4)
imagesc(TRF_hetunion)
colorbar
title('HET -uni - ON')
colormap(redblue)
caxis([-5 5])

subplot(3,2,5)
imagesc(TRF_wtunioff)
colorbar
title('WT - uni - OFF')
colormap(redblue)
caxis([-5 5])

subplot(3,2,6)
imagesc(TRF_hetunioff)
colorbar
title('HET - uni - OFF')
colormap(redblue)
caxis([-5 5])


%%  MEAN + SEM 


TRF_wtbioff = trf2(wtbioff, 1:35); 
TRF_wtunion = trf2(wtunion, 1:35); 
TRF_wtunioff = trf2(wtunioff, 1:35); 

TRF_hetbioff = trf2(hetbioff, 1:35); 
TRF_hetunion = trf2(hetunion, 1:35); 
TRF_hetunioff = trf2(hetunioff, 1:35);     

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
% mean_hetunion = mean(TRF_hetunion);
mean_hetunion = (TRF_hetunion);
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
% col = [255/255 114/255 32/255]; % orange ptchd1
% col = 'm';
col = 'r'; 

subplot(3,3,4)
plot(x, y1b, 'w')
hold on 
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetbioff', 'Color', col, 'LineWidth', 1.3)
title('HET - BI')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

subplot(3,3,5)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunion', 'Color', col, 'LineWidth', 1.3)
title('HET - UNI - ON')
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:')

subplot(3,3,6)
plot(x, y5b, 'w')
hold on 
plot(x, y6b, 'w')
patch([x fliplr(x)], [y5b fliplr(y6b)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunioff', 'Color', col, 'LineWidth', 1.3)
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
plot(mean_hetbioff', 'Color', col, 'LineWidth', 1.3)
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
plot(mean_hetunion', 'Color', col, 'LineWidth', 1.3)
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
plot(mean_hetunioff', 'Color', col, 'LineWidth', 1.3)
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
patch([x fliplr(x)], [y1b fliplr(y2b)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetbioff', 'Color', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:', 'LineWidth', 1)
% title('BI - OFF')
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.5;
box off 


subplot(1,3,2)
plot(x, y3, 'w')
hold on
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunion', 'k', 'LineWidth', 1.3)
% plot(x, y3b, 'w')
% hold on 
% plot(x, y4b, 'w')
% patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunion', 'Color', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:', 'LineWidth', 1)
% title('UNI - ON')
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.5;
box off 

subplot(1,3,3)
plot(x, y5, 'w')
hold on
plot(x, y6, 'w')
patch([x fliplr(x)], [y5 fliplr(y6)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtunioff', 'k', 'LineWidth', 1.3)
plot(x, y5b, 'w')
hold on 
plot(x, y6b, 'w')
patch([x fliplr(x)], [y5b fliplr(y6b)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetunioff', 'Color', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:', 'LineWidth', 1)
% title('UNI - OFF')
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.5;
box off 




figure
imagesc(TRF_hetbion)
colorbar
colormap(redblue)
% m2 = mean(mean(TRF_wton));
% caxis([m2-(m2/2), m2+(m2/2)]) 
caxis([-5 5])


%% JUST BI OFF!!


figure
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_wtbioff', 'k', 'LineWidth', 1.3)
plot(x, y1b, 'w')
hold on 
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_hetbioff', 'Color', col, 'LineWidth', 1.3)
axis([0 35 -4 4])
plot([30 30], [-4 4], 'k:', 'LineWidth', 1.3)
% title('BI - OFF')
box off
xlabel('Time before spike - ms')
ylabel('Sensitivity (z-score)')
xticks([0, 15, 30])
xticklabels({'-500', '-250', '0'})
ax = gca;
ax.TickDir = 'out';
ax.FontSize = 16;
ax.LineWidth = 1.2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% PER ANIMAL PLOT
% 
% all_animals = unique(gtrf_table.Animal);
% ani = all_animals(1);
% 
% GTRF_sorted = sortrows(all_GTRF, 51);
% 
% %% Heatmaps of TempRF - sorted by spatial RF grouping. 
% 
% anion = find(gtrf_table.OnOff ==1 & gtrf_table.Animal == ani); 
% anioff = find(gtrf_table.OnOff ==0  & gtrf_table.Animal == ani); 
% %%  MEAN + SEM 
% 
% TRF = GTRF_sorted(:, 1:45); 
% 
% %zscore
% % for  j= 1:226
% %     TRF(j, :) = zscore(TRF(j,:));
% % end   
% 
% TRF_anion = TRF(anion, 1:45); 
% TRF_anioff = TRF(anioff, 1:45); 
% 
% nWTon = numel(anion); 
% nWToff = numel(anioff); 
% 
% mean_WTon = (mean(TRF_anion)); 
% mean_WToff = (mean(TRF_anioff)); 
% 
% x = (1:1:45);
% 
% semWTon = std(TRF_anion)/sqrt(nWTon); 
% y1 = mean_WTon+semWTon;
% y2 = mean_WTon-semWTon;
% 
% semWToff = std(TRF_anioff)/sqrt(nWToff); 
% y1b = mean_WToff+semWToff;
% y2b = mean_WToff-semWToff;
% 
% % Plot 
% figure
% subplot(3,2,[1,3])
% imagesc(TRF_anion)
% colorbar
% title('ON')
% colormap(redblue)
% % caxis([-3 3])
% 
% subplot(3,2,[2,4])
% imagesc(TRF_anioff)
% colorbar
% title('OFF')
% colormap(redblue)
% % caxis([-3 3])
% 
% subplot(3,2,5)
% plot(x, y1, 'w')
% hold on
% plot(x, y2, 'w')
% patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
% plot(mean_WTon', 'k', 'LineWidth', 1.3)
% title('ON')
% % axis([0 30 -3 3])
% 
% subplot(3,2,6)
% plot(x, y1b, 'w')
% hold on
% plot(x, y2b, 'w')
% patch([x fliplr(x)], [y1b fliplr(y2b)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
% plot(mean_WToff', 'k', 'LineWidth', 1.3)
% title('OFF')
% % axis([0 30 -3 3])
% 
% sgtitle(string(ani))
% 
% 
% 
% %% ALL ANI 
% i=3;
% 
% ani = all_animals(i);
% 
% % allani = find(GTRF_sorted(:,50)==ani); 
% % TRF = GTRF_sorted(:, 1:45); 
% 
% allani = find(all_GTRF(:,50)==ani); 
% TRF = all_GTRF(:, 1:45); 
% 
% depthvals = num2str(all_GTRF(allani, 51)); 
% TRF_ani = TRF(allani, 1:45);     
% 
% n_ANI = numel(allani); 
% 
% mani = mean(TRF_ani); 
% 
% x = (1:1:45);
% 
% semANI = std(TRF_ani)/sqrt(n_ANI); 
% y1b = mani+semANI;
% y2b = mani-semANI;
% 
% 
% 
% % Plot 
% figure
% subplot(4,1,[1,2,3])
% imagesc(TRF_ani)
% colorbar
% colormap(redblue)
% % caxis([-3 3])
% caxis([-60 60])
% yticks([1:1:n_ANI])
% yticklabels({depthvals})
% 
% 
% subplot(4,1,4)
% plot(x, y1b, 'w')
% hold on
% plot(x, y2b, 'w')
% patch([x fliplr(x)], [y1b fliplr(y2b)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
% plot(mani', 'k', 'LineWidth', 1.3)
% box off
% 
% % axis([0 30 -3 3])
% 
% sgtitle(string(ani))
% 
% 
%






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SPATIAL ARRAY - 1D Spatial Receptive Field - average of dissection of 2D SpatRF through x, y, and both diagonals. 

% col = [255/255 114/255 32/255]; % orange ptchd1
% col = 'm'; 
col = 'r';

% Open all files '*GRF.mat'
files = dir('*GRF.mat'); 
num_files = numel(files);

% 1 - Make array with ALL The SPATIAL ARRAYS: 
all_GSRF_FULL = [];

for i = 1:num_files
    filename = files(i).name;
    load(filename)
    
    % Normalise each cell.
    num = numel(GRF_full(:,1, 1));
    GRF_FULL_r2 = zeros(num, 161, 161); 
        
    for j = 1:num
        GRF_FULL_r2(j, :, :) = zscore(GRF_full(j,:, :));
    end
    
    all_GSRF_FULL = vertcat(all_GSRF_FULL, GRF_FULL_r2); 
end


%% 2 - 1D average across middle vertical, horizontal and both diagonals. 

n_all = numel(all_GSRF_FULL(:,1));

all_gsrf_1D = zeros(n_all, 161); 

for i = 1:n_all
        data = squeeze(all_GSRF_FULL(i, :,:));
        a1 = data(80, :);
        a2 = data(:, 80);
        a3 = diag(data);
        a4 = diag(flip(data));
        a5 = mean(horzcat(a1',a2,a3,a4)');
        all_gsrf_1D(i, :) = a5;
end 

%% 3 - Make array with info about all of the spatial arrays. 
gsrf_info = [];

for i = 1:num_files
    filename = files(i).name;
    load(filename)
    
    %Remove extra depth column and column for 'grf'
    GRF_info(:,4) = [];
    GRF_info(:,7) = [];
    
    % Normalise each cell.
    num = numel(GRF_full(:,1, 1));
    
    for j = 1:num
        
        data = squeeze(zscore(GRF_full(j,:, :))); %2D
        a1 = data(80, :);
        a2 = data(:, 80);
        a3 = diag(data);
        a4 = diag(flip(data));
        DATA = mean(horzcat(a1',a2,a3,a4)'); % 1D
        
        average_val = mean(DATA);
        [max_val, maxi] = max(max(data));
        [min_val, mini] = min(min(data));
        
        GRF_info(j, 9) = average_val;
        GRF_info(j, 10) = max_val;
        GRF_info(j, 11) = maxi;
        GRF_info(j, 12) = min_val;
        GRF_info(j, 13) = mini;
        
        % Find the value in the centre of the image - should be area with
        % highest variance - find out if on/off centre. 
        average60to100 = mean(DATA(70:90));
%         averageall = mean(GRF_r2(j, :));
        if average60to100 >= 0
            onoff = 1;
            GRF_info(j, 14) = onoff;
        elseif average60to100 < 0
            onoff = 0;
            GRF_info(j, 14) = onoff;
        end
        
        GRF_info(j, 15) = average60to100; 
        
        if onoff == 1
            maxon = max(DATA(70:90));
            d = DATA;
            d = d-(maxon/2);
            vals = sign(d);
            vals2 = diff(vals);
            vall = find(vals2 ~=0);
            if numel(vall) == 2
                difval = diff(vall);
            elseif numel(vall)~=2
                difval = NaN;
            end
            
        elseif onoff == 0
            minoff = min(DATA);
            d = DATA;
            d = d+(abs(minoff)/2);
            vals = sign(d);
            vals2 = diff(vals);
            vall = find(vals2 ~=0);
            if numel(vall) == 2
                difval = diff(vall);
            elseif numel(vall)~=2
                difval = NaN;
            end
        end
        
        GRF_info(j, 16) = difval;
        
    end
    
    gsrf_info = vertcat(gsrf_info, GRF_info);
    
end


%%

% Add WT/HET
for j = 1:n_all
    if gsrf_info(j,7)==1389 || gsrf_info(j,7)==2710 || gsrf_info(j,7)==4366 || gsrf_info(j,7)==2869 ||gsrf_info(j,7)==3558 ||gsrf_info(j,7)==4123 ||gsrf_info(j,7)==7270 ||gsrf_info(j,7)==7616 || gsrf_info(j,7)==7788
        gsrf_info(j,17) = 1; %WT
    else 
        gsrf_info(j,17) = 2; %HET
    end 
end 

% % Remove cells which were removed from TEMP RF above. 
% all_GSRF_FULL(vals_rm, :, :) = []; 
% gsrf_info(vals_rm, :) = [];

for k = 1:n_all
  if gsrf_info(k,17)==1 && gsrf_info(k,14)==1
      gsrf_info(k,18) = 1; 
  elseif gsrf_info(k,17)==1 && gsrf_info(k,14)==0
      gsrf_info(k,18) = 2; 
  elseif gsrf_info(k,17)==2 && gsrf_info(k,14)==1
      gsrf_info(k,18) = 3; 
  elseif gsrf_info(k,17)==2 && gsrf_info(k,14)==0
      gsrf_info(k,18) = 4; 
  end 
end 


allWT = find(gsrf_info(:,17)==1);
allHET = find(gsrf_info(:,17)==2);

% all_less500 = find(gsrf_info(:,4)<500); 
% all_GSRF_FULL(all_less500, :) = [];
% gsrf_info(all_less500, :) = [];

allWT_on = find(gsrf_info(:,17)==1 & gsrf_info(:,14)==1);
allHET_on = find(gsrf_info(:,17)==0 & gsrf_info(:,14)==1);
allWT_off = find(gsrf_info(:,17)==1 & gsrf_info(:,14)==2);
allHET_off = find(gsrf_info(:,17)==0 & gsrf_info(:,14)==2);

all_GSRF_WT = all_gsrf_1D(allWT,:, :);
all_GSRF_HET = all_gsrf_1D(allHET,:, :);

gsrf_info_WT = gsrf_info(allWT,:);
gsrf_info_HET = gsrf_info(allHET,:);


%%

figure; imagesc(all_GSRF_HET)



%% Sort by DEPTH

% gsrf_info
% Col 16 = whp
% 17 = geno
% 18 = group
% 14 = onoff
% 5 = depth
% 7 = animal






%% Heatmap plot split by genotype and On/OFF. 

% col = [255/255 114/255 32/255]; % orange ptchd1

% PLOT 
figure
subplot(2,2,1)
imagesc(all_gsrf_1D(allWT_on,:))
colorbar
title('WT - ON')
colormap(redblue)
caxis([-4 4])

subplot(2,2,2)
imagesc(all_gsrf_1D(allWT_off,:))
colorbar
title('WT - OFF')
colormap(redblue)
caxis([-4 4])

subplot(2,2,3)
imagesc(all_gsrf_1D(allHET_on,:))
colorbar
title('HET - ON')
colormap(redblue)
caxis([-4 4])

subplot(2,2,4)
imagesc(all_gsrf_1D(allHET_off,:))
colorbar
title('HET - OFF')
colormap(redblue)
caxis([-4 4])

%% Plot of MEAN + SEM for Geno + ON/OFF

WTon = all_gsrf_1D(allWT_on,1:161); 
WToff = all_gsrf_1D(allWT_off,1:161);
HETon = all_gsrf_1D(allHET_on,1:161);
HEToff = all_gsrf_1D(allHET_off,1:161); 

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
axis([0 161 -1 3])
box off

subplot(3,2,2)
plot(x, y1b, 'w')
hold on
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WToff', 'k', 'LineWidth', 1.3)
title('WT - OFF')
axis([0 161 -3 1])
box off

subplot(3,2,3)
plot(x, y3, 'w')
hold on 
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HETon', 'Color', col, 'LineWidth', 1.3)
title('HET - ON')
axis([0 161 -1 3])
box off

subplot(3,2,4)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HEToff', 'Color', col, 'LineWidth', 1.3)
title('HET - OFF')
axis([0 161 -3 1])
box off

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
plot(mean_HETon', 'Color', col, 'LineWidth', 1.3)
title('ON')
axis([0 161 -1 3])
box off

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
plot(mean_HEToff', 'Color', col, 'LineWidth', 1.3)
title('OFF')
axis([0 161 -3 1])
box off

%% JUST COMBINED WT/HET - LEFt = ON, RIGHT = OFF

figure
subplot(1,2,1)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WTon', 'k', 'LineWidth', 1.3)
plot(x, y3, 'w')
hold on 
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HETon', 'Color', col, 'LineWidth', 1.3)
title('ON')
axis([0 161 -1 3])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')
ax = gca;
ax.TickDir = 'out';
ax.FontSize = 16;
ax.LineWidth = 1.2;

subplot(1,2,2)
plot(x, y1b, 'w')
hold on
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WToff', 'k', 'LineWidth', 1.3)
plot(x, y3b, 'w')
hold on 
plot(x, y4b, 'w')
patch([x fliplr(x)], [y3b fliplr(y4b)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HEToff', 'Color', col, 'LineWidth', 1.3)
title('OFF')
axis([0 161 -3 1])
box off
xlabel('Pixels')
ylabel('Sensitivity (z-score)')
ax = gca;
ax.TickDir = 'out';
ax.FontSize = 16;
ax.LineWidth = 1.2;


% [p,h] = ranksum(mean_HEToff, mean_WToff)

%% Boxplot  - whp for both  + Indiv cells plotted on top. 

var = gsrf_info(:,16);
group = gsrf_info(:,18);
col = 'r'; 

figure
hold on
n = numel(gsrf_info(:,1));
for i = 1:n
    jit = rand(1)/3-0.15; 
    
    if gsrf_info(i, 18) == 1
        xval = 1; 
        yval = gsrf_info(i, 16); 
        plot(xval+jit, yval, 'Marker', 'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 10)
       %  hold on 
    elseif gsrf_info(i, 18) == 2
        xval = 2; 
        yval = gsrf_info(i, 16); 
        plot(xval+jit, yval, 'Marker', 'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 10)
%         hold on 
    elseif  gsrf_info(i, 18) == 3
        xval = 3; 
        yval = gsrf_info(i, 16); 
        plot(xval+jit, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 10)
%         hold on 
    elseif gsrf_info(i, 18) == 4
        xval = 4; 
        yval = gsrf_info(i, 16); 
        plot(xval+jit, yval, 'Marker', 'o',  'Color', col, 'MarkerSize', 10)
%         hold on 
    end 
end 

boxplot(var, group, 'Color', 'k')
set(findobj(gca,'type','line'),'linew',1.5)
hold on 
xticks([1,2,3,4])
xticklabels({'WT-ON', 'WT-OFF', 'HET-ON', 'HET-OFF'})
xtickangle(45)
ylabel('Width (Pixels)')
% title('Width at Half Maximum')

box off
axis([0 5 0 50])
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out';
ax.LineWidth = 1.1;

[p,t,stats] = anova1(var, group)
[c,m,h] = multcompare(stats, 'CType', 'bonferroni', 'Display', 'off')

%% STATS
% 
% nanmean([whp_WTon, whp_WToff])
% nanmean([whp_HETon,whp_HEToff])
% 
% range([whp_WToff])
% range([whp_HEToff])
% 
% [p,h] =ranksum(gsrf_info(allWT_off,16), gsrf_info(allHET_off, 16))
% [h,p] = kstest2(gsrf_info(allWT_off,16), gsrf_info(allHET_off, 16))
% 
% [p,h] =ranksum(gsrf_info(allWT_on,16), gsrf_info(allHET_on, 16))
% [h,p] = kstest2(gsrf_info(allWT_on,16), gsrf_info(allHET_on, 16))
% 

%[p,h] =ranksum(mean_WToff, mean_HEToff)

%% 

whp_WTon = [];
whp_WToff = []; 
whp_HETon = []; 
whp_HEToff = []; 

for i = 1:n
    if gsrf_info(i, 18) == 1
        val = gsrf_info(i, 16); 
        whp_WTon = [whp_WTon, val]; 
    elseif gsrf_info(i, 18) == 2
        val = gsrf_info(i, 16); 
        whp_WToff = [whp_WToff, val]; 
    elseif gsrf_info(i, 18) == 3
        val = gsrf_info(i, 16); 
        whp_HETon = [whp_HETon, val]; 
    elseif gsrf_info(i, 18) == 4
        val = gsrf_info(i, 16); 
        whp_HEToff = [whp_HEToff, val]; 
    end 
end

nanmean(whp_WToff)

[p,h] = ranksum(whp_WToff, whp_HEToff)


%%

% GRF info 
% Col 1 = ID
% COl 2 = Amp
% Col 3 = Channel
% Col 4 = spikes 
% Col 5 = depth 
% Col 6 = Date
% Col 7 = Animal
% Col 8 = Trial 

%%  ADD DEPTH!!! 

% ADD column of depth to values.
all_GSRF(:, 162) = gsrf_info(:,5); 

% ADD column of 'group' to values.
all_GSRF(:, 163) = gsrf_info(:,18); 

% SORT by depth! 
all_GSRF_sorted = sortrows(all_GSRF, 162);


%% PLOTS BY DEPTH

allWT_on = find(all_GSRF_sorted(:,163)==1);
allWT_off = find(all_GSRF_sorted(:,163)==2);
allHET_on = find(all_GSRF_sorted(:,163)==3);
allHET_off = find(all_GSRF_sorted(:,163)==4);

wton_depth = num2str(all_GSRF_sorted(allWT_on,162)); 
wtoff_depth = num2str(all_GSRF_sorted(allWT_off,162)); 
heton_depth = num2str(all_GSRF_sorted(allHET_on,162)); 
hetoff_depth = num2str(all_GSRF_sorted(allHET_off,162)); 

% PLOT 
figure
subplot(1,2,1)
imagesc(all_GSRF_sorted(allWT_on,1:161))
colorbar
title('WT - ON')
colormap(redblue)
yticks(1:1:numel(allWT_on))
yticklabels({(wton_depth)})

subplot(1,2,2)
imagesc(all_GSRF_sorted(allWT_off,1:161))
colorbar
title('WT - OFF')
colormap(redblue)
yticks(1:1:numel(allWT_off))
yticklabels({(wtoff_depth)})

sgtitle('WT - by depth')


figure
subplot(1,2,1)
imagesc(all_GSRF_sorted(allHET_on,1:161))
colorbar
title('HET - ON')
colormap(redblue)
yticks(1:1:numel(allHET_on))
yticklabels({(heton_depth)})

subplot(1,2,2)
imagesc(all_GSRF_sorted(allHET_off,1:161))
colorbar
title('HET - OFF')
colormap(redblue)
yticks(1:1:numel(allHET_off))
yticklabels({(hetoff_depth)})

sgtitle('HET - by depth')


%% WHP vs Depth - cool plot! 

figure
for j = 1:n_all
    
    if gsrf_info(j,18) == 1 || gsrf_info(j,18) == 3
        subplot(1,2,1)
        xval = gsrf_info(j, 16); % whp
        yval = gsrf_info(j, 5); %depth
        
        if gsrf_info(j,18)==1
            marker = 'ko'; 
        elseif gsrf_info(j,18)==3
            marker = 'ro';
        end 
        plot(xval, -yval, marker, 'MarkerSize', 9)
        hold on 
%          axis([0 50000 -1600 -500])
        title('ON')
        ylabel('Depth from surface - um')
        xlabel('Width - pixels')
        ylim([-2000 -500])
        xlim([0 50])
        box off
        ax1 = gca;
        ax1.LineWidth = 1.2;
        ax1.TickDir = 'out';
        ax1.FontSize = 18; 
       
        
    elseif gsrf_info(j,18) == 2 || gsrf_info(j,18) == 4
        subplot(1,2,2)
        xval = gsrf_info(j, 16); % whp
        yval = gsrf_info(j, 5); %depth
        
        if gsrf_info(j,18)==2
            marker = 'ko'; 
        elseif gsrf_info(j,18)==4
            marker = 'ro';
        end 
        
        plot(xval, -yval, marker, 'MarkerSize', 9)
        hold on 
%          axis([0 50000 -1600 -500])
        title('OFF')
        ylabel('Depth from surface - um')
        xlabel('Width - pixels')
        ylim([-2000 -500])
        box off
         xlim([0 50])
        ax2 = gca;
        ax2.LineWidth = 1.2;
        ax2.TickDir = 'out';
        ax2.FontSize = 18; 
    end 
end
sgtitle('Width at half-peak by Depth')










%% PER ANIMAL!!!!

% WT - 7270, 7788, 7475, 7616
% HET - 7269, 7790, 7476, 7614

ani = 4123; 
all_ANI = find(gsrf_info(:,7)==ani); 
% allANI_on = find(gsrf_info(:,7)== ani & gsrf_info(:,14)==1);
% allANI_off = find(gsrf_info(:,7)== ani & gsrf_info(:,14)==0);

all_GSRF_ANI = all_gsrf_1D(all_ANI, :); 
all_GSRF_ANI_sorted = sortrows(all_GSRF_ANI, 162); 

allANI_on = find(all_GSRF_ANI_sorted(:,163)== 1 | all_GSRF_ANI_sorted(:,163)== 3);
allANI_off = find(all_GSRF_ANI_sorted(:,163)== 2 | all_GSRF_ANI_sorted(:,163)== 4);

anion_depth = num2str(all_GSRF_ANI_sorted(allANI_on,162)); 
anioff_depth = num2str(all_GSRF_ANI_sorted(allANI_off,162)); 

%% Heatmap plot split by genotype and On/OFF - sorted by DEPTH - with depth as ylabel. 

% PLOT 
% figure
subplot(1,2,1)
imagesc(all_GSRF_ANI_sorted(allANI_on,1:161))
colorbar
title('ON')
colormap(redblue)
caxis([-3 3])
yticks(1:1:numel(allANI_on))
yticklabels({(anion_depth)})
ylabel('Depth')
hold off

subplot(1,2,2)
imagesc(all_GSRF_ANI_sorted(allANI_off,1:161))
colorbar
title('OFF')
colormap(redblue)
caxis([-3 3])
yticks(1:1:numel(allANI_off))
yticklabels({(anioff_depth)})
ylabel('Depth')

sgtitle(string(ani))
hold off


%% Plot of MEAN + SEM for Geno + ON/OFF

ANIon = all_GSRF_ANI_sorted(allANI_on,1:161); 
ANIoff = all_GSRF_ANI_sorted(allANI_off,1:161);

nANIon = numel(allANI_on); 
nANIoff = numel(allANI_off); 

mean_ANIon = mean(ANIon); 
mean_ANIoff = mean(ANIoff); 

x = (1:1:161);

semANIon = std(ANIon)/sqrt(nANIon); 
y1 = mean_ANIon+semANIon;
y2 = mean_ANIon-semANIon;

semANIoff = std(ANIoff)/sqrt(nANIoff); 
y1b = mean_ANIoff+semANIoff;
y2b = mean_ANIoff-semANIoff;
     
% PLOT 
figure
subplot(1,2,1)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_ANIon', 'k', 'LineWidth', 1.3)
title('ON')
axis([0 161 -1 3])
axis square

subplot(1,2,2)
plot(x, y1b, 'w')
hold on
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_ANIoff', 'k', 'LineWidth', 1.3)
title('OFF')
axis([0 161 -3 1])
axis square

%% COMBINED SUBPLOT WITH HEATMAPS AND MEAN + SEM!! Good plot! 
% Heatmap sorted by depth. 

% PLOT 
figure
subplot(3,2,[1,3])
imagesc(all_GSRF_ANI_sorted(allANI_on,1:161))
% colorbar
title('ON')
colormap(redblue)
caxis([-3 3])
yticks(1:1:numel(allANI_on))
yticklabels({(anion_depth)})
ylabel('Depth')

subplot(3,2,[2,4])
imagesc(all_GSRF_ANI_sorted(allANI_off,1:161))
% colorbar
title('OFF')
colormap(redblue)
caxis([-3 3])
yticks(1:1:numel(allANI_off))
yticklabels({(anioff_depth)})
ylabel('Depth')

subplot(3,2,5)
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_ANIon', 'r', 'LineWidth', 1.3)
% plot(ANIon', 'r', 'LineWidth', 1.3)
title('ON')
axis([0 161 -1 3])
% axis square

subplot(3,2,6)
plot(x, y1b, 'w')
hold on
plot(x, y2b, 'w')
patch([x fliplr(x)], [y1b fliplr(y2b)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_ANIoff', 'r', 'LineWidth', 1.3)
title('OFF')
axis([0 161 -3 1])
% axis square

sgtitle(string(ani))

%% Dot plot - WHP by Depth - for a single animal - each dot is a cell. 

% WT - 7270, 7788, 7616
% HET - 7269, 7476, 7614

all_animals = unique(gsrf_info(:,7));

for k = 1:6 
ani = all_animals(k); 

all_ANI = find(gsrf_info(:,7)==ani); 

ani_info = gsrf_info(all_ANI, :); 
num_ani = numel(ani_info(:,1)); 

allANI_on = find(ani_info(:,18)== 1 | ani_info(:,18)== 3);
allANI_off = find(ani_info(:,18)== 2 | ani_info(:,18)== 4);

subplot(1,6,k) 

for j = 1:num_ani
    
    if ani_info(j,18) == 1 
        marker = 'k.';
        sz = 30;
    elseif ani_info(j,18) == 2
        marker = 'ko';
        sz = 15; 
     elseif ani_info(j,18) == 3
        marker = 'r.';
        sz = 30; 
    elseif  ani_info(j,18) == 4
        marker = 'ro';
        sz = 15; 
    end 
        
    xval = ani_info(j, 16); % whp
    yval = ani_info(j, 5); %depth
        
    plot(xval, -yval, marker, 'MarkerSize', sz)
    hold on 
    axis([0 50 -1600 -500])
    if k ==1
        ylabel('Depth from surface - um')
        xlabel('Width')
    end
    
end
% title('Width at Half Max by depth')
ax = gca;
ax.FontSize = 12; 
ax.TickDir = 'out';
box off
ax.LineWidth = 1.2;
ylim([-2000 -500])

end 


%% STATS - Per Animal
% 1 - ANIMAL PAIRS!!!! 


% whp_WTon = [];
% whp_WToff = []; 
% whp_HETon = []; 
% whp_HEToff = []; 
% 
% for i = 1:226
%     if gsrf_info(i, 7) == 7616 || gsrf_info(i, 7) == 7614
%         if gsrf_info(i, 18) == 1
%             val = gsrf_info(i, 4); %16 is whp
%             whp_WTon = [whp_WTon, val];
%         elseif gsrf_info(i, 18) == 2
%             val = gsrf_info(i, 4);
%             whp_WToff = [whp_WToff, val];
%         elseif gsrf_info(i, 18) == 3
%             val = gsrf_info(i, 4);
%             whp_HETon = [whp_HETon, val];
%         elseif gsrf_info(i, 18) == 4
%             val = gsrf_info(i, 4);
%             whp_HEToff = [whp_HEToff, val];
%         end
%     end
% end
% 
% % ON
% stats_data = zeros(2,5); 
% 
% stats_data(1,1) = nanmean(whp_WTon); 
% stats_data(1,2) = numel(whp_WTon); 
% stats_data(1,3) = numel(find(isnan(whp_WTon))); 
% 
% stats_data(2,1) = nanmean(whp_HETon); 
% stats_data(2,2) = numel(whp_HETon); 
% stats_data(2,3) = numel(find(isnan(whp_HETon))); 
% 
% [p,h] =ranksum(whp_WTon, whp_HETon);
% [h2,p2] = kstest2(whp_WTon, whp_HETon);
% 
% stats_data(1,4) = p; 
% stats_data(1,5) = p2; 
% 
% 
% 
% % OFF
% stats_data = zeros(2,3); 
% 
% stats_data(1,1) = nanmean(whp_WToff); 
% stats_data(1,2) = numel(whp_WToff); 
% stats_data(1,3) = numel(find(isnan(whp_WToff))); 
% 
% stats_data(2,1) = nanmean(whp_HEToff); 
% stats_data(2,2) = numel(whp_HEToff); 
% stats_data(2,3) = numel(find(isnan(whp_HEToff))); 
% 
% [p,h] =ranksum(whp_WToff, whp_HEToff)
% [h2,p2] = kstest2(whp_WToff, whp_HEToff)
% 
% stats_data(1,4) = p; 
% stats_data(1,5) = p2; 
% 
% 
% %
% 
% a1 = zeros(1,3);
% a2 = zeros(1,3);
% 
% mean(a1)
% mean(a2)
% [h,p] =ranksum(a1,a2)
% [h,p] =kstest2(a1,a2)

%% % of ON/ OFF cells split into 3 depths - for GSRF only. 

gsrf_info = array2table(gsrf_info, 'VariableNames', {'ID', 'Amp', 'Channel', 'Spikes', 'Depth', 'Good', 'Animal', 'Trial', 'ONOFF', 'Geno', 'Group'});

data = zeros(6,3); 

% Depth 1. < 800
group1 = find(gsrf_info.Depth<800);
depth1 = gsrf_info(group1, :);

% Depth 2 >= 800 til <1250
group2 = find(gsrf_info.Depth>=800 & gsrf_info.Depth<1250);
depth2 = gsrf_info(group2, :);

% Depth 3 >= 1250
group3 = find(gsrf_info.Depth>=1250);
depth3 = gsrf_info(group3, :);


%% Group 1 

tbl = zeros(5,2); 

% total cells Wt/HET
allWT = find(depth1.Geno == 1); 
allHET = find(depth1.Geno == 0); 

tbl(1,1) =numel(allWT);
tbl(1,2) =numel(allHET);

% num cells on 
% num cells off 
allWTon = find(depth1.Group == 1); 
allWToff = find(depth1.Group == 2); 
allHETon = find(depth1.Group == 3); 
allHEToff = find(depth1.Group == 4); 

tbl(2,1) = numel(allWTon);
tbl(3,1) = numel(allWToff);
tbl(2,2) = numel(allHETon);
tbl(3,2) = numel(allHEToff);

% %on/off
tbl(4,1) = (numel(allWTon)/ numel(allWT))*100; 
tbl(4,2) = (numel(allHETon)/ numel(allHET))*100; 
tbl(5,1) = (numel(allWToff)/ numel(allWT))*100; 
tbl(5,2) = (numel(allHEToff)/ numel(allHET))*100; 






%% Test to see what RF looks like if squeeze in x/y

% i = 11; 
% figure
% data = squeeze(all_GSRF_FULL(i, :,:));
% subplot(4,1,1)
% imagesc(data)
% axis off
% colormap(redblue)
% 
% subplot(4,1,2) % % % % % % THIS SEEMS THE 'BEST' VIEW - CAPTURES RF. 
% a1 = mean(squeeze(all_GSRF_FULL(i, :,:)));
% imagesc(a1)
% colormap(redblue)
% axis off
% 
% subplot(4,1,3)
% a2 = mean(squeeze(all_GSRF_FULL(i, :,:))');
% imagesc(a2)
% axis off
% colormap(redblue)
% 
% subplot(4,1,4)
% a3 = a1.*a2; 
% imagesc(a3)
% axis off
% colormap(redblue)


%%
% 
% i = 11; 
% figure
% data = squeeze(all_GSRF_FULL(i, :,:));
% subplot(4,1,1)
% imagesc(data)
% axis square
% axis off
% colormap(redblue)
% 
% subplot(4,1,2) % % % % % % THIS SEEMS THE 'BEST' VIEW - CAPTURES RF. 
% a1 = data(50:110, :);
% imagesc(a1)
% colormap(redblue)
% axis square
% axis off
% 
% subplot(4,1,3)
% a2 = data(:, 50:110);
% imagesc(a2)
% axis square
% axis off
% colormap(redblue)
% 
% subplot(4,1,4)
% a3 = diag(data');
% imagesc(a3)
% axis square
% axis off
% colormap(redblue)

%%
%%% 
i = 5; 
figure
data = squeeze(all_GSRF_FULL(i, :,:));
subplot(5,1,1)
a1 = data(80, :);
imagesc(a1)
axis square
axis off
colormap(redblue)

subplot(5,1,2)
a2 = data(:, 80);
imagesc(a2');
axis square
axis off
colormap(redblue)

subplot(5,1,3)
a3 = diag(data);
imagesc(a3')
axis square
axis off
colormap(redblue)

subplot(5,1,4)
a4 = diag(flip(data));
imagesc(a4')
axis square
axis off
colormap(redblue)

subplot(5,1,5)
a5 = mean(horzcat(a1',a2,a3,a4)');
imagesc(a5)
axis square
axis off
colormap(redblue)

% a5 is what you want or Spatial RF!!



%% Plot individual RFs

% 
% for i = 1:50
% data = squeeze(all_GSRF_FULL(i,:,:));
% data2 = (all_GSRF_FULL(i,1:160));
% imagesc(data)
% xticks([])
% yticks([])
% colormap(redblue)
% m2 = mean(data2);
% caxis([m2-(m2/2), m2+(m2/2)])
% depth = gsrf_info(i,5);
% title(string(depth))
% box off
% axis square
% axis off
% % x = input('Type');
% hold off
% end 










































%% PLOT  - Make a plot of the 2D spatial receptive fields PER ANIMAL. 

% PER ANIMAL
ani = 3557; 
all_ANI = find(gsrf_info(:,7)==ani); 

num_RF = numel(all_ANI); 

all_GSRF_ANI = all_GSRF_FULL(all_ANI, :, :); 
info_ANI = gsrf_info(all_ANI, :); 

figure
for i = 1:num_RF
    subplot(11,8,i)
    data = squeeze(all_GSRF_ANI(i,:,:));
    data2 = (all_GSRF_ANI(i,1:160)); 
    imagesc(data)
    xticks([])
    yticks([])
%     colorbar
    colormap(redblue)
    m2 = mean(data2);
    caxis([m2-(m2/2), m2+(m2/2)]) 
    depth = info_ANI(i,5); 
    title(string(depth))
end 
sgtitle(string(ani))




%%
% Depth 1. < 800
group1 = find(info_table.Depth<800);
depth1 = info_table(group1, :);

% Depth 2 >= 800 til <1250
group2 = find(info_table.Depth>=800 & info_table.Depth<1250);
depth2 = info_table(group2, :);

% Depth 3 >= 1250
group3 = find(info_table.Depth>=1250);
depth3 = info_table(group3, :);


%% For Info_table

tbl = zeros(14,1); 

allWT = find(info_table.Geno == 1); 
allHET = find(info_table.Geno == 0); 

allWTgood = find(info_table.Geno == 1 & info_table.Good ==1); 
allHETgood = find(info_table.Geno == 0 & info_table.Good ==1); 

tbl(1,1) = numel(allWT);
tbl(2,1) = numel(allHET);

tbl(3,1) = numel(allWTgood);
tbl(4,1) = numel(allHETgood);

tbl(5,1) = (numel(allWTgood)/ numel(allWT))*100; 
tbl(6,1) = (numel(allHETgood)/ numel(allHET))*100; 

tbl(7,1) = mean(info_table.Spikes(allWT));
tbl(8,1) = mean(info_table.Spikes(allHET));

tbl(9,1) = range(info_table.Spikes(allWT));
tbl(10,1) = range(info_table.Spikes(allHET));

% Depth

tbl(11,1) = mean(info_table.Depth(allWT));
tbl(12,1) = mean(info_table.Depth(allHET));

tbl(13,1) =range(info_table.Depth(allWT));
tbl(14,1) =range(info_table.Depth(allHET));

%% Group 1 

tbl = zeros(14,1); 

allWT = find(depth3.Geno == 1); 
allHET = find(depth3.Geno == 0); 

allWTgood = find(depth3.Geno == 1 & depth3.Good ==1); 
allHETgood = find(depth3.Geno == 0 & depth3.Good ==1); 

tbl(1,1) = numel(allWT);
tbl(2,1) = numel(allHET);

tbl(3,1) = numel(allWTgood);
tbl(4,1) = numel(allHETgood);

tbl(5,1) = (numel(allWTgood)/ numel(allWT))*100; 
tbl(6,1) = (numel(allHETgood)/ numel(allHET))*100; 

tbl(7,1) = mean(depth3.Spikes(allWT));
tbl(8,1) = mean(depth3.Spikes(allHET));

tbl(9,1) = range(depth3.Spikes(allWT));
tbl(10,1) = range(depth3.Spikes(allHET));

% Depth

tbl(11,1) = mean(depth3.Depth(allWT));
tbl(12,1) = mean(depth3.Depth(allHET));

tbl(13,1) =range(depth3.Depth(allWT));
tbl(14,1) =range(depth3.Depth(allHET));


%% Per ANiMAL

group3 = find(info_table.Animal==7790);
depth3 = info_table(group3, :);

tbl = zeros(14,1); 

allWT = numel(depth3(:,1)); 
allWTgood = find(depth3.Good ==1); 

tbl(1,1) = (allWT);

tbl(3,1) = numel(allWTgood);

tbl(5,1) = (numel(allWTgood)/ (allWT))*100; 


tbl(7,1) = mean(depth3.Spikes);

tbl(9,1) = range(depth3.Spikes);

% Depth

tbl(11,1) = mean(depth3.Depth);

tbl(13,1) =range(depth3.Depth);

 
% for i = 1:2005
%     if info_table.Geno(i) == 1
%     marker = 'k.';
%     elseif info_table.Geno(i) == 0 
%         marker = 'r.';
%     end 
%     
%     xval = info_table.Spikes(i);
%     yval = -info_table.Depth(i);
%     
%     plot(xval, yval, marker, 'MarkerSize', 20)
%     hold on 
% end 


%% Channel Info
tbl = zeros(32,2); 

tbl(:,1) = unique(info_table.Channel);

for i = 1:32
    ch = tbl(i,1); 
    allch = find(info_table.Channel == ch);
    val = numel(allch); 
    tbl(i,2) = val;
    
    valsgood = find(info_table.Good(allch)==1); 
    val2 = numel(valsgood);
    tbl(i,3) = val2; 
    
    chspikes = mean(info_table.Spikes(allch)); 
    tbl(i,4) = chspikes; 
    
    champ = mean(info_table.Amp(allch)); 
    tbl(i,5) = champ; 
    
end 
    
%     
% a1 = [11.4, 13.9, 11.5];
% a2 = [7.2, 11.0, 6.91]; 
% [p,h] = ttest(a1,a2);
% p = 1;
% h = 0.0168;

%% BOXPLOT
% Animal average WHP. Manually entered data from number sheet. 

% data = zeros(12,2);
% var = data(:,1);
% group = data(:,2);
% 
% boxplot(var, group, 'Color', 'k')
% hold on 
% for i = 1:12
%     xval = data(i,2);
%     yval = data(i,1); 
%     
%     if xval ==1 || xval ==2 
%         marker = 'ko';
%     else 
%         marker = 'ro';
%     end 
%     
%     plot(xval, yval, marker, 'MarkerSize', 8)
% end 
% 
% [p,h, stats] = anova1(var, group);
% multcompare(stats)
% 
% wton = find(data(:,2)==1);
% wtoff = find(data(:,2)==2);
% heton = find(data(:,2)==3);
% hetoff = find(data(:,2)==4);
% 
% vals1 = data(wton, 1);
% vals2 = data(heton,1); 
% 
% [p,h] =ranksum(vals1, vals2)
% 
% mean(vals1)
% mean(vals2)





%% OLD 

% % Open all files '*GRF.mat'
% files = dir('*GRF.mat'); 
% num_files = numel(files);
% 
% % Make array with ALL The SPATIAL ARRAYS: 
% 
% all_GSRF = [];
% gsrf_info = [];
% 
% for i = 1:num_files
%     filename = files(i).name;
%     load(filename)
%     
%     %Remove extra depth column and column for 'grf'
%     GRF_info(:,4) = [];
%     GRF_info(:,7) = [];
%     
%     % Normalise each cell.
%     num = numel(GRF_r(:,1));
%     GRF_r2 = []; 
%         
%     for j = 1:num
%         GRF_r2(j, :) = zscore(GRF_r(j,:));
%         
%         average_val = mean(GRF_r2(j,:));
%         [max_val, maxi] = max(GRF_r2(j,:));
%         [min_val, mini] = min(GRF_r2(j,:));
%         
%         GRF_info(j, 9) = average_val;
%         GRF_info(j, 10) = max_val;
%         GRF_info(j, 11) = maxi;
%         GRF_info(j, 12) = min_val;
%         GRF_info(j, 13) = mini;
%         
%         average60to100 = mean(GRF_r2(j, 70:90));
% %         averageall = mean(GRF_r2(j, :));
%         if average60to100 >= 0
%             onoff = 1;
%             GRF_info(j, 14) = 1;
%         elseif average60to100 < 0
%             onoff = 0;
%             GRF_info(j, 14) = 0;
%         end
%         
%         GRF_info(j, 15) = average60to100; 
%         
%         if onoff == 1
%             maxon = max(GRF_r2(j,70:90));
%             data = GRF_r2(j,40:120);
%             data = data-(maxon/2);
%             vals = sign(data);
%             vals2 = diff(vals);
%             vall = find(vals2 ~=0);
%             if numel(vall) == 2
%                 difval = diff(vall);
%             elseif numel(vall)~=2
%                 difval = NaN;
%             end
%             
%         elseif onoff == 0
%             minoff = min(GRF_r2(j,70:90));
%             data = GRF_r2(j,40:120);
%             data = data+(abs(minoff)/2);
%             vals = sign(data);
%             vals2 = diff(vals);
%             vall = find(vals2 ~=0);
%             if numel(vall) == 2
%                 difval = diff(vall);
%             elseif numel(vall)~=2
%                 difval = NaN;
%             end
%         end
%         
%         GRF_info(j, 16) = difval;
%         
%     end
%     
%     all_GSRF = vertcat(all_GSRF, GRF_r2);
%     gsrf_info = vertcat(gsrf_info, GRF_info);
%     
% end
% 
% 
% n_all = numel(all_GSRF(:,1));
% 
% % Add WT/HET
% % for j = 1:n_all
% %     if gsrf_info(j,7)==7270 || gsrf_info(j,7)==7788 || gsrf_info(j,7)==7475 ||gsrf_info(j,7)==7616 ||gsrf_info(j,7)==2832 || gsrf_info(j,7)==2830
% %         gsrf_info(j,17) = 1; %WT
% %     else 
% %         gsrf_info(j,17) = 0; %HET
% %     end 
% % end 
% 
% for j = 1:n_all
% %         if all_GTRF(j,50)==1389 || all_GTRF(j,50)==2710 || all_GTRF(j,50)==4366 ||all_GTRF(j,50)==2869 ||all_GTRF(j,50)==3558 || all_GTRF(j,50)==4123 || all_GTRF(j,50)==7270 || all_GTRF(j,50)==7616 || all_GTRF(j,50)==7788    
%     if gsrf_info(j,7)==1389 || gsrf_info(j,7)==2710 || gsrf_info(j,7)==4366 || gsrf_info(j,7)==2869 ||gsrf_info(j,7)==3558 ||gsrf_info(j,7)==4123 ||gsrf_info(j,7)==7270 ||gsrf_info(j,7)==7616 || gsrf_info(j,7)==7788
%         gsrf_info(j,17) = 1; %WT
%     else 
%         gsrf_info(j,17) = 0; %HET
%     end 
% end
% 
% 
% for k = 1:n_all
%   if gsrf_info(k,17)==1 && gsrf_info(k,14)==1
%       gsrf_info(k,18) = 1; 
%   elseif gsrf_info(k,17)==1 && gsrf_info(k,14)==0
%       gsrf_info(k,18) = 2; 
%   elseif gsrf_info(k,17)==0 && gsrf_info(k,14)==1
%       gsrf_info(k,18) = 3; 
%   elseif gsrf_info(k,17)==0 && gsrf_info(k,14)==0
%       gsrf_info(k,18) = 4; 
%   end 
% end 
% 
% save('211210_ALL_GSRF_Cul3_N6.mat', 'all_GSRF', 'gsrf_info');
% 
% allWT = find(gsrf_info(:,17)==1);
% allHET = find(gsrf_info(:,17)==0);
% 
% % all_less500 = find(gsrf_info(:,4)<500); 
% % all_GSRF(all_less500, :) = [];
% % gsrf_info(all_less500, :) = [];
% 
% allWT_on = find(gsrf_info(:,17)==1 & gsrf_info(:,14)==1);
% allHET_on = find(gsrf_info(:,17)==0 & gsrf_info(:,14)==1);
% allWT_off = find(gsrf_info(:,17)==1 & gsrf_info(:,14)==0);
% allHET_off = find(gsrf_info(:,17)==0 & gsrf_info(:,14)==0);
% 
% % all_GSRF_WT = all_GSRF(allWT,:);
% % all_GSRF_HET = all_GSRF(allHET,:);
% % 
% % gsrf_info_WT = gsrf_info(allWT,:);
% % gsrf_info_HET = gsrf_info(allHET,:);
% 
% 
% 











