% Combine LOOM analysis arrays:
% Cretaed by Burnett - 28/02/22
clear
close all

% LOAD table with information about the depth of the sSC in different animals.
% load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Setd5_Rec_Depth_Info.mat', 'animal_depth_info_table')

% ANALYSIS RESULTS FILES STORED HERE:
% "/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/SPIKESS/5Looms/Loom_Spike_Res"

%% For all the analysis files want to:

het_animals = [7269, 7476, 7614];

Loom_Table = table();

files = dir('*LOOM*');
nfiles =numel(files);

for i = 1:nfiles
    
    fname = files(i).name;
    load(fname)
    
    % Ammend the depth of each cell!
    animal_num = fname(8:11);
    row = find(animal_depth_info_table.Animal == str2num(animal_num));
%     sSC_depth = animal_depth_info_table.sSC_Depth(row); % depth of beginning of sSC
        sSC_depth = animal_depth_info_table.SC_DEPTH_NEW(row); %     

    probe = probe_depths(1); %probe depth
    n_gcl = numel(ids_depth(:,1)); %number of units in this recording
    
    for jj = 1:n_gcl
        d = probe - (800 - ids_depth(jj, 2));
%         d2 = sSC_depth - d; 
        d2 = (sSC_depth*-1)-d;
        real_depth(jj,1) = d2;
        
        if ismember(str2num(animal_num), het_animals)
            genoo(:, 1) = ones(n_gcl, 1)*2;
        else
            genoo(:, 1) = ones(n_gcl, 1);
        end
        
        % Total number of spikes over the stimulus  - used to find low
        % firing cells.
        total_sp(jj, 1) = sum(L251(jj, :)) + sum(L252(jj,:));
        
    end
    
    % Make TABLE
    Date = dates;
    Ani = anis;
    Exp = tests;
    Geno = genoo;
    ProbeDepth = probe_depths;
    
    Depth = real_depth;
    TotalSpikes = total_sp;
    P_VALUE = pvals_responsive_alllooms;
    
    tbl = table(Date, Ani, Exp, Geno, P_VALUE, ProbeDepth, TotalSpikes, Depth, L1REP, L1REP_norm, L25_ST1, L25_ST2, L251, L252, L25_norm1, L25_norm2);
    
    Loom_Table = vertcat(Loom_Table, tbl);
    
    clearvars tbl Date Ani Exp Geno ProbeDepth TotalSpikes PVALS Depth real_depth total_sp genoo ids_depth
    
end

% save('220729_Setd5_Looom_table_PVAL.mat', 'Loom_Table', 'animal_depth_info_table');


%% Generate 'Loom_Variables' table to quantify max / average firing and time to max

wid = 5; 

% Loom on at 60 - 

for i = 1:height(Loom_Table)
    
    % Only the first loom presentation - not averaged across all loom
    % presentations.
    
    l11 = Loom_Table.L251(i,20:350);
    l12 = Loom_Table.L252(i,20:350);
    
    l1av = mean(vertcat(l11, l12));
    maxduringL1 = max(l1av(40:80));
    row_max = find(l1av(40:80)==max(l1av(40:80)));
    if numel(row_max)>1
        row_max = row_max(1);
    end
    
    % mean max - 1 frame either side of max.
    mmax = nanmean(l1av((row_max-1)+40:(row_max+1)+40));
    
    %Time to peak
    t2max = (row_max)/60; 
    
    Loom_Table.L1ONLY_rep1(i, :) = l11; % 20:109 = First loom!
    Loom_Table.L1ONLY_rep2(i, :) = l12;
    Loom_Table.L1ONLY_AV(i,:) = l1av;
%     Loom_Table.L1ONLY_rep1N(i, :) = Loom_Table.L25_norm1(i,20:350); % 20:109 = First loom!
%     Loom_Table.L1ONLY_rep2N(i, :) = Loom_Table.L25_norm2(i,20:350);
%     
    % Find average firing BEFORE LOOMS to then look at the change in firing.
    b41 = [Loom_Table.L251(i, 1:60),Loom_Table.L251(i, 300:360),Loom_Table.L251(i, 610:680),Loom_Table.L251(i, 925:1000),Loom_Table.L251(i, 1245:1325)];
    b42 = [Loom_Table.L252(i, 1:60),Loom_Table.L252(i, 300:360),Loom_Table.L252(i, 610:680),Loom_Table.L252(i, 925:1000),Loom_Table.L252(i, 1245:1325)];
    
    b4 = vertcat(b41, b42);
    avb4 = mean(b4);
    AVB4 = mean(avb4);
    
    Loom_Table.B4LOOMS1(i,:) = b41;
    Loom_Table.B4LOOMS2(i,:) = b42;
    Loom_Table.MeanBefore(i, :) = avb4;
    
    Loom_Table.AVB4(i) = AVB4;
    Loom_Table.MaxL1(i) = mmax;
    Loom_Table.T2Max(i) = t2max; % Time in s
    Loom_Table.DeltaLoom(i) = maxduringL1/AVB4;
    
    % Average over 5 loom reps: 
    rep1  =Loom_Table.L251(i, 1:1590);
    rep2 = Loom_Table.L252(i, 1:1590);
    reps = vertcat(rep1, rep2);
    av25 = mean(reps);
    Loom_Table.L25AV(i, 1:1590) = av25;

end


Loom_Variables = table();

wid = 2; 
vv = [51, 93, 136, 179, 222];
ls = [40, 84,  128, 172, 216]; % 44 loom 
lr = 44;

for i = 1:height(Loom_Table)
    
    % Average of 2 * 5 Looms
    data = Loom_Table.L1ONLY_AV(i,:);
    
    % Max1 
    m1 = nanmean(data(vv(1)-wid:vv(1)+wid));
    av1 = nanmean(data(ls(1):ls(1)+lr));
    
    % Max2
    m2 = nanmean(data(vv(2)-wid:vv(2)+wid));
    av2 = nanmean(data(ls(2):ls(2)+lr));
    
    % Max3
    m3 = nanmean(data(vv(3)-wid:vv(3)+wid));
    av3 = nanmean(data(ls(3):ls(3)+lr));
    
    %Max4
    m4 = nanmean(data(vv(4)-wid:vv(4)+wid));
    av4 = nanmean(data(ls(4):ls(4)+lr));
    
    %Max5
    m5 = nanmean(data(vv(5)-wid:vv(5)+wid));
    av5 = nanmean(data(ls(5):ls(5)+lr));
    
    Loom_Variables.m1(i) = m1;
    Loom_Variables.m2(i) = m2;
    Loom_Variables.m3(i) = m3;
    Loom_Variables.m4(i) = m4;
    Loom_Variables.m5(i) = m5;
    Loom_Variables.maxAV(i) = nanmean([m1, m2, m3, m4,m5]);
    
    Loom_Variables.av1(i) = av1;
    Loom_Variables.av2(i) = av2;
    Loom_Variables.av3(i) = av3;
    Loom_Variables.av4(i) = av4;
    Loom_Variables.av5(i) = av5;
    Loom_Variables.avAV(i) = nanmean([av1, av2, av3, av4, av5]);
    
end 


%% Max over the 5 looms - Fig. S4p

allWT = find(Loom_Table.Geno == 1 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<0 & Loom_Table.Depth> -1000);
allHET = find(Loom_Table.Geno == 2 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<0 & Loom_Table.Depth> -1000);

vals_WT = zeros(1,5);
vals_HET = zeros(1,5);
sem_WT = zeros(1,5);
sem_HET = zeros(1,5);

for k = 1:height(Loom_Variables)
    
    vals_WT(1) = nanmean(Loom_Variables.av1(allWT)*60);
    vals_WT(2) = nanmean(Loom_Variables.av2(allWT)*60);
    vals_WT(3) = nanmean(Loom_Variables.av3(allWT)*60);
    vals_WT(4) = nanmean(Loom_Variables.av4(allWT)*60);
    vals_WT(5) = nanmean(Loom_Variables.av5(allWT)*60);
    
    vals_HET(1) = nanmean(Loom_Variables.av1(allHET)*60);
    vals_HET(2) = nanmean(Loom_Variables.av2(allHET)*60);
    vals_HET(3) = nanmean(Loom_Variables.av3(allHET)*60);
    vals_HET(4) = nanmean(Loom_Variables.av4(allHET)*60);
    vals_HET(5) = nanmean(Loom_Variables.av5(allHET)*60);
    
    sem_WT(1) = nanstd(Loom_Variables.av1(allWT)*60)/sqrt(numel(allWT));
    sem_WT(2) = nanstd(Loom_Variables.av2(allWT)*60)/sqrt(numel(allWT));
    sem_WT(3) = nanstd(Loom_Variables.av3(allWT)*60)/sqrt(numel(allWT));
    sem_WT(4) = nanstd(Loom_Variables.av4(allWT)*60)/sqrt(numel(allWT));
    sem_WT(5) = nanstd(Loom_Variables.av5(allWT)*60)/sqrt(numel(allWT));
    
    sem_HET(1) = nanstd(Loom_Variables.av1(allHET)*60)/sqrt(numel(allHET));
    sem_HET(2) = nanstd(Loom_Variables.av2(allHET)*60)/sqrt(numel(allHET));
    sem_HET(3) = nanstd(Loom_Variables.av3(allHET)*60)/sqrt(numel(allHET));
    sem_HET(4) = nanstd(Loom_Variables.av4(allHET)*60)/sqrt(numel(allHET));
    sem_HET(5) = nanstd(Loom_Variables.av5(allHET)*60)/sqrt(numel(allHET));
end 

figure
errorbar(1:1:5, vals_WT, sem_WT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:5, vals_HET, sem_HET, 'o', 'CapSize', 0, 'Color', [1 0.8 0.8], 'MarkerFaceColor', [1 0.8 0.8], 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75)
errorbar(1:1:5, vals_WT, sem_WT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:5, vals_HET, sem_HET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.4],'Marker', 'none', 'LineWidth', 1.75)

xlim([0.5 5.5])
ylim([0 10])
box off
ax = gca;
ax.LineWidth = 1.2;
ax.TickDir = 'out'; 
ax.TickLength = [0.02 0.02];
f = gcf;
f.Position = [680   867   268   231]; 

% STATS - repeated measures ANOVA
nWT = numel(allWT);
nHET = numel(allHET);

WT1 = (Loom_Variables.av1(allWT)*60);
WT2 = (Loom_Variables.av2(allWT)*60);
WT3 = (Loom_Variables.av3(allWT)*60);
WT4 = (Loom_Variables.av4(allWT)*60);
WT5 = (Loom_Variables.av5(allWT)*60);

HET1 = (Loom_Variables.av1(allHET)*60);
HET2 = (Loom_Variables.av2(allHET)*60);
HET3 = (Loom_Variables.av3(allHET)*60);
HET4 = (Loom_Variables.av4(allHET)*60);
HET5 = (Loom_Variables.av5(allHET)*60);

y = vertcat(WT1, WT2, WT3, WT4, WT5, HET1, HET2, HET3, HET4, HET5);
gp1 = horzcat(ones(1, nWT*5), ones(1, nHET*5)*2);
gp2 = horzcat(ones(1, nWT), ones(1, nWT)*2,  ones(1, nWT)*3,  ones(1, nWT)*4,  ones(1, nWT)*5,  ones(1, nHET),  ones(1, nHET)*2, ones(1, nHET)*3, ones(1, nHET)*4, ones(1, nHET)*5);


% n way ANOVA
[p, tbl, stats] = anovan(y, {gp1', gp2'}, 'model', 'interaction', 'varnames', {'Geno', 'Loom'})
results = multcompare(stats, 'Dimension', [1,2], 'CType', 'bonferroni')

%
[p, h] = ranksum(WT1, HET1)


%% BOXPLOT for Max and Average firing - Figs. S4k,l

val = 12; % 12 = Average firing (Hz), % 6 = max firing (Hz)
dataWT = Loom_Variables{allWT, val}*60;
dataHET = Loom_Variables{allHET, val}*60;

nWT = numel(dataWT)
nHET = numel(dataHET)

mWT = nanmean(dataWT)
mHET = nanmean(dataHET)
[p, h] = ranksum(dataWT, dataHET)

% PLOT

xWT = ones(nWT, 1);
xHET = ones(nHET, 1)*2;

d = vertcat(dataWT, dataHET);
gp = vertcat(xWT, xHET);

b = boxplot(d, gp, 'Colors', [0 0 0; 0 0 0], 'Symbol','w.');
set(b, 'linew', 1.25);
ylim([-2 20])
box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
f = gcf;
f.Position = [1053  427  172  257]; 

%% Time to Max Sp

val = 20; 
dataWT = Loom_Table{allWT, val};
dataHET = Loom_Table{allHET, val};

nWT = numel(dataWT)
nHET = numel(dataHET)

mWT = nanmean(dataWT)
mHET = nanmean(dataHET)
[p, h] = ranksum(dataWT, dataHET)

% PLOT

xWT = ones(nWT, 1);
xHET = ones(nHET, 1)*2;

d = vertcat(dataWT, dataHET);
gp = vertcat(xWT, xHET);

b = boxplot(d, gp, 'Colors', [0 0 0; 0 0 0], 'Symbol','w.');
set(b, 'linew', 1.25);

box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
f = gcf;
f.Position = [1053  427  172  257]; 

%% Average of just ONE loom - for plotting Mean/SEM - avreage over the 5 looms - only for the first loom bout. 

ls = [40, 84,  128, 172, 216];

for i = 1:height(Loom_Table)
    
    d1 = Loom_Table.L1ONLY_rep1(i,:);
    d2 = Loom_Table.L1ONLY_rep2(i,:);
    
    data = [];
    for j = 1:5
        d = d1(ls(j)-3:ls(j)+46);
        dd = d2(ls(j)-3:ls(j)+46);
        if j ==5
            d = d1(ls(j)-6:ls(j)+43);
            dd = d2(ls(j)-6:ls(j)+43);
        end
     data = vertcat(data, d, dd);
    end
    
    mean_1L = nanmean(data);
    Loom_Table.REPS10_1L(i) = {data};
    Loom_Table.Av10_1L(i) = {mean_1L};
end 

% save('220801_LoomTable_OldT_PVal.mat', 'Loom_Table', 'Loom_Variables');

%% Find the most responsive cells for WT/ HET - Plot raster plots / heatmaps of firing over 10 reps

wt_resp = find(Loom_Table.Geno == 1 & Loom_Table.P_VALUE<0.0005);
het_resp = find(Loom_Table.Geno == 2 & Loom_Table.P_VALUE<0.0005);

%% RM ANOVA

i = 1; 

WT1 = table2array(Loom_Variables(allWT, i));
WT2 = table2array(Loom_Variables(allWT, i+1));
WT3 = table2array(Loom_Variables(allWT, i+2));
WT4 = table2array(Loom_Variables(allWT, i+3));
WT5 = table2array(Loom_Variables(allWT, i+4));

HET1 = table2array(Loom_Variables(allHET, i));
HET2 = table2array(Loom_Variables(allHET, i+1));
HET3 = table2array(Loom_Variables(allHET, i+2));
HET4 = table2array(Loom_Variables(allHET, i+3));
HET5 = table2array(Loom_Variables(allHET, i+4));

y = vertcat(WT1, WT2, WT3, WT4, WT5, HET1, HET2, HET3, HET4, HET5); 
gp1 = vertcat(ones(1,numel(WT1))', (ones(1,numel(WT2))*2)', (ones(1,numel(WT3))*3)', (ones(1,numel(WT4))*4)', (ones(1,numel(WT5))*5)',ones(1,numel(HET1))',(ones(1,numel(HET2))*2)', (ones(1,numel(HET3))*3)', (ones(1,numel(HET4))*4)', (ones(1,numel(HET5))*5)'); % CONTRAST
gp2 = vertcat(ones(1,numel(allWT)*5)', (ones(1,numel(allHET)*5)*2)'); % GENO


% n way ANOVA
[p, tbl, stats] = anovan(y, {gp1, gp2}, 'model', 'interaction', 'varnames', {'Loom', 'Geno'})
results = multcompare(stats, 'Dimension', [1,2], 'CType', 'bonferroni')


%% Add 'p-value' from ZETA test from FLASH stimulus to loom table.
% Loom_Table.P_VALUE_FLASH = Flash_Table.P_VALUE;

%% Save this new table
% save('220530_Setd5_Loom_Table.mat', 'Loom_Table', 'all_responsive_cells');

% b = Loom_Table.P_VALUE;
% figure; histogram(a, 'Normalization', 'pdf'); hold on; histogram(b, 'Normalization', 'pdf')

%% WT/HET cells - set variables for how to define WT/ HET

% % All cells
allWT = find(Loom_Table.Geno == 1 & Loom_Table.Depth<0 & Loom_Table.Depth>-350 & Loom_Table.TotalSpikes>10); 
allHET = find(Loom_Table.Geno == 2 & Loom_Table.Depth<0 & Loom_Table.Depth>-350 & Loom_Table.TotalSpikes>10); 

% Resp cells - all depths
allWT = find(Loom_Table.Geno == 1 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);
allHET = find(Loom_Table.Geno == 2 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);

% % Resp cells - SUPERFICIAL
% allWT = find(Loom_Table.Geno == 1 & Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.P_VALUE<0.05); % & Loom_Table.Depth<0 & Loom_Table.Depth>-2000 ); %& Loom_Table.Depth>-2000 . % & Loom_Table.Depth>-1000
% allHET = find(Loom_Table.Geno == 2 & Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.P_VALUE<0.05); % & Loom_Table.Depth<0 & Loom_Table.Depth>-2000); % & Loom_Table.Depth>-2000
% 
% % Resp cells - DEEP
% allWT = find(Loom_Table.Geno == 1 & Loom_Table.Depth<=-400 & Loom_Table.Depth>-1000 & Loom_Table.P_VALUE<0.05); % & Loom_Table.Depth<0 & Loom_Table.Depth>-2000 ); %& Loom_Table.Depth>-2000 . % & Loom_Table.Depth>-1000
% allHET = find(Loom_Table.Geno == 2 & Loom_Table.Depth<=-400 & Loom_Table.Depth>-1000 & Loom_Table.P_VALUE<0.05); % & Loom_Table.Depth<0 & Loom_Table.Depth>-2000); % & Loom_Table.Depth>-2000

% [p,h] = ranksum(dWT, dHET)

%% Plot MEAN + SEM  - Fig. S4n, o

% Col 14 = L251
% Col 16 = P-VALUE
% Col 19 = L1ONLYAV
% Col 22 = before AV
% Col 23 = AVB4
% Col 24 = Max L1
% Col 25 = T2M
% Col 26 = Delta Loom
% Col 27 = L25AV
% Col 29 = 1 Loom - average over 10 (2 x 5 - first loom bouts)
close

allWT = find(Loom_Table.Geno == 1 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<0 & Loom_Table.Depth> -1000);
allHET = find(Loom_Table.Geno == 2 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<0 & Loom_Table.Depth> -1000);


val = 25; %24; % 29 for Setd5, 24 for Ptchd1

data = cell2mat(Loom_Table{1, val});
maxval = numel(data(1,:));
   
% close all

x = 1:1:maxval;

figure

% WT
data2 = cell2mat(Loom_Table{allWT, val})*60; % L1REP norm
% Variables for Mean/SEM
n_trials = numel(allWT);
mWT = smoothdata(nanmean(data2));
semWT = nanstd(data2(:, 1:maxval))/sqrt(n_trials);
y1 = mWT+semWT;
y2 = mWT-semWT;

% Plot SEM
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mWT, 'k', 'LineWidth', 1.2);

% HET
data3 = cell2mat(Loom_Table{allHET, val})*60; % L1REP norm
% Variables for Mean/SEM
n_trials = numel(allHET);
mHET = smoothdata(nanmean(data3));
semHET = nanstd(data3(:, 1:maxval))/sqrt(n_trials);
y3 = mHET+semHET;
y4 = mHET-semHET;

% Plot SEM
plot(x, y3, 'w')
hold on
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mHET, 'r', 'LineWidth', 1.2);
box off

xlim([0 maxval])


ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.025 0.025];
ax.LineWidth = 1.5;
ax.XAxis.Visible = 'on';
% xticks([0:12:43])
% xticks([40:60:300]) %5 LOOMS
xticks([3:15:50])% 1 LOOM
ylim([0 38])
xticklabels('')

% each bin = 1/60 second.

f = gcf;
% f.Position = [391   420   504   214]; % 5 LOOMS
f.Position = [680   568   255   237]; % ONE LOOM

% [p,h] = ranksum(mWT, mHET)

% plot([40 40], [-0.1 0.5], 'k:', 'LineWidth', 1.5)
% plot([85 85], [-0.1 0.5], 'k:', 'LineWidth', 1.5)
% 
% % xticks(4:12:88)
% % xticklabels({'-0.6','-0.4','-0.2','0', '0.2', '0.4', '0.6','0.8'})
% xticks([10, 40, 85])
% xticklabels({'-0.5', '0', '0.75'})
% % title('All Cells > 10SP - -1000 to -2000')
% % ylim([-0.01 0.23])
% ylabel('Normalised Firing')
% xlabel('Time (s)')
% ax.FontSize = 18;
% ax.TickLength = [0.01 0.01];
% 
% f = gcf;
% f.Position = [ 634   372   360   307];
% 
% % xlim([10 330])
% % ylim([-0.05 0.5])
% % xlim([37 87])


%% Boxplot - variables.

val = 24;
dataWT = Loom_Table{allWT, val};
dataHET = Loom_Table{allHET, val};

nWT = numel(dataWT);
nHET = numel(dataHET);

xWT = ones(nWT, 1);
xHET = ones(nHET, 1)*2;

d = vertcat(dataWT, dataHET);
gp = vertcat(xWT, xHET);

b = boxplot(d, gp, 'Color', 'k', 'Symbol','w.');
set(b, 'linew', 1.25);

mWT = nanmean(dataWT)
mHET = nanmean(dataHET)
[p, h] = ranksum(dataWT, dataHET)
[h, p] = kstest2(dataWT, dataHET)

% Test for normality
[h, p] = kstest(dataWT)
[h, p] = kstest(dataHET)



%%  Run through Loom_Table - Find MAX and TOTAL FIRING

nrows= height(Loom_Table);

for i2 = 1:nrows
    
    % L1REP - norm
    data = Loom_Table{i2, 9};
    
    maxval = max(data);
    totalsp = sum(data);
    
    Loom_Table.MAXVAL(i2) = maxval;
    Loom_Table.TOTSP(i2) = totalsp;
    
end


% FROM RAW
for i2 = 1:nrows
    
    % L1REP - norm
    data = Loom_Table{i2, 8};
    
    maxval = max(data);
    totalsp = sum(data);
    
    Loom_Table.MAXRAW(i2) = maxval;
    Loom_Table.ALLSP(i2) = totalsp;
    
end


%% HISTOGRAMS

allWT = find(Loom_Table.Geno == 1 & Loom_Table.TotalSpikes>10);
allHET = find(Loom_Table.Geno == 2  & Loom_Table.TotalSpikes>10);

% TOTSP
data1 = Loom_Table{allWT, 24};
data2 = Loom_Table{allHET, 24};

figure; histogram(data1, -10:0.5:20, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf'); hold on
histogram(data2, -10:0.5:20, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
box off
ax = gca;
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.FontSize = 15;
xlabel('Delta Firing')
ylabel('PDF')

% MAX FIRING
data1 = Loom_Table{allWT, 16};
data2 = Loom_Table{allHET, 16};

figure; histogram(data1, -0.05:0.1:4.8, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf'); hold on
histogram(data2, -0.05:0.1:4.8, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
box off
ax = gca;
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.FontSize = 15;
xlabel('Delta Firing')
ylabel('PDF')


% MAX RAW 
data1 = Loom_Table{allWT, 18};
data2 = Loom_Table{allHET, 18};

figure; histogram(data1, -10:0.5:20, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf'); hold on
histogram(data2, -10:0.5:20, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
box off
ax = gca;
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.FontSize = 15;
xlabel('Delta Firing')
ylabel('PDF')

% MAX FIRING
data1 = Loom_Table{allWT, 19};
data2 = Loom_Table{allHET, 19};

figure; histogram(data1, 0:1:55, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf'); hold on
histogram(data2,0:1:55, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
box off
ax = gca;
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.FontSize = 15;
xlabel('Delta Firing')
ylabel('PDF')


[p,h] = ranksum(data1, data2)
[h, p] = kstest2(data1, data2)




%% ONLY THE VERY FIRST LOOM PRESENTATION - 
% Loom_Table.L1ONLY_rep1 = []; 
% Loom_Table.L1ONLY_rep2 = []; 

for i = 1:height(Loom_Table)
 Loom_Table.L1ONLY_rep1(i, :) = Loom_Table.L25_norm(i,20:350); % 20:109 = First loom! 
%  Loom_Table.L1ONLY_rep2(i, :) = Loom_Table.L25_norm2(i,20:350);
end 

















%% EXTRA CODE - NOT USED IN PAPER

% for i = 1:height(Loom_Table)
%     
%     d1 = Loom_Table.L251(i, :);
%     d2 = Loom_Table.L252(i, :); 
%     
%     dtotal = sum(d1)+sum(d2);
%     avfiring = dtotal/53;
%     
%     Loom_Table.AVFIRING(i) = avfiring;
% end 
% 
% dataa = Loom_Table.AVFIRING;
% figure; histogram(dataa)

%%
% alld = LT2.AVV;
% figure; imagesc(alld)
% 
% db = LT2.AV_RAW{237}*60;
% figure; imagesc(db); colormap(gray); caxis([0 100])
% f = gcf; 
% f.Position = [251   499   831   127];
% 
% 
% % RASTER PLOTS 
% db2 = Loom_Table.L25_ST1{9};
% 
% figure
% n_spikes = numel(db2);
% for k = 1:n_spikes
%     xval = db2(k);
%     yval = 1;
%     plot(xval, yval, 'k.', 'MarkerSize', 10);
%     hold on 
% end 


%%

close all
x = 1:1:280;

val = 7; 

figure
% WT

data2 = (LT2{allWT, val})*60; % L1REP norm
% Variables for Mean/SEM
n_trials = numel(allWT);
mWT = smooth(nanmean(data2(:, 1:280)))';
semWT = nanstd(data2(:, 1:280))/sqrt(n_trials);
y1 = mWT+semWT;
y2 = mWT-semWT;

% Plot SEM
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mWT, 'k', 'LineWidth', 1.2);

% HET

data3 = (LT2{allHET, val}*60); % L1REP norm
% Variables for Mean/SEM
n_trials = numel(allHET);
mHET = smooth(nanmean(data3(:, 1:280)))';
semHET = nanstd(data3(:, 1:280))/sqrt(n_trials);
y3 = mHET+semHET;
y4 = mHET-semHET;

% Plot SEM
plot(x, y3, 'w')
hold on
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mHET, 'r', 'LineWidth', 1.2);
box off

f = gcf;
f.Position = [257   425   825   317];

xlim([0 280])
ylim([0 10])
ax = gca;
ax.TickDir = 'out'; 
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
[p,h] = ranksum(mWT, mHET);

%% KMEANS - SORT WT and HET together - KEY PLOTS - HEATMAP AND MEAN+SEM

% data = [Loom_Table.L1REP, Loom_Table.Geno];
data = [LT2.AVV, Loom_Table.Geno];

% Only for cells between 0 and -400. 
% all_depths = find(Loom_Table.Depth<=0 & Loom_Table.P_VALUE_FLASH<0.05 & Loom_Table.Depth>-1000); % & Loom_Table.Depth>-1500
all_depths = find(Loom_Table.Ani==7270 & Loom_Table.P_VALUE<0.05 | Loom_Table.Ani == 7269 & Loom_Table.P_VALUE<0.05); % & Loom_Table.Depth>-1500
% all_depths = find(Flash_Table.Depth>0 | Flash_Table.P_VALUE>=0.05);

data = data(all_depths, :);
data(:, 1:281) = data(:, 1:281)*60;
% data(all_depths, :) = []; 

% % % % % Evaluate the optimal number of clusters:
eva = evalclusters(data, 'kmeans','CalinskiHarabasz', 'KList', [1:50]); % 'silhouette',
figure; plot(eva)
 
%% Set the elbow value as the number of clusters. 
n_k = 6;

data(:, 283) = kmeans(data(:, 1:281), n_k);
data = sortrows(data, 283);

allWT = find(data(:, 282)==1);
allHET = find(data(:, 282)==2);

nWT = numel(allWT);
nHET = numel(allHET);

%% 

% % % % % % % PLOT HEATMAPS

spn = 15; 

figure
ax = subplot(1,spn,1:spn-1);
imagesc(data(allWT, 1:281)); caxis([0 100])
hold on
% plot([30 30], [0 nWT], 'w:', 'LineWidth', 1.2)
% plot([90 90], [0 nWT], 'w:', 'LineWidth', 1.2)
box off
ax.XTick = []; 
colormap(ax(1), redblue)
% caxis(ax, [0 60])
ax.YTick = [];

ax2 = subplot(1,spn,spn);
imagesc(data(allWT, 283))
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


%
figure
ax = subplot(1,spn,1:spn-1);
imagesc(data(allHET,1:281)); caxis([0 100])
hold on
% plot([30 30], [0 nHET], 'w:', 'LineWidth', 1.2)
% plot([90 90], [0 nHET], 'w:', 'LineWidth', 1.2)
box off
ax.XTick = []; 
colormap(ax(1), redblue)
% caxis(ax, [0 60])
ax.YTick = [];

ax2 = subplot(1,spn,spn);
imagesc(data(allHET, 283))
colormap(ax2, gray)
f = gcf;
f.Position = [699 1 356   804];
box on
ax2.LineWidth = 0.5;
ax2.Color = 'k';
ax2.TickLength = [0 0];
ax2.XTick = [];
ax2.YTick = []; 
    
%% % % % % % % Plots of the Mean/SEM of each cluster:

x = 1:1:281;

% WT % % % % %
figure
for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 283)== j & data(:,282) == 1);
    if numel(all_type)>=1
        % Separate these rows into a different array
        data2 = data(all_type, 1:281);
        
        % Variables for Mean/SEM
        
        n_trials = numel(all_type);
        if numel(all_type)>1
            mWT = nanmean(data2);
        else
            mWT = (data2);
        end
        
        semWT = nanstd(data2)/sqrt(n_trials);
        y1 = mWT+semWT;
        y2 = mWT-semWT;
        
        % Make the plot
        subplot(n_k,1,j); hold on
        
        % Find Min/Max
            maxx = max(nanmean(data2))+0.075;
            minn = min(nanmean(data2))-0.075;
        
%         maxx = 3;
%         minn = 0;
        
        % Rectangle of light stim
%         rectangle('Position', [0 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
%         %     rectangle('Position', [30 minn 60 diff([minn maxx])], 'FaceColor', [1 1 0 0.1], 'EdgeColor', 'none')
%         rectangle('Position', [90 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
%         
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
        ax.LineWidth = 1;
%         xlim([ 0 42])
    end
end
f2 = gcf;
f2.Position = [365  1   256   804];


% HET % % % % %
hold on 
for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 283)== j & data(:,282) == 2);
    
    if numel(all_type)>=1
    % Separate these rows into a different array
    data2 = data(all_type, 1:281);
    
    % Variables for Mean/SEM
    n_trials = numel(all_type);
    
     if numel(all_type)>1
        mHET = nanmean(data2);
     else
        mHET = (data2);
     end
    
%     mHET = nanmean(data2);
    semHET = nanstd(data2)/sqrt(n_trials);
    y1 = mHET+semHET;
    y2 = mHET-semHET;
    
    % Make the plot
    subplot(n_k,1,j); hold on
    
    % Find Min/Max
    maxx = max(nanmean(data2))+0.075;
    minn = min(nanmean(data2))-0.075;
%      maxx = 3;
%     minn = 0;
    
    % Rectangle of light stim
%     rectangle('Position', [0 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
% %     rectangle('Position', [30 minn 60 diff([minn maxx])], 'FaceColor', [1 1 0 0.1], 'EdgeColor', 'none')
%     rectangle('Position', [90 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
%     
    % Plot SEM
    plot(x, y1, 'w')
    plot(x, y2, 'w')
    patch([x fliplr(x)], [y1 fliplr(y2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    % Plot mean
    plot(mHET, 'r', 'LineWidth', 1.2);
    box off
    ylim([minn maxx])
    ax = gca;
    ax.TickDir = 'out';
    ax.XAxis.Visible = 'off';
    ax.LineWidth = 1;
%     xlim([0 42])
    end 
    
end
f3 = gcf;
f3.Position = [1056  1 256 804];


%% PIE CHART

expl = ones(1,8);

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


%% Loom responsive cells


all_ani = unique(Loom_Table.Ani);

WTVALS = [];
HETVALS = []; 

het_animals = [7269, 7476, 7614];
        
for j = 1:numel(all_ani)
    ani= all_ani(j);
    
    ani_all =find(Loom_Table.Ani == ani & Loom_Table.Depth<=0 & Loom_Table.Depth>-1000); % Flash_Table.Depth<-500 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500
    n_all = numel(ani_all);
    ani_rows = find(Loom_Table.Ani == ani & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<=0 & Loom_Table.Depth>-1000); % Flash_Table.Depth<-500 & Flash_Table.Depth>-1000
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


% Look at distribution of ALL p-values
% allwt = find(Loom_Table.Geno == 1 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);
% allhet = find(Loom_Table.Geno ==2 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);
% 
% w1 = Loom_Table.P_VALUE(allwt);
% h1 = Loom_Table.P_VALUE(allhet);
% 
% [p,h] = kstest2(w1, h1)


%% % of Loom REsponsive cells - pooled across all animals
% Number of units in each condition: 

    all_wt = numel(find(Loom_Table.Geno == 1 & Loom_Table.Depth<=-400 & Loom_Table.Depth>-1000 & Loom_Table.TotalSpikes<3500)) 
    all_het = numel(find(Loom_Table.Geno == 2 & Loom_Table.Depth<=-400 & Loom_Table.Depth>-1000 & Loom_Table.TotalSpikes<3500)) 

% Responsive cells: 

    allwt_resp = numel(find(Loom_Table.Geno == 1 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<=-400 & Loom_Table.Depth>-1000 & Loom_Table.TotalSpikes<3500))
    allhet_resp = numel(find(Loom_Table.Geno == 2 & Loom_Table.P_VALUE<0.05 & Loom_Table.Depth<=-400 & Loom_Table.Depth>-1000 & Loom_Table.TotalSpikes<3500))
    
    wt_ratio = allwt_resp/all_wt
    het_ratio = allhet_resp/all_het
    
    figure
    bar(1, (wt_ratio), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on
    bar(2, (het_ratio),'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    axis([0 3 0 1])
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.TickDir = 'out';
    ax.TickLength = [0.02, 0.02];
    ax.LineWidth= 1.5;
    f = gcf;
    f.Position = [680   750   230   348];
    
    

%% BAR CHART + ANIMAL POINTS - % OF LOOM RESPONSIVE CELLS
% Plus errorbar

close
semwt = nanstd(WTVALS)/sqrt(numel(WTVALS)); 
semhet = nanstd(HETVALS)/sqrt(numel(HETVALS)); 

figure
bar(1, nanmean(WTVALS), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
bar(2, nanmean(HETVALS),'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

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
   
    marker = 'r.';
    scatter(xval, yval, 1000,  marker, 'jitter', 'on');
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
errorbar(2, nanmean(HETVALS), semhet, 'r', 'LineWidth', 1.5);

f = gcf;
f.Position = [680   750   230   348]; 

%% 
%% Plot adaptation - 25 Looms

allWT = find(Loom_Table.Geno == 1 & Loom_Table.TotalSpikes>10 & Loom_Table.Depth<0); %& Loom_Table.Depth>-2000
allHET = find(Loom_Table.Geno == 2  & Loom_Table.TotalSpikes>10 & Loom_Table.Depth< 0); % & Loom_Table.Depth>-2000

% allWT = find(Loom_Table.Geno == 1 & Loom_Table.TotalSpikes>10 );
% allHET = find(Loom_Table.Geno == 2  & Loom_Table.TotalSpikes>10);

x = 1:1:1590;

figure

% WT

data2 = Loom_Table{allWT, 15}; % L1REP norm
% Variables for Mean/SEM
n_trials = numel(allWT);
mWT = nanmean(data2);
semWT = nanstd(data2)/sqrt(n_trials);
y1 = mWT+semWT;
y2 = mWT-semWT;

% Plot SEM
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mWT, 'k', 'LineWidth', 1.2);

% HET

data3 = Loom_Table{allHET, 15}; % L1REP norm
% Variables for Mean/SEM
n_trials = numel(allHET);
mHET = nanmean(data3);
semHET = nanstd(data3)/sqrt(n_trials);
y3 = mHET+semHET;
y4 = mHET-semHET;

% Plot SEM
plot(x, y3, 'w')
hold on
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mHET, 'r', 'LineWidth', 1.2);


box off
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1;

ylabel('Normalised Firing')
xlabel('Time (s)')
ax.FontSize = 18;
ax.TickLength = [0.02 0.02];

