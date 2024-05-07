%% PATCH DATA ANALYSIS
% Author: Burnett, 2022. 

% Analyse Data from Peter Koppensteiner from the patch clamp experiments. 
% Folders of each animal with .abf file. 
% Each .abf file contains the data from one single 'run' 
% Need to assess which .abf files correspond to which type of experiment
% type - e.g. 10pA rectangular steps, mini EPSC..

% fn = '2021_11_22_0008.abf';

%% Make cell array 'd_all' with the firing data from all experiments in a new row. 
% For experiments with multiple sweeps e.g. 10pA steps - each step is
% stored in a different column. 
clear
close all

files = dir('*.abf');
nfiles = length(files);

d_all = []; 

for i = 1:nfiles
fn = files(i).name;
[d,si,h]=abfload(fn);
d = squeeze(d);
n = numel(d(1,:));

%% Make array of data from current injection experiments. - add real data to 'data_info'
% Make cell array 'd_all'
% Each row is a cell - for info see 'data_info' 
% Each column is a current injection step. - starts from col3. 
[nn, n_current_steps]  = size(d);

for jj = 1:n_current_steps
d_all{i,jj} = d(:,jj);
end 
 
end 


%% 
% figure;
% for k = 1:13
%     plot(d(:, k))
% end 

%% FIRING versus CURRENT STEPS

% Which rows correspond to the 10pA Rect Step experiments? 
% C3 - M1: [2,3,4,5,8,10,12];
% C3 - M2: 

%Find rows with 'sweeps'
cell_rows = []; 

for kk = 1:numel(d_all(:,1))
if ~isempty(d_all{kk,15})
    cell_rows = [cell_rows, kk];
end
end
% cell_rows = [4,6,8,9,10,12];

n_cells = numel(cell_rows);

% % Number of steps
nsteps = size(d_all, 2);

% create array where each column is a different cell and each row is a
% 10pA current step. Contains values of number of spikes for each current
% step. Find the number of times the mV > 0mV

spikedata = NaN(nsteps, n_cells); 

for i2 = 1:n_cells

% Which row to examine. 
row = cell_rows(i2); 

% 1 - Plot the data to observe roughly what is going on: 
% figure
% for ji = 1:nsteps 
%     d = d_all(row, :);
%     plot(d{:, ji}')
%     hold on 
% end 

% 2 - Calculate the number of spikes within each step. 
lim1 = 8000;
lim2 = 62500; 

for ji = 1:nsteps 
d = d_all(row, :);
spkdata = d{:, ji};

if ~isempty(spkdata)
    spkdata = spkdata(lim1:lim2-1);% clip data to step time.
    signdata = sign(spkdata);
    sdata2 = find(signdata>0);
    nsd = numel(sdata2);

    if nsd>1
        for k = 2:nsd
            sdata2(k,2) = sdata2(k,1) -sdata2(k-1,1);
        end
        n_spikes =  numel(find(sdata2(:,2)>1))+1;
    else
        n_spikes = 0 ;
    end
    spikedata(ji, i2) = n_spikes;
else
    spikedata(ji, i2) = NaN;
end

end 

end 

% save('/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/C3/CurrentInjection_Data/Spike_Data_Setd5_C3_M5_SC.mat', 'spikedata');

save_folder = '/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/Ptchd1/Analysis/10pA_Steps_Spike_Data';
% save_folder = '/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/Cul3/10pA_CurrentInj_Data'; 
fname = 'Spike_Data_Ptchd1_M11_DTX_PAG.mat'; 

% save(fname, 'spikedata', 'd_all');
save(fullfile(save_folder, fname),  'spikedata', 'd_all');

%% Plot all current injections on same figure

close
cell_i = 3; 

% figure
offset = 0; 
for i = 1:numel(d_all(1,:))
curr_i = i; 
dd = cell2mat(d_all(cell_i, curr_i));
plot(dd(1:75000)-offset, 'k')
hold on 
curr_i = curr_i +1; 
offset = offset+120;
end 
box off
ax = gca; 
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.TickDir  = 'out';
hold off
f = gcf;
f.Position = [440   516   232   282];


%% Combine spiking data from the different animals to make a table.

% folder with data per animal: /Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/Ptchd1/Analysis/10pA_Steps_Spike_Data

SPIKE_DATA = [];

files = dir('*.mat');
het_animals = [2,3,6,8,9,10,11];

for i = 1 %numel(files)
    fname = files(i).name;
    load(fname)
    
    n_cells = numel(spikedata(1,:));
    data = spikedata(1:24, :)'; 
    
%     ani_num = str2num(fname(end-8));
    
    anii = repmat(ani_num, n_cells, 1);
    if ismember(ani_num, het_animals)
        genoo = ones(n_cells, 1)*2;
    else
        genoo = ones(n_cells, 1);
    end
    celll = [1:1:n_cells]';
        
    dataa = [anii, celll, genoo, data];

    SPIKE_DATA = vertcat(SPIKE_DATA, dataa);
    
end 

SPIKE_DATA = array2table(SPIKE_DATA, 'VariableNames', {'Ani', 'Cell', 'Geno', 'ST1', 'ST2', 'ST3', 'ST4', 'ST5', 'ST6', 'ST7', 'ST8', 'ST9', 'ST10', 'ST11', 'ST12', 'ST13', 'ST14', 'ST15', 'ST16', 'ST17', 'ST18', 'ST19', 'ST20', 'ST21', 'ST22', 'ST23', 'ST24'});

save('221107_SPIKE_DATA_N11_Ptchd1_10pA_Steps_wDTX.mat', 'SPIKE_DATA')



%% Errorbar

% 34 current inj steps. 
load('/Users/lauraburnett/Data_Analysis_Mac/PATCH/curr_inj_vals.mat', 'curr_inj')

inj_steps = curr_inj(:, 2);
n_steps = numel(inj_steps);

data_per_day = zeros(2,n_steps); 
sem_per_day = zeros(2,n_steps); 
n_per_day  =zeros(2,n_steps); 
% 
allWT = find((SPIKE_DATA.Geno) ==1);
allHET = find((SPIKE_DATA.Geno) ==2);

for j = 4:26 %1:n_steps
    data_WT = SPIKE_DATA{allWT, j}; 
    data_HET = SPIKE_DATA{allHET, j}; 
    data_per_day(1,j-3) = nanmean(data_WT); 
    data_per_day(2,j-3) = nanmean(data_HET); 
    n_per_day(1,j-3) = numel(data_WT); 
    n_per_day(2,j-3) = numel(data_HET); 
    sem_per_day(1,j-3) = nanstd(data_WT)/sqrt(numel(data_WT));
    sem_per_day(2,j-3) = nanstd(data_HET)/sqrt(numel(data_HET));
end 


figure
errorbar(inj_steps(1:23), data_per_day(1,1:23), sem_per_day(1,1:23),'o', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'k') %28 
hold on 
errorbar(inj_steps(1:23), data_per_day(2,1:23), sem_per_day(2,1:23), 'o', 'Color', col,'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'w')
% xticks(1:1:16)
% xticklabels(string(inj_steps(1:16)))
ylabel('AP Frequency (Hz)')
box off
set(gca, 'FontSize', 18)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.5;
xlim([-10 150])
ylim([-5 80])

%% MEAN - SEM - PATCH 

allWT = find(SPIKE_DATA.Geno == 1 );
allHET = find(SPIKE_DATA.Geno == 2 );

valsWT2 = table2array(SPIKE_DATA(allWT, 4:27))';
valsHET2 = table2array(SPIKE_DATA(allHET, 4:27))';

% col = 'm';
v = 21;
inj_steps = curr_inj(3:3+v, 2);

nWT = numel(valsWT2(1,:)); 
nHET = numel(valsHET2(1,:));

mean_WT = nanmean(valsWT2'); 
mean_HET = nanmean(valsHET2'); 

x = (1:1:v);

semWT = nanstd(valsWT2')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(valsHET2')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col, 'LineWidth', 1.3)

% xticks(1:1:15)
% xticklabels(string(inj_steps(1:15)))

xticks([1:4:24])
xticklabels(string(inj_steps(1:4:24)))

xlabel('Current (pA)')
ylabel('Firing Rate (Hz)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.02 0.02];
xlim([0.5 24])
ylim([0 100])

f = gcf;
f.Position = [742   492   304   258]; 

%% MEAN SEM - WITH DENDROTOXIN 

allWT = find(SPIKE_DATA.Geno == 1 & SPIKE_DATA.DTX ==0);
allHET = find(SPIKE_DATA.Geno == 2 & SPIKE_DATA.DTX ==0);
allHETDTX = find(SPIKE_DATA.Geno == 2 & SPIKE_DATA.DTX ==1);

valsWT2 = table2array(SPIKE_DATA(allWT, 4:27))';
valsHET2 = table2array(SPIKE_DATA(allHET, 4:27))';
valsHETDTX = table2array(SPIKE_DATA(allHETDTX, 4:27))';

% col = 'm';
v = 24;
inj_steps = curr_inj(3:3+v, 2);

nWT = numel(valsWT2(1,:)); 
nHET = numel(valsHET2(1,:));
nHETD = numel(valsHETDTX(1,:));

mean_WT = nanmean(valsWT2'); 
mean_HET = nanmean(valsHET2'); 
mean_HETD = nanmean(valsHETDTX'); 

x = (1:1:v);

semWT = nanstd(valsWT2')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(valsHET2')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

semHETD = nanstd(valsHETDTX')/sqrt(nHETD); 
y5 = mean_HETD+semHETD;
y6 = mean_HETD-semHETD;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col, 'LineWidth', 1.3)

plot(x, y5(1:v), 'w')
hold on 
plot(x, y6(1:v), 'w')
patch([x fliplr(x)], [y5(1:v) fliplr(y6(1:v))],  [0.5 0.15 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HETD(1:v)', 'Color', [0.5 0.15 1], 'LineWidth', 1.3)

% xticks(1:1:15)
% xticklabels(string(inj_steps(1:15)))

xticks([1:4:24])
xticklabels(string(inj_steps(1:4:24)))

% xlabel('Current (pA)')
% ylabel('Firing Rate (Hz)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.02 0.02];
xlim([0.5 24])
ylim([0 100])

f = gcf;
f.Position = [742   492   304   258]; 


%% STATS for with DTX





%% Each ani

inj_steps = curr_inj(:, 2);
n_steps = numel(inj_steps);

data_per_day = zeros(2,n_steps); 
sem_per_day = zeros(2,n_steps); 
n_per_day  =zeros(2,n_steps); 
% 
allWT = find((SPIKE_DATA.Ani) ==1);

for j = 4:26 %1:n_steps
    data_WT = SPIKE_DATA{allWT, j}; 
    data_per_day(1,j-3) = nanmean(data_WT); 
    n_per_day(1,j-3) = numel(data_WT); 
    sem_per_day(1,j-3) = nanstd(data_WT)/sqrt(numel(data_WT));
end 


% figure
errorbar(inj_steps(1:23), data_per_day(1,1:23), sem_per_day(1,1:23),'o-', 'Color', 'r', 'MarkerSize', 5, 'LineWidth', 1.4, 'MarkerFaceColor', 'w') %28 
hold on 


%% Each CELL

figure
for i = 1:50

        g = SPIKE_DATA.Geno(i);
        d = SPIKE_DATA.DTX(i);
        
        if g == 1 && d == 0 
            col = 'k';
        elseif g == 2 && d == 0 
            col = [1 0.55 0.2];
        elseif g == 2 && d == 1
            col = [0.5 0.15 1]; %'m';
        end 
        
        x = inj_steps(1:20);
        data = SPIKE_DATA{i, 4:23};
        
        plot(x, data, '-', 'Color', col, 'LineWidth', 1)
        hold on
end 


%% Plot trace of spiking to 50pA injection from d_all
close all

i = 17;
data = d_all{i,4};
figure; plot(data(1:75000, :), 'Color', 'k'); %[1 0.5 0.15]); % [0.5 0.15 1]
xlim([0 75000])
ylim([-80 60])
box off
ax = gca;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
f = gcf;
f.Position = [355   439   514   201];

f = gcf;
f.Position = [219   456   560   191];

% Cul3 HETs = 2,3,6

%%
saveas(gcf, '/Users/lauraburnett/Data_Analysis_Mac/PATCH/FIGS/Cul3/50pA_Traces/M1C1_50pA_WT.svg');
% saveas(gcf, '/Users/lauraburnett/Data_Analysis_Mac/PATCH/FIGS/Ptchd1/50pA_Traces/M6C3_50pA_WT.svg');

close all


%%
% %% 100pA 
% 
% % i = 1;
% cell_i = i; 
% curr_i = 6; % 6 = 50pA
% % 
% % if i == 4 || i ==5 || i ==6 || i == 14 || i == 15 || i ==16 || i ==17 || i == 18 || i == 19
% %     col = 'k';
% % else 
% %     col = 'r';
% % end 
% 
% % figure
% dd = cell2mat(d_all(cell_i, curr_i));
% plot(dd(1:75000), 'Color', col)
% box off
% ax = gca; 
% ax.XAxis.Visible = 'off';
% ax.TickDir  = 'out';
% ax.LineWidth = 1.5;
% ax.FontSize = 20;
% hold off
% i = i+1; 
% 
% %% Make an array where 
% % col 1 = the MAX number of spikes elicited by the cell
% % Col 2 = COL (i.e. curr inj) where APs are highest. 
% % col 3 =  current injection that led to this value. 
% % from 'data'
% % COl 4 = geno. 
% 
% max_data = zeros(21,4);
% 
% for jj = 1:21
%     max_data(jj,1) = max(data(jj, :));
%     inj_val = find(data(jj, :) == max(data(jj,:)));
%     if numel(inj_val)>1
%         inj_val = inj_val(1);
%     end 
%     max_data(jj, 2) = inj_val;
%     inj_val2 = inj_steps(inj_val);
%     max_data(jj, 3) = inj_val2;
%     max_data(jj,4) = data_info.geno(jj);
% end 


%% Plots of firing at 3 different current injections. 

% col = 'r'; 
% 
% subplot(3,1,1)
% plot(d(1:70000, 4), 'Color', col)
% hold on
% box off
% xticks([])
% % yticks([])
% ylim([-70 70])
% hold off
% ax1 = gca;
% ax1.XAxis.Visible = 'off';
% xlabel('10pA')
%     
% 
% subplot(3,1,2)
% if n >=10
% plot(d(1:70000, 10), 'Color', col)
% else 
%     disp(n)
% end 
% hold on
% box off
% xticks([])
% % yticks([])
% xlabel('70pA')
% ylim([-70 70])
% hold off
% ax2 = gca;
% ax2.XAxis.Visible = 'off';
% 
% subplot(3,1,3)
% if n >=18
% plot(d(1:70000, 18), 'Color', col)
% end 
% hold on
% box off
% xticks([])
% % yticks([])
% xlabel('150pA')
% ylim([-70 70])
%     hold off
%     ax3 = gca;
% ax3.XAxis.Visible = 'off';
%     
% sgtitle('M8-C2-59')
% hold off


%%


% val = 0; 
% for i = 1:n
%     plot(d(1:70000, i)-val, 'Color', [0.3 0.3 0.3])
% %     subplot(18,1,i)
%     hold on 
%     box off
%     val = val +120; 
%     xticks([])
%     yticks([])
% end 
% title('M8-C2-59')
% hold off


%%

% figure
% plot(d(1:70000, 1), 'Color', 'k')
% ax = gca;
% ax.XAxis.Visible = 'off'

%%

% APanalysis(fn)
% [N ,out1] = spike_times2(d(:,10),10)



%%

inj_steps = data(:,1);
data = data(:, [2:end]);
data2 = data2';
data = data';

animal = data2(:,1);
cell = data2(:,2);
filen = data2(:,3);
geno = data2(:,4);

data_info = table(animal, cell, filen, geno);

figure
for i = 1:21
    if data_info.geno(i) == 1
        col = 'k';
    else 
        col = 'r'; 
    end 
plot(inj_steps, data(i, :), '-','Color', col)
hold on 
input('')
end

%%

allWT = find(data_info.geno == 1);
allHET = find(data_info.geno == 2);

figure
subplot(1,2,1)
plot(inj_steps, data(allWT, :), 'k-', 'LineWidth', 1.5)
axis([0 150 0 100])
box off
ax1 = gca;
ax1.TickDir  = 'out';
ax1.LineWidth = 1.75;
ax1.FontSize = 18;
hold on
plot([0 150], [0 80], 'k:')

subplot(1,2,2)
plot(inj_steps, data(allHET, :), 'r-', 'LineWidth', 1.5)
axis([0 150 0 100])
box off
ax2 = gca;
ax2.TickDir  = 'out';
ax2.LineWidth = 1.75;
ax2.FontSize = 18;
hold on
plot([0 150], [0 80], 'k:')

%% Errorbar

% 34 current inj steps. 

n_steps = numel(inj_steps);

data_per_day = zeros(2,n_steps); 
sem_per_day = zeros(2,n_steps); 
n_per_day  =zeros(2,n_steps); 
% 
allWT = find((data_info.geno) ==1);
allHET = find((data_info.geno) ==2);

for j = 1:n_steps
    data_WT = data(allWT, j); 
    data_HET = data(allHET, j); 
    data_per_day(1,j) = nanmean(data_WT); 
    data_per_day(2,j) = nanmean(data_HET); 
    n_per_day(1,j) = numel(data_WT); 
    n_per_day(2,j) = numel(data_HET); 
    sem_per_day(1,j) = nanstd(data_WT)/sqrt(numel(data_WT));
    sem_per_day(2,j) = nanstd(data_HET)/sqrt(numel(data_HET));
end 


figure
errorbar(inj_steps(1:16), data_per_day(1,1:16), sem_per_day(1,1:16),'o', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'k') %28 
hold on 
errorbar(inj_steps(1:16), data_per_day(2,1:16), sem_per_day(2,1:16), 'o', 'Color', col,'MarkerSize', 10, 'LineWidth', 1.4, 'MarkerFaceColor', 'w')
% xticks(1:1:16)
% xticklabels(string(inj_steps(1:16)))
ylabel('AP Frequency (Hz)')
box off
set(gca, 'FontSize', 18)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.5;
xlim([-10 150])
ylim([-5 80])

% xtickangle = 45;


%% Plot histogram of time to first spike at 3 conditions. 150/75. 



%% Mean SEM - WT/ HET - Red/Black

% v = 23;
%  
% speed_WT = [];
% speed_HET = []; 
% 
% for i = 1:21
%     if data_info.geno(i) == 1   
%         G = data(i,:);
%         speed_WT = vertcat(speed_WT, G); 
%     elseif data_info.geno(i) == 2   
%         F = data(i, :);
%         speed_HET = vertcat(speed_HET, F);
%     end
% end

col = 'r';
v = 20;
inj_steps = curr_inj(3:3+v, 2);

nWT = numel(valsWT2(1,:)); 
nHET = numel(valsHET2(1,:));

mean_WT = nanmean(valsWT2'); 
mean_HET = nanmean(valsHET2'); 

x = (1:1:v);

semWT = nanstd(valsWT2')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(valsHET2')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col, 'LineWidth', 1.3)

% xticks(1:1:15)
% xticklabels(string(inj_steps(1:15)))

xticks([1:4:15])
xticklabels(string(inj_steps(1:4:15)))

xlabel('Current (pA)')
ylabel('Firing Rate (Hz)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.02 0.02];
xlim([0.5 15])
ylim([0 80])

f = gcf;
f.Position = [742   492   304   258]; 

% 
% save('211014_CurrentInj.mat', 'data', 'data_info', 'inj_steps', 'max_data', 'd_all');

%% Mean SEM - WT/ HET - Dendrotoxin - Blue/Magenta

col1 = 'b';
col2 = [0.7 0 1]; %'m';

v = 21;
inj_steps = curr_inj(3:3+v, 2);

nWT = numel(DTX_WT(1,:)); 
nHET = numel(DTX_HET(1,:));

mean_WT = nanmean(DTX_WT'); 
mean_HET = nanmean(DTX_HET'); 

x = (1:1:v);

semWT = nanstd(DTX_WT')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(DTX_HET')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

% figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], col1, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', col1, 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col2, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col2, 'LineWidth', 1.3)

% xticks(1:1:15)
% xticklabels(string(inj_steps(1:15)))

xticks([1:4:15])
xticklabels(string(inj_steps(1:4:15)))

xlabel('Current (pA)')
ylabel('Firing Rate (Hz)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.02 0.02];
xlim([0.5 15])
ylim([0 90])

f2 = gcf;
f2.Position = [742   492   304   258]; 




%% SC cells 

v = 21;
inj_steps = curr_inj(3:3+v, 2);

WTVALS = []; 
HETVALS = []; 

for i = 1:15
    valsW = find(SC_WT(:,i)>0);
    val = valsW(1);
    inj_valW = inj_steps(val);
    WTVALS = [WTVALS, inj_valW];

    valsH = find(SC_HET(:,i)>0);
    val = valsH(1);
    inj_valH = inj_steps(val);
    HETVALS = [HETVALS, inj_valH];
end 

[p,h] = ranksum(WTVALS, HETVALS)
nanmean(WTVALS)
nanmean(HETVALS)

%% PAG cells - rheobase

WTVALS = []; 
HETVALS = []; 

for i = 1:18
    valsW = find(valsWT2(:,i)>0);
    val = valsW(1);
    inj_valW = inj_steps(val);
    WTVALS = [WTVALS, inj_valW];
end 

for i = 1:24
    valsH = find(valsHET2(:,i)>0);
    val = valsH(1);
    inj_valH = inj_steps(val);
    HETVALS = [HETVALS, inj_valH];
end 

[p,h] = ranksum(WTVALS, HETVALS)
[h,p] = kstest2(WTVALS, HETVALS)
nanmean(WTVALS)
nanmean(HETVALS)



%% Current needed to first spike. 
WTVALS = []; 
HETVALS = []; 

for i = 1:21
    vals = find(data(i, :)>0);
    val = vals(1);
    inj_val = inj_steps(val);
    geno = data_info.geno(i);
    if geno ==1 
        WTVALS = [WTVALS, inj_val];
    elseif geno == 2
        HETVALS = [HETVALS, inj_val];
    end 
end 

%% Current to >1 Spike
WTVALS = []; 
HETVALS = []; 

for i = 1:21
    vals = find(data(i, :)>1);
    val = vals(1);
    inj_val = inj_steps(val)
    geno = data_info.geno(i);
    if geno ==1 
        WTVALS = [WTVALS, inj_val];
    elseif geno == 2
        HETVALS = [HETVALS, inj_val];
    end 
end 

wt = WTVALS;
het = HETVALS;

curr2spike = table(wt, het);

% Violin Plot. 
violinplot(curr2spike,[], 'ViolinColor', [0.7 0.7 0.7])
xlabel('Genotype')
ylabel('Current Inj to first spike')
ax = gca; ax.FontSize = 18;

% Box plot + DOTS
boxplot([WTVALS', HETVALS'], [ones(1,19)', ones(1,24)'*2])

mWT = nanmean(WTVALS);
mHET = nanmean(HETVALS);
semWT = nanstd(WTVALS)/ sqrt(numel(WTVALS));
semHET = nanstd(HETVALS)/ sqrt(numel(HETVALS));

x1 = ones(1, 19);
x2 = ones(1, 24)*2;

figure
scatter(x1, WTVALS,'SizeData', 200, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, HETVALS, 'SizeData', 200,'MarkerEdgeColor', [1 0.4 0.4], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
errorbar(1, mWT, semWT, 'o','Color', [0.2 0.2 0.2], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(2, mHET, semHET, 'o', 'Color', [1 0.2 0.2],'Marker', 'none', 'LineWidth', 1.75)
plot([0.75 1.25], [mWT mWT], 'k', 'LineWidth', 2)
plot([1.75 2.25], [mHET mHET], 'r', 'LineWidth', 2)
ylim([-10 100])
xlim([0.5 2.5])
ax= gca;
ax.TickDir = 'out';
box off
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickLength = [0.02 0.02];

[p,h] = ranksum(WTVALS, HETVALS)


%% Current needed to MAX SPIKES
WTVALS = []; 
HETVALS = []; 

for i = 1:21
    vals = find(data(i, :) == max(data(i, :)));
    val = vals(1);
    inj_val = inj_steps(val);
    geno = data_info.geno(i);
    if geno ==1 
        WTVALS = [WTVALS, inj_val];
    elseif geno == 2
        HETVALS = [HETVALS, inj_val];
    end 
end 

wt = WTVALS;
het = HETVALS;

curr2max = table(wt, het);

violinplot(curr2max,[], 'ViolinColor', [0.7 0.7 0.7])
xlabel('Genotype')
ylabel('Current Inj to max spiking')
ax = gca; ax.FontSize = 18;



%%

for i = 1:21
    vals = find(data(i, :)>0);
    val = vals(1);
    inj_val = inj_steps(val);
    data_info.curr2spike(i) = inj_val;
    
    vals2 = find(data(i, :) == max(data(i, :)));
    val2 = vals2(1);
    inj_val2 = inj_steps(val2);
    data_info.curr2max(i) = inj_val2;
    
    a1 = inj_val2 - inj_val;
    data_info.range(i) = a1;
end 

figure
for i = 1:21
    if data_info.geno(i) == 1
        col = 'k';
    else 
        col = 'r';
    end 
    plot(data_info.curr2spike(i), data_info.curr2max(i), 'o', 'Color', col, 'MarkerSize', 10);
    hold on 
end 
box off


[p,h] = ranksum(WTVALS, HETVALS)




%% SOMA SIZE


WPAG = reshape(allWT_PAG, [24,1]);
HPAG = reshape(allHET_PAG, [60,1]);

nanmean(WPAG)
nanmean(HPAG)
[p,h] = ranksum(WPAG, HPAG)
[h, p] = kstest2(WPAG, HPAG)

WSC = reshape(allWT_SC, [9,1]);
HSC = reshape(allHET_SC, [9,1]);

nanmean(WSC)
nanmean(HSC)
[p,h] = ranksum(WSC, HSC)
[h, p] = kstest2(WSC, HSC)

% all cells individuallly
subplot(1,2,1)
boxplot(allWT_PAG)
ylim([0 250])
box off
xlabel('Mouse')
ylabel('Soma Size (um^2)')
title('WT- SC')
subplot(1,2,2)
boxplot(allHET_PAG)
ylim([0 250])
box off
xlabel('Mouse')
ylabel('Soma Size (um^2)')
title('HET - SC')

% mean across cells
subplot(1,2,1)
boxplot(WPAG)
ylim([0 250])
box off
xlabel('WT')
ylabel('Soma Size (um^2)')
subplot(1,2,2)
boxplot(HPAG)
ylim([0 250])
box off
xlabel('HET')
ylabel('Soma Size (um^2)')

save('211201_PATCHDATA_SOMASIZE.mat', 'allWT_PAG', 'allHET_PAG', 'allHET_SC', 'allWT_SC', 'WSC', 'WPAG', 'HSC', 'HPAG');


%% SOMA SIZE _ PLOT. 
allNAN = find(isnan(WPAG));
WPAG(allNAN) = []; 

allNAN = find(isnan(WPAG));
WPAG(allNAN) = []; 


% Box plot + DOTS
% boxplot([WPAG', HPAG'], [ones(1,19), ones(1,32)*2])

% Errorbar + points
nWT = numel(WPAG);
nHET = numel(HPAG);

mWT = nanmean(WPAG);
mHET = nanmean(HPAG);

semWT = nanstd(WPAG)/sqrt(nWT);
semHET = nanstd(HPAG)/sqrt(nHET);

x1 = ones(1, nWT);
x2 = ones(1, nHET)*2;

figure
scatter(x1, WPAG,'SizeData', 100, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, HPAG, 'SizeData', 100,'MarkerEdgeColor', [1 0.4 0.4], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
errorbar(1, mWT, semWT, 'o','Color', [0.2 0.2 0.2], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(2, mHET, semHET, 'o', 'Color', [1 0.2 0.2],'Marker', 'none', 'LineWidth', 1.75)
plot([0.75 1.25], [mWT mWT], 'k', 'LineWidth', 2)
plot([1.75 2.25], [mHET mHET], 'r', 'LineWidth', 2)

ylim([0 400])
xlim([0.5 2.5])
ax= gca;
ax.TickDir = 'out';
box off
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickLength = [0.02 0.02];
ax.XTickLabels = {''};

f=gcf;
f.Position = [680   734   268   364]; 

[p,h] = ranksum(WPAG, HPAG)
% p = 0.8991

%% SC soma size plot

%% SOMA SIZE _ PLOT. 
allNAN = find(isnan(WSC));
WSC(allNAN) = []; 

allNAN = find(isnan(HSC));
HSC(allNAN) = []; 


% Box plot + DOTS
% boxplot([WPAG', HPAG'], [ones(1,19), ones(1,32)*2])

% Errorbar + points
nWT = numel(WSC);
nHET = numel(HSC);

mWT = nanmean(WSC);
mHET = nanmean(HSC);

semWT = nanstd(WSC)/sqrt(nWT);
semHET = nanstd(HSC)/sqrt(nHET);

x1 = ones(1, nWT);
x2 = ones(1, nHET)*2;

figure
scatter(x1, WSC,'SizeData', 100, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, HSC, 'SizeData', 100,'MarkerEdgeColor', [1 0.4 0.4], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
errorbar(1, mWT, semWT, 'o','Color', [0.2 0.2 0.2], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(2, mHET, semHET, 'o', 'Color', [1 0.2 0.2],'Marker', 'none', 'LineWidth', 1.75)
plot([0.75 1.25], [mWT mWT], 'k', 'LineWidth', 2)
plot([1.75 2.25], [mHET mHET], 'r', 'LineWidth', 2)

ylim([0 250])
xlim([0.5 2.5])
ax= gca;
ax.TickDir = 'out';
box off
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickLength = [0.02 0.02];
ax.XTickLabels = {''};

f=gcf;
f.Position = [680   734   268   364]; 

[p,h] = ranksum(WSC, HSC)
[h, p] = ttest2(WSC, HSC)
% 0.6688
nanmean(WSC) % 139.8750
nanmean(HSC) % 129.5714


mWT = nanmean(dataWT)
mHET =  nanmean(dataHET)
[p, h] = ranksum(mWT, mHET)

nanmean(mWT)
nanmean(mHET)


wtvals = dataWT(:,1); 
hetvals = dataHET(:,1); 

mWT = nanmean(wtvals)
semWT = nanstd(wtvals)/sqrt(6)

mWT = nanmean(hetvals)
semWT = nanstd(hetvals)/sqrt(10)


%% OPTOPATCH

% For 10Hz stimulation 
% 4 HET mice, 3 WT mice. 
% 04/10/21 - 07/10/21


% Copied and pasted data from Peter's analysis - from the excel
% spreadsheets. I then copied this data into the RESULTS numbers file found
% in the Burnett_et_al folder

%%

dataWT = [];
dataHET = []; 

save('EPSC_Relative_Amplitude_10Hz_C1_Setd5_OPTO.mat', 'dataWT', 'dataHET');

%% EPSC Amp

% dataWT/dataHET are arrays where each column is a different cell and each
% row is the EPSC max amplitude (pA) for each stimulation from 10Hz opto
% stimulation. 

col = 'r';
v = 10;

% MEAN/ SEM plot: 

% Make the data absolute because Peter's plot has the axis going from 0 at
% origin to -500 at top. i.e. flipped. Make abs then change axis labels
% later. 
dataWT = abs(dataWT);
dataHET = abs(dataHET);

% 
nWT = numel(dataWT(1,:)); 
nHET = numel(dataHET(1,:));

mean_WT = nanmean(dataWT'); 
mean_HET = nanmean(dataHET'); 

x = (1:1:v);

semWT = nanstd(dataWT')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(dataHET')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col, 'LineWidth', 1.3)


xticks([0,5,10])
xlabel('Stimulus #')
ylabel('EPSC Amplitude (pA)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.025 0.025];
xlim([0 11])
ylim([0 500])

yticks([0:100:500])
yticklabels({'0', '-100', '-200', '-300', '-400', '-500'})

f = gcf;
f.Position = [742   492   304   258]; 


%% STATS FOR OPTOPATCH
% 
% y = vertcat(dataWT', dataHET'); % remove last 3 cols from pre-wt
% gp1 = vertcat(ones(1,6)', (ones(1,10)*2)'); % GENO
% Geno = table([1,2], 'VariableNames', {'Geno'});

y = vertcat(dataWT, dataHET); % remove last 3 cols from pre-wt
gp1 = vertcat(ones(1,6)', (ones(1,10)*2)'); % GENO
Geno = table([1,2], 'VariableNames', {'Geno'});

t = array2table(y);
t2 = array2table(gp1);

t3 = [t2, t];

rm = fitrm(t3, 'y1-y10~gp1')
ranovatbl = ranova(rm)
     
rm.Coefficients
% [p, tbl, stats] = anovan(y, {gp1}, 'model', 1, 'varnames', {'Geno'})
% tbl = multcompare(stats, 'Dimension', [1,2], 'CType', 'bonferroni')


% var = exit_analysis.NumOutTrig;
% grp = exit_analysis.Group;
% grp2 = exit_analysis.geno;
%
% [p,t,stats] = anovan(var, {grp2, grp}, 'model', 'interaction')
% [c,m,h] = multcompare(stats, 'CType', 'bonferroni')


%% Cumulative

col = 'r';
v = 10;

% MEAN/ SEM plot: 

% Make the data absolute because Peter's plot has the axis going from 0 at
% origin to -500 at top. i.e. flipped. Make abs then change axis labels
% later. 
dataWT = abs(dataWT);
dataHET = abs(dataHET);

% 
nWT = numel(dataWT(1,:)); 
nHET = numel(dataHET(1,:));

mean_WT = nanmean(dataWT'); 
mean_HET = nanmean(dataHET'); 

x = (1:1:v);

semWT = nanstd(dataWT')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(dataHET')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col, 'LineWidth', 1.3)


xticks([0,5,10])
xlabel('Stimulus #')
ylabel('Cum. EPSC Amplitude (pA)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.025 0.025];
xlim([0 11])
ylim([0 4000])

yticks([0:1000:4000])
yticklabels({'0', '-1000', '-2000', '-3000', '-4000'})

f = gcf;
f.Position = [742   492   304   258]; 


%% Relative

col = 'r';
v = 10;

% MEAN/ SEM plot: 

nWT = numel(dataWT(1,:)); 
nHET = numel(dataHET(1,:));

mean_WT = nanmean(dataWT'); 
mean_HET = nanmean(dataHET'); 

x = (1:1:v);

semWT = nanstd(dataWT')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(dataHET')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col, 'LineWidth', 1.3)


xticks([0,5,10])
xlabel('Stimulus #')
ylabel('Rel. EPSC Amplitude')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.025 0.025];
xlim([0 11])
ylim([0 1.1])

yticks([0:0.5:1])

f = gcf;
f.Position = [742   492   304   258]; 









%% Traces!!! 

mousee = 2;

het_animals = [2,4,5,8];
if ismember(mousee, het_animals)
    col = 'r';
else 
    col = 'k';
end 

files = dir('*.abf');
nfiles = length(files);

%%

% Which file! 
i = 11; 
% M4 - WC - 3/19/31 (whole cell 10Hz)
% M4 - Firing to Opto = 7/8/9/23/24/34/35
% M7 - WC - 3/13/14/31
% M7 - Firing to Opto - 7/8/9/18/19/20/35/36/37

fn = files(i).name;
[d,si,h]=abfload(fn);
d = squeeze(d);
% n = numel(d(1,:));

n = 3;

%WT
colvals = linspace(0.8, 0.2, n);

figure
for j = 1:n
    subplot(n,1,j)
    plot(d(:, j), 'Color', [colvals(j), colvals(j), colvals(j)])
    hold on
end

%HET
% n=10;
colvals = linspace(0.8, 0.2, n);

figure
for j = 1:n
 subplot(n,1,j)   
%  plot(d(:, j), 'Color', [1, colvals(j), colvals(j)])
 plot(d(:, j), 'Color', col)
 hold on 
 box off
 ylim([-80 50])
 xlim([20000 80000])
 ax = gca;
 ax.XAxis.Visible = 'off'; 
 yticklabels({''})
end 
f = gcf;
f.Position = [216 42  378  1048]; 

%% PLOT1 - VC - OPTO 

xlim([20000 80000])
ylim([-300 50])
box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
f.Position = [59   681   786   353];

%% PLOT 2 - CC - OPTO
% xlim([20000 80000]) %10 HZ
xlim([20000 100000]) %20 Hz
ylim([-80 40])
box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
f.Position = [111   605   599   412]; 
  
  
  
  
  
  
  
  
%% Plot just one trace. 

% If just one. 
figure;
plot(d)

%% Cell attached  - 10Hz

figure 
plot(d(:,5), col)
xlim([20000 80000])

ylim([-50 40])

box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
% f.Position = [111  650  1463   367]; %large
 f.Position = [111   738   543   279]; %small

title('Cell Attached- 10Hz')


%% ****** - - - - - - - - Whole cell - 10Hz

% JUST ONE TRACE
figure 
plot(d(:,1), col)
xlim([20000 80000])
ylim([-150 100])
box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
% f.Position = [111  650  1463   367]; %large
%  f.Position = [111   738   543   279]; %small
 f.Position = [59   681   786   353];

title('Whole Cell - 10Hz')


% 5 TRACES
%HET 
n=5;
colvals = linspace(0.1, 0.7, n);

figure
for j = 1:n
 plot(d(:, j), 'Color', [1, colvals(j), colvals(j)])
 hold on 
end 
xlim([20000 80000])
ylim([-1000 100])
box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
 f.Position = [696   675   786   353];

%%
n_steps = 10; 
colvals = linspace(0.8, 0.2, n_steps);

figure
for j = 1:n_steps
 plot(d(:, j), 'Color', [colvals(j), colvals(j), colvals(j)])
 hold on 
end 
xlim([500 3000])
ylim([-130 10])
box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
% f.Position = [111  650  1463   367]; %large
 f.Position = [111   738   543   279]; %small
title('0.05 Ramp')

%% ******* - - - -  - - 0.2ms 10Hz opto 

%WT
colvals = linspace(0.8, 0.2, n);

figure
for j = 1:n
 plot(d(:, j), 'Color', [colvals(j), colvals(j), colvals(j)])
 hold on 
end 


%HET
colvals = linspace(0.8, 0.2, n);

figure
for j = 1:n
 plot(d(:, j), 'Color', [1, colvals(j), colvals(j)])
 hold on 
end 

xlim([20000 80000])
ylim([-80 40])
box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
% f.Position = [111  650  1463   367]; %large
%  f.Position = [111   738   543   279]; %small
  f.Position = [111   605   599   412]; 
% title('0.2ms 10Hz')

%% 10pA Current steps

colvals = linspace(0.8, 0.2, n);

figure
for j = 1:n
 plot(d(:, j), 'Color', [colvals(j), colvals(j), colvals(j)])
 hold on 
end 

xlim([0 80000])
ylim([-120 40])
box off
ax =gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out'; 
ax.FontSize = 22;
ax.TickLength = [0.015 0.015];
f = gcf;
% f.Position = [111  650  1463   367]; %large
 f.Position = [111   738   543   279]; %small
title('10pA Current Steps')




%% PAIRD PULSE RATIO - Ratio of PEAK AMPLUTIDE between the 2nd and 1st EPSCs in train - OPTOPATCH 

PPR_WT = []; 
PPR_HET = []; 

PPR_WT = array2table(PPR_WT, 'VariableNames', {'MouseID', 'CellID', 'PPR'});
PPR_HET = array2table(PPR_HET, 'VariableNames', {'MouseID', 'CellID', 'PPR'});

save('PPR_2ndto1st_Setd5.mat', 'PPR_WT', 'PPR_HET')

WPAG = PPR_WT.PPR;
HPAG = PPR_HET.PPR;

[p,h] = ranksum(WPAG, HPAG)

% Errorbar + points
nWT = numel(WPAG);
nHET = numel(HPAG);

mWT = nanmean(WPAG);
mHET = nanmean(HPAG);

semWT = nanstd(WPAG)/sqrt(nWT);
semHET = nanstd(HPAG)/sqrt(nHET);

x1 = ones(1, nWT);
x2 = ones(1, nHET)*2;

figure
scatter(x1, WPAG,'SizeData', 100, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, HPAG, 'SizeData', 100,'MarkerEdgeColor', [1 0.4 0.4], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
errorbar(1, mWT, semWT, 'o','Color', [0.2 0.2 0.2], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(2, mHET, semHET, 'o', 'Color', [1 0.2 0.2],'Marker', 'none', 'LineWidth', 1.75)
plot([0.75 1.25], [mWT mWT], 'k', 'LineWidth', 2)
plot([1.75 2.25], [mHET mHET], 'r', 'LineWidth', 2)

ylim([0.4 1.4])
xlim([0.5 2.5])
ax= gca;
ax.TickDir = 'out';
box off
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickLength = [0.02 0.02];
ax.XTickLabels = {''};

f=gcf;
f.Position = [680   734   268   364]; 















%% Probability of spiking per pulse - 10 pulses. 

% 1 - Define which animal the recording is from: 
% This will not change for several files. 
mousee = 8;
het_animals = [2,4,5,8];
if ismember(mousee, het_animals)
    col = 'r';
    genoo = 2;
else 
    col = 'k';
    genoo = 1;
end 

% Create directory of files in folder. 
files = dir('*.abf');
nfiles = length(files);

%% Set variables that will change with each file! % % % % % THIS NEEDS TO BE CHANGED EACH TIME. 

i = 24; 
celll = 4;
pulsewidth = '5';
pw = 5;
lightintensity = 5; 
freq = 10; 

% Open the file. 
fn = files(i).name;
[d,si,h]=abfload(fn);
d = squeeze(d);
n = numel(d(1,:));

% Plot the traces for each sweep. 
colvals = linspace(0.8, 0.2, n);

figure
for j = 1:n
    subplot(n,1,j)
    plot(d(:, j), 'Color', [colvals(j), colvals(j), colvals(j)])
    hold on
    box off
    axis off
end

% Plot all traces on top of each other. 
% figure
% for j = 1:n
%     plot(d(:, j), 'Color', [colvals(j), colvals(j), colvals(j)])
% hold on 
% end

%% Run through each sweep and determine if a spike is elicited after 

% LOAD THE SUMMARY TABLE
% load(strcat('/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/OPTOPATCH/PROB_SPIKING_TABLE.mat'), 'PROB_SPIKE')

% 
if freq == 10
    
    pulse1 = 26954;
    pulse2 = 32128;
    pulse3 = 37163;
    pulse4 = 42024;
    pulse5 = 47180;
    pulse6 = 52200;
    pulse7 = 57300;
    pulse8 = 62400;
    pulse9 = 67575;
    pulse10 = 72175;
    
    probSPIKE = zeros(n, 18);
    rngval = 1000;
    
    for jj = 1:n
        
        sweep = d(:, jj);
        ss = sign(sweep);
        sd = diff(ss);
        valsSPIKE = find(sd == 2);
        
        probSPIKE(jj, 1) = i;
        probSPIKE(jj, 2) = mousee;
        probSPIKE(jj, 3) = celll;
        probSPIKE(jj, 4) = jj;
        probSPIKE(jj, 5) = genoo;
        probSPIKE(jj, 6) = pw;
        probSPIKE(jj, 7) = lightintensity;
        probSPIKE(jj, 8) = freq;
        
        a1 = find(valsSPIKE > pulse1-rngval & valsSPIKE < pulse1+rngval);
        if ~isempty(a1)
            probSPIKE(jj, 9) = 1;
        end
        
        a2 = find(valsSPIKE > pulse2-rngval & valsSPIKE < pulse2+rngval);
        if ~isempty(a2)
            probSPIKE(jj, 10) = 1;
        end
        
        a3 = find(valsSPIKE > pulse3-rngval & valsSPIKE < pulse3+rngval);
        if ~isempty(a3)
            probSPIKE(jj, 11) = 1;
        end
        
        a4 = find(valsSPIKE > pulse4-rngval & valsSPIKE < pulse4+rngval);
        if ~isempty(a4)
            probSPIKE(jj, 12) = 1;
        end
        
        a5 = find(valsSPIKE > pulse5-rngval & valsSPIKE < pulse5+rngval);
        if ~isempty(a5)
            probSPIKE(jj, 13) = 1;
        end
        
        a6 = find(valsSPIKE > pulse6-rngval & valsSPIKE < pulse6+rngval);
        if ~isempty(a6)
            probSPIKE(jj, 14) = 1;
        end
        
        a7 = find(valsSPIKE > pulse7-rngval & valsSPIKE < pulse7+rngval);
        if ~isempty(a7)
            probSPIKE(jj, 15) = 1;
        end
        
        a8 = find(valsSPIKE > pulse8-rngval & valsSPIKE < pulse8+rngval);
        if ~isempty(a8)
            probSPIKE(jj, 16) = 1;
        end
        
        a9 = find(valsSPIKE > pulse9-rngval & valsSPIKE < pulse9+rngval);
        if ~isempty(a9)
            probSPIKE(jj, 17) = 1;
        end
        
        a10 = find(valsSPIKE > pulse10-rngval & valsSPIKE < pulse10+rngval);
        if ~isempty(a10)
            probSPIKE(jj, 18) = 1;
        end
        
        
    end
    
elseif freq == 20
    
    pulse1 = 27200;
    pulse2 = 29600;
    pulse3 = 32100;
    pulse4 = 34700;
    pulse5 = 37200;
    pulse6 = 39700;
    pulse7 = 42300;
    pulse8 = 44600;
    pulse9 = 47200;
    pulse10 = 49800;
    
    probSPIKE = zeros(n, 18);
    rngval = 1000;
    
    for jj = 1:n
        
        sweep = d(:, jj);
        ss = sign(sweep);
        sd = diff(ss);
        valsSPIKE = find(sd == 2);
        
        probSPIKE(jj, 1) = i;
        probSPIKE(jj, 2) = mousee;
        probSPIKE(jj, 3) = celll;
        probSPIKE(jj, 4) = jj;
        probSPIKE(jj, 5) = genoo;
        probSPIKE(jj, 6) = pw;
        probSPIKE(jj, 7) = lightintensity;
        probSPIKE(jj, 8) = freq;
        
        a1 = find(valsSPIKE > pulse1-rngval & valsSPIKE < pulse1+rngval);
        if ~isempty(a1)
            probSPIKE(jj, 9) = 1;
        end
        
        a2 = find(valsSPIKE > pulse2-rngval & valsSPIKE < pulse2+rngval);
        if ~isempty(a2)
            probSPIKE(jj, 10) = 1;
        end
        
        a3 = find(valsSPIKE > pulse3-rngval & valsSPIKE < pulse3+rngval);
        if ~isempty(a3)
            probSPIKE(jj, 11) = 1;
        end
        
        a4 = find(valsSPIKE > pulse4-rngval & valsSPIKE < pulse4+rngval);
        if ~isempty(a4)
            probSPIKE(jj, 12) = 1;
        end
        
        a5 = find(valsSPIKE > pulse5-rngval & valsSPIKE < pulse5+rngval);
        if ~isempty(a5)
            probSPIKE(jj, 13) = 1;
        end
        
        a6 = find(valsSPIKE > pulse6-rngval & valsSPIKE < pulse6+rngval);
        if ~isempty(a6)
            probSPIKE(jj, 14) = 1;
        end
        
        a7 = find(valsSPIKE > pulse7-rngval & valsSPIKE < pulse7+rngval);
        if ~isempty(a7)
            probSPIKE(jj, 15) = 1;
        end
        
        a8 = find(valsSPIKE > pulse8-rngval & valsSPIKE < pulse8+rngval);
        if ~isempty(a8)
            probSPIKE(jj, 16) = 1;
        end
        
        a9 = find(valsSPIKE > pulse9-rngval & valsSPIKE < pulse9+rngval);
        if ~isempty(a9)
            probSPIKE(jj, 17) = 1;
        end
        
        a10 = find(valsSPIKE > pulse10-rngval & valsSPIKE < pulse10+rngval);
        if ~isempty(a10)
            probSPIKE(jj, 18) = 1;
        end
        
        
    end
    
    
end

pSP_Table = array2table(probSPIKE, 'VariableNames', {'FileID', 'Mouse', 'Cell', 'Sweep', 'Geno', 'PulseWidth', 'Light', 'Frequency', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10'});


av_table = pSP_Table(1, :);

for k= 1:10
   av_table{1,9+k-1} = mean(table2array(pSP_Table(:,9+k-1))); 
end 

save(strcat('/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/OPTOPATCH/ProbSpiking_Tables/ProbSpiking_M', string(mousee), '_C', string(celll), '_', pulsewidth, 'ms_', string(lightintensity), '_', string(freq), 'Hz.mat'), 'pSP_Table', 'av_table')

% COMBINE
% PROB_SPIKE = table();
PROB_SPIKE = vertcat(PROB_SPIKE, av_table);

save(strcat('/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/OPTOPATCH/PROB_SPIKING_TABLE.mat'), 'PROB_SPIKE')


%%

psp_results = zeros(16,8);

for i = 1:16
    
    pw = psp_results(i, 1);
    lp = psp_results(i, 2);
    fr = psp_results(i, 3);
    
    allWT = find(PROB_SPIKE.Geno ==1 & PROB_SPIKE.PulseWidth == pw & PROB_SPIKE.Light == lp & PROB_SPIKE.Frequency == fr);
    allHET = find(PROB_SPIKE.Geno ==2 & PROB_SPIKE.PulseWidth == pw & PROB_SPIKE.Light == lp & PROB_SPIKE.Frequency == fr);
    
    nWT = numel(allWT);
    nHET = numel(allHET);
    
    psp_results(i, 4) = nWT;
    psp_results(i, 5) = nHET;
    
    %
    if nWT>1
        valsWT = PROB_SPIKE.P1(allWT);
        mWT = mean(valsWT);
    elseif nWT == 1
        mWT = PROB_SPIKE.P1(allWT);
    else
        mWT = NaN;
    end
    psp_results(i, 6) = mWT;
    
    
    if nHET>1
        valsHET = PROB_SPIKE.P1(allHET);
        mHET = mean(valsHET);
    elseif nHET == 1
        mHET = PROB_SPIKE.P1(allHET);
    else
        mHET = NaN;
    end
    psp_results(i, 7) = mHET;
    
    if nWT >=1 && nHET >= 1
        [p,h] = ranksum(valsWT, valsHET);
        psp_results(i, 8) = p;
    else
        psp_results(i, 8) = NaN;
    end
end

psp_results_table = array2table(psp_results, 'VariableNames', {'PW', 'LP', 'FR', 'nWT', 'nHET', 'mWT', 'mHET', 'p_RS'});

save(strcat('/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/OPTOPATCH/PROB_SPIKING_RESULTS_P_AV.mat'), 'psp_results_table')
%% ALL CELLS 


    allWT = find(PROB_SPIKE.Geno ==1  & PROB_SPIKE.Frequency == 10 & PROB_SPIKE.Light == 25); %& PROB_SPIKE.PulseWidth == 5
    allHET = find(PROB_SPIKE.Geno ==2 & PROB_SPIKE.Frequency == 10  & PROB_SPIKE.Light == 25); % & PROB_SPIKE.PulseWidth == 5
    
    valsWT = PROB_SPIKE.P1(allWT);
    valsHET = PROB_SPIKE.P1(allHET);
    
    mWT = mean(valsWT)
    mHET = mean(valsHET)
    
    [p,h] = ranksum(valsWT, valsHET)
    
    
    
    %%
    
    for jj = 1:35 
        av_p  = mean(PROB_SPIKE{jj, 9:18});
        PROB_SPIKE.P_AV(jj) = av_p;
    end
    
    
    
    
%% Compare Wt/HEt/ with and without DTX - spiking at 120pA 
    
    % DATA:
    % '/Users/lauraburnett/Data_Analysis_Mac/PATCH/220507_PATCH_DATA_3Cs_PAG_Setd5_DTX.mat

    % row 13 = 120pA
    
    vWT = valsWT2(13,:);
    n1 = numel(vWT);
    mWT = nanmean(vWT);
    
    vHET = valsHET2(13,:);
    n3 = numel(vHET);
    mHET = nanmean(vHET);
      
%     vWT_DTX = DTX_WT(13,:);
%     n2 = numel(vWT_DTX);
%     mWT_DTX = nanmean(vWT_DTX);
    
    vHET_DTX = DTX_HET(13,:);
    n4 = numel(vHET_DTX);
    mHET_DTX = nanmean(vHET_DTX);  
    
    y = vertcat(vWT', vWT_DTX', vHET', vHET_DTX');
    
    gp1 = vertcat(ones(1,n1)', ones(1,n2)', (ones(1,n3)*2)', (ones(1,n4)*2)'); % GENO
    gp2 = vertcat(zeros(1,n1)', ones(1,n2)', (zeros(1,n3))', (ones(1,n4))'); % DTX
    
    gp3 = vertcat(ones(1,n1)', (ones(1,n2)')*2, (ones(1,n3)')*3, ((ones(1,n4))')*4); % DTX

    % STATISTICS - use n-way analysis of variance. Not repeated measures
    % since only 1 time point. 
    
    % Test for equal variance in each group
    % Bartlett's test
  [p, stats] = vartestn(y, gp3)
  
  % Test for normally distributed data
  [h,p] = kstest(vWT) % H = 1 - not normal
    [h,p] = kstest(vHET) % H = 1 - not normal
      [h,p] = kstest(vWT_DTX) % H = 1 - not normal
        [h,p] = kstest(vHET_DTX) % H = 1 - not normal
  
[p, tbl, stats] = kruskalwallis(y, gp3)  
%     [p, tbl, stats] = anovan(y, {gp1, gp2}, 'model', 'interaction', 'varnames', {'Geno', 'DTX'})
    results = multcompare(stats, 'Dimension', [1,2], 'CType', 'bonferroni')
    
    % RESULTS
   % 1.0000    2.0000  -19.9912   -2.3056   15.3801    1.0000
   % 1.0000    3.0000    1.9338   19.6194   37.3051    0.0206
   % 1.0000    4.0000  -27.6379   -7.8248   11.9883    1.0000
   % 2.0000    3.0000    4.7111   21.9250   39.1389    0.0047
   % 2.0000    4.0000  -24.9125   -5.5192   13.8740    1.0000
   % 3.0000    4.0000  -46.8375  -27.4442   -8.0510    0.0011
   
%% THIS IS TOO MUCH ^^^^^^ 

% Just want to compare 

% 1 - Setd5++ without DTX to  Setd5++ WITH DTX
[p,h] = ranksum(vWT, vWT_DTX)
% p = 0.4824

% 2 - Setd5+- without DTX to  Setd5+- WITH DTX
[p,h] = ranksum(vHET, vHET_DTX)
% p = 6.5020e-04

% 3 - Setd5++ vs Setd5 +- BEFORE DTX
[p,h] = ranksum(vWT, vHET)
% p  = 7.2778e-04

% 4 - Setd5++ vs Setd5+- AFTER DTX
[p,h] = ranksum(vWT_DTX, vHET_DTX)
% p = 0.3863

% Multiple ranksums with correction for multiple testing?? 

    
    
    %% Repeated measues anova - Statistics

% Use this two-way repeated measure ANOVA test to test for the effect of
% a-DTX on the firing of PAG cells - looking at ALL amounts of current
% injected!! Not just 120pA. 

% repeated because same cells over time. Is it repeated measures because
% these are different cells....

% But non-parametric?? Which test to use?? 

% DATA of interest are:

% ValsWT2
% valsHET2
% DTX_WT
% DTX_HET

PRE_WT = valsWT2';
PRE_HET = valsHET2';

% POST_WT = DTX_WT'; 
POST_HET = valsHETDTX'; 

% 1 - CHECK WT vs WT POST

y = vertcat(PRE_WT, POST_WT); % remove last 3 cols from pre-wt
gp1 = vertcat(ones(1,n1)', (ones(1,n2)*2)'); % GENO

t = array2table(y);
t2 = array2table(gp1);

t3 = [t,t2];

rm = fitrm(t3, 'y1-y20~gp1') %'WithinDesign',1:1:7, WithinModel', 'separatemeans'
ranovatbl = ranova(rm)

% 2 - CHECK HET vs HET POST

y = vertcat(PRE_HET, POST_HET); % remove last 3 cols from pre-wt
gp1 = vertcat(ones(1,nHET)', (ones(1,nHETDTX)*2)'); % GENO

t = array2table(y);
t2 = array2table(gp1);

t3 = [t,t2];

rm = fitrm(t3, 'y1-y20~gp1') %'WithinDesign',1:1:7, WithinModel', 'separatemeans'
ranovatbl = ranova(rm)

 
% 3 - CHECK WT vs HET PRE

y = vertcat(PRE_WT, PRE_HET); % remove last 3 cols from pre-wt
gp1 = vertcat(ones(1,nWT)', (ones(1,nHET)*2)'); % GENO

t = array2table(y);
t2 = array2table(gp1);

t3 = [t,t2];

rm = fitrm(t3, 'y1-y24~gp1') %'WithinDesign',1:1:7, WithinModel', 'separatemeans'
ranovatbl = ranova(rm)
    
% OTHER TESTS 

[p, tbl, stats] = friedman()

[p, tbl, stats] = kruskalwallis(y', gp1')  
results = multcompare(stats, 'Dimension', [1,2], 'CType', 'bonferroni')
    
    

%% Adding new column to SPIKE_DATA - group - includes DTX/ Geno


allWT = find(SPIKE_DATA.Geno == 1 & SPIKE_DATA.DTX ==0);
allHET = find(SPIKE_DATA.Geno == 2 & SPIKE_DATA.DTX ==0);
allHETDTX = find(SPIKE_DATA.Geno == 2 & SPIKE_DATA.DTX ==1);

valsWT2 = table2array(SPIKE_DATA(allWT, 4:27));
valsHET2 = table2array(SPIKE_DATA(allHET, 4:27));
valsHETDTX = table2array(SPIKE_DATA(allHETDTX, 4:27));

nWT = numel(valsWT2(:,1)); 
nHET = numel(valsHET2(:,1));
nHETD = numel(valsHETDTX(:,1));


y = vertcat(valsWT2, valsHET2, valsHETDTX); % remove last 3 cols from pre-wt
Geno = vertcat(ones(1,nWT)', (ones(1,nHET)*2)', (ones(1,nHETD)*2)'); % GENO
DTX = vertcat(ones(1,nWT)', ones(1, nHET)', (ones(1,nHETD)*2)');
GP = t3.Var27;

[p, tbl, stats] = kruskalwallis(y', GP')  
results = multcompare(stats, 'Dimension', [1,2], 'CType', 'Bonferroni')
   

t = array2table(y);
G = array2table(Geno);
D = array2table(DTX);

t3 = [t,G,D];

rm = fitrm(t3, 'y1-y24~Geno+DTX');
ranovatbl = ranova(rm)

%%

% FOR REPEATED MEASURES CAN WORK FROm SPIKE_DATA TABLE!!!! 
rm = fitrm(SPIKE_DATA, 'ST1-ST24~Geno'); % Cannot do interaction between geno and DTX at the moment because don't have WT DTX. Insufficient. 
ranovatbl = ranova(rm)
% test for sphericity
tbl = mauchly(rm)

% WT versus HET after DTX
rm = fitrm(SPIKE_DATA3, 'ST1-ST24~DTX'); % Cannot do interaction between geno and DTX at the moment because don't have WT DTX. Insufficient. 
ranovatbl = ranova(rm)

% HET before versus HET after DTX
rm = fitrm(SPIKE_DATA4, 'ST1-ST24~DTX'); % Cannot do interaction between geno and DTX at the moment because don't have WT DTX. Insufficient. 
ranovatbl = ranova(rm)


%% Test at highest current inj
y = SPIKE_DATA.ST24;
GP = SPIKE_DATA.Geno;
D = SPIKE_DATA.DTX;
GROUP = t3.Var27;

% 
% [p, tbl, stats] = anovan(Y2, {Geno, DTX}, 'varnames', {'Geno', 'DTX'})
% results = multcompare(stats, 'Dimension', [1,2], 'CType', 'tukey-kramer')










    %% BOXPLOT
    
    
%%%%% Combine arrays %%%% 
speed_ALL = horzcat(vWT, vWT_DTX, vHET, vHET_DTX);

x1 = ones(1, n1);
x2 = ones(1, n2)*2;
x3 = ones(1, n3)*3;
x4 = ones(1, n4)*4;

% PLOT 
figure
scatter(x1, vWT,'SizeData', 200, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, vWT_DTX, 'SizeData', 200,'MarkerEdgeColor', [0 0 0.9], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
scatter(x3, vHET, 'SizeData', 200,'MarkerEdgeColor', 'r', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
scatter(x4, vHET_DTX, 'SizeData', 200,'MarkerEdgeColor', [0.8 0 1], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot(y, {gp1, gp2}, 'Colors', 'k');
set(b, 'linew', 1.25);

xlim([0.5 4.5])
ylim([-10 120])

box off
xticks([1,2,3,4])
xticklabels({''})
ax = gca;
ax.FontSize = 25;
ax.TickDir = 'out'; 
ax.LineWidth = 2;

f = gcf;
f.Position = [704   207   355   572]; 

    
    
    
    
    %% EXAMPLE TRACES
    
% For DTX PLOTS 

% Row 8 = 50pA 
% Row 13 = 100pA

% Wt - with DTX = C3 - M1 - Cell4 - File 011 - 50pA
% WT - Without DTX = C2 - M6- Cell1 - File 007 - 50pA 

% Create directory of files in folder. 
files = dir('*.abf');
nfiles = length(files);

%% Set variables that will change with each file! % % % % % THIS NEEDS TO BE CHANGED EACH TIME. 

i = 8; 

% Open the file. 
fn = files(i).name;
[d,si,h]=abfload(fn);
d = squeeze(d);
n = numel(d(1,:));

subplot(2,1,1)
plot(d(1:75000, 8), 'Color', 'b')
box off
ylim([-80 60])
xlim([0 75000])





%% 

y = vertcat(SC_WT', SC_HET'); % remove last 3 cols from pre-wt
gp1 = vertcat(ones(1,15)', (ones(1,15)*2)'); % GENO

t = array2table(y);
t2 = array2table(gp1);

t3 = [t,t2];

% tbl = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
[tbl,rm] = simple_mixed_anova(y, gp1, {'CurrentInput'}, {'Genotype'})

tbl = multcompare(rm, 'Genotype', 'By', 'CurrentInput', 'ComparisonType', 'bonferroni')
% M= margmean(rm, {'CurrentInput', 'Genotype'})


 
    %% STATS!!!
   
    % Current injection data here: '/Users/lauraburnett/Data_Analysis_Mac/PATCH/220507_PATCH_DATA_3Cs_PAG_Setd5_DTX.mat'
    
%% 1-  one way - Repeated measures ANOVA for Current Injections!     
% % 19 WT cells, 24 HET cells
% Within factors = Current Input - 23 levels
% Between factors = Genotype - 2 levels

nWT = numel(valsWT2(1,:)); % Number of cells - WT
nHET = numel(valsHET2(1,:)); % Number of cells - WT

between_factors = vertcat(ones(nWT, 1), ones(nHET,1)*2);
arrayy = vertcat(valsWT2', valsHET2');

% tbl = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
[tbl,rm] = simple_mixed_anova(arrayy, between_factors, {'CurrentInput'}, {'Genotype'})

tbl = multcompare(rm, 'Genotype', 'By', 'CurrentInput', 'ComparisonType', 'bonferroni')
% M= margmean(rm, {'CurrentInput', 'Genotype'})

%% 2 - Repeated measures ANOVA model - 2 way! - Genotype and DTX

% % 19 + 20WT cells, 24 +13 HET cells

% Within factor 1 = Current Input - 20 levels

% Between factors = Genotype - 2 levels
   % adding DTX as 2nd between factor - 2 levels - 1 or 0 
% - - - - - - 
    
% Need to 'chop' the number of current inpput levels because DTX-WT only
% has 20 steps. 

nWT = numel(valsWT2(1,:)); % Number of cells - WT
nHET = numel(valsHET2(1,:)); % Number of cells - WT

% nWTDTX = numel(DTX_WT(1,:));
nHETDTX = numel(valsHETDTX(1,:));

arrayy = vertcat(valsWT2(1:20, :)', valsHET2(1:20, :)', valsHETDTX(1:20, :)');

between_factors(:, 1) = vertcat(ones(nWT, 1),  ones(nHET,1)*2,  ones(nHETDTX,1)*2);
between_factors(:, 2) = vertcat(zeros(nWT, 1),  zeros(nHET,1),  ones(nHETDTX,1));
    

[tbl,rm] = simple_mixed_anova(arrayy, between_factors, {'CurrentInput'}, {'Genotype', 'DTX'})  

tbl = multcompare(rm, 'Genotype', 'By', 'DTX', 'ComparisonType', 'bonferroni')

tbl = multcompare(rm, 'DTX', 'By', 'Genotype', 'ComparisonType', 'bonferroni')

%%

nWT = numel(valsWT2(1,:)); % Number of cells - WT

nHET = numel(valsHET2(1,:)); % Number of cells - WT
nHETDTX = numel(valsHETDTX(1,:));


  y = vertcat(valsWT2(13, :)', valsHET2(13, :)', valsHETDTX(13, :)');
  gp1 = vertcat(ones(nWT, 1), zeros(nHET,1), zeros(nHETDTX,1)); % genotype
  gp2 = vertcat(zeros(nWT, 1), zeros(nHET,1), ones(nHETDTX,1));  %dtx 
    
[p, tbl, stats] = anovan(y, {gp1, gp2}, 'model', 2, 'varnames', {'Geno', 'DTX'})
tbl = multcompare(stats, 'Dimension', [1,2], 'CType', 'bonferroni')


    
 
 %% 3 -  Repeated measures ANOVA for Current Injections - with / without DTx 
 
 %  A - WT 
 
% % 19 WT cells, 24 HET cells
% Within factors = Current Input - 20 levels
% Between factors = DTX - 2 levels

nWT = numel(valsWT2(1,:)); % Number of cells - WT
nWTDTX = numel(DTX_WT(1,:));

between_factors = vertcat(ones(nWT, 1), zeros(nWTDTX,1));
arrayy = vertcat(valsWT2(1:20, :)', DTX_WT');

% tbl = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
[tbl,rm] = simple_mixed_anova(arrayy, between_factors, {'CurrentInput'}, {'DTX'})

 %  B - HET 
 
%  24 HET cells, 13 HET_DTX cells
% Within factors = Current Input - 23 levels
% Between factors = DTX - 2 levels

nHET = numel(valsHET2(1,:)); % Number of cells - WT
nHETDTX = numel(DTX_HET(1,:));

between_factors = vertcat(ones(nHET, 1), zeros(nHETDTX,1));
arrayy = vertcat(valsHET2', DTX_HET');

% tbl = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
[tbl,rm] = simple_mixed_anova(arrayy, between_factors, {'CurrentInput'}, {'DTX'})

tbl = multcompare(rm, 'DTX', 'By', 'CurrentInput', 'ComparisonType', 'bonferroni')


% C - WT without a-DTX versus HET with a-DTX

nWT = numel(valsWT2(1,:)); 
nHETDTX = numel(DTX_HET(1,:));
    
between_factors = vertcat(ones(nWT, 1), zeros(nHETDTX,1));
arrayy = vertcat(valsWT2', DTX_HET'); 
    
[tbl,rm] = simple_mixed_anova(arrayy, between_factors, {'CurrentInput'}, {'DTX'})
tbl = multcompare(rm, 'DTX', 'By', 'CurrentInput', 'ComparisonType', 'bonferroni')
 
  

%%
    
nWT = numel(valsWT2(1,:)); % Number of cells - WT
nWTDTX = numel(DTX_WT(1,:));
nHET = numel(valsHET2(1,:)); % Number of cells - WT
nHETDTX = numel(DTX_HET(1,:));


  y = vertcat(valsWT2(13, :)', valsHET2(13, :)', DTX_WT(13, :)', DTX_HET(13, :)');
  gp1 = vertcat(ones(nWT, 1), zeros(nHET,1), ones(nWTDTX,1), zeros(nHETDTX,1)); % genotype
  gp2 = vertcat(zeros(nWT, 1), zeros(nHET,1), ones(nWTDTX,1), ones(nHETDTX,1));  %dtx 
    
[p, tbl, stats] = anovan(y, {gp1, gp2}, 'model', 2, 'varnames', {'Geno', 'DTX'})
tbl = multcompare(stats, 'Dimension', [1,2], 'CType', 'bonferroni')



%% MEAN SEM - Firing to Current Input - SUPERIOR cOlICULUS CELLS
% Mean SEM - WT/ HET - Red/Black

col = 'r';
v = 20;
inj_steps = curr_inj(3:3+v, 2);

nWT = numel(SC_WT(1,:)); 
nHET = numel(SC_HET(1,:));

mean_WT = nanmean(SC_WT'); 
mean_HET = nanmean(SC_HET'); 

x = (1:1:v);

semWT = nanstd(SC_WT'); %/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(SC_HET'); %/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', col, 'LineWidth', 1.3)

% xticks(1:1:15)
% xticklabels(string(inj_steps(1:15)))

xticks([1:4:15])
xticklabels(string(inj_steps(1:4:15)))

xlabel('Current (pA)')
ylabel('Firing Rate (Hz)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.02 0.02];
xlim([0.5 15])
ylim([0 100])

f = gcf;
f.Position = [742   492   304   258]; 



%% STATS

between_factors = vertcat(ones(nWT, 1), ones(nHET,1)*2);
arrayy = vertcat(SC_WT', SC_HET');

% tbl = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
[tbl,rm] = simple_mixed_anova(arrayy, between_factors, {'CurrentInput'}, {'Genotype'})

tbl = multcompare(rm, 'Genotype', 'By', 'CurrentInput', 'ComparisonType', 'bonferroni')
% M= margmean(rm, {'CurrentInput', 'Genotype'})



%% 

figure

for i = 1:6
    xvals = 1:1:10;
    yvals = dataWT(i, :);
    plot(xvals, yvals, 'k-');
    hold on
end

for i = 1:10
    xvals = 1:1:10;
    yvals = dataHET(i, :);
    plot(xvals, yvals, 'r-');
    hold on
end

nanmean(nanmean(dataWT))




