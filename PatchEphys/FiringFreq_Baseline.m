% Firing frequency at baseline
% Burnett - 28/02/23


% Open the cell-attached file. 
% Remove at least 1 minute off the beginning of the recording. ( 1200 000 )
% Find the firing frequency 
% Save figs of 5s and 60s time. 

%% Path to save the figs. 
fig_save_path = '/Users/lauraburnett/Data_Analysis_Mac/PATCH/ANALYSIS/FiringFreq-Baseline/FIGS_5s';

%% Load file. 

files = dir('*.abf');
nfiles = length(files);

cohort = 2;
mouse = 8;
celll = 1; 
glu = 0;

%% LOAD AND INVESTIGATE

fn = files(i).name;
[d,si,h]=abfload(fn);

figure; plot(d)

% Total length of the recording
total_length = length(d)/20000;

%% PLOT - 5s period of time. 

figure; plot(d(end-100000:end))
% figure; plot(d(500000:600000))
% figure; plot(d(end-350000:end-250000))
% figure; plot(d(350000:450000))

f = gcf;
f.Position = [300   381   443   304];
box off
ylim([-100 50])
xlim([0 100000])

% Save the figure
% savefig(f, fullfile(fig_save_path, 'Setd5_Cohort1_M3_C1_5s_Baseline.fig'))
fname = strcat('Setd5_Cohort', string(cohort), '_M', string(mouse), '_C', string(celll), '_5s_Baseline.svg');
saveas(f, fullfile(fig_save_path, fname))


%% Quantify 

d3 = d(end-1200000:end);
% d3 = d(end-1400000:end-200000);
% d3 = d(end-600000:end);
% d3 = d(1:300000);

% Find peaks in the cropped recording. 
[pks, loc] = findpeaks(d3, 'MinPeakHeight', 10, 'MinPeakDistance', 15);

% Plot the PEAKS - check that they match with the spikes. 
figure; plot(1:1:length(d3), d3,loc,pks,'o') 

% Quantify the time of the recording and the number of spikes in that
% recording to find the firing frequency in Hz. 
time_rec = length(d3)/20000;
spikes_rec = numel(pks);

firing_freq_base = (spikes_rec/time_rec)*60; 

%% Plot 60s period of time. 

figure; plot(d3);
f = gcf;
f.Position = [554   383   725   220];
box off
ylim([-100 50])
xlim([0 1200000])

% Save the figure
fname = strcat('Setd5_Cohort', string(cohort), '_M', string(mouse), '_C', string(celll), '_60s_Baseline.svg');
saveas(f, fullfile(fig_save_path, fname))


%% Enter data into a table. 

% Generate table. 
% base_firing = table('Size', [45, 8], 'VariableTypes', repmat({'double'}, 1,8), 'VariableNames', {'Cohort', 'Mouse',  'Cell', 'Glu', 'T_tot', 'T_rec', 'S_rec', 'Freq'});

% Enter data into table. 
base_firing.Cohort(kk) = cohort; 
base_firing.Mouse(kk) = mouse; 
base_firing.Cell(kk) = celll; 
base_firing.Glu(kk) = glu; 
base_firing.T_tot(kk) = total_length; 
base_firing.T_rec(kk) = time_rec; 
base_firing.S_rec(kk) = spikes_rec; 
base_firing.Freq(kk) = firing_freq_base; 

close all
clearvars -except kk base_firing fig_save_path files
kk= kk+1; 

%%

save('230228_FreqFiring_Baseline.mat', 'base_firing');




%% PLOT - FIRING PER GROUP - BOX And SCATTER

[C, ia, ic] = unique(base_firing(:, [4,9]));

for i = 1:4
    
    rows = find(ic == i);
    vals = base_firing.Freq(rows);
    n_vals(i) = numel(vals);
    valuess{i} = vals;
    
end 

x1 = ones(1, n_vals(1))*4;
x2 = ones(1, n_vals(2))*3;
x3 = ones(1, n_vals(3))*2;
x4 = ones(1, n_vals(4))*1;

figure
scatter(x1, valuess{1},'SizeData', 100, 'MarkerEdgeColor', 'm', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, valuess{2}, 'SizeData', 100,'MarkerEdgeColor', 'b', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
scatter(x3, valuess{3},'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
scatter(x4, valuess{4}, 'SizeData', 100,'MarkerEdgeColor', 'r', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([valuess{1}', valuess{2}', valuess{3}', valuess{4}'], [ones(1,n_vals(1))*4, ones(1,n_vals(2))*3, ones(1, n_vals(3))*2, ones(1, n_vals(4))], 'Colors', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
box off
xlim([0.5 4.5])
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ax.TickLength = [0.025 0.025];
ax.LineWidth = 2;


%%

data2 = []; 

[p, tbl, stats] = kruskalwallis(data)

g = vertcat(ones(17,1), ones(17,1)*2, ones(17,1)*3, ones(17,1)*4);
dunn(data2(:,1)', data2(:,2)')

save('FREQ_FIRING_DATA.mat', 'data', 'data2')












