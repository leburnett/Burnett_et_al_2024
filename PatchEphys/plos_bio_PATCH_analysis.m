%%  Analyse Patch data for PLOS BIOLogy Revisions
% Burnett - 19 Feb 2024

clear
close all

animal = 'M8';

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

% FIRING versus CURRENT STEPS

% Which rows correspond to the 10pA Rect Step experiments? 
% C3 - M1: [2,3,4,5,8,10,12];
% C3 - M2: 

%Find rows with 'sweeps'
cell_rows = []; 

for kk = 1:numel(d_all(:,1))
    if ~isempty(d_all{kk,8})
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

spikedata = zeros(nsteps, n_cells); 

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
            spikedata(ji, i2) = 0;
        end
    end

end 

save_folder = '/Users/burnettl/Documents/ISTA/PLOS_BIOL_REVISIONS/PATCH/DATA_C1/CURRENT_STEP_TRACES/10pA_steps_spike_data';
fname = ['Spike_Data_Setd5_C1_', animal, '.mat']; 
save(fullfile(save_folder, fname),  'spikedata', 'd_all');








%% Combine spiking data from the different animals to make a table.
clear
SPIKE_DATA = [];

files = dir('*.mat');
het_animals = [2,4,5,8];

for i = 1:numel(files)
    fname = files(i).name;
    load(fname)
    
    n_cells = numel(spikedata(1,:));
    data = spikedata'; 

    if length(data)<36
        current_len = length(data);
        array_to_add = NaN(n_cells, 36-current_len); 
        data = [data, array_to_add];
    end 

    ani_num = str2double(fname(end-4));
    
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

SPIKE_DATA = array2table(SPIKE_DATA, 'VariableNames', {'Ani', 'Cell', 'Geno', 'ST1', 'ST2', 'ST3', 'ST4', 'ST5', 'ST6', 'ST7', 'ST8', 'ST9', 'ST10', 'ST11', 'ST12', 'ST13', 'ST14', 'ST15', 'ST16', 'ST17', 'ST18', 'ST19', 'ST20', 'ST21', 'ST22', 'ST23', 'ST24', 'ST25', 'ST26', 'ST27', 'ST28', 'ST29', 'ST30', 'ST31', 'ST32', 'ST33', 'ST34', 'ST35','ST36'});

save('240219_SPIKE_DATA_Setd5_10pA_Steps.mat', 'SPIKE_DATA')
% Use this to plot C-I trace per cell! 


%% MEAN - SEM - PATCH 

col = 'r';
allWT = find(SPIKE_DATA.Geno == 1 );
allHET = find(SPIKE_DATA.Geno == 2 );

valsWT2 = table2array(SPIKE_DATA(allWT, 6:26))';
valsHET2 = table2array(SPIKE_DATA(allHET, 6:26))';

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

xticks([1:4:v])
xticklabels(string(inj_steps(1:4:v)))

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


%% Plot for the individual cells. 

save_path = '/Users/burnettl/Documents/ISTA/PLOS_BIOL_REVISIONS/PATCH/PLOTS/CI_10pA';

for idx = 1:height(SPIKE_DATA)

    close
    animal = num2str(SPIKE_DATA.Ani(idx));
    cell_id = num2str(SPIKE_DATA.Cell(idx)); 
    genoo = SPIKE_DATA.Geno(idx);
    data = table2array(SPIKE_DATA(idx, 6:26))';

    if genoo == 1
        col = 'k';
    elseif genoo == 2
        col = 'r';
    end

    v = 21;
    inj_steps = curr_inj(3:3+v, 2);
    figure
    plot(data(1:v), col, 'LineWidth', 1.3)
    hold on 
    plot([1, 21], [0, 100], 'k:', 'LineWidth', 1)
    xticks([1:4:v])
    xticklabels(string(inj_steps(1:4:v)))
    xlabel('Current (pA)')
    ylabel('Firing Rate (Hz)')
    box off
    set(gca, 'FontSize', 16)
    ax = gca;
    ax.TickDir  = 'out';
    ax.LineWidth = 1.75;
    ax.TickLength = [0.02 0.02];
    xlim([0.5 24])
    ylim([0 100])
    title(strcat('M' , animal, '- C', cell_id))
    f = gcf;
    f.Position = [742   492   304   258]; 
    fname = strcat("CI_plot_10pA_M", animal, "_C", cell_id, ".svg");
    saveas(f, fullfile(save_path, fname));

end 





nWT = numel(find(SPIKE_DATA.Geno == 1)); % Number of cells - WT
nHET = numel(find(SPIKE_DATA.Geno == 2));  % Number of cells - HET

y = SPIKE_DATA{:, 5:end};

gp1 = vertcat(ones(nWT, 1), zeros(nHET, 1));


datatbl = array2table(SPIKE_DATA);
% rm = fitrm(datatbl, 'ST1-ST23~Geno');
% ranovatbl = ranova(rm)
% tbl = multcompare(rm, 'Group')


[tbl, rm] = simple_mixed_anova(y, gp1, {'CurrentInput'}, {'Geno'})
tbl = multcompare(rm, 'Geno', 'By', 'CurrentInput')























































%% OPTOGENETICS PATCH DATA 

%% % % % % % % % % % - Current Clamp  % % % % % % % % % %

clear all
close all

CC_path = '/Users/burnettl/Documents/ISTA/PLOS_BIOL_REVISIONS/PATCH/PLOTS/CC_10Hz'; 

animal = 'M8';
cell_id = 'C4';
file_to_open = 2; 

%%
files = dir('*.abf');
nfiles = length(files);

fn = files(file_to_open).name;
[d,si,h]=abfload(fn);
d = squeeze(d);
n = numel(d(1,:));

%% #1 - plot all reps: 

offset = 0;
figure;
for k = 1:n
    plot(d(:, k)-offset, 'k')
    hold on
    offset = offset+120;
end 
title([animal, '-', cell_id])
xlim([0, 100000])
box off
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.XAxis.Visible = 'off';
ylabel('Voltage (mV)');
f = gcf; 
f.Position = [601   359   387   559];

filename  = ['CC_10Hz_AllReps_', animal, '_', cell_id, '.svg'];
saveas(f, fullfile(CC_path, filename))

%%


%% % % % % % % % % % - Voltage Clamp  % % % % % % % % % %
clear all
close all

VC_path = '/Users/burnettl/Documents/ISTA/PLOS_BIOL_REVISIONS/PATCH/PLOTS/VC_10Hz';

animal = 'M8';
cell_id = 'C4';
file_to_open = 2; 

%%

% Change into folder with VC data: 

files = dir('*.abf');
nfiles = length(files);

fn = files(file_to_open).name;
[d,si,h]=abfload(fn);
d = squeeze(d);
n = numel(d(1,:));

if n<10
    n_reps = n;
else 
    n_reps = 10;
end 

%  #2 - Plot all reps: 
close 
offset = 0;
figure;
for k = 1:n_reps
    plot(d(:, k)-offset, 'k')
    hold on
    offset = offset+250;
end 
title([animal, '-', cell_id])
xlim([0, 100000])
box off
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.XAxis.Visible = 'off';
ylabel('Current (pA)');
f = gcf; 
f.Position = [601   310   377   608];

filename  = ['VC_10Hz_Reps_', animal, '_', cell_id, '.svg'];
saveas(f, fullfile(VC_path, filename))


%% # 3 - Plot average 

close
figure
aa = d(:, 1:n_reps);
a = nanmean(aa');
plot(a, 'k')
xlim([0, 100000])
ylim([-250, 0])
box off
title([animal, '-', cell_id])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.FontSize = 14;
ax.XAxis.Visible = 'off';
f2 = gcf; 
f2.Position = [405   688   460   173];

filename  = ['VC_10Hz_AVERAGE_', animal, '_', cell_id, '.svg'];
saveas(f2, fullfile(VC_path, filename))



%%



















