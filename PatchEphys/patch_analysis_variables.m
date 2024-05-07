%% PATCH - making box plot with data points for PATCH data. 
% 17/03/22
% Burnett

% Take data from Peter's analysis: 
% Want to make plots for 
% - Input resistance
% - Tau 
% - Membrane capacitance
% - Resting Membrane Potential

% Which animals are WT/ HET

% COHORT 1 -  all PAG
% M1 - 1185 - W

% M2 - 1186 - H
% M3 - 1348 - W
% M4 - 1347 - H
% M5 - 1415 - H
% M6 - 1414 - W
% M7 - 1416 - W
% M8 - 1417 - H

% COHORT 2 - PAG/SC/DTX
% M1 - 1992 - H - C1-C3 = PAG, C4 - C6 = SC
% M2 - 1990 - W - 1-4 = PAG, 5-7 = SC
% M3 - 1938 - H - all PAG, TDX...
% M4 - 1942 - H - 1-4 = PAG, 5-6 = SC
% M5 - 1943 - W - 1-3 PAG, 4-6 = SC
% M6 - 1936 - W - 1-3 = PAG, 4-6 = SC
% M7 - 1937 - H - 1-3 = PAG, 4-6 = SC
% M8 - 2139 - H - 4 PAG, DTX
% M9 - 4666 - H - 5 PAG, DTX

%% Manually entered the data from the excel spreadsheet Peter generated. Only had the data for PAG cell.s 
patch_variables = zeros(50, 12);
patch_variables_table = array2table(patch_variables, 'VariableNames', {'Animal', 'Geno', 'Cell', 'PAG', 'Rin', 'T', 'Cm', 'RMP', 'SpFr', 'SpAmp', 'SpRise', 'SpDecay'});

save('220327_Setd5_PatchVariables_TYPE.mat', 'patch_variables', 'patch_variables_table');


%% Make plots about variables from patch data:
allWT = find(patch_variables_table.Geno == 1);
allHET = find(patch_variables_table.Geno ==2);

% Column val
i = 4; 

% BOXPLOT WITH SCATTER POINTS

figure
wtvals = patch_variables_table{allWT, 4+i};
hetvals = patch_variables_table{allHET, 4+i};

var = patch_variables_table{:, 4+i};
grp = patch_variables_table{:,2};

b = boxplot(var, grp, 'Color', 'k');
set(b , 'LineWidth', 1.25)
hold on
scatter(ones(1,19), wtvals, 150, 'o','jitter', 'on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0.5 0.5 0.5]);
scatter(ones(1,29)*2, hetvals, 150, 'ro', 'jitter', 'on', 'jitterAmount', 0.1);

box off
ylim([-75 -45])
xlim([0.5 2.5])
ax = gca;
% ax.XAxis.Visible = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.TickLength = [0.03 0.03];
ax.FontSize = 18;
% ax.YTick = [0,1,2,3,4];
f = gcf;
f.Position = [704   207   355   572]; 



ax.YTick = [-75, -65, -55, -45];

    



%%

figure
% wtvals = patch_variables_table{allWT, 4+i};
% hetvals = patch_variables_table{allHET, 4+i};

% wtvals = [];
% hetvals = [];

mWT = nanmean(wtvals);
mHET = nanmean(hetvals);

figure
scatter(ones(1,15), wtvals, 150, 'o','jitter', 'on', 'jitterAmount', 0.25, 'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on 
scatter(ones(1,15)*2, hetvals, 150, 'o', 'jitter', 'on', 'jitterAmount', 0.25,  'MarkerEdgeColor', [1 0.5 0.5]);
% Errorbar
errorbar(1, mWT, nanstd(wtvals)/sqrt(15), 'k', 'LineWidth', 1.5);
errorbar(2, mHET, nanstd(hetvals)/sqrt(15), 'r', 'LineWidth', 1.5);
% bar 
plot([0.75 1.25], [mWT mWT], 'k', 'LineWidth', 3.5)
plot([1.75 2.25], [mHET mHET], 'r', 'LineWidth', 3.5)

box off
xlim([0.25 2.75])
ax = gca;
ax.XAxis.Visible = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.TickLength = [0.025 0.025];
ax.FontSize = 18;
% ax.YTick = [0,1,2,3,4];

f = gcf;
f.Position = [680   634   241   464]; %[680   670   234   428]; 

ylim([0 140])

% [p,h] = ranksum(wtvals, hetvals)
[p, h] = ttest2(wtvals, hetvals)
nanmean(wtvals)
nanmean(hetvals)

[H, pValue] = swtest(wtvals)
[H, pValue] = swtest(hetvals)

% save('220401_SC_IntrinsicProp_RMP_N6N7.mat', 'wtvals', 'hetvals');






%% Need to analyse these properties for the SC cells. 

% Need to current current input - spikes output data.
% Need the HYPERPOLARISING sweeps (-10pA, -20A)

% RMP = average MP prior to injection. 
% Input R = slope of linear fit to hyperpolarising current voltage
% relationship
% T  = exp fit of mem pot change at -10pA or -20pA
% cm = t/Rin






%% Make plots about variables from patch data - WITH CELL TYPE! 

% allWT = find(patch_variables_table.Geno == 1 & cell2mat(patch_variables_table.data) == 1);
% allHET = find(patch_variables_table.Geno ==2 & cell2mat(patch_variables_table.data) == 1);

% GROUP - GENO + TYPE. 
allWTG = find(cell2mat(patch_variables_table.data1) == 1);
allHETG = find(cell2mat(patch_variables_table.data1) == 3);

allWTB = find(cell2mat(patch_variables_table.data1) == 2);
allHETB = find(cell2mat(patch_variables_table.data1) == 4);

%%

% Column val
i = 8; 

% BOXPLOT WITH SCATTER POINTS
close

figure
wtvalsG = patch_variables_table{allWTG, 6+i};
hetvalsG = patch_variables_table{allHETG, 6+i};
wtvalsB = patch_variables_table{allWTB, 6+i};
hetvalsB = patch_variables_table{allHETB, 6+i};


var = patch_variables_table{:, 6+i};
grp = patch_variables_table{:,6};

b = boxplot(var, grp, 'Color', 'k');
set(b , 'LineWidth', 1.25)
hold on
scatter(ones(1,numel(allWTG)), wtvalsG, 150, 'o','jitter', 'on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0.5 0.5 0.5]);
scatter(ones(1,numel(allHETG))*3, hetvalsG, 150, 'ro', 'jitter', 'on', 'jitterAmount', 0.1);
scatter(ones(1,numel(allWTB))*2, wtvalsB, 150, 'o','jitter', 'on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0.5 0.5 1]);
scatter(ones(1,numel(allHETB))*4, hetvalsB, 150, 'mo', 'jitter', 'on', 'jitterAmount', 0.1);

box off
ax = gca;
% ax.XAxis.Visible = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.TickLength = [0.03 0.03];
ax.FontSize = 18;
% ax.YTick = [0,1,2,3,4];
f = gcf;
f.Position = [704   394   326   385]; 


%%
[p,h] = ranksum(wtvals, hetvals)
[h, p] = ttest2(wtvals, hetvals)





















