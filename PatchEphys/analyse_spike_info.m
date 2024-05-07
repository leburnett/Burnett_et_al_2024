% Script for analysis of spike info 
% Burnett - 09/08/2022

%% Analysis of data from table: 'spike_table' generated in 'make_spike_info_table.m'

% ALL SPIKES FROM ALL CELLS FROM ALL SWEEPS
allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0); 
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0); 

allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1);

%% REMOVE SPIKES THAT HAPPENED WTIH 0pA!! 

all_03 = find(spike_table.Sweep == 3);
spike_table(all_03, :) = []; 

%% ALL SPIKES FROM ALL CELLS FROM SWEEPS 4:23 (10pA - 200pA) - - - - -DEFAULT 
allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24);  % sweep 4 = 10pA - sweep 23 = 200pA
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24); 

allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24);


%% SPIKES 1-25 FROM ALL CELLS FROM SWEEPS 4:23 (10pA - 200pA)
allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN>2 & spike_table.SpikeN<30);  % sweep 4 = 10pA - sweep 23 = 200pA
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN>2 & spike_table.SpikeN<30); 

allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN>2 & spike_table.SpikeN<30);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN>2 & spike_table.SpikeN<30);

%% All spikes from 180pA sweep 
allWT = find(spike_table.Geno == 1 & spike_table.Sweep == 13); 
allHET = find(spike_table.Geno == 2 & spike_table.Sweep == 13); 



%% ALL SPIKES - ALL cells - WITHIN one sweep. 
% Looking at the average 1st spike/ average of all spikes in one sweep. 
% 
% sw = 21;
% 
% allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep ==sw); %
% allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep ==sw); %& spike_table.Sweep ==sw
% 
% allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep ==sw);
% allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep ==sw);

%% TEST - check where on plot dvdt>10mv/ms is.

% vWT = spike_table.V(1, st:nd);
% dvWT = spike_table.dvdt(1,  st:nd);
% 
% plot(vWT, dvWT, 'k'); hold on 
% 
% % dv2WT = diff(dvWT);
% % figure; plot(dv2WT, 'k')
% % find(dvWT>10)
% plot(vWT(463), dvWT(463), 'b.', 'MarkerSize', 10)

%% 1 - DATA - Average AP - WT/ HET
st = 500;
nd = 1500;

vWT = spike_table.V(allWT, st:nd);
vHET = spike_table.V(allHET,  st:nd);

dvWT = spike_table.dvdt(allWT,  st:nd);
dvHET = spike_table.dvdt(allHET,  st:nd);

vWTD = spike_table.V(allWTDTX,  st:nd);
vHETD = spike_table.V(allHETDTX,  st:nd);
% 
dvWTD = spike_table.dvdt(allWTDTX,  st:nd);
dvHETD = spike_table.dvdt(allHETDTX,  st:nd);

%% PLOT - WT versus HET - NO DTX - JUST THE MEAN 
col = [1 0.55 0.2];


if sw == 21
    st = 750;
    nd = 1250;
elseif sw == 1
    st = 1;
    nd = 2000;
else 
    st = 700;
    nd = 1300; 
end 

% AP average
figure; plot(st:1:nd, nanmean(vWT(:, st:nd)), 'k', 'LineWidth', 1.2)
hold on; plot(st:1:nd, nanmean(vHET(:, st:nd)), 'Color', col, 'LineWidth', 1.2)
% hold on; plot(st:1:nd, nanmean(vWTD(:, st:nd)),  'Color', [0.5 0.7 0.9], 'LineWidth', 1.2)
% hold on; plot(st:1:nd, nanmean(vHETD(:, st:nd)), 'Color', [1 0.6 0.7], 'LineWidth', 1.2)
% xlim([0 2000])
xlim([st nd])
ylim([-70 50])
f = gcf;
f.Position = [440   535   315   263];
box off
ax = gca; 
ax.TickDir = 'out'; 
ax.LineWidth = 1;
ax.TickLength = [0.02 0.02];
grid on

% Phase plane plot
figure; plot(nanmean(vWT(:, st:nd)), nanmean(dvWT(:, st:nd)), 'k', 'LineWidth', 1.2);
hold on; plot(nanmean(vHET(:, st:nd)), nanmean(dvHET(:, st:nd)), 'Color', col, 'LineWidth', 1.2);
% hold on; plot(nanmean(vWTD(:, st:nd)), nanmean(dvWTD(:, st:nd)), 'Color', [0.5 0.7 0.9], 'LineWidth', 1.2)
% hold on; plot(nanmean(vHETD(:, st:nd)), nanmean(dvHETD(:, st:nd)), 'Color', [1 0.6 0.7], 'LineWidth', 1.2)
xlim([-70 45])
ylim([-150 170])
% xlim([-70 45])
% ylim([-200 225])
f = gcf;
f.Position = [440   535   315   263];
box off
ax = gca; 
ax.TickDir = 'out'; 
ax.LineWidth = 1;
ax.TickLength = [0.02 0.02];
grid on

%% MEAN + SEM - average action potential - (t x v)

% WT - HET - overlaid
st = 1;
nd = 1001;
% x = (st:1:nd);
x = (1:1:nd-st+1);

nWT = 19;%numel(vWT(:,1)); 
nHET = 21; % numel(vHET(:,1));

mean_WT = nanmean(vWT); 
mean_HET = nanmean(vHET); 


semWT = nanstd(vWT)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

semHET = nanstd(vHET)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT', 'k-', 'LineWidth', 1.3)
plot(x, y3, 'w-')
hold on 
plot(x, y4, 'w-')
patch([x fliplr(x)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET', '-','Color', col, 'LineWidth', 1.3)


box off
ax = gca; 
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];
grid off
% xlim([250 750]) 
% xlim([800 1200]) 
% xlim([150 500])
f = gcf;
f.Position = [440   535   315   263];


xlim([400 700])
ylim([-70 50])
% WTD - HETD

nWT = 12; %numel(vWTD(:,1)); 
nHET = 13; %numel(vHETD(:,1));
% 
mean_WT = nanmean(vWTD); 
mean_HET = nanmean(vHETD); 
% 
% 
semWT = nanstd(vWTD)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
% 
semHET = nanstd(vHETD)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;
% 
% figure
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT', 'b-', 'LineWidth', 1.3)
plot(x, y3, 'w-')
hold on 
plot(x, y4, 'w-')
patch([x fliplr(x)], [y3 fliplr(y4)], 'm' , 'FaceAlpha', 0.2, 'EdgeColor', 'none') %[0.5 0.15 1]
plot(mean_HET', '-','Color', 'm', 'LineWidth', 1.3)


%% MEAN / SEM - phase plane plots - (dvdt x v) 
% Need to find the average per cell - otherwise cells with many spikes with have a disproportionate effect on the final shape. 

st = 250;
nd = 750;

x = nanmean(vWT(:, st:nd));
x2 = nanmean(vHET(:, st:nd));

nWT = 16; %numel(dvWT(:,1)); 
nHET = 20; %numel(dvHET(:,1));

mean_WT = nanmean(dvWT(:, st:nd)); 
mean_HET = nanmean(dvHET(:, st:nd)); 

semWT = nanstd(dvWT(:, st:nd))/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

semHET = nanstd(dvHET(:, st:nd))/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

% WT - HET
figure
% WT
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x, mean_WT', 'k-', 'LineWidth', 1.3)
% HET
plot(x2, y3, 'w-')
hold on 
plot(x2, y4, 'w-')
patch([x2 fliplr(x2)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x2, mean_HET', '-','Color', col, 'LineWidth', 1.3)

% % % % % % % % % % % % % % % 

% DTX - WT - HET
x = nanmean(vWTD(:, :));
x2 = nanmean(vHETD(:, :));
% 
nWT = 12; %numel(dvWTD(:,1)); 
nHET = 13;% numel(dvHETD(:,1));
% 
mean_WT = nanmean(dvWT); 
mean_HET = nanmean(dvHET); 
% 
semWT = nanstd(dvWTD)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
% 
semHET = nanstd(dvHETD)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

% 
% 
% % WT
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x, mean_WT', 'b-', 'LineWidth', 1.3)
% % HET
plot(x2, y3, 'w-')
hold on 
plot(x2, y4, 'w-')
patch([x2 fliplr(x2)], [y3 fliplr(y4)],  'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x2, mean_HET', '-','Color', 'm', 'LineWidth', 1.3)


box off
ax = gca; 
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];
grid off
f = gcf;
f.Position = [440   535   315   263];

xlim([-70 50])
ylim([-150 175])

%% MEAN SEM - dvdt


% WT - HET - overlaid
st = 1;
nd = 501;
x = (st:1:nd);

nWT = 19 ;%numel(vWT(:,1)); 
nHET = 21; % numel(vHET(:,1));

mean_WT = nanmean(dvWT); 
mean_HET = nanmean(dvHET); 


semWT = nanstd(dvWT)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

semHET = nanstd(dvHET)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT', 'k-', 'LineWidth', 1.3)
plot(x, y3, 'w-')
hold on 
plot(x, y4, 'w-')
patch([x fliplr(x)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET', '-','Color', col, 'LineWidth', 1.3)

box off
ax = gca; 
ax.TickDir = 'out'; 
ax.LineWidth = 1;
ax.TickLength = [0.02 0.02];
grid on
xlim([100 400])
f = gcf;
f.Position = [440   535   315   263];

% DTX

nWT = 12; %numel(vWTD(:,1)); 
nHET = 13; %numel(vHETD(:,1));

mean_WT = nanmean(dvWTD); 
mean_HET = nanmean(dvHETD); 


semWT = nanstd(dvWTD)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

semHET = nanstd(dvHETD)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;


plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT', 'b-', 'LineWidth', 1.3)
plot(x, y3, 'w-')
hold on 
plot(x, y4, 'w-')
patch([x fliplr(x)], [y3 fliplr(y4)],  'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET', '-','Color', 'm', 'LineWidth', 1.3)





%% STATS %% %% %% %% 
wtvals = [];
hetvals = [];

for kk = 1:numel(cellids)

allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.UniCellID ==kk & spike_table.Sweep == 19); 
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.UniCellID ==kk & spike_table.Sweep == 19); 

allWT2 = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.UniCellID ==kk & spike_table.Sweep == 4); 
allHET2 = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.UniCellID ==kk & spike_table.Sweep == 4); 

if ~isempty(allWT)
    wtv = nanmean(spike_table.MaxAmp(allWT));
    wtv2 = nanmean(spike_table.MaxAmp(allWT2));
    ww = wtv2-wtv;
    wtvals = vertcat(wtvals, ww);
elseif ~isempty(allHET)

    htv = nanmean(spike_table.MaxAmp(allHET));
    htv2 = nanmean(spike_table.MaxAmp(allHET2));
    hh = htv2-htv;
    hetvals = vertcat(hetvals, hh);
end

end 

%%  
% wtvals = spike_table.MaxAmp(allWT);
% hetvals = spike_table.MaxAmp(allHET);
nanmean(wtvals)
nanmean(hetvals)
numel(wtvals)
numel(hetvals)
[h, p] = kstest2(wtvals, hetvals)
[p,h] = ranksum(wtvals, hetvals)


[cellids,ac,bc] = unique(spike_table(:, [1,2,3,4])); % Recs with / wihtout DTX are together! 
spike_table.UniCellID = bc;

%%
% REMOVE NANS
aaa = find(isnan(spike_table.APThresh));
spike_table(aaa,:) = [];

save('220811_Spike_Table_wHyperpol.mat', 'spike_table');


%% PLOT DTX
%% AP shape
% WT-DTx vs WT+DTX
figure; plot(1:1:2001, nanmean(vWT), 'k')
hold on; plot(1:1:2001, nanmean(vWTD), 'b')
xlim([0 2000])
ylim([-70 50])
f = gcf;
f.Position = [440   535   315   263];

% HET-DTx vs HET+DTX
figure; plot(1:1:2001, nanmean(vHET), 'r')
hold on; plot(1:1:2001, nanmean(vHETD), 'm')
xlim([0 2000])
ylim([-70 50])
f = gcf;
f.Position = [440   535   315   263];

% WT+DTx versus HET+DTX
figure; plot(1:1:2001, nanmean(vWTD), 'b')
hold on; plot(1:1:2001, nanmean(vHETD), 'm')
xlim([0 2000])
ylim([-70 50])
f = gcf;
f.Position = [440   535   315   263];

%% Phase plane plot
% WT-DTx vs WT+DTX
figure; plot(nanmean(vWT(:, st:nd)), nanmean(dvWT(:, st:nd)), 'k');
hold on; plot(nanmean(vWTD(:, st:nd)), nanmean(dvWTD(:, st:nd)), 'b');
xlim([-70 40])
ylim([-200 200])
f = gcf;
f.Position = [440   535   315   263];


% HET-DTx vs HET+DTX
figure; plot(nanmean(vHET(:, st:nd)), nanmean(dvHET(:, st:nd)), 'r');
hold on; plot(nanmean(vHETD(:, st:nd)), nanmean(dvHETD(:, st:nd)), 'm');
xlim([-70 50])
ylim([-200 200])
f = gcf;
f.Position = [440   535   315   263];


% WT+DTx versus HET+DTX
figure; plot(nanmean(vWTD(:, st:nd)), nanmean(dvWTD(:, st:nd)), 'b');
hold on; plot(nanmean(vHETD(:, st:nd)), nanmean(dvHETD(:, st:nd)), 'm');
xlim([-70 50])
ylim([-200 200])
f = gcf;
f.Position = [440   535   315   263];

%% ALL 4 overlaid

figure; plot(nanmean(vWT(:, st:nd)), nanmean(dvWT(:, st:nd)), 'k', 'LineWidth', 1.2); hold on;
% hold on; plot(nanmean(vWTD(:, st:nd)), nanmean(dvWTD(:, st:nd)), 'Color', [0.5 0.7 0.9], 'LineWidth', 1.2);
plot(nanmean(vHET(:, st:nd)), nanmean(dvHET(:, st:nd)), 'r', 'LineWidth', 1.2);
plot(nanmean(vHETD(:, st:nd)), nanmean(dvHETD(:, st:nd)), 'Color', [1 0.7 0.8], 'LineWidth', 1.2);
xlim([-70 40])
ylim([-200 200])
f = gcf;
f.Position = [440   535   315   263];



%% PHASE PLANE PLOT - change in the average pp ACROSS sweeps (average of all spikes in that sweep). 

figure
a = 0.8;
b = 1;

% Lims for plots
ymax = 225;
ymin = -200;
xmin = -70;
xmax = 45;

% Which spike number to check?
%  spn = 51;

% Start and end of data
st = 800;
nd = 1200;

for i = [5,7,9,11,13,15,17,19] % 10/20, 60, 120, 180pA steps [4,6,8,10,12,14,16,18,20] % 19=180pA - 23=200pA
sw = i; 
    
% if sw > 8
%     st = 750;
%     nd = 1250;
% else
%     st = 500;
%     nd = 1500;
% end 

allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep ==sw); % & spike_table.SpikeN ==spn
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep ==sw);

allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep ==sw);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep ==sw);

vWT = spike_table.V(allWT, :);
vHET = spike_table.V(allHET, :);

dvWT = spike_table.dvdt(allWT, :);
dvHET = spike_table.dvdt(allHET, :);

vWTD = spike_table.V(allWTDTX, :);
vHETD = spike_table.V(allHETDTX, :);

dvWTD = spike_table.dvdt(allWTDTX, :);
dvHETD = spike_table.dvdt(allHETDTX, :);


subplot(2,1,1); 
% plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
plot(nanmean(vWT(:, st:nd)), nanmean(dvWT(:, st:nd)), 'Color', [a a a], 'LineWidth', 1.2); hold on
xlim([xmin xmax])
ylim([ymin ymax])
grid on

% % % % % subplot(2,2,3); 
% % % % % % plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% % % % % % plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% % % % % % plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % % % plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % % % plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % % % plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % % plot(nanmean(vWTD(:, st:nd)), nanmean(dvWTD(:, st:nd)), 'Color', [a a b], 'LineWidth', 1.2); hold on
% % % % % xlim([xmin xmax])
% % % % % ylim([ymin ymax])
% % % % % grid on

subplot(2,1,2); 
% plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
plot(nanmean(vHET(:, st:nd)), nanmean(dvHET(:, st:nd)), 'Color', [b a a], 'LineWidth', 1.2); hold on 
xlim([xmin xmax])
ylim([ymin ymax])
grid on

% % % % subplot(2,2,4); 
% % % % % plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% % % % % plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% % % % % plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % % plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % % plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % % plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % % % plot(nanmean(vHETD(:, st:nd)), nanmean(dvHETD(:, st:nd)), 'Color', [b a b], 'LineWidth', 1.2); hold on 
% % % % xlim([xmin xmax])
% % % % ylim([ymin ymax])
% % % % grid on

a = a-0.08;
b = b-0.025;
end 

f = gcf;
f.Position = [172   136   708   603];

%% AP SHAPE - change in the average AP shape ACROSS sweeps (average of all spikes in that sweep). 

figure
a = 0.8;
b = 1;

ymax = 45;
ymin = -70;

st = 801; 
nd = 1200;
xvals = st:1:nd;

xmin = st;
xmax = nd;

spn = 51; 

for i =  [5,7,9,11,13,15,17,19,21,23] % 10/20, 60, 120, 180pA steps [4,6,8,10,12,14,16,18,20]
sw = i; 
    
allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep ==sw & spike_table.SpikeN <spn);
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep ==sw & spike_table.SpikeN <spn);

% allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep ==sw & spike_table.SpikeN < spn);
% allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep ==sw & spike_table.SpikeN < spn);

vWT = spike_table.V(allWT, :);
vHET = spike_table.V(allHET, :);

% vWTD = spike_table.V(allWTDTX, :);
% vHETD = spike_table.V(allHETDTX, :);


subplot(1,2,1); 
plot(xvals, nanmean(vWT(:, st:nd)), 'Color', [a a a], 'LineWidth', 1.2); hold on;
xlim([xmin xmax])
ylim([ymin ymax])
grid on

% subplot(2,2,3); 
% plot(xvals, nanmean(vWTD(:, st:nd)), 'Color', [a a b], 'LineWidth', 1.2);hold on;
% xlim([xmin xmax])
% ylim([ymin ymax])
% grid on

subplot(1,2,2); 
plot(xvals, nanmean(vHET(:, st:nd)), 'Color', [b a a], 'LineWidth', 1.2);hold on;
xlim([xmin xmax])
ylim([ymin ymax])
grid on

% subplot(2,2,4); 
% plot(xvals, nanmean(vHETD(:, st:nd)), 'Color', [b a b], 'LineWidth', 1.2);hold on;
% xlim([xmin xmax])
% ylim([ymin ymax])
% grid on

a = a-0.08;
b = b-0.025;
end 

f = gcf;
f.Position = [172   136   708   603];










%% WITHIN SWEEPS!!! 







%% PHASE PLANE - change in phase plane across spikes WITHIN the same sweep. 

figure
a = 0.8;
b = 1;

% Lims for plots
ymax = 220;
ymin = -200;
xmin = -70;
xmax = 45;

% Which SWEEP to choose
sw = 15;
% 4 = 10pA
% 5 = 20pA
% 9 = 60pA
% 13 = 100pA
% 15 = 120pA
% 18 = 150pA
% 21 = 180pA
% 23 = 200pA

% Start and end of data
if sw == 5
    st = 600;
    nd = 1400;
    maxsp = 7;
elseif sw==9
    st = 750;
    nd = 1250;
    maxsp = 20;
elseif sw == 15
    st = 750;
    nd = 1250;
    maxsp = 40;
elseif sw == 21 || sw == 23
    st = 800;
    nd = 1200;
    maxsp = 40;
end

if maxsp == 20
    am = 0.04;
    bm = 0.009;  
elseif maxsp == 7
    am = 0.09;
    bm = 0.025; 
elseif maxsp == 40
    am = 0.018;
    bm = 0.005; 
elseif maxsp == 50 
    am = 0.01;
    bm = 0.001; 
end 

for i = 1:1:maxsp %[1,2,3,4,5,8,10,12,15,20,25,30,35,40] %[5,7,9,11,13,15,17,19] % 10/20, 60, 120, 180pA steps [4,6,8,10,12,14,16,18,20]
spn = i; 

allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep ==sw & spike_table.SpikeN ==spn); % & spike_table.SpikeN ==spn
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep ==sw & spike_table.SpikeN ==spn);

allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep ==sw & spike_table.SpikeN ==spn);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep ==sw & spike_table.SpikeN ==spn);

vWT = spike_table.V(allWT, :);
vHET = spike_table.V(allHET, :);

dvWT = spike_table.dvdt(allWT, :);
dvHET = spike_table.dvdt(allHET, :);

vWTD = spike_table.V(allWTDTX, :);
vHETD = spike_table.V(allHETDTX, :);

dvWTD = spike_table.dvdt(allWTDTX, :);
dvHETD = spike_table.dvdt(allHETDTX, :);

subplot(1,2,1); 
% plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
plot(nanmean(vWT(:, st:nd)), nanmean(dvWT(:, st:nd)), 'Color', [a a a], 'LineWidth', 1.2); hold on
ylim([ymin ymax])
xlim([xmin xmax])
grid on

% subplot(2,2,3); 
% % plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% % plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% % plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot(nanmean(vWTD(:, st:nd)), nanmean(dvWTD(:, st:nd)), 'Color', [a a b], 'LineWidth', 1.2); hold on
% xlim([xmin xmax])
% grid on
% ylim([ymin ymax])

subplot(1,2,2); 
% plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
plot(nanmean(vHET(:, st:nd)), nanmean(dvHET(:, st:nd)), 'Color', [b a a], 'LineWidth', 1.2); hold on
xlim([xmin xmax])
grid on
ylim([ymin ymax])
% 
% subplot(2,2,4); 
% % plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% % plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% % plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot(nanmean(vHETD(:, st:nd)), nanmean(dvHETD(:, st:nd)), 'Color', [b a b], 'LineWidth', 1.2); hold on
% xlim([xmin xmax])
% grid on
% ylim([ymin ymax])

a = a-am;
b = b-bm;
end 

f = gcf;
f.Position = [172   136   708   603];

%% AP shape - change over subsequent spikes - subplot with all 4 conditions. 

figure
a = 0.8;
b = 1;

ymax = 45;
ymin = -70;

st = 801; 
nd = 1200;
xvals = st:1:nd;

xmin = st;
xmax = nd;

% Which SWEEP to choose
sw = 23;
% 4 = 10pA
% 5 = 20pA
% 9 = 60pA
% 13 = 100pA
% 15 = 120pA
% 18 = 150pA
% 21 = 180pA
% 23 = 200pA

% Start and end of data
if sw == 5
    maxsp = 7;
elseif sw==9
    maxsp = 20;
elseif sw == 15 || sw == 21 || sw == 23
    maxsp = 40;
end

if maxsp == 20
    am = 0.04;
    bm = 0.009;  
elseif maxsp == 7
    am = 0.09;
    bm = 0.025; 
elseif maxsp == 40
    am = 0.018;
    bm = 0.005; 
elseif maxsp == 50 
    am = 0.01;
    bm = 0.001; 
end 


for i =  1:1:maxsp %[1,2,3,4,5,8,10,12,15,20,25,30,35,40] 
spn = i; 
    
allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep ==sw & spike_table.SpikeN == spn);
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep ==sw & spike_table.SpikeN == spn);

allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep ==sw & spike_table.SpikeN == spn);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep ==sw & spike_table.SpikeN == spn);

vWT = spike_table.V(allWT, :);
vHET = spike_table.V(allHET, :);

vWTD = spike_table.V(allWTDTX, :);
vHETD = spike_table.V(allHETDTX, :);


subplot(2,2,1); 
plot(xvals, nanmean(vWT(:, st:nd)), 'Color', [a a a], 'LineWidth', 1.2); hold on;
xlim([xmin xmax])
ylim([ymin ymax])
grid on

subplot(2,2,3); 
plot(xvals, nanmean(vWTD(:, st:nd)), 'Color', [a a b], 'LineWidth', 1.2);hold on;
xlim([xmin xmax])
ylim([ymin ymax])
grid on

subplot(2,2,2); 
plot(xvals, nanmean(vHET(:, st:nd)), 'Color', [b a a], 'LineWidth', 1.2);hold on;
xlim([xmin xmax])
ylim([ymin ymax])
grid on

subplot(2,2,4); 
plot(xvals, nanmean(vHETD(:, st:nd)), 'Color', [b a b], 'LineWidth', 1.2);hold on;
xlim([xmin xmax])
ylim([ymin ymax])
grid on

a = a-am; %a-0.05;
b = b-bm; %b-0.01;
end 

f = gcf;
f.Position = [172   136   708   603];





%% PHASE PLANE PLOT - SUBPLOT - WT/HET overlaid : WTD/HTD overlaid - change in the average pp ACROSS sweeps (average of all spikes in that sweep). 

figure
a = 0.8;
b = 1;

% Start and end of data
st = 800;
nd = 1200;

% Lims for plots - PP 
% ymax = 220;
% ymin = -200;
% xmin = -70;
% xmax = 45;

% AP SHAPE 
ymax = 45;
ymin = -70;
xmin = st;
xmax = nd;

% Which spike number to check?
spn = 1;

for i = [5,7,9,11,13,15,17,19] % 10/20, 60, 120, 180pA steps [4,6,8,10,12,14,16,18,20]
sw = i; 
    
% if sw > 8
%     st = 750;
%     nd = 1250;
% else
%     st = 500;
%     nd = 1500;
% end 

allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep ==sw ); % & spike_table.SpikeN ==spn
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep ==sw);

% allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep ==sw);
% allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep ==sw);

vWT = spike_table.V(allWT, :);
vHET = spike_table.V(allHET, :);

dvWT = spike_table.dvdt(allWT, :);
dvHET = spike_table.dvdt(allHET, :);

% vWTD = spike_table.V(allWTDTX, :);
% vHETD = spike_table.V(allHETDTX, :);
% 
% dvWTD = spike_table.dvdt(allWTDTX, :);
% dvHETD = spike_table.dvdt(allHETDTX, :);


s1 = subplot(1,2,1); 
% plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)

% PHASE P 
% plot(nanmean(vWT(:, st:nd)), nanmean(dvWT(:, st:nd)), 'Color', [a a a], 'LineWidth', 1.2); hold on
% plot(nanmean(vHET(:, st:nd)), nanmean(dvHET(:, st:nd)), 'Color', [b a a], 'LineWidth', 1.2); 
% AP SHAPE
plot(st:1:nd, nanmean(vWT(:, st:nd)),'Color', [a a a], 'LineWidth', 1.2); hold on
plot(st:1:nd, nanmean(vHET(:, st:nd)),'Color', [b a a], 'LineWidth', 1.2); 

xlim([xmin xmax])
ylim([ymin ymax])
grid on

% s2 = subplot(1,2,2); 
% % plot([20 20], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  hold on;
% % plot([30 30], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);  
% % plot([-60 -60], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([-40 -40], [ymin ymax], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([xmin xmax], [-100 -100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% % plot([xmin xmax], [100 100], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
% 
% %PHASE PLANE
% plot(nanmean(vWTD(:, st:nd)), nanmean(dvWTD(:, st:nd)), 'Color', [a a b], 'LineWidth', 1.2); hold on
% plot(nanmean(vHETD(:, st:nd)), nanmean(dvHETD(:, st:nd)), 'Color', [b a b], 'LineWidth', 1.2); 

% AP SHAPE
% plot(st:1:nd, nanmean(vWTD(:, st:nd)), 'Color', [a a b], 'LineWidth', 1.2); hold on
% plot(st:1:nd, nanmean(vHETD(:, st:nd)), 'Color', [b a b], 'LineWidth', 1.2); 

% xlim([xmin xmax])
% ylim([ymin ymax])
% grid on

% linkaxes([s1, s2], 'xy')

a = a-0.1;
b = b-0.025;
% a = a-0.1;
% b = b-0.05;
end 

f = gcf;
% PP
f.Position = [109   395   856   284]; %[109         316        1146         363]; 
% AP Shape
% f.Position = [109   335   815   344]; 







%% STATS AND MORE





%% STATS 

% Check what columns are in the table. 
% spike_table(1,:)

% % 10 - MaxAmp
% 14 - AP thresh
% 15 = Rise T
% 16 - whp
% 17 - decay T
% 19 - min amp

%% ALL SPIKES FROM ALL CELLS FROM SWEEPS 4:23 (10pA - 200pA) - - - - -DEFAULT 
allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24);  % sweep 4 = 10pA - sweep 23 = 200pA
allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24); 

allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24);


%% SPIKES 1-25 FROM ALL CELLS FROM SWEEPS 4:23 (10pA - 200pA)
% allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN==1);  % sweep 4 = 10pA - sweep 23 = 200pA
% allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN==1); 
% 
% allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN==1);
% allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep>3 & spike_table.Sweep<24 & spike_table.SpikeN==1);

%% ALL SPIKES - ALL cells - WITHIN one sweep. 
% Looking at the average 1st spike/ average of all spikes in one sweep. 
% 
% sw = 21;
% 
% allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep ==sw); %
% allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep ==sw); %& spike_table.Sweep ==sw
% 
% allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep ==sw);
% allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep ==sw);

%% Variable values - WT/ HET

stats_values = zeros(6, 20);
roww = 1;

for k = [10, 14, 15, 16, 17, 21, 22, 23] % Set k as column in spike table to be assessed.
    
    vWT = spike_table{allWT, k};
    vHET = spike_table{allHET, k};
    vWTD = spike_table{allWTDTX, k};
    vHETD = spike_table{allHETDTX, k};
    
    % Number of spikes in total
    stats_values(roww, 1) = numel(vWT); % Num WT
    stats_values(roww, 2) = numel(vHET);
    stats_values(roww, 3) = numel(vWTD);
    stats_values(roww, 4) = numel(vHETD);
    
    % MEAN VALUES
    stats_values(roww, 5) = nanmean(vWT); 
    stats_values(roww, 6) = nanmean(vHET);
    stats_values(roww, 7) = nanmean(vWTD);
    stats_values(roww, 8) = nanmean(vHETD);
    
    % RANGE 
    stats_values(roww, 9) = range(vWT); 
    stats_values(roww, 10) = range(vHET);
    stats_values(roww, 11) = range(vWTD);
    stats_values(roww, 12) = range(vHETD);
    
    % Wilcoxon rank sum - WT vs het
    [p, h] = ranksum(vWT, vHET);
    stats_values(roww, 13) = p; 
    
    % Wilcoxon rank sum - WT vs WTD
    [p, h] = ranksum(vWT, vWTD);
    stats_values(roww, 14) = p; 
    
    % Wilcoxon rank sum - HET vs HETD
    [p, h] = ranksum(vHETD, vHET);
    stats_values(roww, 15) = p; 
    
    % Wilcoxon rank sum - WT vs HETD
    [p, h] = ranksum(vWT, vHETD);
    stats_values(roww, 16) = p;
    
    % Kolmogorov Smirnoff - WT vs het
    [h, p] = kstest2(vWT, vHET);
    stats_values(roww, 17) = p; 
    
    % Kolmogorov Smirnoff - WT vs WTD
    [h, p] = kstest2(vWT, vWTD);
    stats_values(roww, 18) = p; 
    
    % Kolmogorov Smirnoff  - HET vs HETD
    [h, p] = kstest2(vHETD, vHET);
    stats_values(roww, 19) = p; 
    
    % Kolmogorov Smirnoff  - WT vs HETD
    [h, p] = kstest2(vWT, vHETD);
    stats_values(roww, 20) = p;
    
    roww = roww+1;
end 


%% stats - average per sweep 

var_cols = [10, 14, 15, 16, 17,  21, 22, 23];

for k = 1:numel(var_cols)
    cc = var_cols(k); % Run through the 6 variables of interest.
    
    % Empty array for results.
    vw = zeros(4, 20);
    vh = zeros(4, 20);
    vwd = zeros(4, 20);
    vhd = zeros(4, 20);
    
    for jj = 4:1:23 % 20 steps
        
        allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep==jj );
        allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep==jj );
        allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep==jj );
        allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep==jj);
        
        vWT = spike_table{allWT, cc};
        vHET = spike_table{allHET, cc};
        vWTD = spike_table{allWTDTX, cc};
        vHETD = spike_table{allHETDTX, cc};
        
        vw(1,jj-3) = nanmean(vWT);
        vh(1,jj-3) = nanmean(vHET);
        vwd(1,jj-3) = nanmean(vWTD);
        vhd(1,jj-3) = nanmean(vHETD);
        
        vw(2,jj-3) = nanstd(vWT)/ sqrt(numel(vWT));
        vh(2,jj-3) = nanstd(vHET)/ sqrt(numel(vHET));
        vwd(2,jj-3) = nanstd(vWTD)/ sqrt(numel(vWTD));
        vhd(2,jj-3) = nanstd(vHETD)/ sqrt(numel(vHETD));
        
        vw(3,jj-3) = numel(vWT);
        vh(3,jj-3) = numel(vHET);
        vwd(3,jj-3) = numel(vWTD);
        vhd(3,jj-3) = numel(vHETD);
        
        vw(4,jj-3) = nanvar(vWT);
        vh(4,jj-3) = nanvar(vHET);
        vwd(4,jj-3) = nanvar(vWTD);
        vhd(4,jj-3) = nanvar(vHETD);
        
    end
    
    if k ==1
        stats.PeakAmp.WT.av = vw(1, :);
        stats.PeakAmp.WT.sem = vw(2, :);
        stats.PeakAmp.WT.num = vw(3, :);
        stats.PeakAmp.WT.var = vw(4, :);
        stats.PeakAmp.HET.av = vh(1, :);
        stats.PeakAmp.HET.sem = vh(2, :);
        stats.PeakAmp.HET.num = vh(3, :);
        stats.PeakAmp.HET.var = vh(4, :);
        stats.PeakAmp.WTD.av = vwd(1, :);
        stats.PeakAmp.WTD.sem = vwd(2, :);
        stats.PeakAmp.WTD.num = vwd(3, :);
        stats.PeakAmp.WTD.var = vwd(4, :);
        stats.PeakAmp.HETD.av = vhd(1, :);
        stats.PeakAmp.HETD.sem = vhd(2, :);
        stats.PeakAmp.HETD.num = vhd(3, :);
        stats.PeakAmp.HETD.var = vhd(4, :);
        
        [p1,h] = ranksum(vw(1,:), vh(1,:));
        stats.PeakAmp.RS.wthet = p1;
        [p2,h] = ranksum(vw(1,:), vwd(1,:));
        stats.PeakAmp.RS.wtwtd = p2;
        [p3,h] = ranksum(vh(1,:), vhd(1,:));
        stats.PeakAmp.RS.hethetd = p3;
        [p4,h] = ranksum(vw(1,:), vhd(1,:));
        stats.PeakAmp.RS.wthetd = p4;
        
        
    elseif k ==2
        stats.APT.WT.av = vw(1, :);
        stats.APT.WT.sem = vw(2, :);
        stats.APT.WT.num = vw(3, :);
        stats.APT.WT.var = vw(4, :);
        stats.APT.HET.av = vh(1, :);
        stats.APT.HET.sem = vh(2, :);
        stats.APT.HET.num = vh(3, :);
        stats.APT.HET.var = vh(4, :);
        stats.APT.WTD.av = vwd(1, :);
        stats.APT.WTD.sem = vwd(2, :);
        stats.APT.WTD.num = vwd(3, :);
        stats.APT.WTD.var = vwd(4, :);
        stats.APT.HETD.av = vhd(1, :);
        stats.APT.HETD.sem = vhd(2, :);
        stats.APT.HETD.num = vhd(3, :);
        stats.APT.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.APT.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.APT.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.APT.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.APT.RS.wthetd = p;
        
    elseif k == 3
        stats.Rise.WT.av = vw(1, :);
        stats.Rise.WT.sem = vw(2, :);
        stats.Rise.WT.num = vw(3, :);
        stats.Rise.WT.var = vw(4, :);
        stats.Rise.HET.av = vh(1, :);
        stats.Rise.HET.sem = vh(2, :);
        stats.Rise.HET.num = vh(3, :);
        stats.Rise.HET.var = vh(4, :);
        stats.Rise.WTD.av = vwd(1, :);
        stats.Rise.WTD.sem = vwd(2, :);
        stats.Rise.WTD.num = vwd(3, :);
        stats.Rise.WTD.var = vwd(4, :);
        stats.Rise.HETD.av = vhd(1, :);
        stats.Rise.HETD.sem = vhd(2, :);
        stats.Rise.HETD.num = vhd(3, :);
        stats.Rise.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.Rise.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.Rise.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.Rise.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.Rise.RS.wthetd = p;
        
    elseif k ==4
        stats.WHP.WT.av = vw(1, :);
        stats.WHP.WT.sem = vw(2, :);
        stats.WHP.WT.num = vw(3, :);
        stats.WHP.WT.var = vw(4, :);
        stats.WHP.HET.av = vh(1, :);
        stats.WHP.HET.sem = vh(2, :);
        stats.WHP.HET.num = vh(3, :);
        stats.WHP.HET.var = vh(4, :);
        stats.WHP.WTD.av = vwd(1, :);
        stats.WHP.WTD.sem = vwd(2, :);
        stats.WHP.WTD.num = vwd(3, :);
        stats.WHP.WTD.var = vwd(4, :);
        stats.WHP.HETD.av = vhd(1, :);
        stats.WHP.HETD.sem = vhd(2, :);
        stats.WHP.HETD.num = vhd(3, :);
        stats.WHP.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.WHP.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.WHP.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.WHP.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.WHP.RS.wthetd = p;
        
    elseif k ==5
        stats.Decay.WT.av = vw(1, :);
        stats.Decay.WT.sem = vw(2, :);
        stats.Decay.WT.num = vw(3, :);
        stats.Decay.WT.var = vw(4, :);
        stats.Decay.HET.av = vh(1, :);
        stats.Decay.HET.sem = vh(2, :);
        stats.Decay.HET.num = vh(3, :);
        stats.Decay.HET.var = vh(4, :);
        stats.Decay.WTD.av = vwd(1, :);
        stats.Decay.WTD.sem = vwd(2, :);
        stats.Decay.WTD.num = vwd(3, :);
        stats.Decay.WTD.var = vwd(4, :);
        stats.Decay.HETD.av = vhd(1, :);
        stats.Decay.HETD.sem = vhd(2, :);
        stats.Decay.HETD.num = vhd(3, :);
        stats.Decay.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.Decay.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.Decay.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.Decay.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.Decay.RS.wthetd = p;
    
      elseif k == 6
        stats.MaxDvDt.WT.av = vw(1, :);
        stats.MaxDvDt.WT.sem = vw(2, :);
        stats.MaxDvDt.WT.num = vw(3, :);
        stats.MaxDvDt.WT.var = vw(4, :);
        stats.MaxDvDt.HET.av = vh(1, :);
        stats.MaxDvDt.HET.sem = vh(2, :);
        stats.MaxDvDt.HET.num = vh(3, :);
        stats.MaxDvDt.HET.var = vh(4, :);
        stats.MaxDvDt.WTD.av = vwd(1, :);
        stats.MaxDvDt.WTD.sem = vwd(2, :);
        stats.MaxDvDt.WTD.num = vwd(3, :);
        stats.MaxDvDt.WTD.var = vwd(4, :);
        stats.MaxDvDt.HETD.av = vhd(1, :);
        stats.MaxDvDt.HETD.sem = vhd(2, :);
        stats.MaxDvDt.HETD.num = vhd(3, :);
        stats.MaxDvDt.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.MaxDvDt.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.MaxDvDt.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.MaxDvDt.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.MaxDvDt.RS.wthetd = p;
        
     elseif k == 7
        stats.MinDvDt.WT.av = vw(1, :);
        stats.MinDvDt.WT.sem = vw(2, :);
        stats.MinDvDt.WT.num = vw(3, :);
        stats.MinDvDt.WT.var = vw(4, :);
        stats.MinDvDt.HET.av = vh(1, :);
        stats.MinDvDt.HET.sem = vh(2, :);
        stats.MinDvDt.HET.num = vh(3, :);
        stats.MinDvDt.HET.var = vh(4, :);
        stats.MinDvDt.WTD.av = vwd(1, :);
        stats.MinDvDt.WTD.sem = vwd(2, :);
        stats.MinDvDt.WTD.num = vwd(3, :);
        stats.MinDvDt.WTD.var = vwd(4, :);
        stats.MinDvDt.HETD.av = vhd(1, :);
        stats.MinDvDt.HETD.sem = vhd(2, :);
        stats.MinDvDt.HETD.num = vhd(3, :);
        stats.MinDvDt.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.MinDvDt.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.MinDvDt.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.MinDvDt.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.MinDvDt.RS.wthetd = p;
        
        
    elseif k == 8
        stats.MinAmp.WT.av = vw(1, :);
        stats.MinAmp.WT.sem = vw(2, :);
        stats.MinAmp.WT.num = vw(3, :);
        stats.MinAmp.WT.var = vw(4, :);
        stats.MinAmp.HET.av = vh(1, :);
        stats.MinAmp.HET.sem = vh(2, :);
        stats.MinAmp.HET.num = vh(3, :);
        stats.MinAmp.HET.var = vh(4, :);
        stats.MinAmp.WTD.av = vwd(1, :);
        stats.MinAmp.WTD.sem = vwd(2, :);
        stats.MinAmp.WTD.num = vwd(3, :);
        stats.MinAmp.WTD.var = vwd(4, :);
        stats.MinAmp.HETD.av = vhd(1, :);
        stats.MinAmp.HETD.sem = vhd(2, :);
        stats.MinAmp.HETD.num = vhd(3, :);
        stats.MinAmp.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.MinAmp.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.MinAmp.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.MinAmp.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.MinAmp.RS.wthetd = p;
        
    end
    
    
end


%%

figure; plot(stats.PeakAmp.WT.av, 'k'); hold on; plot(stats.PeakAmp.HET.av, 'r'); plot(stats.PeakAmp.WTD.av, 'b'); plot(stats.PeakAmp.HETD.av, 'm')
figure; plot(stats.APT.WT.av, 'k'); hold on; plot(stats.APT.HET.av, 'r'); plot(stats.APT.WTD.av, 'b'); plot(stats.APT.HETD.av, 'm')
figure; plot(stats.WHP.WT.av, 'k'); hold on; plot(stats.WHP.HET.av, 'r'); plot(stats.WHP.WTD.av, 'b'); plot(stats.WHP.HETD.av, 'm')
figure; plot(stats.Decay.WT.av, 'k'); hold on; plot(stats.Decay.HET.av, 'r'); plot(stats.Decay.WTD.av, 'b'); plot(stats.Decay.HETD.av, 'm')
figure; plot(stats.MinAmp.WT.av, 'k'); hold on; plot(stats.MinAmp.HET.av, 'r'); plot(stats.MinAmp.WTD.av, 'b'); plot(stats.MinAmp.HETD.av, 'm')
figure; plot(stats.Rise.WT.av, 'k'); hold on; plot(stats.Rise.HET.av, 'r'); plot(stats.Rise.WTD.av, 'b'); plot(stats.Rise.HETD.av, 'm')
figure; plot(stats.MaxDvDt.WT.av, 'k'); hold on; plot(stats.MaxDvDt.HET.av, 'r'); plot(stats.MaxDvDt.WTD.av, 'b'); plot(stats.MaxDvDt.HETD.av, 'm')
figure; plot(stats.MinDvDt.WT.av, 'k'); hold on; plot(stats.MinDvDt.HET.av, 'r'); plot(stats.MinDvDt.WTD.av, 'b'); plot(stats.MinDvDt.HETD.av, 'm')


%% STATS PER SPIKE 

var_cols = [10, 14, 15, 16, 17, 21, 22, 23];

sw = 19; 

for k = 1:numel(var_cols)
    cc = var_cols(k); % Run through the 6 variables of interest.
    
    % Empty array for results.
    vw = zeros(4, 15);
    vh = zeros(4, 15);
    vwd = zeros(4, 15);
    vhd = zeros(4, 15);
    
    for jj = 1:1:15 %4:1:23 % 20 steps
        
        allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep==sw & spike_table.SpikeN == jj);
        allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep==sw & spike_table.SpikeN == jj);
        allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep==sw  & spike_table.SpikeN == jj);
        allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep==sw & spike_table.SpikeN == jj);

%         allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0 & spike_table.Sweep>10 & spike_table.Sweep<23 & spike_table.SpikeN == jj);
%         allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0 & spike_table.Sweep>10 & spike_table.Sweep<23 & spike_table.SpikeN == jj);
%         allWTDTX =  find(spike_table.Geno == 1 & spike_table.DTX ==1 & spike_table.Sweep>10 & spike_table.Sweep<23 & spike_table.SpikeN == jj);
%         allHETDTX =  find(spike_table.Geno == 2 & spike_table.DTX ==1 & spike_table.Sweep>10 & spike_table.Sweep<23 & spike_table.SpikeN == jj);

        vWT = spike_table{allWT, cc};
        vHET = spike_table{allHET, cc};
        vWTD = spike_table{allWTDTX, cc};
        vHETD = spike_table{allHETDTX, cc};
        
        vw(1,jj) = nanmean(vWT);
        vh(1,jj) = nanmean(vHET);
        vwd(1,jj) = nanmean(vWTD);
        vhd(1,jj) = nanmean(vHETD);
        
        vw(2,jj) = nanstd(vWT)/ sqrt(numel(vWT));
        vh(2,jj) = nanstd(vHET)/ sqrt(numel(vHET));
        vwd(2,jj) = nanstd(vWTD)/ sqrt(numel(vWTD));
        vhd(2,jj) = nanstd(vHETD)/ sqrt(numel(vHETD));
        
        vw(3,jj) = numel(vWT);
        vh(3,jj) = numel(vHET);
        vwd(3,jj) = numel(vWTD);
        vhd(3,jj) = numel(vHETD);
        
        vw(4,jj) = nanvar(vWT);
        vh(4,jj) = nanvar(vHET);
        vwd(4,jj) = nanvar(vWTD);
        vhd(4,jj) = nanvar(vHETD);
        
    end
    
    if k ==1
        stats.PeakAmp.WT.av = vw(1, :);
        stats.PeakAmp.WT.sem = vw(2, :);
        stats.PeakAmp.WT.num = vw(3, :);
        stats.PeakAmp.WT.var = vw(4, :);
        stats.PeakAmp.HET.av = vh(1, :);
        stats.PeakAmp.HET.sem = vh(2, :);
        stats.PeakAmp.HET.num = vh(3, :);
        stats.PeakAmp.HET.var = vh(4, :);
        stats.PeakAmp.WTD.av = vwd(1, :);
        stats.PeakAmp.WTD.sem = vwd(2, :);
        stats.PeakAmp.WTD.num = vwd(3, :);
        stats.PeakAmp.WTD.var = vwd(4, :);
        stats.PeakAmp.HETD.av = vhd(1, :);
        stats.PeakAmp.HETD.sem = vhd(2, :);
        stats.PeakAmp.HETD.num = vhd(3, :);
        stats.PeakAmp.HETD.var = vhd(4, :);
        
        [p1,h] = ranksum(vw(1,:), vh(1,:));
        stats.PeakAmp.RS.wthet = p1;
        [p2,h] = ranksum(vw(1,:), vwd(1,:));
        stats.PeakAmp.RS.wtwtd = p2;
        [p3,h] = ranksum(vh(1,:), vhd(1,:));
        stats.PeakAmp.RS.hethetd = p3;
        [p4,h] = ranksum(vw(1,:), vhd(1,:));
        stats.PeakAmp.RS.wthetd = p4;
        
        
    elseif k ==2
        stats.APT.WT.av = vw(1, :);
        stats.APT.WT.sem = vw(2, :);
        stats.APT.WT.num = vw(3, :);
        stats.APT.WT.var = vw(4, :);
        stats.APT.HET.av = vh(1, :);
        stats.APT.HET.sem = vh(2, :);
        stats.APT.HET.num = vh(3, :);
        stats.APT.HET.var = vh(4, :);
        stats.APT.WTD.av = vwd(1, :);
        stats.APT.WTD.sem = vwd(2, :);
        stats.APT.WTD.num = vwd(3, :);
        stats.APT.WTD.var = vwd(4, :);
        stats.APT.HETD.av = vhd(1, :);
        stats.APT.HETD.sem = vhd(2, :);
        stats.APT.HETD.num = vhd(3, :);
        stats.APT.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.APT.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.APT.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.APT.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.APT.RS.wthetd = p;
        
    elseif k == 3
        stats.Rise.WT.av = vw(1, :);
        stats.Rise.WT.sem = vw(2, :);
        stats.Rise.WT.num = vw(3, :);
        stats.Rise.WT.var = vw(4, :);
        stats.Rise.HET.av = vh(1, :);
        stats.Rise.HET.sem = vh(2, :);
        stats.Rise.HET.num = vh(3, :);
        stats.Rise.HET.var = vh(4, :);
        stats.Rise.WTD.av = vwd(1, :);
        stats.Rise.WTD.sem = vwd(2, :);
        stats.Rise.WTD.num = vwd(3, :);
        stats.Rise.WTD.var = vwd(4, :);
        stats.Rise.HETD.av = vhd(1, :);
        stats.Rise.HETD.sem = vhd(2, :);
        stats.Rise.HETD.num = vhd(3, :);
        stats.Rise.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.Rise.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.Rise.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.Rise.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.Rise.RS.wthetd = p;
        
    elseif k ==4
        stats.WHP.WT.av = vw(1, :);
        stats.WHP.WT.sem = vw(2, :);
        stats.WHP.WT.num = vw(3, :);
        stats.WHP.WT.var = vw(4, :);
        stats.WHP.HET.av = vh(1, :);
        stats.WHP.HET.sem = vh(2, :);
        stats.WHP.HET.num = vh(3, :);
        stats.WHP.HET.var = vh(4, :);
        stats.WHP.WTD.av = vwd(1, :);
        stats.WHP.WTD.sem = vwd(2, :);
        stats.WHP.WTD.num = vwd(3, :);
        stats.WHP.WTD.var = vwd(4, :);
        stats.WHP.HETD.av = vhd(1, :);
        stats.WHP.HETD.sem = vhd(2, :);
        stats.WHP.HETD.num = vhd(3, :);
        stats.WHP.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.WHP.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.WHP.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.WHP.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.WHP.RS.wthetd = p;
        
    elseif k ==5
        stats.Decay.WT.av = vw(1, :);
        stats.Decay.WT.sem = vw(2, :);
        stats.Decay.WT.num = vw(3, :);
        stats.Decay.WT.var = vw(4, :);
        stats.Decay.HET.av = vh(1, :);
        stats.Decay.HET.sem = vh(2, :);
        stats.Decay.HET.num = vh(3, :);
        stats.Decay.HET.var = vh(4, :);
        stats.Decay.WTD.av = vwd(1, :);
        stats.Decay.WTD.sem = vwd(2, :);
        stats.Decay.WTD.num = vwd(3, :);
        stats.Decay.WTD.var = vwd(4, :);
        stats.Decay.HETD.av = vhd(1, :);
        stats.Decay.HETD.sem = vhd(2, :);
        stats.Decay.HETD.num = vhd(3, :);
        stats.Decay.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.Decay.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.Decay.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.Decay.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.Decay.RS.wthetd = p;
        
    elseif k == 6
        stats.MaxDvDt.WT.av = vw(1, :);
        stats.MaxDvDt.WT.sem = vw(2, :);
        stats.MaxDvDt.WT.num = vw(3, :);
        stats.MaxDvDt.WT.var = vw(4, :);
        stats.MaxDvDt.HET.av = vh(1, :);
        stats.MaxDvDt.HET.sem = vh(2, :);
        stats.MaxDvDt.HET.num = vh(3, :);
        stats.MaxDvDt.HET.var = vh(4, :);
        stats.MaxDvDt.WTD.av = vwd(1, :);
        stats.MaxDvDt.WTD.sem = vwd(2, :);
        stats.MaxDvDt.WTD.num = vwd(3, :);
        stats.MaxDvDt.WTD.var = vwd(4, :);
        stats.MaxDvDt.HETD.av = vhd(1, :);
        stats.MaxDvDt.HETD.sem = vhd(2, :);
        stats.MaxDvDt.HETD.num = vhd(3, :);
        stats.MaxDvDt.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.MaxDvDt.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.MaxDvDt.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.MaxDvDt.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.MaxDvDt.RS.wthetd = p;
        
     elseif k == 7
        stats.MinDvDt.WT.av = vw(1, :);
        stats.MinDvDt.WT.sem = vw(2, :);
        stats.MinDvDt.WT.num = vw(3, :);
        stats.MinDvDt.WT.var = vw(4, :);
        stats.MinDvDt.HET.av = vh(1, :);
        stats.MinDvDt.HET.sem = vh(2, :);
        stats.MinDvDt.HET.num = vh(3, :);
        stats.MinDvDt.HET.var = vh(4, :);
        stats.MinDvDt.WTD.av = vwd(1, :);
        stats.MinDvDt.WTD.sem = vwd(2, :);
        stats.MinDvDt.WTD.num = vwd(3, :);
        stats.MinDvDt.WTD.var = vwd(4, :);
        stats.MinDvDt.HETD.av = vhd(1, :);
        stats.MinDvDt.HETD.sem = vhd(2, :);
        stats.MinDvDt.HETD.num = vhd(3, :);
        stats.MinDvDt.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.MinDvDt.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.MinDvDt.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.MinDvDt.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.MinDvDt.RS.wthetd = p;
        
        
    elseif k == 8
        stats.MinAmp.WT.av = vw(1, :);
        stats.MinAmp.WT.sem = vw(2, :);
        stats.MinAmp.WT.num = vw(3, :);
        stats.MinAmp.WT.var = vw(4, :);
        stats.MinAmp.HET.av = vh(1, :);
        stats.MinAmp.HET.sem = vh(2, :);
        stats.MinAmp.HET.num = vh(3, :);
        stats.MinAmp.HET.var = vh(4, :);
        stats.MinAmp.WTD.av = vwd(1, :);
        stats.MinAmp.WTD.sem = vwd(2, :);
        stats.MinAmp.WTD.num = vwd(3, :);
        stats.MinAmp.WTD.var = vwd(4, :);
        stats.MinAmp.HETD.av = vhd(1, :);
        stats.MinAmp.HETD.sem = vhd(2, :);
        stats.MinAmp.HETD.num = vhd(3, :);
        stats.MinAmp.HETD.var = vhd(4, :);
        
        [p,h] = ranksum(vw(1,:), vh(1,:));
        stats.MinAmp.RS.wthet = p;
        [p,h] = ranksum(vw(1,:), vwd(1,:));
        stats.MinAmp.RS.wtwtd = p;
        [p,h] = ranksum(vh(1,:), vhd(1,:));
        stats.MinAmp.RS.hethetd = p;
        [p,h] = ranksum(vw(1,:), vhd(1,:));
        stats.MinAmp.RS.wthetd = p;
        
    end
    
    
end


%%
%         [p,h] = ranksum(vWT, vHET);
%         stats.PeakAmp.RS.wthet = p; 
%         [p,h] = ranksum(vWT, vWTD);
%         stats.PeakAmp.RS.wtwtd = p; 
%         [p,h] = ranksum(vHET, vHETD);
%         stats.PeakAmp.RS.hethetd = p; 
%         [p,h] = ranksum(vWT, vHETD);
%         stats.PeakAmp.RS.wthetd = p; 





%%

% Number of cells - unique by mouse, cell and cohort. 
[univals, ia, ic] = unique(spike_table(:, [3,5, 18]));
MaxAmpTable = groupsummary(spike_table, {'Geno', 'Sweep', 'DTX'}, {'mean'}, {'MaxAmp', 'MinAmp'});

% univals = the unique vals in each group
% ia = the first row in spike_table that belongs to that group
% ic = a vector the height of the table x 1- which contains the group # for each row in the table. 

[G, ID] = findgroups(spike_table(:, [1,2,4]));
% G contains the 'group' id for each row in spike_table.
% ID contains the details of the different groups. 

%
MaxAmpTable = groupsummary(spike_table, {'Geno', 'Cell', 'Ani', 'Cohort'}, {'mean', 'std', 'range', 'max', 'min'}, 'MaxAmp');

allWT = find(MaxAmpTable.Geno ==1);
allHET = find(MaxAmpTable.Geno ==2);

wtvals = (MaxAmpTable.mean_MaxAmp(allWT));
hetvals = (MaxAmpTable.mean_MaxAmp(allHET));

[p,h] = ranksum(wtvals, hetvals)

%% All spikes - grouped by geno - across sweeps
MaxAmpTable = groupsummary(spike_table, {'Geno', 'Sweep'}, 'mean', {'MaxAmp', 'RiseT', 'DecayT', 'whp', 'APThresh'});

allWT = find(MaxAmpTable.Geno ==1 & MaxAmpTable.Sweep>3 & MaxAmpTable.Sweep<22);
allHET = find(MaxAmpTable.Geno ==2 & MaxAmpTable.Sweep>3 & MaxAmpTable.Sweep<22);

figure;
plot(MaxAmpTable.Sweep(allWT), MaxAmpTable.mean_(allWT), 'k', 'LineWidth', 1.2);
hold on 
plot(MaxAmpTable.Sweep(allHET), MaxAmpTable.mean_whp(allHET), 'r', 'LineWidth', 1.2);
title('whp')

%% All spikes - grouped by geno - across spikes
MaxAmpTable = groupsummary(spike_table, {'Geno', 'Sweep', 'SpikeN'}, 'mean', {'MaxAmp', 'RiseT', 'DecayT', 'whp', 'APThresh'});

figure
for i = 1:10

allWT = find(MaxAmpTable.Geno ==1 & MaxAmpTable.SpikeN==i & MaxAmpTable.Sweep>3 & MaxAmpTable.Sweep<22);
allHET = find(MaxAmpTable.Geno ==2 & MaxAmpTable.SpikeN==i & MaxAmpTable.Sweep>3 & MaxAmpTable.Sweep<22);

% allWT = find(MaxAmpTable.Geno ==1 & MaxAmpTable.SpikeN<76 & MaxAmpTable.Sweep == 9);
% allHET = find(MaxAmpTable.Geno ==2 & MaxAmpTable.SpikeN<76 & MaxAmpTable.Sweep == 9);

hold on 
plot(MaxAmpTable.Sweep(allWT), MaxAmpTable.mean_RiseT(allWT), 'k', 'LineWidth', 1.2);
hold on 
plot(MaxAmpTable.Sweep(allHET), MaxAmpTable.mean_RiseT(allHET), 'r', 'LineWidth', 1.2);
title('peak')

x = input('');

end 


figure
for i = 4:21

allWT = find(MaxAmpTable.Geno ==1 & MaxAmpTable.Sweep==i & MaxAmpTable.SpikeN<60);
allHET = find(MaxAmpTable.Geno ==2 & MaxAmpTable.Sweep==i &  MaxAmpTable.SpikeN<60);

% allWT = find(MaxAmpTable.Geno ==1 & MaxAmpTable.SpikeN<76 & MaxAmpTable.Sweep == 9);
% allHET = find(MaxAmpTable.Geno ==2 & MaxAmpTable.SpikeN<76 & MaxAmpTable.Sweep == 9);

hold on 
plot(MaxAmpTable.SpikeN(allWT), MaxAmpTable.mean_MaxAmp(allWT), 'k', 'LineWidth', 1.2);
hold on 
plot(MaxAmpTable.SpikeN(allHET), MaxAmpTable.mean_MaxAmp(allHET), 'r', 'LineWidth', 1.2);
title('Max')

x = input('');

end 

%%

% Grouped by geno:
[genoid, genogps] = findgroups(spike_table(:, 3));

% Grouped by cell:
[cellid, cellgps] = findgroups(spike_table(:, [3,1,2,4, 18]));

% Grouped geno + all spikes
[spikeid, spinfo] = findgroups(spike_table(:, [3,6]));

splitapply(@mean,spike_table.MaxAmp,genoid)
mvals = splitapply(@nanmean,spike_table.whp,cellid)

[p,h] = ranksum(mvals(1:19), mvals(20:end))


allDTX = find(spike_table.DTX ==1);
sp_table2 = spike_table; % keep original
spike_table(allDTX, :) = [];


%% Generate subtable without DTX!

% CREATE COPY WITH ALL DATA in memory
spike_table2 = spike_table;

% find rows with spikes that are from dtx experiments
allDTx = find(spike_table.DTX == 1);
% remove them 
spike_table(allDTx, :) = []; 

%% Generate subtable ONLY with the VERY first spike ever generated for each cell! 

% Group by cell
[cellid, cellgps] = findgroups(spike_table(:, [3,1,2,4,18]));

rows_firstsweep = [];

for k = 1:numel(cellgps(:,1))
    
    allcell = find(cellid(:,1)==k);
    firstsweep = min(spike_table.Sweep(allcell));
    cellgps.FirstSweep(k) = firstsweep;
    
    allrows = find(cellid(:,1)==k & spike_table.Sweep==firstsweep);
    
    rows_firstsweep = vertcat(rows_firstsweep, allrows);
end 

spike_table = spike_table(rows_firstsweep, :);

% save('220811_Spike_Table_OnlySpikesin1stSweep_PerCell.mat', 'spike_table');

save('221104_Cul3_Spikes_from_Rheobase.mat', 'spike_table');

%% Errobar - Peak amplitude - across sweeps

dWT = stats.PeakAmp.WT.av;
dHET = stats.PeakAmp.HET.av;

sWT = stats.PeakAmp.WT.sem;
sHET = stats.PeakAmp.HET.sem;

figure
errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0], 'MarkerFaceColor', [1 0.7 0.7] ,'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0],'Marker', 'none', 'LineWidth', 1.75)

box off
xlim([0 21])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;

f = gcf;
f.Position = [469   470   327   261]; % [ 469   443   334   288]; 


   
%% MEAN/SEM shaded - Peak amplitude - across sweeps - WTHET in background - DTX in foreground

dWT = stats.Decay.WT.av;
dHET = stats.Decay.HET.av;

sWT = stats.Decay.WT.sem;
sHET = stats.Decay.HET.sem;

% MEANSEM Line
x = 1:1:20;
mean_WT = dWT; 
mean_HET = dHET; 
y1 = mean_WT+sWT;
y2 = mean_WT-sWT;
y3 = mean_HET+sHET;
y4 = mean_HET-sHET;

figure
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
%     plot(mean_WT', '-', 'Color', 'k', 'LineWidth', 1.3)
plot(x, y3, 'w-')
hold on 
plot(x, y4, 'w-')
patch([x fliplr(x)], [y3 fliplr(y4)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
%     plot(mean_HET', '-','Color', 'r', 'LineWidth', 1.3)


% ADD DTX 

dWT = stats.Decay.WTD.av;
dHET = stats.Decay.HETD.av;

sWT = stats.Decay.WTD.sem;
sHET = stats.Decay.HETD.sem;
% 

x = 1:1:20;
mean_WT = dWT; 
mean_HET = dHET; 
y1 = mean_WT+sWT;
y2 = mean_WT-sWT;
y3 = mean_HET+sHET;
y4 = mean_HET-sHET;


plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(mean_WT', '-', 'Color', 'b', 'LineWidth', 1.3)
plot(x, y3, 'w-')
hold on 
plot(x, y4, 'w-')
patch([x fliplr(x)], [y3 fliplr(y4)], [0.8 0 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(mean_HET', '-','Color', [0.8 0 0.8], 'LineWidth', 1.3)

box off
xlim([0 17])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;

f = gcf;
f.Position = [334   445   355   205]; % [469   470   327   261]; % [ 469   443   334   288]; 




% errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 1],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
% hold on 
% errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0], 'MarkerFaceColor', [1 0.7 1] ,'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
% errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 1], 'Marker', 'none', 'LineWidth', 1.75)
% errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 1],'Marker', 'none', 'LineWidth', 1.75)



%%  ACROSS SWEEPS

dWT = stats.Decay.WT.av;
dHET = stats.Decay.HET.av;

sWT = stats.Decay.WT.sem;
sHET = stats.Decay.HET.sem;

figure
errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0], 'MarkerFaceColor', [1 0.7 0.7] ,'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0],'Marker', 'none', 'LineWidth', 1.75)

box off
xlim([0 21])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;

f = gcf;
f.Position = [469   470   327   261]; % [ 469   443   334   288]; 


dWT = stats.PeakAmp.WTD.av;
dHET = stats.PeakAmp.HETD.av;

sWT = stats.PeakAmp.WTD.sem;
sHET = stats.PeakAmp.HETD.sem;

errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 1],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0], 'MarkerFaceColor', [1 0.7 1] ,'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 1], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 1],'Marker', 'none', 'LineWidth', 1.75)





%% WITHIN SWEEP - ACROSS SPIKES

dWT = stats.PeakAmp.WT.av;
dHET = stats.PeakAmp.HET.av;

sWT = stats.PeakAmp.WT.sem;
sHET = stats.PeakAmp.HET.sem;

if sw == 5
    modv = 0.85;
elseif sw == 9 
    modv = 0.65;
elseif sw == 15
    modv = 0.55;
elseif sw == 19
    modv = 0.4;
end 
    
c1 = [modv modv modv];
c2 = [1 modv modv];

figure
errorbar(1:1:15, dWT, sWT, 'o', 'CapSize', 0, 'Color', c1, 'MarkerFaceColor', c1,'MarkerEdgeColor', 'none',  'MarkerSize', 6, 'LineWidth', 1.75)
hold on 
errorbar(1:1:15, dHET, sHET, 'o', 'CapSize', 0, 'Color', c2, 'MarkerFaceColor', c2 ,'MarkerEdgeColor','none', 'MarkerSize', 6, 'LineWidth', 1.75) %[1 0.8 0.8]

% figure
% hold on
% scatter(1:1:15, dWT,100, 'o', 'MarkerFaceColor', c1, 'MarkerEdgeColor', c1);
% scatter(1:1:15, dHET,100, 'o', 'MarkerFaceColor', c2, 'MarkerEdgeColor', c2);


% errorbar(1:1:15, dWT, sWT, 'o', 'CapSize', 0, 'Color', c1, 'Marker', 'none', 'LineWidth', 1.75)
% errorbar(1:1:15, dHET, sHET, 'o', 'CapSize', 0, 'Color', c2,'Marker', 'none', 'LineWidth', 1.75)

xlim([0 16])
ylim([0 50])
box off
% ADD DTX ON TOP
% dWT = stats.PeakAmp.WTD.av;
% dHET = stats.PeakAmp.HETD.av;
% 
% sWT = stats.PeakAmp.WTD.sem;
% sHET = stats.PeakAmp.HETD.sem;
% 
% errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 1],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
% hold on 
% errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0], 'MarkerFaceColor', [1 0.7 1] ,'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
% errorbar(1:1:20, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 1], 'Marker', 'none', 'LineWidth', 1.75)
% errorbar(1:1:20, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 1],'Marker', 'none', 'LineWidth', 1.75)




%% Add MAX dvdt (700:1000), MIN dvdt(1000:1400) , minnnn V (1000-1250)

figure; plot(spike_table.V(1,:))

for i = 1:height(spike_table)
    
    V = spike_table.V(i,:);
    dvdt = spike_table.dvdt(i,:);
    
    % MAXDVDT
    maxdvdt = max(dvdt(700:1000));
    spike_table.MaxDVDT(i) = maxdvdt;
    
    % MINDVDT 
    mindvdt = min(dvdt(1000:1400));
     spike_table.MinDVDT(i) = mindvdt;
     
    % MinAmp 
    minamp = min(V(1000:1400));
    spike_table.MinAmp2(i) = minamp;
end 

%% Add AVERAGE dvdt (700:1000), AVERAGE dvdt(1000:1400)

for i = 1:height(spike_table)
   
    dvdt = spike_table.dvdt(i,:);
    
    % avRISEdvdt
    avRISEdvdt = nanmean(dvdt(700:1000));
    spike_table.AvDVDTrise(i) = avRISEdvdt;
   
    % avDECAYdvdt 
    avDECAYdvdt = nanmean(dvdt(1000:1400));
    spike_table.AvDVDTdecay(i) = avDECAYdvdt;
     
end 

%% For each sweep - determine the % of WT/ HET cells firing

prop_firing = zeros(4, 20);

for kk = 1:20 
    
    allrows = find(spike_table.Sweep == kk & spike_table.DTX == 0); % NO DTX
    subtable = spike_table(allrows, :);
    
    gp = unique(subtable(:, [1,2,3,4]));
    
    w = find(gp.Geno ==1);
    h = find(gp.Geno ==2);
    
    prop_firing(1,kk) = (numel(w)/19);
    prop_firing(2,kk) = (numel(h)/22);
    prop_firing(3,kk) = (numel(w));
    prop_firing(4,kk) = (numel(h));
end 

figure; plot(prop_firing(1,:), 'k'); hold on; plot(prop_firing(2,:), 'r'); 
figure; plot(prop_firing(3,:), 'k'); hold on; plot(prop_firing(4,:), 'r'); 

bar(1:1:20, prop_firing(1,:)); hold on; bar(1:1:20, prop_firing(2,:));




%% Analysing variables at RHEOBASE (Current to first reach AP thresh) - all 4 conditions. 

% 1 - FIRST - create subtable with only the spikes that were generated during the first current injection sweep to elicit APS. 

% CREATE COPY WITH ALL DATA in memory
spike_table2 = spike_table;

% Group by cell
[cellid, cellgps] = findgroups(spike_table(:, [1,2,4]));

rows_firstsweep = [];

for k = 1:numel(cellgps(:,1))
    
    allcell = find(cellid(:,1)==k);
    firstsweep = min(spike_table.Sweep(allcell));
    cellgps.FirstSweep(k) = firstsweep;
    
    allrows = find(cellid(:,1)==k & spike_table.Sweep==firstsweep);
    
    rows_firstsweep = vertcat(rows_firstsweep, allrows);
end 

spike_table = spike_table(rows_firstsweep, :);

% save('221007_Ptchd1_Spikes_ONLY_Rheobase_N8.mat', 'spike_table');

% 2 - FIND - CURRENT TO FIRST SPIKE. 
wt_firstspike = find(spike_table.Geno == 1 & spike_table.DTX == 0 & spike_table.SpikeN ==1);
het_firstspike = find(spike_table.Geno == 2 & spike_table.DTX == 0 & spike_table.SpikeN ==1);
% wtD_firstspike = find(spike_table.Geno == 1 & spike_table.DTX == 1 & spike_table.SpikeN ==1);
% hetD_firstspike = find(spike_table.Geno == 2 & spike_table.DTX == 1 & spike_table.SpikeN ==1);

wtvals = (spike_table.Sweep(wt_firstspike)-3)*10;
hetvals = (spike_table.Sweep(het_firstspike)-3)*10;
% wtDvals = (spike_table.Sweep(wtD_firstspike)-3)*10;
% hetDvals = (spike_table.Sweep(hetD_firstspike)-3)*10;

% STATS
[p,h] = ranksum(wtvals, hetvals)
% [p,h] = ranksum(wtvals, wtDvals)
% [p,h] = ranksum(wtvals, hetDvals)
% [p,h] = ranksum(hetvals, hetDvals)


nanmean(wtvals)
nanmean(hetvals)

%% BOXPLOT - WT/HET

n_wt = numel(wtvals);
n_het = numel(hetvals);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

figure
scatter(x1, wtvals,'SizeData', 150, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, hetvals, 'SizeData', 150,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvals', hetvals'], [ones(1,n_wt), ones(1,n_het)*2], 'Colors', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
box off
axis([0.5 2.5 -10 100])
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ax.TickLength = [0.025 0.025];
ax.LineWidth = 2;

f = gcf;
f.Position = [704   207   355   572]; 
f.Renderer = 'painters';

% histogram(wtvals, 0:10:90, 'Normalization', 'cdf'); hold on; histogram(hetvals, 0:10:90, 'Normalization', 'cdf');

%% 4 - FIND - average dvdt/ peak amp / min amp. 


wt_firstspike = find(spike_table.Geno == 1 & spike_table.DTX == 0); %& spike_table.SpikeN ==1
het_firstspike = find(spike_table.Geno == 2 & spike_table.DTX == 0);
% wtD_firstspike = find(spike_table.Geno == 1 & spike_table.DTX == 1);
% hetD_firstspike = find(spike_table.Geno == 2 & spike_table.DTX == 1);

wtvals = spike_table.ExpDecayTau(wt_firstspike);
hetvals = spike_table.ExpDecayTau(het_firstspike);
% wtDvals = spike_table.RiseT(wtD_firstspike);
% hetDvals = spike_table.RiseT(hetD_firstspike);

% STATS
% [p,h] = ranksum(wtvals, wtDvals)
% [p,h] = ranksum(wtvals, hetDvals)
% [p,h] = ranksum(hetvals, hetDvals)

%% PLOT - WT versus HET 
n_wt = numel(wtvals);
n_het = numel(hetvals);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

figure
% scatter(x1, wtvals,'SizeData', 50, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
% scatter(x2, hetvals, 'SizeData', 50,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvals', hetvals'], [ones(1,n_wt), ones(1,n_het)*2], 'Colors', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
box off
xlim([0.5 2.5])
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ax.TickLength = [0.025 0.025];
ax.LineWidth = 2;

% ylim([-0.01 0.15])
ylim([-75 -50])
ylim([0 4])
ylim([0 1.75])
ylim([0 1.5])
ylim([-50 -10])
ylim([10 70])

f = gcf;
f.Position = [ 705   451   216   327]; %[704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

nanmean(wtvals)
nanmean(hetvals)
[p,h] = ranksum(wtvals, hetvals)
[h, p] = kstest2(wtvals, hetvals)





%% PLOT WTD versus HETD

n_wt = numel(wtDvals);
n_het = numel(hetDvals);

figure
b = boxplot([wtDvals', hetDvals'], [ones(1,n_wt), ones(1,n_het)*2], 'Colors', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
box off
xlim([0.5 2.5]) 
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ax.TickLength = [0.025 0.025];
ax.LineWidth = 2;

ylim([0 8])

f = gcf;
f.Position = [704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

nanmean(wtDvals)
nanmean(hetDvals)
[p,h] = ranksum(wtDvals, hetDvals)
[h,p] = kstest2(wtDvals, hetDvals)


%% 
[p,h] = ranksum(wtvals, wtDvals)
[p,h] = ranksum(hetvals, hetDvals)



%% 

st = 1050;
nd = 1075;

x = 1050:1:1075;
vWT = spike_table.V(allWT, st:nd);
vHET = spike_table.V(allHET,  st:nd);

% figure; plot(nanmean(vWT), 'k'); hold on; plot(nanmean(vHET), 'r');

mWT = nanmean(vWT);
mHET = nanmean(vHET);

figure; plot(nanmean(vWT),  'k'); hold on; 
f = fit(mWT',x', 'exp2');
plot(f, 'b')




%% TAU - EXP DECAY 

data = spike_table.V(i, :);
figure; plot(data)
y = data(1015:1100)';
x = [1:1:numel(y)]';

f0 = fit(x,y,g, 'StartPoint', [[ones(size(x)), -exp(-x)]\y; 1]);


xx = linspace(1,100,100);
plot(x,y,'o'); hold on; plot(xx,f0(xx),'r-')

%%

g = fittype('a-b*exp(-c*x)');


for i = 1:height(spike_table)
    data = spike_table.V(i, :);
    y = data(1015:1100)';
    x = [1:1:numel(y)]';
    
    f0 = fit(x,y,g, 'StartPoint', [[ones(size(x)), -exp(-x)]\y; 1]);

    spike_table.ExpDecayTau(i) = f0.c;
    
end 



%% PLOTS FOR EACH CELL!!!
% Group by cell
[cellid, cellgps] = findgroups(spike_table(:, [1,2,3, 4]));

st = 500;
nd = 1500;
% col = [1 0.65 0];

figure
for k = 1:numel(cellgps(:,1))
    
    subplot(4,8,k)
    all_cell = find(cellid == k);
    
    vWT = spike_table.V(all_cell, st:nd);
    dvWT = spike_table.dvdt(all_cell,  st:nd);

    gg = spike_table.Geno(all_cell(1));
    if gg == 1
        col = 'k';
    else
        col = 'r';%[1 0.65 0];
    end 

    st2 = 1;
    nd2 = 1000;
    
% AP average
% subplot(2,1,1);
%  plot(st2:1:nd2, nanmean(vWT(:, st2:nd2)), col, 'LineWidth', 1.2)
% hold on; 
% xlim([st2 nd2])
% ylim([-80 60])
% f = gcf;
% f.Position = [440   535   315   263];
% box off
% ax = gca; 
% ax.TickDir = 'out'; 
% ax.LineWidth = 1;
% ax.TickLength = [0.02 0.02];

% Phase plane plot
% subplot(2,1,2);
% aaa = nanmean(vWT(:, st2:nd2));
% bbb = nanmean(dvWT(:, st2:nd2));
% 
% vvv = vertcat(vvv, aaa);
% dvvv = vertcat(dvvv, bbb);

plot(nanmean(vWT(:, st2:nd2)), nanmean(dvWT(:, st2:nd2)), 'Color', col, 'LineWidth', 1.2);
hold on; 
xlim([-75 55])
ylim([-400 400])
f = gcf;
f.Position = [440   535   315   263];
box off
ax = gca; 
ax.TickDir = 'out'; 
ax.LineWidth = 1;
ax.TickLength = [0.02 0.02];
title(string(k))
end 

f = gcf;
f.Position = [260 20 749  1030]; 


%% KMEANS - cluster cells

[idx,C,sumd,D] = kmeans(dvvv, 6);
% figure; imagesc(D)

dataa = [dvvv, idx];
dataa = sortrows(dataa, 1001);

imagesc(dataa)

[cellid, cellgps] = findgroups(spike_table(:, [1,2,3, 4]));
cellgps.idx = idx;
cellgps = sortrows(cellgps, 5);

allW = find(cellgps.Geno ==1);
allH = find(cellgps.Geno ==2);

figure; subplot(1,5,1:4); imagesc(dataa(allW, 250:750)); subplot(1,5,5); imagesc(dataa(allW, 1001))
figure; subplot(1,5,1:4); imagesc(dataa(allH, 250:750)); subplot(1,5,5); imagesc(dataa(allH, 1001))



%% STATS PER CELL 

values_per_cell = []; 

for k = 1:numel(cellgps(:,1))
    
    ddd = [];
    all_cell = find(cellid == k);
    
    mousee = spike_table.Ani(all_cell(1));
    celll = spike_table.Cell(all_cell(1));
    genoo = spike_table.Geno(all_cell(1));
    
    ddd(1) = mousee;
    ddd(2) = celll;
    ddd(3) = genoo;
    
    idn = 4;
    
    for p = [10, 14, 15, 16, 17, 19, 25] %23,25] % 22
       
        vals = spike_table{all_cell, p};
        ddd(idn) = nanmean(vals);
        idn = idn+1;
    end 

    values_per_cell = vertcat(values_per_cell, ddd);

end 

% 10 - max amp
% 14 - APT
% 15 - Rise T
% 16 - WHP
% 17 - Decay T
% 19 - Min Amp

% 22 - Decay 

% 22 - Max DVDT
% 23 - MinDVDT
% 25 - Tau 

save('221114_Setd5_SpikeShapeProps_Values_Per_Cell.mat', 'values_per_cell')

%% PLOTS - PER CELL - DOT BOX

allWT = find(values_per_cell(:,3)==1);
allHET = find(values_per_cell(:,3)==2);

collumn = 4;

wtvals = values_per_cell(allWT, collumn);
hetvals = values_per_cell(allHET, collumn);


n_wt = numel(wtvals);
n_het = numel(hetvals);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

figure
scatter(x1, wtvals,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, hetvals, 'SizeData', 100,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvals', hetvals'], [ones(1,n_wt), ones(1,n_het)*2], 'Colors', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
box off
xlim([0.5 2.5])
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ax.TickLength = [0.025 0.025];
ax.LineWidth = 2;

ylim([0 0.3])

f = gcf;
f.Position = [1017 351  226  322]; %[1373  50  226 322]; %[ 705   451   216   327]; %[704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

nanmean(wtvals)
nanmean(hetvals)
[p,h] = ranksum(wtvals, hetvals)
% [p,h] = ttest2(wtvals, hetvals)
[h, p] = kstest2(wtvals, hetvals)



% 16, 23, 25
% Ptchd1 - significantly different in WHP, Min DVDT and Exp Tau. 

%% Intrinsic properties - DOT BOX

% col = [1 0.55 0.2];
col = 'm';

allWT = find(tbl.Geno ==1); 
allHET = find(tbl.Geno ==2); 

wtvals = tbl.Rin(allWT);
hetvals = tbl.Rin(allHET);

n_wt = numel(wtvals);
n_het = numel(hetvals);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

figure
scatter(x1, wtvals,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, hetvals, 'SizeData', 100,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvals', hetvals'], [ones(1,n_wt), ones(1,n_het)*2], 'Colors', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
box off
xlim([0.5 2.5])
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ax.TickLength = [0.025 0.025];
ax.LineWidth = 2;

ylim([-75 -50])

f = gcf;
f.Position = [1045  313  226  322]; %[1373  50  226 322]; %[ 705   451   216   327]; %[704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

nanmean(wtvals)
nanmean(hetvals)
[p,h] = ranksum(wtvals, hetvals)
[p,h] = ttest2(wtvals, hetvals)
[h, p] = kstest2(wtvals, hetvals)

