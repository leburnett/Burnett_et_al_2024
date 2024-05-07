% COMBINE FLASH SPIKING DATA - PUPIL DIAMETER - RUNNING OF ANIMAL.
% Created by Burnett - 27/06/2022

% Location of FLASH / PUPIL DIAMETER DATA from Olga. 
% /Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Behaviour/Flash

% Location of pupil size data:
% /Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/SP_DATA/Eye_DLC/eye_files

% Dynamic range of WT/ HET across recording? 


%% ANALYSIS OF PUPIL DILATION - WT versus HET

% DATA = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Setd5_Motion_Pupil_Spikes_Combined/Olga_Flash_Eye_NEW'; 

% 1 - no knowledge of running. 
% Combine values of PUPIL over 1 rep (A, av_rep_3) and 10 reps (B, r_thru)
% Combine values of SPIKING - average over all 'units' to 1 rep (C) and to 10 reps (D)
% - average over units with more than 2 spikes during stimulus - to 1 rep
% (E) and 10 reps (F). 

FLASH_PUPIL_SPIKES = table();

files = dir('*.mat');

for i = 1:numel(files)
    fname = files(i).name;
    rrun = fname(end-10);
    load(fname)
    
%     % meta info
%     FLASH_PUPIL_SPIKES.Animal{i} = stim_info.meta.animal;
%     FLASH_PUPIL_SPIKES.Date{i} = stim_info.meta.date;
%     FLASH_PUPIL_SPIKES.Test{i} = stim_info.meta.test_number;
%     FLASH_PUPIL_SPIKES.Run{i} = rrun;
%     FLASH_PUPIL_SPIKES.Depth{i} = stim_info.meta.depth;
%     FLASH_PUPIL_SPIKES.Geno{i} = stim_info.meta.isWT;
%     
%     FLASH_PUPIL_SPIKES.AVI_F1{i} = stim_info.timing.avi_first_frame;
%     FLASH_PUPIL_SPIKES.AVI_FN{i} = stim_info.timing.avi_last_frame;
%     
%     FLASH_PUPIL_SPIKES.SpikesPerUnit{i} = stim_info.spiking_info.clusters_total_spikes;
%     
%     gunits = find(stim_info.spiking_info.clusters_total_spikes >5);
%     
%     % Flash - Spiking
%     FLASH_PUPIL_SPIKES.Spiking_1REP{i} = stim_info.spiking_info.av_all_fr;
%     FLASH_PUPIL_SPIKES.Spiking_1REP_GU{i} = nanmean(stim_info.spiking_info.clusters_reps_mean(gunits, :));
%     FLASH_PUPIL_SPIKES.Spiking_1REP_STD{i} = stim_info.spiking_info.std_all_fr;
%     
%     FLASH_PUPIL_SPIKES.Spiking_10REPS{i} = nanmean(stim_info.spiking_info.clusters_thru); 
%     FLASH_PUPIL_SPIKES.Spiking_10REPS_GU{i} =  nanmean(stim_info.spiking_info.clusters_thru(gunits, :)); 
%     
%     % Pupil 
%     FLASH_PUPIL_SPIKES.Pupil_1REP{i} = stim_info.eye_info.av_rep_3;
%     FLASH_PUPIL_SPIKES.Pupil_10REP{i} = stim_info.eye_info.r_thru;
%     FLASH_PUPIL_SPIKES.Pupil_STD{i} = stim_info.eye_info.std_rep_3;
%     
    FLASH_PUPIL_SPIKES.Pupil_ALLREPS{i} = stim_info.eye_info.all_rep;

end 


%% Add MOTION INFORMATION

for i = 1:54
    
    ani = FLASH_PUPIL_SPIKES.Animal{i};
    date = FLASH_PUPIL_SPIKES.Date{i};
    test = FLASH_PUPIL_SPIKES.Test{i};
    depth = FLASH_PUPIL_SPIKES.Depth{i};
    
    fname = strcat('MOTION_YN_',ani,'_', date, '_', ani, '_Test', test, '_', depth, 'um.mat');
    load(fname, 'motion_yn')
    
    data_eye = FLASH_PUPIL_SPIKES.Pupil_10REP{i};
    n_sp_eye = numel(data_eye);
    
    tstart = FLASH_PUPIL_SPIKES.AVI_F1{i};
    tend = FLASH_PUPIL_SPIKES.AVI_FN{i};
    
    motion_data_stim = motion_yn(tstart:tend);
    n_mot = numel(motion_data_stim);
    
    % Interpolate the motion YN values from 30FPS to 20kHZ
    val = 1/(n_sp_eye/n_mot);
    movement_interp = interp1(1:1:n_mot, motion_data_stim, 1:val:n_mot);
    
    FLASH_PUPIL_SPIKES.MOVING{i} = movement_interp;
end 

save('220628_FLASH_PUPIL_SPIKES_TABLE_Setd5.mat', 'FLASH_PUPIL_SPIKES')


%% ADD STIMULUS - 1 /0 for On/Off - from stim_info.Timing.all_reps_time. - future.

%% 

% 1 - Plot ONE rep of pupil dilation!
t = numel(stim_info.eye_info.av_rep_3);
figure
rectangle('Position', [0 4 t/4 2], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
hold on 
rectangle('Position', [3*t/4 4 t/4 2], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(stim_info.eye_info.av_rep_3, 'k')

% 2 - Plot 10 reps of pupil dilation. 
t = numel(stim_info.eye_info.r_thru);
figure
rectangle('Position', [0 4 t/4 2], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
hold on 
rectangle('Position', [3*t/4 4 t/4 2], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(stim_info.eye_info.r_thru, 'k')

% 3 - One rep - spiking activity 

% t = numel(stim_info.spiking_info.av_all_fr);
% 8638
% 27966
figure
rectangle('Position', [0 0 9000 0.3], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
hold on 
rectangle('Position', [28000 0 t/4 0.3], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(stim_info.spiking_info.av_all_fr)

% 4 - 10 reps - spiking activity - all units
figure 
plot(nanmean(stim_info.spiking_info.clusters_thru), 'k')

% 5 - Only units with spiking > 5 spikes  - 10 reps - spiking !!!
gunits = find(stim_info.spiking_info.clusters_total_spikes >5);

hold on 
plot(nanmean(stim_info.spiking_info.clusters_thru(gunits, :)), 'r')

% 6 - 1 rep - spiking - only gunits. !!!
hold on 
plot(nanmean(stim_info.spiking_info.clusters_reps_mean(gunits, :)), 'r')

%%
FLASH_PUPIL_SPIKES = a1;

allrows = find(FLASH_PUPIL_SPIKES.Test == "4");
FLASH_PUPIL_SPIKES = FLASH_PUPIL_SPIKES(allrows, :);
%% AVERAGE PER ANIMAL

all_animals = unique(FLASH_PUPIL_SPIKES.Animal);

allWT = [];
allHET = []; 

% figure
for kk = 1:6
    
    anii = string(all_animals{kk});
    all_ani = find(FLASH_PUPIL_SPIKES.Animal == anii & string(FLASH_PUPIL_SPIKES.Test)=="1" | FLASH_PUPIL_SPIKES.Animal == anii & string(FLASH_PUPIL_SPIKES.Test)=="2");

for j = 1:numel(all_ani)
    
    % ALL PUPIL REPS
    data = cell2mat(FLASH_PUPIL_SPIKES.Pupil_ALLREPS(all_ani(j)));
    
    % FIND WHICH REPS THE MICE ARE MOVING> 
    move = cell2mat(FLASH_PUPIL_SPIKES.MOVING(all_ani(j)));
    sz = numel(move);
    binsize = ceil(sz/10);
    edges = [1:binsize:sz];
    edges = [edges, sz];

    % SUM 1/0 DURING REP TIME - count as running if > 0.5 
    move_YN = zeros(1,10);
    for p = 1:10
        move_YN(1,p) = mean(move(edges(p):edges(p+1)));
    end 
    move_YN(move_YN>0.5) = 1;
    
%     % PLOT 
    reps_running = find(move_YN==1);
    reps_stat = find(move_YN == 0);
%    
%     for k = reps_running
%         plot(data(k, :), 'm')
%         hold on 
%     end
%     for k = reps_stat
%         plot(data(k, :), 'k')
%         hold on 
%     end
%    
%     title('GN7788')
    
        g = cell2mat(FLASH_PUPIL_SPIKES.Geno(all_ani(j)));
        
        if g == 1
%             subplot(1,2,1)
%             for k = reps_running
%                 plot(data(k, :), 'k')
%                 hold on
%             end
               allWT = vertcat(allWT, (data(reps_running, 1:36500))); 

        elseif g ==0 
%             subplot(1,2,2)
%             for k = reps_running
%                 plot(data(k, :), 'r')
%                 hold on
%             end

            allHET = vertcat(allHET, (data(reps_running, 1:36500))); 
        end 
        
end

end 

%% MEAN/ SEM - WT / HET - PUPIL - Only while running! 


x = (1:1:36500);

mean_WT = nanmean(allWT);
mean_HET = nanmean(allHET);

semWT = nanstd(allWT)/sqrt(numel(allWT(:,1))) ; %/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

semHET = nanstd(allHET)/sqrt(numel(allHET(:,1))); %/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;


figure
% yyaxis right
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
ax=gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
xlim([-2500 39000])
xticks([38000/4 (3*38000/4)])
% ylim([10 14])
f= gcf;
f.Position = [784   589   226   209];

% WHEN SAVING AS FLAT IMAGE NOT LAYER
set(f, 'Renderer', 'painters');


[p, h] = ranksum(mean(allWT), mean(allHET))
[h, p] = kstest2(mean(allWT), mean(allHET))

%% Max/ Min / Range of pupil size. 

nWT = numel(allWT(:,1))
nHET = numel(allHET(:,1))

%

Wmax = [];
Wmin = [];
Wrng = [];

% WT 
for i1 = 1:nWT
    maxx = max(allWT(i1,:));
    minn = min(allWT(i1,:));
    rangee = range(allWT(i1,:));
    

    Wmax = [Wmax, maxx];
    Wmin = [Wmin, minn];
    Wrng = [Wrng, rangee];
end 

% HET 
Hmax = [];
Hmin = [];
Hrng = [];

for i1 = 1:nHET
    maxx = max(allHET(i1,:));
    minn = min(allHET(i1,:));
    rangee = range(allHET(i1,:));
    

    Hmax = [Hmax, maxx];
    Hmin = [Hmin, minn];
    Hrng = [Hrng, rangee];
end 

% STATS - max
mWT = mean(Wmax)
mHET = mean(Hmax)
[p,h] = ranksum(Wmax, Hmax)

mWT = mean(Wmin)
mHET = mean(Hmin)
[p,h] = ranksum(Wmin, Hmin)

mWT = mean(Wrng)
mHET = mean(Hrng)
[p,h] = ranksum(Wrng, Hrng)

%% 

n_wt = numel(Wrng(1,:));
n_het = numel(Hrng(1,:));

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

col = 'r';
scatter(x1, Wrng,'SizeData', 200, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, Hrng, 'SizeData', 200,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([Wrng, Hrng], [ones(1,23), ones(1,63)*2], 'Color', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
ax.FontSize = 30;
box off
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ax.TickLength = [0.03 0.03];
ax.LineWidth = 2;

f = gcf;
f.Position = [704   207   355   572]; 
ylim([0 5])


%% Spiking to FLASH

allWT = [];
allHET = []; 

for k = 1:height(FLASH_PUPIL_SPIKES)
    data = cell2mat(FLASH_PUPIL_SPIKES.Spiking_1REP_GU(k));
    sz_data = numel(data);
    d = NaN(1, 40000);
    d(1:sz_data) = data;
    
    g = FLASH_PUPIL_SPIKES.Geno{k};
    if g ==1 
        allWT = vertcat(allWT, d);
    else 
        allHET = vertcat(allHET, d);
    end 
end 

x = (1:1:36500);

allWT = allWT(:, 1:36500);
allHET = allHET(:, 1:36500);

mean_WT = nanmean(allWT);
mean_HET = nanmean(allHET);

semWT = nanstd(allWT)/sqrt(numel(allWT(:,1))) ; %/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

semHET = nanstd(allHET)/sqrt(numel(allHET(:,1))); %/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;


figure
% yyaxis right
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
ax=gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
xlim([-2500 39000])
xticks([38000/4 (3*38000/4)])
ylim([10 14])
f= gcf;
f.Position = [784   589   226   209];



