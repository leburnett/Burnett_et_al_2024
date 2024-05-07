% Analysing LFP data - oscillations.
% Created 20/10/14 by Burnett with help from Nardin. 
load('t_amplifier.mat')

files = dir('amplifier*');
numcol = length(t_amplifier);
lfp = zeros(1, numcol);

for i = 1:32
    file_name = files(i).name;
    load(file_name)
    lfp = lfp + amplifier_data_SUB; 
end

lfp = lfp/32;

spectrogram(lfp, 5000, 3000,[5:5:100], 20000, 'yaxis')
spectrogram(lfp, 5000, 3000,[5:2:50], 20000, 'yaxis')





%% Look at CSD using CSD script from Timothy Olsen - MATLAB. 


%%
% Spiking = spikesz

% Pupil diamter
diam = data.dist;

% all_high_spike = find(info_array.n_spikes>50);

% PLOT
figure
subplot(3,1,1)
plot(nanmean(spikesz(:, 1:25000)))

subplot(3,1,2)
plot(diam(1:25000))

subplot(3,1,3)
plot(motion_yn(1:25000)')
ylim([-0.1 1.1])

title('7270-1686')



%% MOTION

files = dir('*.mat');
prop_running = zeros(numel(files), 2);
het_animals = [7269, 7476, 7614, 7790];

files = dir('*.mat');

for j = 1:numel(files)
    
    filee = files(j).name;
    load(filee)
    
    totalf = numel(motion_yn);
    totalmotion = sum(motion_yn);
    
    prop_motion = totalmotion/totalf;
    
    prop_running(j) = prop_motion; 
    
    ani = str2num(filee(end-20:end-17));
    trial = str2num(filee(end-11));
    
    if ismember(ani, het_animals)
        prop_running(j, 2) = 2;
        prop_running(j, 3) = trial;
    else 
        prop_running(j, 2) = 1;
        prop_running(j, 3) = trial;
    end 
    
end 


%% Animal averages - % running over experiments. 

files = dir('*.mat');
prop_running = zeros(numel(files), 2);
het_animals = [7269, 7476, 7614, 7790];
files = dir('*.mat');
all_animals = [7269, 7270, 7476, 7475, 7614, 7616, 7788, 7790];

WTVALS = []; 
HETVALS = []; 

for k = 1:8
    % Set animal to analyse
    ani_target = all_animals(k);
    ani_vals = [];
%     tf = 0;
    
    for j=1:numel(files)
        filee = files(j).name;
        ani = str2num(filee(end-20:end-17));
        trial = str2num(filee(end-11));
        % If file correspond to the target animal and the trial is 1-4. 
        if ani == ani_target && trial<5%tf==0;
            load(filee)
            totalf = numel(motion_yn);
            totalmotion = sum(motion_yn);
            prop_motion = totalmotion/totalf;
            ani_vals = [ani_vals, prop_motion];
            tf =1;
        end
    end
    
    ani_mean = nanmean(ani_vals);
    
    if ismember(ani_target, het_animals)
        HETVALS = [HETVALS, ani_mean];
    else 
        WTVALS = [WTVALS, ani_mean];
    end 
%     tf =0; 
    
end

%% BAR CHART + ANIMAL POINTS - % OF TIME SPENT RUNNING
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
% 

mWT = nanmean(WTVALS)
mHET = nanmean(HETVALS)
% [p,h] = ranksum(WTVALS, HETVALS)
% [h,p] = ttest2(WTVALS, HETVALS)




%%

wtvals = zeros(2,4);
hetvals = zeros(2,4);

for k = 1:4
    
    v = find(prop_running(:,3)==k & prop_running(:,2)==1);
    v2 = find(prop_running(:,3)==k & prop_running(:,2)==2);
    
    vals = prop_running(v,1);
    vals2 = prop_running(v2,1);
    
    wtvals(1,k) = nanmean(vals);
    hetvals(1,k) =  nanmean(vals2);
    
    wtvals(2,k) = nanstd(vals)/sqrt(numel(vals));
    hetvals(2,k) = nanstd(vals2)/sqrt(numel(vals2));
    
end 


% errorbar(wtvals(1,:), wtvals(2,:), 'k')
% hold on 
% errorbar(hetvals(1,:), hetvals(2,:), 'r')
% 

%% Mean SEM plot across trials. 
mean_WT = wtvals(1,:);
mean_HET = hetvals(1,:);

semWT = wtvals(2,:);
semHET = hetvals(2,:);

y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

y3 = mean_HET+semHET;
y4 = mean_HET-semHET;


figure
plot(1:1:4, y1, 'w-')
hold on
plot(1:1:4, y2, 'w-')
patch([1:1:4 fliplr(1:1:4)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT', 'k-', 'LineWidth', 1.3)

plot(1:1:4, y3, 'w-')
hold on 
plot(1:1:4, y4, 'w-')
patch([1:1:4 fliplr(1:1:4)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET', '-','Color', col, 'LineWidth', 1.3)

xlim([0.5 4.5])
ylim([0 0.5])
box off
ax=gca;
ax.TickDir = 'out'; 
yticks([0:0.1:0.5])
xticks([1:1:4])
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;


%%  Plot the % running for each trial as separate lines on plot. 

files = dir('*.mat');
prop_running = zeros(numel(files), 2);
het_animals = [7269, 7476, 7614, 7790];
files = dir('*.mat');
all_animals = [7269, 7270, 7476, 7475, 7614, 7616, 7788, 7790];

WTVALS = []; 
HETVALS = []; 

for k = 7:8
    % Set animal to analyse
    ani_target = all_animals(k);
    ani_vals = [];
 
    for j=1:numel(files)
        filee = files(j).name;
        ani = str2num(filee(end-20:end-17));
        trial = str2num(filee(end-11));
        % If file correspond to the target animal and the trial is 1-4. 
        if ani == ani_target && trial<5
            load(filee)
            totalf = numel(motion_yn);
            totalmotion = sum(motion_yn);
            prop_motion = totalmotion/totalf;
            ani_vals = [ani_vals, prop_motion];
        end
    end
    
%     ani_mean = nanmean(ani_vals);
    
    if ismember(ani_target, het_animals)
        HETVALS = vertcat(HETVALS, ani_vals);
    else 
        WTVALS = vertcat(WTVALS, ani_vals);
    end 
    
end

v = 0.6;

SEMWT = nanstd(WTVALS)/sqrt(numel(3));
SEMHET = nanstd(HETVALS)/sqrt(numel(3));

mean_WT = nanmean(WTVALS);
mean_HET = nanmean(HETVALS);


y1 = mean_WT+SEMWT;
y2 = mean_WT-SEMWT;

y3 = mean_HET+SEMHET;
y4 = mean_HET-SEMHET;


figure
plot(1:1:4, y1, 'w-')
hold on
plot(1:1:4, y2, 'w-')
patch([1:1:4 fliplr(1:1:4)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
plot(mean_WT', 'k-', 'LineWidth', 1.3)

plot(1:1:4, y3, 'w-')
hold on 
plot(1:1:4, y4, 'w-')
patch([1:1:4 fliplr(1:1:4)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')
plot(mean_HET', '-','Color', col, 'LineWidth', 1.3)


plot(1:1:4, WTVALS(1,:), 'Color', [v v v], 'LineWidth', 1)
hold on 
plot(1:1:4, WTVALS(2,:), 'Color', [v v v], 'LineWidth', 1)
plot(1:1:4, WTVALS(3,:), 'Color', [v v v], 'LineWidth', 1)
plot(1:1:4, HETVALS(1,:), 'Color', [1 v v],  'LineWidth', 1)
plot(1:1:4, HETVALS(2,:), 'Color', [1 v v], 'LineWidth', 1)
plot(1:1:4, HETVALS(3,:), 'Color', [1 v v], 'LineWidth', 1)

plot(1:1:4, nanmean(WTVALS), 'k', 'LineWidth', 2.5)
plot(1:1:4, nanmean(HETVALS), 'r', 'LineWidth', 2.5)

ylim([0 0.6])
xlim([0.5 4.5])
box off
ax=gca;
ax.TickDir = 'out'; 
yticks([0:0.1:0.5])
xticks([1:1:4])
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;






%% Using the data from Olga - 'stim_info' structure. 

% data = spikes to flash - all reps
data_spikes = mean(stim_info.spiking_info.clusters_thru); 

data_eye = stim_info.eye_info.r_thru;

% Start and end AVI frames for motion 
tstart = stim_info.timing.avi_first_frame; 
tend = stim_info.timing.avi_last_frame; 

figure;
subplot(2,1,1)
yyaxis left
plot(data_eye)
% subplot(3,1,2)
yyaxis right
plot(data_spikes)
subplot(2,1,2)
plot(motion_yn(tstart:tend))
ylim([-0.1 1.1])
% sgtitle('7616-test2-run2')
sgtitle('new')
f=gcf;
f.Position =[86   698   378   367]; 

% figure; plot(mean(spikesz(:, tstart/(100/60):tend/(100/60))))

%%
% 
% figure
% plot(data_spikes)
% 
% rep_time = stim_info.timing.all_reps_time;
% 
% hold on 
% for k = 1:10
%     xval =rep_time(k,2);
%     plot([xval xval], [0 0.5], 'k:')
%     
%     xval2 =rep_time(k,3);
%     plot([xval2 xval2], [0 0.5], 'k:')
%     hold on 
% end 

%%

% Camera = 30FPS
% Spiking = 20Hz
% Eye was interpolated from 30fps into 20Hz.

% data = spikes to flash - all reps
data_spikes = mean(stim_info.spiking_info.clusters_thru); 

data_eye = stim_info.eye_info.r_thru;

% Start and end AVI frames for motion 
tstart = stim_info.timing.avi_first_frame; 
tend = stim_info.timing.avi_last_frame; 


motion_data_stim = motion_yn(tstart:tend); 
n_mot = numel(motion_data_stim);
n_sp_eye = numel(data_eye);
val = 1/(n_sp_eye/n_mot);

movement_interp = interp1(1:1:n_mot, motion_data_stim, 1:val:n_mot);

% Find where mouse is moving and put shaded area in background of plot!
M1 = max(data_eye);
M2 = max(data_spikes);

if M1>M2
    maxval = M1;
else
    maxval = M2;
end

moving = find(movement_interp ==1);

if ~isempty(moving)
    
    if numel(moving(1,:)) == numel(movement_interp(1,:))
        figure
        rectangle('Position', [1, 0, numel(moving(1,:)), maxval], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
        hold on
        yyaxis left
        plot(data_eye)
        yyaxis right
        plot(data_spikes)
        box off
    else
        
        moving(2,2:end) = diff(moving(1,:));
        blocks = find(moving(2,:)>1);
        
        n_blocks = numel(blocks);
        block_vals = ones(n_blocks+1, 2);
        
        for k = 1:n_blocks
            
            if k ==1
                block_vals(1,2) = moving(1, blocks(k)-1);
            end
            
            block_vals(k+1, 1) = moving(1, blocks(k));
            
            if k ~= n_blocks
                block_vals(k+1, 2) = moving(1, blocks(k)-1);
            else
                block_vals(k+1, 2) = moving(1, end);
            end
        end
        
        figure
        for i = 1:n_blocks+1
            rectangle('Position', [block_vals(i,1), 0, (block_vals(i,2)-block_vals(i,1)), maxval], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none')
            hold on
        end
        yyaxis left
        plot(data_eye)
        % subplot(3,1,2)
        yyaxis right
        plot(data_spikes)
        box off
        
    end
    
elseif isempty(moving) 
    
    figure
    yyaxis left
    plot(data_eye)
    % subplot(3,1,2)
    yyaxis right
    plot(data_spikes)
    box off
    
end

f = gcf;
f.Position =  [922   478   623   206];
title('7269-T3-R2')









% data = spikes to flash - all reps
data_spikes = mean(stim_info.spiking_info.clusters_thru); 
data_eye = stim_info.eye_info.r_thru;
% Start and end AVI frames for motion 
tstart = stim_info.timing.avi_first_frame; 
tend = stim_info.timing.avi_last_frame; 
   figure
    yyaxis left
    plot(data_eye)
    % subplot(3,1,2)
    yyaxis right
    plot(data_spikes)
    box off























%%  CSD - Average over 9 flashes
% 40000
% dd = CSDoutput'; 

data = zeros(30, 40001, 9);

for i = 1:9
    x1 = 135000+(40000*(i-1));
    st = x1-20000;
    nd = x1+20000;
    
    data(:, :, i) = dd(2:31, st:nd); 
end 

q = mean(data, 3);
figure; imagesc(q); caxis([-100000000 100000000]); colormap(redblue)
box off
ax = gca;
ax.XAxis.Visible = 'on';



