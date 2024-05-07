% Combine GRATINGS analysis arrays:
% Cretaed by Burnett - 22/02/22

% LOAD table with information about the depth of the sSC in different animals. 
load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Setd5_Rec_Depth_Info.mat', 'animal_depth_info_table')

% ANALYSIS RESULTS FILES STORED HERE:
% "/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/SPIKESS/Gratings/Grating_Per_Rec"

%% For all the analysis files want to:

% - Combine the ds_idx arrays - with info about the DS ofeach unit recorded. 
%  - - - Need to ammend 'depth' in column 10 of these neurons 
%       a = probe depth - (800 - depth) 
%       b = a - sSC depth. 

%  sp_pt = spikes per time - across all different gratings.  - for example
%  traces. 

% Also make combined table with Col 1 = rep1 , col2 = rep2, col3 = rep3;
% Used for individual trials. 
het_animals = [7269, 7476, 7614];

Grating_Table = table(); 

files = dir('*GRATING*');
nfiles =numel(files);

for i = 1:nfiles
    
    fname = files(i).name;
    load(fname)
    
    % Ammend the depth of each cell! 
    animal_num = fname(8:11);
    depth = fname(15:18);
    datee = fname(1:6);

%     row = find(animal_depth_info_table.Animal == str2num(animal_num));
    row = find(CSDtable.Ani == str2num(animal_num) &  CSDtable.RecDepth == str2num(depth) &  CSDtable.Date == str2num(datee));
    
%     sSC_depth = animal_depth_info_table.sSC_Depth(row); % depth of beginning of sSC
      sSC_depth = CSDtable.SGSDEPTH(row); % depth of beginning of sSC
      
    probe = ds_ind(1,4); %probe depth
    n_gcl = numel(ds_ind(:,1)); %number of units in this recording
    
    for jj = 1:n_gcl
        d = probe - (800 - ds_ind(jj, 10));
        d2 = sSC_depth - d;
        ds_ind(jj,11) = d2; 
        
        if ismember(str2num(animal_num), het_animals)
            ds_ind(:, 12) = ones(n_gcl, 1)*2;
        else
            ds_ind(:, 12) = ones(n_gcl, 1); 
        end 
        % Total number of spikes over the stimulus  - used to find low
        % firing cells. 
        ds_ind(jj, 13) = sum(spikes_per_time(jj, :));

    end 
   
    % Make TABLE
    
    Date = ds_ind(:,1);
    Ani = ds_ind(:,2);
    Exp = ds_ind(:,3);
    Geno = ds_ind(:,12);
    ProbeDepth = ds_ind(:,4);
    P_Ang = ds_ind(:,5);
    P_Sp = ds_ind(:,6);
    NP_Ang = ds_ind(:,7);
    NP_Sp = ds_ind(:,8);
    DS_Ind = ds_ind(:,9);
    Depth = ds_ind(:,11);
    TotalSpikes = ds_ind(:,13);
    
    SpPT = sp_pt;
    SpPT_RAW = sp_pt_raw;
    SpPS = sp_ps; 
    SpPS_RAW = sp_ps_raw;
    
    SpPT1 = norm_spikes_per_time1;
    SpPS1 = norm_spikes_per_stim1;
    
    if exist('spikes_per_time2', 'var')
        SpPT2 = norm_spikes_per_time2;
        SpPT3 = norm_spikes_per_time;
        SpPS2 = norm_spikes_per_stim2;
        SpPS3 = norm_spikes_per_stim;
        
    elseif ~exist('spikes_per_time2', 'var')
        SpPT2 = norm_spikes_per_time;
        SpPT3 = NaN(n_gcl, 2400);
        SpPS2 = norm_spikes_per_stim;
        SpPS3 = NaN(n_gcl, 8);
    end
    
    tbl = table(Date, Ani, Exp, Geno, ProbeDepth, TotalSpikes, P_Ang, P_Sp, NP_Ang, NP_Sp, DS_Ind, Depth, SpPT, SpPT_RAW, SpPT1, SpPT2, SpPT3, SpPS, SpPS_RAW, SpPS1, SpPS2, SpPS3);
    
    Grating_Table = vertcat(Grating_Table, tbl);
    
end 

% save('220830_Setd5_GRATINGS_table_newCSD_depth.mat', 'Grating_Table', 'CSDtable');

%% GENERAL OVERVIEW - WT/HET

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth<400 & Grating_Table.Depth>-100 & Grating_Table.TotalSpikes>2);  % Flash_Table.P_VALUE <0.05);
allHET = find(Grating_Table.Geno == 2  & Grating_Table.Depth<400 & Grating_Table.Depth>-100 & Grating_Table.TotalSpikes>2); %  & Flash_Table.P_VALUE <0.05);

% % DS Index
% 
mWT = nanmean(Grating_Table.DS_Ind(allWT))
mHET = nanmean(Grating_Table.DS_Ind(allHET))

figure
histogram(Grating_Table.DS_Ind(allWT), 0:0.1:1, 'Normalization', 'pdf')
hold on 
histogram(Grating_Table.DS_Ind(allHET), 0:0.1:1, 'Normalization', 'pdf')

[h, p] = kstest2(Grating_Table.DS_Ind(allWT), Grating_Table.DS_Ind(allHET))


% P_Ang

mWT = nanmean(Grating_Table.P_Ang(allWT));
mHET = nanmean(Grating_Table.P_Ang(allHET));

figure
histogram(Grating_Table.P_Ang(allWT), 0:45:360, 'Normalization', 'pdf')
hold on 
histogram(Grating_Table.P_Ang(allHET), 0:45:360, 'Normalization', 'pdf')

[h, p] = kstest2(Grating_Table.P_Ang(allWT), Grating_Table.P_Ang(allHET))


% P_Sp

mWT = nanmean(Grating_Table.P_Sp(allWT));
mHET = nanmean(Grating_Table.P_Sp(allHET));

figure
histogram(Grating_Table.P_Sp(allWT), 0:10:700, 'Normalization', 'pdf')
hold on 
histogram(Grating_Table.P_Sp(allHET), 0:10:700, 'Normalization', 'pdf')

[h, p] = kstest2(Grating_Table.P_Sp(allWT), Grating_Table.P_Sp(allHET))


%% At correct depth. 

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>-400 & Grating_Table.Depth <0 & Grating_Table.TotalSpikes >5  & Grating_Table.DS_Ind >0.3);
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>-400 & Grating_Table.Depth <0 & Grating_Table.TotalSpikes >5 & Grating_Table.DS_Ind >0.3);


allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>-100 & Grating_Table.Depth <300 & Grating_Table.TotalSpikes >20);
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>-100 & Grating_Table.Depth <300 & Grating_Table.TotalSpikes >20);


%% SORT BY KMEANS

data_WT = Grating_Table{allWT, 13};
data_HET = Grating_Table{allHET, 13};

n_k = 20; 
data_WT(:,2401) = kmeans(data_WT(:, 1:2400), n_k);
data_WT = sortrows(data_WT, 2401);
% figure; imagesc(data_WT(:, 1:2400))

% % Plots of the means of each cluster; 
%  figure
% for j = 1:n_k 
% all_type = find(data_WT(:, 2401)== j);
% subplot(n_k,1,j); plot(nanmean(data_WT(all_type, 1:2400)), 'k'); box off
% end 

%% Plot the mean
figure
plot(smooth(nanmean(Grating_Table{allWT, 13})))
hold on 
plot(smooth(nanmean(Grating_Table{allHET, 13})))



%% Summary plots of spiking over time - per trial and polar plot of the normalised spiking per stim - WT

n_gWT = numel(allWT);

xvls = 0:5*60:40*60;

angls = 0:45:315;
a1 = angls;
angls(9) = angls(1);
angls = deg2rad(angls); % NEEDS TO BE IN RADIANS FOR POLAR PLOT!!

for i2 = 1:n_gWT
    s1 = Grating_Table{allWT(i2), 15};
    s2 = Grating_Table{allWT(i2), 16};
    s3 = Grating_Table{allWT(i2), 17};
    avs = Grating_Table{allWT(i2), 13};
    
    figure
    subplot(1,5,1:4)
    plot(smooth(s1), 'Color', [0 0 1]);
    hold on 
    plot(smooth(s2), 'Color', [0.7 0.7 1]);
    plot(smooth(s3), 'Color', [0.5 0.5 1]);
%     plot(smooth(avs), 'k');
    
%     mv(1) = max((s1));
%     mv(2) = max((s2));
%     mv(3) = max((s3));
%     mvv = max(mv);
    
    for i = 2:9
        plot([xvls(i) xvls(i)], [0 1], 'r:', 'LineWidth', 1.3)
        hold on
    end
    xticks(2.5*60:5*60:37.5*60)
    xticklabels({a1})
    box off
    ax1 = gca;
    ax1.TickDir = 'out';
    ylim([0 1])
    
    subplot(1,5,5)
    
    st1 = Grating_Table{allWT(i2), 20};
    st2 = Grating_Table{allWT(i2), 21};
    st3 = Grating_Table{allWT(i2), 22};
    avst = Grating_Table{allWT(i2), 18};
    
    rho = st1;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.8 0.8 0.8])
    % REP2
    hold on
    rho = st2;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    % REP3
    hold on
    rho = st3;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    
    % MEAN 
    rho = avst;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', 'k', 'LineWidth', 1.3)
    ax = gca;
    ax.RTickLabel = {''};
    ax.ThetaTickLabel = {''};
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    
    sgtitle(string(i2));
    
    f= gcf;
    f.Position = [440   629 1085  169];
    
    x = input(''); 
    close
end 



%% Summary plots of spiking over time - per trial and polar plot of the normalised spiking per stim - HET

n_gWT = numel(allHET);

xvls = 0:5*60:40*60;

angls = 0:45:315;
a1 = angls;
angls(9) = angls(1);
angls = deg2rad(angls); % NEEDS TO BE IN RADIANS FOR POLAR PLOT!!

for i2 = 1:n_gWT
    s1 = Grating_Table{allHET(i2), 15};
    s2 = Grating_Table{allHET(i2), 16};
    s3 = Grating_Table{allHET(i2), 17};
    avs = Grating_Table{allHET(i2), 13};
    
    figure
    subplot(1,5,1:4)
    plot(smooth(s1), 'Color', [0.7 0.7 0.7]);
    hold on 
    plot(smooth(s2), 'Color', [0.7 0.7 0.7]);
    plot(smooth(s3), 'Color', [0.7 0.7 0.7]);
    plot(smooth(avs), 'r');
    
%     mv(1) = max((s1));
%     mv(2) = max((s2));
%     mv(3) = max((s3));
%     mvv = max(mv);
    
    for i = 2:9
        plot([xvls(i) xvls(i)], [0 1], 'k:', 'LineWidth', 1.3)
        hold on
    end
    xticks(2.5*60:5*60:37.5*60)
    xticklabels({a1})
    box off
    ax1 = gca;
    ax1.TickDir = 'out';
    ylim([0 1])
    
    subplot(1,5,5)
    
    st1 = Grating_Table{allHET(i2), 20};
    st2 = Grating_Table{allHET(i2), 21};
    st3 = Grating_Table{allHET(i2), 22};
    avst = Grating_Table{allHET(i2), 18};
    
    rho = st1;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.8 0.8 0.8])
    % REP2
    hold on
    rho = st2;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    % REP3
    hold on
    rho = st3;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    
    % MEAN 
    rho = avst;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', 'r', 'LineWidth', 1.3)
    ax = gca;
    ax.RTickLabel = {''};
    ax.ThetaTickLabel = {''};
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    
    sgtitle(string(i2));
    
    f= gcf;
    f.Position = [440   629 1085  169];
    
    x = input(''); 
    close
end


%% P3 - Only look at cells within a particular depth / organise cells according to depth. 

data_WT = Grating_Table{allWT, 13};
data_WT(:,2401) = Grating_Table{allWT, 12}; %DEPTH
data_WT = sortrows(data_WT, 2401);
figure; imagesc(data_WT(:, 1:2400))


data_HET = Grating_Table{allHET, 13};
data_HET(:,2401) = Grating_Table{allHET, 12}; %DEPTH
data_HET = sortrows(data_HET, 2401);
figure; imagesc(data_HET(:, 1:2400))


%% Plot DS_Idx versus depth. 

figure
for i = 1:height(Grating_Table)
    xval = Grating_Table.DS_Ind(i);
    yval = Grating_Table.Depth(i);
    
    geno = Grating_Table.Geno(i);
    TOTSP = Grating_Table.TotalSpikes(i);
    
    if geno == 1
        col = 'k';
    else 
        col = 'r';
    end 
    
    if TOTSP<30
        marker = '.';
        sz = 75; 
    elseif TOTSP>=30 
        marker = 'o';
        sz = 50;
    end 
    
        scatter(xval, yval,sz, col,  'Marker', marker);
        hold on
    
%     scatter(xval, yval,sz, col,  'Marker', marker);
%     hold on

end 

axis([-0.1 1.1 -1450 600])


%% Distributions!!

% allWT = find(Grating_Table.Geno == 1 & Grating_Table.TotalSpikes >30);
% allHET = find(Grating_Table.Geno == 2 & Grating_Table.TotalSpikes >30);

dstart = 400;
dend = -100;

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart  & Grating_Table.TotalSpikes >2);
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2);

%% PLOTS  - 1 Spiking to Pref Direction

data_WT = Grating_Table.P_Sp(allWT);
data_HET = Grating_Table.P_Sp(allHET);

rng = 0:10:300;

figure; histogram(data_WT, rng, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
hold on 
histogram(data_HET, rng, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
box off

%% PLOT  - 2 Angle of  Pref Direction

data_WT = Grating_Table.P_Ang(allWT);
data_HET = Grating_Table.P_Ang(allHET);

rng = 0:45:360;

figure; histogram(data_WT, rng, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
hold on 
histogram(data_HET, rng, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
box off





%% PLOT  - 3 DS_Index

dstart = 400;
dend = -100;

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2); %& Grating_Table.TotalSpikes >10
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2);

data_WT = Grating_Table.DS_Ind(allWT);
data_HET = Grating_Table.DS_Ind(allHET);

rng = 0:0.075:1;

figure; histogram(data_WT, rng, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
hold on 
histogram(data_HET, rng, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
box off
% xlabel('DSi')
% ylabel('PDF')
ax = gca;
% ax.FontSize = 16;
ax.TickDir = 'out'; 
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;

mWT = nanmean(data_WT);
mHET = nanmean(data_HET);
medWT = nanmedian(data_WT)
medHET = nanmedian(data_HET)

hold on 
plot([medWT medWT], [0 2.2], 'k', 'LineWidth', 1.4)
plot([medHET medHET], [0 2.2], 'r', 'LineWidth', 1.4)

ylim([0 3])

[h,p] = kstest2(data_WT, data_HET)
[h,p] = ranksum(data_WT, data_HET)

% 8.73x10-4; < 10 - less than 0 
% 0.0806 % no specs on firing rate - less than 0 
%  0.9737 : >10 < -350 
% 0.0435 
% 0.0523 : -700 to -1050
% 0.3810 : - 1050 to -2000

f = gcf; 
f.Position = [349   406   359   182]; 
f.Renderer = 'painters';


%% TOP BOTTOM SUBPLOT OF MEDIAN DSI

dstart = 400;
dend = -100;

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2); 
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2);

data_WT = Grating_Table.DS_Ind(allWT);
data_HET = Grating_Table.DS_Ind(allHET);

rng = 0:0.1:1;

mWT = nanmean(data_WT);
mHET = nanmean(data_HET);
medWT = nanmedian(data_WT);
medHET = nanmedian(data_HET);

figure; subplot(2,1,1); histogram(data_WT, rng, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
hold on 
plot([medWT medWT], [0 3], 'k', 'LineWidth', 1.4);
box off
ylabel('PDF')
xlabel('DSi')
ax = gca;
ax.FontSize = 16;
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ylim([0 2.5])


subplot(2,1,2);
histogram(data_HET, rng, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf'); hold on; plot([medHET medHET], [0 4], 'r', 'LineWidth', 1.4);
box off
xlabel('DSi')
ylabel('PDF')
ax = gca;
ax.FontSize = 16;
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ylim([0 2.5])

f = gcf; 
f.Position = [720   587   376   355]; 



[h,p] = kstest2(data_WT, data_HET)


%% DSi versus Total # of Spikes - negative correlation? 

figure
for i = 1:height(Grating_Table)
    xval = Grating_Table.TotalSpikes(i);
    yval = Grating_Table.DS_Ind(i);
    
    geno = Grating_Table.Geno(i);
    tsp = Grating_Table.TotalSpikes(i);
    
    if geno == 1
        col = 'k';
    else 
        col = 'r';
    end 
    
%     if tsp <30
%         marker = '.';
%         sz = 75; 
%     elseif tsp>= 30 
        marker = 'o';
        sz = 50;
        scatter(xval, yval,sz, col,  'Marker', marker);
        hold on
%     end 

end 

xlabel('Total Spikes')
ylabel('DSi')
ax = gca;
box off
ax.FontSize = 16;
ax.LineWidth = 1.2;
ax.TickDir = 'out';

%% PIE CHARTS?
% 1 - Make a pie chart of the total # of WT/HET cells that have a DSi of > 0.3
% 2 - Make a pie chart of the cells that are in the sSC depths - which angle the are selective for. 
 
% pie charts require an array with the numbers corresponding to the size of
% the pie pieces. Want to have a figure with subplots for WT/HET % 
% 
% allWT = find(Grating_Table.Geno == 1);
% allHET = find(Grating_Table.Geno == 2);
% 
% allWT = find(Grating_Table.Geno == 1 & Grating_Table.TotalSpikes >10);
% allHET = find(Grating_Table.Geno == 2 & Grating_Table.TotalSpikes >10);

dstart = 0;
dend =  -2000;

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart  & Grating_Table.TotalSpikes >20);
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >20);

vals(1, 1:10) = histcounts(Grating_Table.DS_Ind(allWT), 0:0.1:1)/numel(allWT);
vals(2, 1:10) = histcounts(Grating_Table.DS_Ind(allHET), 0:0.1:1)/numel(allHET);

nnanwt = isnan(Grating_Table.DS_Ind(allWT));
nnanhet = isnan(Grating_Table.DS_Ind(allHET));

vals(1,11) = numel(find(nnanwt==1))/numel(allWT); % NaNs
vals(2,11) = numel(find(nnanhet==1))/numel(allHET); % NaNs

lb = {'<0.1', '<0.2', '<0.3', '<0.4', '<0.5', '<0.6', '<0.7', '<0.8', '<0.9', '<1', 'NaN'};

X = vals(1,:);
ax1 = subplot(1,2,1);
pie(ax1,X, lb)
title(ax1,'WT');

Y = vals(2,:);
ax2 = subplot(1,2,2);
pie(ax2,Y, lb)
title(ax2,'HET');


%% Polar plot of the AVERAGE direction preferrerd by WTs/ HETs

dstart = 400;
dend =  -100;

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart  & Grating_Table.TotalSpikes >2 & Grating_Table.DS_Ind >0.3);
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2 & Grating_Table.DS_Ind >0.3);

angs = 0:45:315;
angs(9) = angs(1)%
angs = deg2rad(angs);

figure
subplot(1,2,1);
rho = histcounts(Grating_Table.P_Ang(allWT), 0:45:360);
rho(9) = rho(1); 
polarplot(angs,rho, 'Color', 'k', 'LineWidth', 1.2)
ax = gca;
% ax.RTickLabel = {''};
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

subplot(1,2,2);
rho = histcounts(Grating_Table.P_Ang(allHET), 0:45:360);
rho(9) = rho(1); 
polarplot(angs,rho, 'Color', 'r', 'LineWidth', 1.2)
ax2 = gca;
% ax2.RTickLabel = {''};
ax2.ThetaZeroLocation = 'top';
ax2.ThetaDir = 'clockwise';

% sgtitle('0 to -400; >10; DS >0.3')

f = gcf;
f.Position = [323   612   475   290];




%% VIOLIN

DATA = [Grating_Table.DS_Ind, Grating_Table.Geno];

% Remove cells between 0 and -400. % Remove cells with firing < 10
all_depths = find(Grating_Table.Depth > 0 | Grating_Table.TotalSpikes <10);
DATA(all_depths, :) = [];

figure; v = violinplot(DATA(:,1), DATA(:,2), 'BoxColor', [0 0 0], 'ShowData', true);
v(1).ViolinColor = [0.6 0.6 0.6];
v(2).ViolinColor = [1 0 0];

ax = gca;
box off
ax.TickDir = 'out'; 
ax.LineWidth = 1;
ax.FontSize = 16;
xticks([1,2])
xticklabels({'Setd5^+^/^-', 'Setd5^-^/^-'})
xlim([0.5 2.5])
ylabel('DSi')
ax.TickLength = [0.015 0.015];
ylim([-0.1 1.1])



%% Look through DSI polar plots for cells and title = DSI

dstart = 400;
dend =  -100;

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >20);
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >20);

figure
for i2 = 76:100
    
    subplot(5,5,i2-75)
    avst = Grating_Table{allWT(i2), 18};
    rho = avst;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', 'k', 'LineWidth', 2.5)
    ax = gca;
    ax.RTickLabel = {''};
    ax.ThetaTickLabel = {''};
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    dsi = Grating_Table{allWT(i2), 11};
    title(strcat(string(dsi), '-', string(i2)))
    
end 


%% NEW SUBPLOT - WT

n_gWT = numel(allWT);

xvls = 0:5*60:40*60;

angls = 0:45:315;
a1 = angls;
angls(9) = angls(1);
angls = deg2rad(angls); % NEEDS TO BE IN RADIANS FOR POLAR PLOT!!

% col1 = [0 0.5 0]; 
% col2 = [0.6 0  0.6];
% col3 = [1 0.5 0 ];

col1 = [0.7 0.7 0.7];
col2 = col1;
col3 = col1; 

v = 300; 
ymax = 0.7;
lw = 1;
lw2 = 1.3;

% i2 = 12; %158; 
% WT - 130, 158


%% WT EXAMPLES TRACES.

    i2 = 96; %
    
    dsi = Grating_Table.DS_Ind(allWT(i2));

    close
    
    s1 = Grating_Table{allWT(i2), 15};
    s2 = Grating_Table{allWT(i2), 16};
    s3 = Grating_Table{allWT(i2), 17};
    avs = Grating_Table{allWT(i2), 13};
    
    figure
    
     % REP1
    subplot(1,9,1)
    plot(smooth(s1(1:v)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((1:v))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((1:v))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((1:v))), 'Color', 'k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP2
    subplot(1,9,2)
    plot(smooth(s1(v+1:v*2)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v+1:v*2))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v+1:v*2))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v+1:v*2))), 'Color', 'k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    
    % REP3
    subplot(1,9,3)
    plot(smooth(s1(v*2+1:v*3)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*2+1:v*3))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*2+1:v*3))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*2+1:v*3))), 'Color', 'k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP4
    subplot(1,9,4)
    plot(smooth(s1(v*3+1:v*4)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*3+1:v*4))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*3+1:v*4))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*3+1:v*4))), 'Color', 'k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP5
    subplot(1,9,5)
    plot(smooth(s1(v*4+1:v*5)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*4+1:v*5))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*4+1:v*5))), 'Color',col3, 'LineWidth', lw);
    plot(smooth(avs((v*4+1:v*5))), 'Color','k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP6
    subplot(1,9,6)
    plot(smooth(s1(v*5+1:v*6)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*5+1:v*6))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*5+1:v*6))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*5+1:v*6))), 'Color', 'k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP7
    subplot(1,9,7)
    plot(smooth(s1(v*6+1:v*7)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*6+1:v*7))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*6+1:v*7))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*6+1:v*7))), 'Color', 'k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP8
    subplot(1,9,8)
    plot(smooth(s1(v*7+1:v*8)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*7+1:v*8))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*7+1:v*8))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*7+1:v*8))), 'Color', 'k', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    subplot(1,9,9)
    
    st1 = Grating_Table{allWT(i2), 20};
    st2 = Grating_Table{allWT(i2), 21};
    st3 = Grating_Table{allWT(i2), 22};
    avst = Grating_Table{allWT(i2), 18};
    
    rho = st1;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    % REP2
    hold on
    rho = st2;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    % REP3
    hold on
    rho = st3;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    
    % MEAN 
    
    rho = avst;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', 'k', 'LineWidth', 2.5)
    ax = gca;
    ax.RTickLabel = {''};
    ax.ThetaTickLabel = {''};
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
   
    
    f= gcf;
    f.Position =  [131   214   1662   158]; %[1606  184  1860  192];  %[1606  124 2216  252];
    

%% Look through DSI polar plots for cells and title = DSI - HET

dstart = 400;
dend =  -100;

allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >20);

figure
for i2 = 151:177
    
    subplot(5,5,i2-150)
    avst = Grating_Table{allHET(i2), 18};
    rho = avst;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', 'r', 'LineWidth', 2.5)
    ax = gca;
    ax.RTickLabel = {''};
    ax.ThetaTickLabel = {''};
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    dsi = Grating_Table{allHET(i2), 11};
    title(strcat(string(dsi), '-', string(i2)))
    
end 
    
  f = gcf;
  f.Position = [ 89   131   852   661]; 
    
    %%  NEW SUBPLOT - HET

angls = 0:45:315;
a1 = angls;
angls(9) = angls(1);
angls = deg2rad(angls); % NEEDS TO BE IN RADIANS FOR POLAR PLOT!!

% i2 = 5,8,30,34,35,40, 88, 107, 109, 
% NDS - 53, 112,161

v = 300; 
ymax = 0.7;
lw = 1.0;
lw2 = 1.3;

% col1 = [0 0.5 0]; 
% col2 = [0.6 0  0.6];
% col3 = [1 0.5 0 ];

%% HET EXAMPLES TRACES.

    i2 = 174; %
    
    dsi = Grating_Table.DS_Ind(allHET(i2))

    close
    
    s1 = Grating_Table{allHET(i2), 15};
    s2 = Grating_Table{allHET(i2), 16};
    s3 = Grating_Table{allHET(i2), 17};
    avs = Grating_Table{allHET(i2), 13};
    
    
    figure
    
    % REP1
    subplot(1,9,1)
    plot(smooth(s1(1:v)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((1:v))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((1:v))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((1:v))), 'Color', 'r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP2
    subplot(1,9,2)
    plot(smooth(s1(v+1:v*2)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v+1:v*2))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v+1:v*2))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v+1:v*2))), 'Color', 'r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    
    % REP3
    subplot(1,9,3)
    plot(smooth(s1(v*2+1:v*3)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*2+1:v*3))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*2+1:v*3))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*2+1:v*3))), 'Color', 'r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP4
    subplot(1,9,4)
    plot(smooth(s1(v*3+1:v*4)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*3+1:v*4))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*3+1:v*4))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*3+1:v*4))), 'Color', 'r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP5
    subplot(1,9,5)
    plot(smooth(s1(v*4+1:v*5)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*4+1:v*5))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*4+1:v*5))), 'Color',col3, 'LineWidth', lw);
    plot(smooth(avs((v*4+1:v*5))), 'Color','r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP6
    subplot(1,9,6)
    plot(smooth(s1(v*5+1:v*6)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*5+1:v*6))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*5+1:v*6))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*5+1:v*6))), 'Color', 'r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP7
    subplot(1,9,7)
    plot(smooth(s1(v*6+1:v*7)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*6+1:v*7))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*6+1:v*7))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*6+1:v*7))), 'Color', 'r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    % REP8
    subplot(1,9,8)
    plot(smooth(s1(v*7+1:v*8)), 'Color', col1, 'LineWidth', lw);
    hold on 
    plot(smooth(s2((v*7+1:v*8))), 'Color', col2, 'LineWidth', lw);
    plot(smooth(s3((v*7+1:v*8))), 'Color', col3, 'LineWidth', lw);
    plot(smooth(avs((v*7+1:v*8))), 'Color', 'r', 'LineWidth', lw2);
    box off
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.TickDir = 'out';
    ylim([0 ymax])
    
    subplot(1,9,9)
    
    st1 = Grating_Table{allHET(i2), 20};
    st2 = Grating_Table{allHET(i2), 21};
    st3 = Grating_Table{allHET(i2), 22};
    avst = Grating_Table{allHET(i2), 18};
    
    rho = st1;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    % REP2
    hold on
    rho = st2;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    % REP3
    hold on
    rho = st3;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', [0.6 0.6 0.6])
    
    % MEAN 
    rho = avst;
    rho(9) = rho(1);
    polarplot(angls,rho, 'Color', 'r', 'LineWidth', 2.5)
    ax = gca;
    ax.RTickLabel = {''};
    ax.ThetaTickLabel = {''};
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
   
    f= gcf;
    f.Position =  [131   214   1662   158]; 
%     
%     
%     f= gcf;
%     f.Position =  [1606  184  1860  192];  %[1606  124 2216  252];
%     

%%




%% 
close all

dstart = 400;
dend = -100;

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2); 
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2);

data_WT = Grating_Table.DS_Ind(allWT);
data_HET = Grating_Table.DS_Ind(allHET);

ymax = 2.75;

rng = 0:0.075:1;

mWT = nanmean(data_WT);
mHET = nanmean(data_HET);
medWT = nanmedian(data_WT)
medHET = nanmedian(data_HET)

figure; histogram(data_WT, rng, 'FaceColor', 'k', 'FaceAlpha', 0.4, 'Normalization', 'pdf');
hold on 
plot([medWT medWT], [0 ymax], 'k', 'LineWidth', 1.75);
plot([med_ani(2) med_ani(2)], [0 ymax], 'k:', 'LineWidth', 1);
plot([med_ani(5) med_ani(5)], [0 ymax], 'k:', 'LineWidth', 1);
plot([med_ani(6) med_ani(6)], [0 ymax], 'k:', 'LineWidth', 1);
box off
% ylabel('PDF')
% xlabel('DSi')
ax = gca;
% ax.FontSize = 16;
ax.TickDir = 'out'; 
ax.LineWidth = 1.5;
ax.TickLength = [0.025 0.025];
ylim([0 ymax])
ax.XTick = [0, 0.5, 1];
ax.YTick = [0, 1, 2];

f1 = gcf; 
f1.Position = [ 722   650   372   155]; %[722   624   376   225]; 

figure
histogram(data_HET, rng, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'Normalization', 'pdf'); hold on; 
plot([medHET medHET], [0 ymax], 'r', 'LineWidth', 1.75);
plot([med_ani(1) med_ani(1)], [0 ymax], 'r:', 'LineWidth', 1);
plot([med_ani(3) med_ani(3)], [0 ymax], 'r:', 'LineWidth', 1);
plot([med_ani(4) med_ani(4)], [0 ymax], 'r:', 'LineWidth', 1);
box off
% xlabel('DSi')
% ylabel('PDF')
ax = gca;
% ax.FontSize = 16;
ax.TickDir = 'out'; 
ax.LineWidth = 1.5;
ax.TickLength = [0.025 0.025];
ylim([0 ymax])
ax.XTick = [0, 0.5, 1];
ax.YTick = [0, 1, 2];

f = gcf; 
f.Position = [ 722   650   372   155]; %[718   292   376   225]; 



[h,p] = kstest2(data_WT, data_HET)
% 8.73x10-4; < 10 - less than 0 
% 0.0806 % no specs on firing rate - less than 0 
%  0.9737 : >10 < -350 
% 0.0435 
% 0.0523 : -700 to -1050
% 0.3810 : - 1050 to -2000












%% 80th percentile of DSI in sSC 
% For fidning the proportion of DS cells per animal - set a 'DS' cell as a cell with a DSI in the top 20 percentile of all cells.  

allcells = find(Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2); 
da = Grating_Table.DS_Ind(allcells);
figure; histogram(da, 0:0.1:1)
P = prctile(da, [25, 50, 80])
% 80th = 0.65

PRCT80 = 0.65;
%% DSind per animal

all_animals = unique(Grating_Table.Ani);
ani_dsi = zeros(1,6);

het_animals = [7269, 7476, 7614];

for n = 1:6
    anii = all_animals(n); 
    all_ani = find(Grating_Table.Ani == anii & Grating_Table.TotalSpikes >2); %& Grating_Table.Depth <400 & Grating_Table.Depth >-100
    all_ani_dsi = Grating_Table{all_ani, 11}; 
    ani_dsi(1,n) = nanmean(all_ani_dsi);
    
    if ismember(anii, het_animals)
        col = 'r'; 
        ani_geno(n) = 2;
    else
        col = 'k';
        ani_geno(n) = 1;
    end 
    
    med_ani(n) = nanmedian(all_ani_dsi);
    prop_DS(n) = numel(find(all_ani_dsi>PRCT80))/ numel(all_ani);
   
%     % Plot histogram of dsi distribution for each animal
    
    % subplot(6,1,n)
%     histogram(all_ani_dsi, 0:0.1:1, 'FaceColor', col, 'Normalization', 'pdf', 'FaceAlpha', 0.5);
%     hold on 
%     plot([med_ani(n) med_ani(n)], [0 2], col, 'LineWidth', 1.3);
%     box off
    
end 



%% Median DSI for each animal
% Plot Line + SEM for WT/ HET 

close

% 1 - Median 
wti = find(ani_geno == 1);
wtvals = med_ani(wti);
semwt = nanstd(wtvals)/sqrt(3); 

heti = find(ani_geno == 2);
hetvals = med_ani(heti);
semhet = nanstd(hetvals)/sqrt(3); 

figure
bar(1, nanmean(wtvals), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
bar(2, nanmean(hetvals),'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');


for i = 1:6
    xval = ani_geno(i);
    yval = med_ani(i);
    
    if ani_geno(i) == 1
        marker = 'k.';
    else 
        marker = 'r.';
    end 
    
    scatter(xval, yval, 1000,  marker, 'jitter', 'on');
    hold on 
end 
axis([0 3 0 1])
box off
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.02, 0.02];
ax.LineWidth= 1.5;

errorbar(1, nanmean(wtvals), semwt, 'k', 'LineWidth', 1.5);
errorbar(2, nanmean(hetvals), semhet, 'r', 'LineWidth', 1.5);

f = gcf;
f.Position = [680   750   230   348]; 
f.Renderer = 'painters'; 

nanmean(wtvals)
nanmean(hetvals)
[h,p] = ttest(wtvals, hetvals)
[p, h] = ranksum(wtvals, hetvals)


%% Proportion of DSI - per animal

close

wti = find(ani_geno == 1);
wtvals = prop_DS(wti);
semwt = nanstd(wtvals)/sqrt(3); 

heti = find(ani_geno == 2);
hetvals = prop_DS(heti);
semhet = nanstd(hetvals)/sqrt(3); 

figure
bar(1, nanmean(wtvals), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
bar(2, nanmean(hetvals),'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');


for i = 1:6
    xval = ani_geno(i);
    yval = prop_DS(i);
    
    if ani_geno(i) == 1
        marker = 'k.';
    else 
        marker = 'r.';
    end 
    
    scatter(xval, yval, 1000,  marker, 'jitter', 'on');
    hold on 
end 
axis([0 3 0 1])
box off
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.04, 0.04];
ax.LineWidth= 1.5;

errorbar(1, nanmean(wtvals), semwt, 'k', 'LineWidth', 1.5);
errorbar(2, nanmean(hetvals), semhet, 'r', 'LineWidth', 1.5);

f = gcf;
f.Position = [680   750   230   348]; 

nanmean(wtvals)
nanmean(hetvals)
[p,h] = ttest(wtvals, hetvals)
[p, h] = ranksum(wtvals, hetvals)


%% Proportion of DSI - pooled across animals. 

dstart = 400;
dend = -100; 

allWT = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2); 
allHET = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2);

allWTDSI = find(Grating_Table.Geno == 1 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2 & Grating_Table.DS_Ind > 0.6); 
allHETDSI = find(Grating_Table.Geno == 2 & Grating_Table.Depth>dend & Grating_Table.Depth <dstart & Grating_Table.TotalSpikes >2 & Grating_Table.DS_Ind > 0.6);

numel(allWTDSI)/numel(allWT)
numel(allHETDSI)/numel(allHET)

for i = 1:height(Grating_Table)
    val = Grating_Table.DS_Ind(i);
    if val>0.6
        Grating_Table.Selective(i) = 1;
    else
        Grating_Table.Selective(i) = 0;
    end 
end 

wtvals = Grating_Table.Selective(allWT);
hetvals = Grating_Table.Selective(allHET);

    
nanmean(wtvals)
nanmean(hetvals)
[p,h] = ttest(wtvals, hetvals)
[p, h] = ranksum(wtvals, hetvals)
[h, p] = kstest2(wtvals, hetvals)
