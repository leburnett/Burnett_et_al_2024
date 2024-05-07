%% Analyse pre-loom exploration activity
% Created by Burnett - 07/06/22

% Behavioural metrics calculated over the first 5 minutes of exploration in the arena. This period of time is before the mice ever see a looming stimulus. 

% Setd5 Data found here: /Users/lauraburnett/Documents/Burnett_etal/DATA/SETD5/BEHAVIOUR/PreLoom_Exploration
% Ptchd1 data found here: /Users/lauraburnett/Documents/Burnett_etal/DATA/PTCHD1/Behaviour/PreLoom
% Cul3 data found here: /Users/lauraburnett/Documents/Burnett_etal/DATA/CUL3/Behaviour/PreLoom_Exploration

% Script located here: /Users/lauraburnett/Documents/Burnett_etal/Scripts/Behaviour

% Use 'ALL_ACTIVITY_ARRAY' for initial box plots then use 'exit_analysis'
% for trajectory heatmaps and other box plots. 

allWT = find(string(ALL_ACTIVITY_TABLE.Geno) =="wt");
allHET = find(string(ALL_ACTIVITY_TABLE.Geno) =="het");

col = 'm'; 
% col = [255/255 114/255 32/255]; 

%% STATS

% 1 - MEAN VALUES
nanmean(speed_WT)
nanmean(speed_HET)

% 2 - test normality with Ks test. - more than 30 samples.
% [h,p] =kstest(dataWT)
% [h,p] =kstest(dataHET)

% Less than 30 samples
[H, pValue] = swtest(speed_WT)
[H, pValue] = swtest(speed_HET)

% 3 - Levene test for homogenity of variance. 
nWT = numel(speed_WT);
nHET = numel(speed_HET);

col1 = vertcat(speed_WT, speed_HET);
col2 = vertcat(ones(nWT,1), ones(nHET,1)*2);

X = horzcat(col1, col2);

%Levene's test for variance:
Levenetest(X)

% Difference of means - parametric and non-parametric
[p,h] = ttest2(speed_WT, speed_HET)
[p,h] = ranksum(speed_WT, speed_HET)

%% BOXPLOTS - ALL_ACTIVITY_ARRAY

% Col 7 = Speed during the acclim period.
% Col 8 = Max Speed
% Col 9 = Percentage of Time < 2cm/s
% Col 10 = Percentage of time out of the box
% Col 11 = Percentage of time spent in the centre
% Col 12 - Percentage of time spent at the edge of the box. 

val = 7; % Column you want to assess. 

n = height(ALL_ACTIVITY_TABLE);

speed_WT = [];
speed_HET = []; 

for i = 1:n
        if string(ALL_ACTIVITY_TABLE.Geno{i}) == "wt"
            G = (ALL_ACTIVITY_TABLE{i, val});
            speed_WT = vertcat(speed_WT, G);
        elseif string(ALL_ACTIVITY_TABLE.Geno{i}) == "het"
            F = (ALL_ACTIVITY_TABLE{i, val});
            speed_HET = vertcat(speed_HET, F);
        end
end 

% PLOT 
n_wt = numel(speed_WT);
n_het = numel(speed_HET);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

geno = string(ALL_ACTIVITY_TABLE.Geno); 
var = (ALL_ACTIVITY_TABLE{:, val});

figure
b = boxplot(var, geno, 'Colors', 'k');
set(b , 'LineWidth', 1.5)
hold on 
scatter(x1, speed_WT','SizeData', 175, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.05)
scatter(x2, speed_HET', 'SizeData', 175,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05)

xticks([1,2])
xticklabels({''})
ax = gca;
ax.FontSize = 30;
box off
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 2;
% ax.XAxis.Visible = 'off'; 
xlim([0.5 2.5])
% ylabel('Speed (cm s^-1)')
% ylabel('Number of Triggers')
ylim([0 40])

f = gcf;
f.Position = [704   288   275   491]; %[704   207   355   572]; 
% ylabel('Time (s)')
% ylabel('Proportion of Trials')
% ylabel('Delta Speed')
% ylabel('Acceleration cm s^-2')



%% BOXPLOTS - Exit_analysis

% Col 8 = Num OUt
% COl 9 - Number out past trigger
% Col 10 - MEan D
% Col 11 - Max D
% Col 12 - Inter bout interval
% Col 13 - Exit duration

val = 13; % Column you want to assess. 

n = height(exit_analysis);

speed_WT = [];
speed_HET = []; 

for i = 1:n
        if (exit_analysis.geno(i)) == 1
            G = table2array(exit_analysis(i, val));
            speed_WT = vertcat(speed_WT, G);
        elseif (exit_analysis.geno(i)) == 2
            F = table2array(exit_analysis(i, val));
            speed_HET = vertcat(speed_HET, F);
        end
end 

% PLOT 
n_wt = numel(speed_WT);
n_het = numel(speed_HET);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

geno = (exit_analysis.geno); 
var = table2array(exit_analysis(:, val));

figure
b = boxplot(var, geno, 'Colors', 'k');
set(b , 'LineWidth', 1.5)
hold on 
scatter(x1, speed_WT','SizeData', 175, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.05)
scatter(x2, speed_HET', 'SizeData', 175,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05)

xticks([1,2])
xticklabels({''})
ax = gca;
ax.FontSize = 30;
box off
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 2;
% ax.XAxis.Visible = 'off'; 
xlim([0.5 2.5])
% ylabel('Speed (cm s^-1)')
% ylabel('Number of Triggers')
ylim([0 16])

f = gcf;
f.Position = [704   288   275   491]; %[704   207   355   572]; 
% ylabel('Time (s)')
% ylabel('Proportion of Trials')
% ylabel('Delta Speed')
% ylabel('Acceleration cm s^-2')


%% Plot HEATMAP of location in arena where mice spend most of their time. 

positionWT = [];
positionHET = [];

for j = 1:height(exit_analysis)
    
    
%     if j <=14
%         maxD = 512;
%     else
        maxD = 416;
%     end
    
    n_bins = 12;
    bin_size = maxD/n_bins;
    n = numel(cell2mat(exit_analysis.X(j)));
    
    dataa = zeros(n, 2);
    
    xdata = cell2mat(exit_analysis.X(j));
    ydata = cell2mat(exit_analysis.Y(j));
    
    dataa(:,1)=ceil(xdata/bin_size);
    dataa(:,2)=ceil(ydata/bin_size);
    
    NumF = zeros(n_bins,n_bins);
    
    for j2 = 1:n
        NumF(dataa(j2,1), dataa(j2,2)) = NumF(dataa(j2,1),dataa(j2,2)) + 1;
    end
    
    NumF_Norm = NumF/n; % Normalised by the number of frames. 
    
    genoo = exit_analysis.geno(j);
    
    if genoo == 1
        if isempty(positionWT)
            positionWT = NumF_Norm;
        else
            positionWT = (positionWT + NumF_Norm)/2;
        end
    else
        if isempty(positionHET)
            positionHET = NumF_Norm;
        else
            positionHET = (positionHET + NumF_Norm)/2;
        end
    end 
    
end

figure; imagesc(positionWT); colormap(redblue); caxis([0 0.04]); axis off; axis square
figure; imagesc(positionHET); colormap(redblue); caxis([0 0.04]); axis off; axis square

% 1-14 = 512;
% 15-28 = 416; 






%% Find the maximum values of the x,y values.
% This will give us the upper bounds of our data - useful when assigning
% data to bins. 
max_x = max(centroids_array(:,1));
max_y = max(centroids_array(:,2));
upper_bound = 360; 

bin_size = 30; % How large we want to make each bin. 
n_bins = upper_bound/bin_size; % The number of bins of size 'bin_size' that cover the arena.
n = numel(centroids_array(:,1)); % Number of frames. 

%% Assign data to 'bins'

centroids_array(:,3)=ceil(centroids_array(:,1)/bin_size);
centroids_array(:,4)=ceil(centroids_array(:,2)/bin_size);

%% Adding DeltaF/F data into centroids_array for ONE CELL. 
%Then using the 'bin values' that are now in columns 3/4 - use these to
%assign the DF/F value to a certain spatial bin. 

for j = 1 % Only the first column in the calcium imaging file. 15 ROIs/cells in total. use j = 1:15 for all cells.
    
% centroids_array(:,5) = Data_video00(:,j); %Adding data from column 1 of DeltaF/F data. 
% data = centroids_array; %'renaming' to data for ease. 

%% Create array of Spatial Activity

% Create two arrays. Number of values in these arrays == the number of bins.
% Activity will sum up the activity value in each bin. 
% NumF will sum up the number of frames that are being allocated to that
% bin. 

% Activity = zeros(n_bins,n_bins); 
NumF = zeros(n_bins,n_bins); 

for i = 1:n  % Number of rows/frames in the video
    
% Use the values in col3/4 as your indices for Activity/NumF.

% For Activity - add the DF/F value.
% For NumF - add 1 for each row assigned to that bin. 
       
% Sum the values of deltaF/F
% Activity(data(i,3), data(i,4)) = Activity(data(i,3),data(i,4))+(data(i,5)); %This is a cumulative function. 
%Taking the value that is already inside this 'bin' and adding the activity
%value for that row to the current total for that bin. 

%Sum the number of frames 
NumF(data(i,3), data(i,4)) = NumF(data(i,3),data(i,4)) + 1; 
%Again, a cumulative function. 
end

%% Find the activity divided by the number of frames for each bin. 
% Imagine the mouse just sat in one corner for 1/2 of the recording -
% obviously this spatial bin would contain a much higher DF/F value. Need
% to normalise. 

Norm_Act = zeros(12,12);
Norm_Act = Activity./NumF; %Elementwise division of Activity array but Number of Frames array. 

%% Save all the figures
Z = imagesc(Norm_Act);
% colorbar
% caxis([0 50])
% filename = sprintf('FigCell_%d', j);
% saveas(Z,filename);
end 

























