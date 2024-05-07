% Analyse Silicon Probe recordings
% Created by Burnett - 14/06/22

% General analysis about recordings
% Also need to load specific TABLES created for the the different stimuli -
% e.g. FLASH_TABLE from 'comb_flash_analysis.m' for the Flash Stimulus etc.
% 

%% GENERAL STATS - PER ANIMAL

% LOAD table with information about the depth of the sSC in different animals. 
load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Setd5_Rec_Depth_Info.mat', 'animal_depth_info_table')

WTani = find(animal_depth_info_table.Geno == 1);
HETani = find(animal_depth_info_table.Geno == 0);

% het_animals = [7269, 7476, 7614];

%% Average number of spikes over the entire recording - ALL CELLS:
dataWT = animal_depth_info_table.AvSpikesAll(WTani);
dataHET = animal_depth_info_table.AvSpikesAll(HETani);

nanmean(dataWT) % 4914
nanmean(dataHET) % 3954
[p, h] = ranksum(dataWT, dataHET)% p = 0.700

%%  Average number of spikes over the entire recording - VR CELLS:
dataWT = animal_depth_info_table.AvSpikesVR(WTani);
dataHET = animal_depth_info_table.AvSpikesVR(HETani);

nanmean(dataWT) % 7245
nanmean(dataHET) % 7507
[p, h] = ranksum(dataWT, dataHET)% p = 0.700

%% ANALYSIS USING FLASH_TABLE
% data here: '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe'

% 1 - LOAD FLASH_TABLE. 

allWT = find(Flash_Table.Geno == 1 & Flash_Table.Depth<0);
allHET = find(Flash_Table.Geno == 2 & Flash_Table.Depth<0);

allWTRESP = find(Flash_Table.Geno == 1 & Flash_Table.Depth<0 & Flash_Table.P_VALUE<0.05);
allHETRESP = find(Flash_Table.Geno == 2 & Flash_Table.Depth<0 & Flash_Table.P_VALUE<0.05);

%% Plot the depth of the responsive cells: 
n = height(Flash_Table);

figure
for i = 1:n
    
    if Flash_Table.Depth(i) < 0 & Flash_Table.P_VALUE(i) < 0.05
        
        ani= Flash_Table.Ani(i);
        
        % Set colour of dot for individual mouse:
        if ani == 7270
            col = [0 0 0];
            xval = 1;
        elseif ani == 7269
            col = [1 0 0];
            xval = 2;
        elseif ani == 7476
            col = [1 0.4 0.4];
            xval = 6;
        elseif ani == 7614
            col = [1 0.8 0.8];
            xval = 4;
        elseif ani == 7616
            col = [0.4 0.4 0.4];
            xval = 5;
        elseif ani == 7788
            col = [0.8 0.8 0.8];
            xval = 3;
        elseif ani == 7475
            col = [0.8 0.8 0.8];
            xval = 7;
        end
        
        rndnum = (2 * rand - 1)/4;
        xval = xval + rndnum;
        yval = Flash_Table.Depth(i);
        plot(xval, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 15)
        hold on
    end
    
end

box off 
ax= gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
ax.FontSize = 18;
xlim([0 7])
ylim([-1500 100])
f = gcf;
f.Position = [70   250   461   806];


%% DEPTH - STATS - ANIMAL AVERAGES
all_ani = unique(Flash_Table.Ani);

WTVALS = [];
HETVALS = []; 

het_animals = [7269, 7476, 7614];
        
for j = 1:6 
    ani= all_ani(j);
    
    ani_rows = find(Flash_Table.Depth<0 & Flash_Table.Ani == ani & Flash_Table.P_VALUE<0.05); % 
    ani_mean = nanmean(Flash_Table.Depth(ani_rows));
    
    if ismember(ani, het_animals)
        genoo = 2;
    else
        genoo = 1;
    end
    
    if genoo ==1
        WTVALS = [WTVALS, ani_mean];
    elseif genoo ==2
        HETVALS = [HETVALS, ani_mean];
    end
end

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)

%% DEPTH - STATS - POOLED

WTVALS = [];
HETVALS = []; 

for j = 1:n
    
    if Flash_Table.Depth(j)<0 %& Flash_Table.P_VALUE(j)<0.05
        
        val = Flash_Table.Depth(j);
        
        if Flash_Table.Geno(j) ==1
            WTVALS = [WTVALS, val];
        elseif Flash_Table.Geno(j) ==2
            HETVALS = [HETVALS, val];
        end
    end
    
end

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)


%% Find the % of Flash Responsive cells PER ANIMAL

all_ani = unique(Flash_Table.Ani);

WTVALS = [];
HETVALS = []; 

het_animals = [7269, 7476, 7614];
        
for j = 1:numel(all_ani)
    ani= all_ani(j);
    
    ani_all =find(Flash_Table.Ani == ani & Flash_Table.Depth<=-400 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500); % Flash_Table.Depth<-500 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500
    n_all = numel(ani_all);
    ani_rows = find(Flash_Table.Ani == ani & Flash_Table.P_VALUE<0.05 & Flash_Table.Depth<=-400 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500); % Flash_Table.Depth<-500 & Flash_Table.Depth>-1000
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
allwt = find(Flash_Table.Geno == 1 & Flash_Table.Depth<0 & Flash_Table.Depth>-1000);
allhet = find(Flash_Table.Geno ==2 & Flash_Table.Depth<0 & Flash_Table.Depth>-1000);

w1 = Flash_Table.P_VALUE(allwt);
h1 = Flash_Table.P_VALUE(allhet);

[p,h] = kstest2(w1, h1)


%% % of Flash REsponsive cells - pooled across all animals
% Number of units in each condition: 

    all_wt = numel(find(Flash_Table.Geno == 1 & Flash_Table.Depth<=-400 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500)) 
    all_het = numel(find(Flash_Table.Geno == 2 & Flash_Table.Depth<=-400 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500)) 

% Responsive cells: 

    allwt_resp = numel(find(Flash_Table.Geno == 1 & Flash_Table.P_VALUE<0.05 & Flash_Table.Depth<=-400 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500))
    allhet_resp = numel(find(Flash_Table.Geno == 2 & Flash_Table.P_VALUE<0.05 & Flash_Table.Depth<=-400 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500))
    
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
    
    

%% BAR CHART + ANIMAL POINTS - % OF FLASH RESPONSIVE CELLS
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


%%  PLOT P-VAL VERSUS DEPTH

n = height(Flash_Table);

figure
for i = 1:n
    
    if Flash_Table.Depth(i) < 0 & Flash_Table.Depth(i)>-1000% & Flash_Table.P_VALUE(i) < 0.05
        
        ani= Flash_Table.Ani(i);
        
        % Set colour of dot for individual mouse:
        if ani == 7270
            col = [0 0 0];

        elseif ani == 7269
            col = [1 0 0];

        elseif ani == 7476
            col = [1 0.4 0.4];

        elseif ani == 7614
            col = [1 0.8 0.8];
     
        elseif ani == 7616
            col = [0.4 0.4 0.4];
      
        elseif ani == 7788
            col = [0.8 0.8 0.8];
       
        end
        
        xval = (Flash_Table.P_VALUE(i));
        yval = Flash_Table.Depth(i);
        plot(xval, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 10)
        hold on
    end
    
end

box off 
ax= gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
ax.FontSize = 18;
xlim([0 7])
ylim([-1500 100])
f = gcf;
f.Position = [70   250   461   806];

%% FInd MAX and MIN - zscore value - comparison to baseline - LOOKING AT THE DYNAMIC RANGE
% Find these during ON/OFF stages separately:

% ON
for j = 1:height(Flash_Table)
    data = Flash_Table.OFO(j,31:90)*60;
    maxd = max(data);
    mind = min(data);
    
    norm_range = abs(diff([mind ,maxd]));
    Flash_Table.Range_ON(j) = norm_range;
  
    Flash_Table.Min_ON(j) = mind;
    Flash_Table.Max_ON(j) = maxd;
end 

% OFF
for j = 1:height(Flash_Table)
    data = [Flash_Table.OFO(j,1:30), Flash_Table.OFO(j, 91:120)]*60;
    maxd = max(data);
    mind = min(data);
    
    norm_range = abs(diff([mind ,maxd]));
    Flash_Table.Range_OFF(j) = norm_range;
  
    Flash_Table.Min_OFF(j) = mind;
    Flash_Table.Max_OFF(j) = maxd;
end 

% On versus OFF 
for j = 1:height(Flash_Table)
    max1 = Flash_Table.Max_ON(j);
    max2 = Flash_Table.Max_OFF(j);
    maxR = max1/max2;
    Flash_Table.Max(j) = maxR;
    
    min1 = Flash_Table.Min_ON(j);
    min2 = Flash_Table.Min_OFF(j);
    minR = min1/min2;
    Flash_Table.Min(j) = minR;
end 


allWT = find(Flash_Table.Geno ==1  & Flash_Table.P_VALUE<0.05 & Flash_Table.Depth<=0 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500); %& Flash_Table.Max~=Inf
allHET = find(Flash_Table.Geno ==2  & Flash_Table.P_VALUE<0.05 & Flash_Table.Depth<=0 & Flash_Table.Depth>-1000 & Flash_Table.TotalSpikes<3500);

d1 = Flash_Table.Max_ON(allWT);
d2 = Flash_Table.Max_ON(allHET);

d1 = Flash_Table.Min_ON(allWT);
d2 = Flash_Table.Min_ON(allHET);

d1 = Flash_Table.Range_ON(allWT);
d2 = Flash_Table.Range_ON(allHET);

d1 = Flash_Table.Max_OFF(allWT);
d2 = Flash_Table.Max_OFF(allHET);

d1 = Flash_Table.Min_OFF(allWT);
d2 = Flash_Table.Min_OFF(allHET);

d1 = Flash_Table.Range_OFF(allWT);
d2 = Flash_Table.Range_OFF(allHET);

d1 = Flash_Table.Max(allWT);
d2 = Flash_Table.Max(allHET);

nanmean(d1)
nanmean(d2)
[p,h] = kstest2(d1,d2)
[p,h] = ranksum(d1,d2)

figure; histogram(d1, 'Normalization', 'pdf');
hold on;
histogram(d2, 'Normalization', 'pdf')

%% 
% % % % % % % % % % % % % % % % % % % . % % % % % % %
%%











%% %%% ANALYSIS USING LOOM_TABLE
% data here: '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe'

% 1 - LOAD FLASH_TABLE. 

allWT = find(Loom_Table.Geno == 1 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);
allHET = find(Loom_Table.Geno == 2 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);

allWTRESP = find(Loom_Table.Geno == 1 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000 & Loom_Table.PVALS<0.05);
allHETRESP = find(Loom_Table.Geno == 2 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000 & Loom_Table.PVALS<0.05);

% Combine L1ONLY-Rep1 and REP2 

for i = 1:1847
    
    data1 = Loom_Table.L1ONLY_rep1(i,:);
    data2 = Loom_Table.L1ONLY_rep2(i,:);
    
%     figure; imagesc(data1); figure; imagesc(data2)
    d3 = cat(1, data1, data2);
    d3 = nanmean(d3, 1);
    
    Loom_Table.L1ONLY_AV(i, 1:331) = d3;
    
end 


unique(Loom_Table.Ani)

%% Plot the depth of the responsive cells: 
n = height(Loom_Table);

figure
for i = 1:n
    
%     if Loom_Table.Depth(i) < 0 
        
        ani= Loom_Table.Ani(i);
        
        % Set colour of dot for individual mouse:
        if ani == 7270 || ani == 1387
            col = [0 0 0];
            xval = 1;
        elseif ani == 7269 || ani == 1389
            col = [1 0 0];
            xval = 2;
        elseif ani == 7476 || ani == 1394
            col = [1 0.4 0.4];
            xval = 6;
        elseif ani == 7614 || ani == 2709
            col = [1 0.8 0.8];
            xval = 4;
        elseif ani == 7616 || ani == 2710
            col = [0.4 0.4 0.4];
            xval = 5;
        elseif ani == 7788 || ani == 4366
            col = [0.8 0.8 0.8];
            xval = 3;
        elseif ani == 4369
            col = [0.2 0.2 0.2];
            xval = 7;
        end
        
        rndnum = (2 * rand - 1)/5;
        xval = xval + rndnum;
        yval = Loom_Table.Depth(i);
        plot(xval, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 7)
        hold on
%     end
    
end

box off 
ax= gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
ax.FontSize = 18;
xlim([0 7])
ylim([-1500 100])
f = gcf;
f.Position = [70   250   461   806];

%% Plot the depth of all cells - light grey for unresponsive cells and red/black for responsive (p<0.05)
%% Loom TABLE
n = height(Loom_Table);

figure
for i = 1:n
    
%     if Loom_Table.Depth(i) < 0
        ani= Loom_Table.Ani(i);
        
        if Loom_Table.P_VALUE(i)<0.05
            
            % Set colour of dot for individual mouse:
            if ani == 7270
                col = [0 0 0];
                xval = 1;
            elseif ani == 7269
                col = [1 0 0];
                xval = 2;
            elseif ani == 7476
                col = [1 0 0];
                xval = 6;
            elseif ani == 7614
                col = [1 0 0];
                xval = 4;
            elseif ani == 7616
                col = [0 0 0];
                xval = 5;
            elseif ani == 7788
                col = [0 0 0];
                xval = 3;
            elseif ani == 7475
                col = [0 0 0];
                xval = 7;
            end
            
        elseif Loom_Table.P_VALUE(i)>=0.05
            
            % Set colour of dot for individual mouse:
            if ani == 7270
                col = [0.7 0.7 0.7];
                xval = 1;
            elseif ani == 7269
                col = [0.7 0.7 0.7];
                xval = 2;
            elseif ani == 7476
                col = [0.7 0.7 0.7];
                xval = 6;
            elseif ani == 7614
                col = [0.7 0.7 0.7];
                xval = 4;
            elseif ani == 7616
                col = [0.7 0.7 0.7];
                xval = 5;
            elseif ani == 7788
                col = [0.7 0.7 0.7];
                xval = 3;
            elseif ani == 7475
                col = [0.7 0.7 0.7];
                xval = 7;   
            end
            
        end
        
        rndnum = (2 * rand - 1)/4;
        xval = xval + rndnum;
        yval = Loom_Table.Depth(i);
        plot(xval, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 10)
        hold on
%     end
    
end

box off 
ax= gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
ax.FontSize = 18;
xlim([0 7])
ylim([-1500 100])
f = gcf;
f.Position = [70   250   461   806];
% title('Loom')

%% Plot the depth of all cells - light grey for unresponsive cells and red/black for responsive (p<0.05)
%% FLASH TABLE
n = height(Flash_Table);

figure
for i = 1:n
    
%     if Flash_Table.Depth(i) < 0
        ani=  Flash_Table.Ani(i);
        datee = Flash_Table.Date(i);
        
        if Flash_Table.P_VALUE(i)<0.05
            
            % Set colour of dot for individual mouse:
            if ani == 7270
                col = [0 0 0];
                xval = 1;
            elseif ani == 7269
                col = [1 0 0];
                xval = 2;
            elseif ani == 7476 && datee == 200924
                col = [1 0 0];
                xval = 7;
            elseif ani == 7476 && datee == 200925
                col = [1 0 0];
                xval = 8;
            elseif ani == 7614 && datee == 201002
                col = [1 0 0];
                xval = 4;
            elseif ani == 7616 && datee == 201001
                col = [0 0 0];
                xval = 9;
            elseif ani == 7616 && datee == 201005
                col = [0 0 0];
                xval = 10;
            elseif ani == 7788
                col = [0 0 0];
                xval = 3;
            elseif ani == 7475
                col = [0 0 0];
                xval = 6;
            elseif ani == 7614 && datee == 201006
                col = [1 0 0];
                xval = 5;
            end
            
        elseif Flash_Table.P_VALUE(i)>=0.05
            
            if ani == 7270
                col = [0.7 0.7 0.7];
                xval = 1;
            elseif ani == 7269
                col = [0.7 0.7 0.7];
                xval = 2;
            elseif ani == 7476 && datee == 200924
                col = [0.7 0.7 0.7];
                xval = 7;
            elseif ani == 7476 && datee == 200925
                col = [0.7 0.7 0.7];
                xval = 8;
            elseif ani == 7614 && datee == 201002
                col = [0.7 0.7 0.7];
                xval = 4;
            elseif ani == 7616 && datee == 201001
                col = [0.7 0.7 0.7];
                xval = 9;
            elseif ani == 7616 && datee == 201005
                col = [0.7 0.7 0.7];
                xval = 10;
            elseif ani == 7788
                col = [0.7 0.7 0.7];
                xval = 3;
            elseif ani == 7475
                col = [0.7 0.7 0.7];
                xval = 6;
            elseif ani == 7614 && datee == 201006
                col = [0.7 0.7 0.7];
                xval = 5;
            end
            
            % Set colour of dot for individual mouse:
%             if ani == 7270
%                 col = [0.7 0.7 0.7];
%                 xval = 1;
%             elseif ani == 7269
%                 col = [0.7 0.7 0.7];
%                 xval = 2;
%             elseif ani == 7476
%                 col = [0.7 0.7 0.7];
%                 xval = 6;
%             elseif ani == 7614
%                 col = [0.7 0.7 0.7];
%                 xval = 4;
%             elseif ani == 7616
%                 col = [0.7 0.7 0.7];
%                 xval = 5;
%             elseif ani == 7788
%                 col = [0.7 0.7 0.7];
%                 xval = 3;
%             elseif ani == 7475
%                 col = [0.7 0.7 0.7];
%                 xval = 7;   
%             end
            
        end
        
        rndnum = (2 * rand - 1)/4;
        xval = xval + rndnum;
        yval = Flash_Table.Depth(i);
        plot(xval, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 10)
        hold on
%     end
    
end

box off 
ax= gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
ax.FontSize = 18;
xlim([0 7])
ylim([-1500 100])
f = gcf;
f.Position = [70   250   461   806];
% title('Loom')


%% DEPTH - STATS - ANIMAL AVERAGES
all_ani = unique(Loom_Table.Ani);

WTVALS = [];
HETVALS = []; 

het_animals = [7269, 7476, 7614];
        
for j = 1:6 
    ani= all_ani(j);
    
    ani_rows = find(Loom_Table.Depth<0 & Loom_Table.Ani == ani & Loom_Table.PVALS<0.05); % 
    ani_mean = nanmean(Loom_Table.Depth(ani_rows));
    
    if ismember(ani, het_animals)
        genoo = 2;
    else
        genoo = 1;
    end
    
    if genoo ==1
        WTVALS = [WTVALS, ani_mean];
    elseif genoo ==2
        HETVALS = [HETVALS, ani_mean];
    end
end

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)


%% Find the % of Loom Responsive cells PER ANIMAL

all_ani = unique(Loom_Table.Ani);

WTVALS = [];
HETVALS = []; 

het_animals = [7269, 7476, 7614];
        
for j = 1:6
    ani= all_ani(j);
    
    ani_all =find(Loom_Table.Depth<-400 & Loom_Table.Depth>-1000 & Loom_Table.Ani == ani); %5
    n_all = numel(ani_all);
    ani_rows = find(Loom_Table.Depth<-400 & Loom_Table.Depth>-1000 & Loom_Table.Ani == ani & Loom_Table.P_VALUE<0.01); % 
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
% [p,h] = ranksum(WTVALS, HETVALS)
[p,h] = ttest2(WTVALS, HETVALS)



   ani_all =find(Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.Geno == 2); 
   ani_rows = find(Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.Geno == 2 & Loom_Table.P_VALUE<0.01); % 
   numel(ani_all)
   numel(ani_rows)

%% BAR CHART + ANIMAL POINTS - % OF LOOM RESPONSIVE CELLS
% Plus errorbar

% close
semwt = nanstd(WTVALS)/sqrt(3); 
semhet = nanstd(HETVALS)/sqrt(3); 

figure
bar(1, nanmean(WTVALS), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
bar(2, nanmean(HETVALS),'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% WT
for i = 1:3
    xval =  1;
    yval = WTVALS(i);
   
    marker = 'k.';
    scatter(xval, yval, 1000,  marker, 'jitter', 'on');
    hold on 
end 

% HET
for i = 1:3
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


%% Distribution of p-values to loom 

allWT = find(Loom_Table.Geno == 1 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);
allHET = find(Loom_Table.Geno == 2 & Loom_Table.Depth<0 & Loom_Table.Depth>-1000);

PVALSWT = Loom_Table.P_VALUE(allWT);
PVALSHET = Loom_Table.P_VALUE(allHET);

figure; histogram(PVALSWT); hold on; histogram(PVALSHET)
[h,p] = kstest2(PVALSHET, PVALSWT)
% 0.1945 - NS

%% % of LOOM REsponsive cells - pooled across all animals
% Number of units in each condition: 

    all_wt = numel(find(Loom_Table.Geno == 1 & Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.TotalSpikes<3500)) 
    all_het = numel(find(Loom_Table.Geno == 2 & Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.TotalSpikes<3500)) 

% Responsive cells: 

    allwt_resp = numel(find(Loom_Table.Geno == 1 & Loom_Table.P_VALUE<0.01 & Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.TotalSpikes<3500))
    allhet_resp = numel(find(Loom_Table.Geno == 2 & Loom_Table.P_VALUE<0.01 & Loom_Table.Depth<0 & Loom_Table.Depth>-400 & Loom_Table.TotalSpikes<3500))
    
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
    


%% %% FInd MAX and MIN - zscore value - comparison to baseline - LOOKING AT THE DYNAMIC RANGE

for j = 1:1847
    data = Loom_Table.L25_norm1(j,:);
    maxd = max(data);
    mind = min(data);
    
%     norm_range = abs(diff([mind ,maxd]));
%     Flash_Table.Norm_Range(j) = norm_range;
    
    data2 = Loom_Table.L251(j,:);
    maxdSP = max(data2);
    
    Loom_Table.Min(j) = mind;
    Loom_Table.Max(j) = maxd;
    Loom_Table.MaxSpikes(j) = maxdSP;
end 

allWT = find(Loom_Table.Geno ==1  & Loom_Table.P_VALUE_FLASH<0.05 & Loom_Table.Depth<0);
allHET = find(Loom_Table.Geno ==2  & Loom_Table.P_VALUE_FLASH<0.05 & Loom_Table.Depth<0);

d1 = Flash_Table.Max(allWT);
d2 = Flash_Table.Max(allHET);

d1 = Flash_Table.Min(allWT);
d2 = Flash_Table.Min(allHET);

d1 = Flash_Table.MaxSP(allWT);
d2 = Flash_Table.MaxSP(allHET);

d1 = Flash_Table.Norm_Range(allWT);
d2 = Flash_Table.Norm_Range(allHET);

[p,h] = kstest2(d1,d2)
[p,h] = ranksum(d1,d2)

figure; histogram(d1, 'Normalization', 'pdf');
hold on;
histogram(d2, 'Normalization', 'pdf')

figure; histogram(d1);
hold on;
histogram(d2)


%% Correlation between flash responsive and loom responsive cells

allWTF = find(Loom_Table.Geno ==1  & Loom_Table.P_VALUE_FLASH<0.05 & Loom_Table.Depth<0); % 256
allHETF = find(Loom_Table.Geno ==2  & Loom_Table.P_VALUE_FLASH<0.05 & Loom_Table.Depth<0); % 333

allWT = find(Loom_Table.Geno ==1  & Loom_Table.PVALS<0.05 & Loom_Table.Depth<0); % 245
allHET = find(Loom_Table.Geno ==2  & Loom_Table.PVALS<0.05 & Loom_Table.Depth<0); % 507

a = sum(ismember(allWT, allWTF));
a = sum(ismember(allHETF, allHET));

% HET - 276/507 %0.5444
% WT- 149/256 % 0.5820

data = Loom_Table.L251(4,:);


%% KMEANS SORT LOOM RESPONSES:

close all

% data = [(Loom_Table.L1ONLY_AV(:, 37:87))*60, Loom_Table.Geno];
data = [cell2mat(Loom_Table.Av10_1L)*60, Loom_Table.Geno];


% Only for cells between 0 and -400. 
all_depths = find(Loom_Table.Depth<0 & Loom_Table.Depth>-1000 & Loom_Table.P_VALUE<0.01 & Loom_Table.TotalSpikes<3500); % & Flash_Table.Depth>-600 &  Flash_Table.Depth>-1000 
% all_depths = find(Flash_Table.Depth>0 | Flash_Table.P_VALUE>=0.05);

data = data(all_depths, :);
% data(all_depths, :) = []; 

% % % % % Evaluate the optimal number of clusters:
eva = evalclusters(data, 'kmeans','CalinskiHarabasz', 'KList', [1:30]); % 'silhouette',
figure; plot(eva)

%% RUN KMEANS

% Set the elbow value as the number of clusters. 
n_k = 10;

data(:, 52) = kmeans(data(:, 1:50), n_k);
data = sortrows(data, 52);

allWT = find(data(:, 51)==1);
allHET = find(data(:, 51)==2);

nWT = numel(allWT);
nHET = numel(allHET);

x = 1:1:50;

%% HEATMAPS

spn = 15; 

figure
ax = subplot(1,spn,1:spn-1);
imagesc(data(allWT,1:50)); caxis([0 100])
hold on
% plot([30 30], [0 nWT], 'w:', 'LineWidth', 1.2)
% plot([90 90], [0 nWT], 'w:', 'LineWidth', 1.2)
box off
ax.XTick = []; 
colormap(ax(1), redblue)
% caxis(ax, [-0.05 0.35])
% caxis(ax, [-0.2 0.8])
caxis(ax, [-10 90])
ax.YTick = [];

ax2 = subplot(1,spn,spn);
imagesc(data(allWT, 52))
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
f = gcf;
f.Position = [100 600 300 nWT*2];

%
figure
ax = subplot(1,spn,1:spn-1);
imagesc(data(allHET,1:50)); caxis([0 100])
% hold on
% plot([30 30], [0 nHET], 'w:', 'LineWidth', 1.2)
% plot([90 90], [0 nHET], 'w:', 'LineWidth', 1.2)
box off
ax.XTick = []; 
colormap(ax(1), redblue)
% caxis(ax, [-0.05 0.35])
% caxis(ax, [-0.2 0.8])
caxis(ax, [-3 23])
ax.YTick = [];

ax2 = subplot(1,spn,spn);
imagesc(data(allHET, 52))
colormap(ax2, gray)
f = gcf;
f.Position = [699 1 356   804];
box on
ax2.LineWidth = 0.5;
ax2.Color = 'k';
ax2.TickLength = [0 0];
ax2.XTick = [];
ax2.YTick = []; 
f2 = gcf;
f2.Position = [100 600 300 nHET*2];  



%% MEAN / SEM - COMBINED

close all

figure
for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 122)== j & data(:,121) == 1);
    % Separate these rows into a different array
    data2 = data(all_type, 1:120);
    
    % Variables for Mean/SEM
    n_trials = numel(all_type);
    if n_trials==1
        mWT = data2;
    else
        mWT = nanmean(data2);
    end
    semWT = nanstd(data2)/sqrt(n_trials);
    y1 = mWT+semWT;
    y2 = mWT-semWT;
    
    % Make the plot
    subplot(n_k,1,j); hold on
    
    % Find Min/Max
%     maxx = max(nanmean(data2))+0.3;
%     minn = min(nanmean(data2))-0.15;
%     
%     maxx = 42 ;%1.1; % 0.55
%     minn = 0 ; % -0.25; %-0.08

     all_cl= find(data(:, 122)==j); 
     max1 =max(nanmean(data(all_cl, 1:120)));
     maxx = max1 + max1*0.75;
     
     if maxx > 15
         minn = -5;
     else 
         minn = 0;
     end 
% 
%     % Rectangle of light stim
%     rectangle('Position', [0 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
% %     rectangle('Position', [30 minn 60 diff([minn maxx])], 'FaceColor', [1 1 0 0.1], 'EdgeColor', 'none')
%     rectangle('Position', [90 minn 30 diff([minn maxx])], 'FaceColor', [0 0 0 0.05], 'EdgeColor', 'none')
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
%     if j~= n_k
%     ax.YAxis.Visible = 'off';
%     end 
    ax.LineWidth = 1;
    
end

for j = 1:n_k
    % Find the rows belonging to the cluster.
    all_type = find(data(:, 122)== j & data(:,121) == 2);
    % Separate these rows into a different array
    data2 = data(all_type, 1:120);
    
    % Variables for Mean/SEM
    n_trials = numel(all_type);
    if n_trials == 1
        mHET = data2;
    else
        mHET = nanmean(data2);
    end
    semHET = nanstd(data2)/sqrt(n_trials);
    y1 = mHET+semHET;
    y2 = mHET-semHET;
    
    % Make the plot
    subplot(n_k,1,j); hold on
    
    % Find Min/Max
%     maxx = max(nanmean(data2))+0.3;
%     minn = min(nanmean(data2))-0.15;

%      maxx = 40;% 1.1; % 0.55
%     minn = 0; %-0.25; %-0.08

    % Plot SEM
    plot(x, y1, 'w')
    plot(x, y2, 'w')
    patch([x fliplr(x)], [y1 fliplr(y2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    % Plot mean
    plot(mHET, 'r', 'LineWidth', 1.2);
    box off
%     ylim([minn maxx])
    ax = gca;
    ax.TickDir = 'out';
    ax.XAxis.Visible = 'off';
%     if j~= n_k
%     ax.YAxis.Visible = 'off';
%     end 
    ax.LineWidth = 1;
    
end
f4 = gcf;
% f4.Position = [1056  1 256 804];
f4.Position = [ 633    23   263   782]; % [237    41   272   764]; 




%% FULL STIMULUS - firing during grey screen - ACCLIM - at beginning of recording. 

% Full_Table.P_VALUE = Loom_Table.P_VALUE;
% Full_Table.Depth = Loom_Table.Depth;
% 
% allWT = find(Full_Table.Geno == 1 & Full_Table.P_VALUE<0.05 & Full_Table.Depth<0 & Full_Table.Depth> -1000);
% allHET = find(Full_Table.Geno == 2 & Full_Table.P_VALUE<0.05 & Full_Table.Depth<0 & Full_Table.Depth> -1000);

%% Spikes per time = Hz - per second. 

all_animals = unique(Full_Table.Ani);
WTVALS = []; 
HETVALS = [];

for kk = 1:numel(all_animals)
    ani = all_animals(kk);
    
    ani_rows = find(Full_Table.Ani == ani & Full_Table.Depth<0); % & Full_Table.Depth>-300); % & Full_Table.P_VALUE<0.05
    if ~isempty(ani_rows)
        vals = Full_Table.MeanSPT(ani_rows);
        m_ani = nanmean(vals);
        
        g = Full_Table.Geno(ani_rows(1));
        if g == 1
            WTVALS = [WTVALS, m_ani];
        else
            HETVALS = [HETVALS, m_ani];
        end
    end
    
end

% BAR PLOT - POINTS PER ANI

% Ptchd1 
col = [255/255 114/255 32/255]; 

close
semwt = nanstd(WTVALS)/sqrt(numel(WTVALS)); 
semhet = nanstd(HETVALS)/sqrt(numel(HETVALS)); 

figure
bar(1, nanmean(WTVALS), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
bar(2, nanmean(HETVALS),'FaceColor', col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

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
   
    marker = '.';
    scatter(xval, yval, 1000, col, marker, 'jitter', 'on');
    hold on 
end 

% axis([0 3 -1000 0])
axis([0 3 0 2])
box off
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.02, 0.02];
ax.LineWidth= 1.5;

errorbar(1, nanmean(WTVALS), semwt, 'k', 'LineWidth', 1.5);
errorbar(2, nanmean(HETVALS), semhet, 'Color', col, 'LineWidth', 1.5);

f = gcf;
f.Position = [680   750   230   348]; 

nanmean(WTVALS)
nanmean(HETVALS)
[p, h] = ttest2(WTVALS, HETVALS)
[p, h] = ranksum(WTVALS, HETVALS)

% Only for pooled 
% [h, p] = kstest2(WTVALS, HETVALS)

%%  POOLED

wt_rows = find(Full_Table.Geno == 1  & Full_Table.Depth<0 & Full_Table.Depth>-300 & Full_Table.Ani ~= 1387);
het_rows = find(Full_Table.Geno == 2   & Full_Table.Depth<0 & Full_Table.Depth>-300 & Full_Table.Ani ~= 1387);

wtvals = Full_Table.Rng(wt_rows);
m_wt = nanmean(wtvals);

hetvals = Full_Table.Rng(het_rows);
m_het = nanmean(hetvals);


nWT = numel(wtvals);
nHET = numel(hetvals);

mWT = nanmean(wtvals)
mHET = nanmean(hetvals)
[p, h] = ranksum(wtvals, hetvals)
[h, p] = kstest2(wtvals, hetvals)



%% PLOT

xWT = ones(nWT, 1);
xHET = ones(nHET, 1)*2;

d = vertcat(wtvals, hetvals);
gp = vertcat(xWT, xHET);

b = boxplot(d, gp, 'Colors', [0 0 0; 0 0 0], 'Symbol','w.');
set(b, 'linew', 1.2);
ylim([-1 20])
box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
f = gcf;
f.Position = [1053  427  172  257]; 

nWT = numel(wtvals)
nHET = numel(hetvals)
