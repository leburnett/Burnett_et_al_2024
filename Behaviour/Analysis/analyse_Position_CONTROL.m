%% Analyse CONTROL aspects of loom escape response
% Created by Burnett - 07/06/2022

% ALL_XYLOOM_POS_TABLE

% data stored here: /Users/lauraburnett/Documents/Burnett_etal/DATA/SETD5/BEHAVIOUR/Position
% FULL_XY_TABLE

% 1200 values - 
% Loom happens at value 600! 

%% Add new columns with:
% 1 - maximum speed from loom to end
% 2 - Time to maximum speed
% 3 - Distance from centre when loom starts
% 4 - Distance from shelter when loom starts

n = height(FULL_XY_TABLE);
row_loom = 600; 

all_days = unique(FULL_XY_TABLE.Date);


for i = 1:n 
    
    datee = FULL_XY_TABLE.Date{i};
    drow = find(string(all_days) == datee);
    
    % ADD DAY! 
    if sum(ismember([1,6,11,16,21], drow)) == 1
        FULL_XY_TABLE.DAY{i} = 1; 
    elseif sum(ismember([2,7,12,17,22], drow)) == 1
        FULL_XY_TABLE.DAY{i} = 2; 
    elseif sum(ismember([3,8,13,18,23], drow)) == 1
        FULL_XY_TABLE.DAY{i} = 3; 
    elseif sum(ismember([4,9,14,19,24], drow)) == 1
        FULL_XY_TABLE.DAY{i} = 4; 
    elseif sum(ismember([5,10,15,20,25], drow)) == 1
        FULL_XY_TABLE.DAY{i} = 5; 
    end 
        
    
    sp = FULL_XY_TABLE.SPEED{i};
    acc = FULL_XY_TABLE.ACC{i};
    
    FULL_XY_TABLE.SP_AT_LOOM{i} = sp(600); 
    FULL_XY_TABLE.ACC_AT_LOOM{i} = acc(600); 
    
    
    maxsp = max(sp(600:end));
    valmaxsp =  find(sp(600:end) == maxsp);
    if numel(valmaxsp) > 1
        valmaxsp = valmaxsp(1);
    end 
    t2maxsp = valmaxsp/60; 
    
    maxacc = max(acc(600:end));
    valmaxacc =  find(acc(600:end) == maxacc); 
    if numel(valmaxacc) > 1
        valmaxacc = valmaxacc(1);
    end 
    t2maxacc = valmaxacc/60; 
    
    d_cent = FULL_XY_TABLE.DIST_C{i};
    dc_start = d_cent(600);
    
    d_shelt = FULL_XY_TABLE.DIST_SH{i};
    ds_start = d_shelt(600);
    
    FULL_XY_TABLE.MAXSP{i} = maxsp; 
    FULL_XY_TABLE.T2M{i} = t2maxsp; 
    
    FULL_XY_TABLE.MAXACC{i} = maxacc; 
    FULL_XY_TABLE.T2A{i} = t2maxacc; 
    
    FULL_XY_TABLE.DC_START{i} = dc_start; 
    FULL_XY_TABLE.DS_START{i} = ds_start; 
    
    % For finding out the directedness of the trajectory back to the
    % shelter: 
    
    % Find row when mouse back in shelter. 
    rowbackshelt = find(d_shelt(600:end)<0);
    if ~isempty(rowbackshelt)
        rowbackshelt = rowbackshelt(1)+600;
    else
        rowbackshelt = 1200;
    end 
        FULL_XY_TABLE.ROWSHELT{i} = rowbackshelt;
 
     % X, Y position of mouse when loom happened:
     xvals = FULL_XY_TABLE.X{i};
     xloom = xvals(600);
     xshelt = xvals(rowbackshelt);
     yvals = FULL_XY_TABLE.Y{i};
     yloom = yvals(600);
     yshelt = yvals(rowbackshelt);
     
     min_dist_loom_shelt = pdist([xloom, yloom; xshelt, yshelt]); % minimum distance from shelter edge when loom started
     real_dist_loom_shelt = sum(sp(600:rowbackshelt));  % add all the distances between frames from the xy_array between these two rows to find the actual distance travelled.
     dist_ratio = 1 - ((real_dist_loom_shelt-min_dist_loom_shelt) /(real_dist_loom_shelt+min_dist_loom_shelt));
     
     FULL_XY_TABLE.MIN_DIST{i} = min_dist_loom_shelt;
     FULL_XY_TABLE.TRUE_DIST{i} = real_dist_loom_shelt;
     FULL_XY_TABLE.DIST_RATIO{i} = dist_ratio; 
     
     if FULL_XY_TABLE.Geno{i} == "wt"
         FULL_XY_TABLE.GENO{i} = 1;
     else
         FULL_XY_TABLE.GENO{i} = 2;
     end
     
end

 %% Adjust Traj Directedness - invert it! 
n = height(FULL_XY_TABLE);
for i = 1:n 
  
     min_dist_loom_shelt = FULL_XY_TABLE.MIN_DIST{i}; %pdist([xloom, yloom; xshelt, yshelt]); % minimum distance from shelter edge when loom started
     real_dist_loom_shelt = FULL_XY_TABLE.TRUE_DIST{i}  ; % sum(sp(600:rowbackshelt));  % add all the distances between frames from the xy_array between these two rows to find the actual distance travelled.
%      dist_ratio = min_dist_loom_shelt/real_dist_loom_shelt; %((real_dist_loom_shelt-min_dist_loom_shelt)/(real_dist_loom_shelt+min_dist_loom_shelt));
%        dist_ratio = 1 - ((real_dist_loom_shelt-min_dist_loom_shelt)/(real_dist_loom_shelt+min_dist_loom_shelt));
     
       dist_ratio = min_dist_loom_shelt/real_dist_loom_shelt;
       FULL_XY_TABLE.DIST_RATIO2{i} = dist_ratio;   
%        FULL_XY_TABLE.DIST_RATIO{i} = dist_ratio;   
        
end 


%% ANALYSIS:

allWT = find(cell2mat(FULL_XY_TABLE.GENO) ==1 & cell2mat(FULL_XY_TABLE.DS_START)>9);
allHET = find(cell2mat(FULL_XY_TABLE.GENO) ==2 & cell2mat(FULL_XY_TABLE.DS_START)>9);

% Only DAY 5
allWT = find(cell2mat(FULL_XY_TABLE.DAY)==5 &  cell2mat(FULL_XY_TABLE.GENO)==1);
allHET = find(cell2mat(FULL_XY_TABLE.DAY)==5 &  cell2mat(FULL_XY_TABLE.GENO)==2);

% Only 1st Trial
% dataWT = cell2mat(FULL_XY_TABLE.DIST_RATIO(allfirstWT));
% dataHET = cell2mat(FULL_XY_TABLE.DIST_RATIO(allfirstHET));

dataWT = cell2mat(FULL_XY_TABLE.DIST_RATIO(allWT));
dataHET = cell2mat(FULL_XY_TABLE.DIST_RATIO(allHET));

nanmean(dataWT)
nanmean(dataHET)
[p,h] = ranksum(dataWT, dataHET)
[h, p] = kstest2(dataWT, dataHET)

% POSITION OF MOUSE WHEN LOOM STARTS:
% 1 - PLOT DOT - WT / HET

% Boxplot + scatter 
% 2 - Distance from shelter
% 3 - Distance from centre
% 4 - Dshelt - Day 1 only
% 5 - DCentre - D1 only

% 6 - Plot trajectories of first escape from each mouse 
% 7 - Plot trajectoires of the last escape from each mouse. 

% find the row when mouse back in shelter. 
% Calculate distance travelled between loom start and position when mouse
% back in shelter. 
% Calculate minimum distance back to shelter. 
% Calulate difference in trajectory. 

        
%%  PLOT TRAJECTORIES
col = [1 0.4471, 0.1255]; %Ptchd1

% Find only the first ever presentations of the loom:
% allfirstWT = find(cell2mat(FULL_XY_TABLE.DAY)==1 & string(FULL_XY_TABLE.Exp)== "02_Loom" & cell2mat(FULL_XY_TABLE.GENO)==1 & string(FULL_XY_TABLE.Loom)=="01");
% allfirstHET = find(cell2mat(FULL_XY_TABLE.DAY)==1 & string(FULL_XY_TABLE.Exp)== "02_Loom" & cell2mat(FULL_XY_TABLE.GENO)==2 & string(FULL_XY_TABLE.Loom)=="01");

% ALL
allfirstWT = find(cell2mat(FULL_XY_TABLE.GENO)==1);
allfirstHET = find(cell2mat(FULL_XY_TABLE.GENO)==2);

% Only DAY 1 
%  allfirstWT = find(cell2mat(FULL_XY_TABLE.DAY)==1 &  cell2mat(FULL_XY_TABLE.GENO)==1);
% allfirstHET = find(cell2mat(FULL_XY_TABLE.DAY)==1 &  cell2mat(FULL_XY_TABLE.GENO)==2);

% % Only DAY 5
%  allfirstWT = find(cell2mat(FULL_XY_TABLE.DAY)==5 &  cell2mat(FULL_XY_TABLE.GENO)==1);
% allfirstHET = find(cell2mat(FULL_XY_TABLE.DAY)==5 &  cell2mat(FULL_XY_TABLE.GENO)==2);


% Only plot the FIRST LOOM - WT

figure
for v = 1:numel(allfirstWT)
     q = allfirstWT(v);
    
     rowshelt = FULL_XY_TABLE.ROWSHELT{q};
     loom_row = 600; 
    
     if q < 8 %56  %Box size was LARGER - 512
         imsize = 518;
         xvals = FULL_XY_TABLE.X{q};
         xvals = xvals *(416/imsize);
         
         xloom = xvals(600);
         xshelt = xvals(rowshelt);
         
         yvals = FULL_XY_TABLE.Y{q};
         yvals = yvals*(416/imsize);
         
         yloom = yvals(600);
         yshelt = yvals(rowshelt);
         
     elseif q >= 9 && q <14 %q >= 72 && q <89
         
         imsize = 524;
         xvals = FULL_XY_TABLE.X{q};
         xvals = xvals *(416/imsize);
         
         xloom = xvals(600);
         xshelt = xvals(rowshelt);
         
         yvals = FULL_XY_TABLE.Y{q};
         yvals = yvals*(416/imsize);
         
         yloom = yvals(600);
         yshelt = yvals(rowshelt);
         
     elseif q>= 14 %89% Newer experiments - box size = 416
         imsize = 416;
         
         xvals = FULL_XY_TABLE.X{q};
         xloom = xvals(600);
         xshelt = xvals(rowshelt);
         yvals = FULL_XY_TABLE.Y{q};
         yloom = yvals(600);
         yshelt = yvals(rowshelt);
     end
     
     
     imsize2 = 416; 
    for j = 600:rowshelt  
        x = xvals(j);
        y = imsize2 -yvals(j); 
        x2 = xvals(j+1);
        y2 = imsize2 - yvals(j+1);
        
        plot([x, x2],[y,y2],'k', 'LineWidth', 1.2)
        hold on 
    end 
    plot(xvals(600), imsize2-yvals(600), 'Color', col, 'Marker', '.', 'MarkerSize', 25)

end 
axis([0 416 0 416])
axis square
ax = gca;
ax.XTick = []; 
ax.YTick = [];

f= gcf;
f.Renderer = "painters"; 

%% HETS
figure
for v = 1:numel(allfirstHET)
     q = allfirstHET(v);
    
     rowshelt = FULL_XY_TABLE.ROWSHELT{q};
     loom_row = 600; 
    
     if q < 8 %56  %Box size was LARGER - 512
         imsize = 518;
         xvals = FULL_XY_TABLE.X{q};
         xvals = xvals *(416/imsize);
         
         xloom = xvals(600);
         xshelt = xvals(rowshelt);
         
         yvals = FULL_XY_TABLE.Y{q};
         yvals = yvals*(416/imsize);
          
         yloom = yvals(600);
         yshelt = yvals(rowshelt);
         
     elseif q >= 9 && q <14 %q >= 72 && q <89
         
         imsize = 524;
         xvals = FULL_XY_TABLE.X{q};
         xvals = xvals *(416/imsize);
         
         xloom = xvals(600);
         xshelt = xvals(rowshelt);
         
         yvals = FULL_XY_TABLE.Y{q};
         yvals = yvals*(416/imsize);
         
         yloom = yvals(600);
         yshelt = yvals(rowshelt);
         
     elseif q>= 14 %89% Newer experiments - box size = 416
         imsize = 416;
         
         xvals = FULL_XY_TABLE.X{q};
         xloom = xvals(600);
         xshelt = xvals(rowshelt);
         yvals = FULL_XY_TABLE.Y{q};
         yloom = yvals(600);
         yshelt = yvals(rowshelt);
     end
     
     
     imsize2 = 416; 
     if rowshelt<1200
         for j = 600:rowshelt
             x = xvals(j);
             y = imsize2 -yvals(j);
             x2 = xvals(j+1);
             y2 = imsize2 - yvals(j+1);
             
             plot([x, x2],[y,y2],'Color', col, 'LineWidth', 1.2)
             hold on
         end
     elseif rowshelt == 1200
         for j = 600:rowshelt-1
             x = xvals(j);
             y = imsize2 -yvals(j);
             x2 = xvals(j+1);
             y2 = imsize2 - yvals(j+1);
             
             plot([x, x2],[y,y2],'Color', col, 'LineWidth', 1.2)
             hold on
         end
     end 
    plot(xvals(600), imsize2-yvals(600), 'k.', 'MarkerSize', 25)

end 
axis([0 416 0 416])
axis square
ax = gca;
ax.XTick = []; 
ax.YTick = [];

f= gcf;
f.Renderer = "painters"; 

%% PLOT ONLY THE POSITION OF THE MOUSE WHEN THE LOOM STARTED

loom_row = 600; 

figure
for q = 1:n
    
    genooo = FULL_XY_TABLE.GENO{q};
    
    if genooo == 1
        col = 'k';
    else
        col = 'm';
    end 
    
    
%      if q < 56  %Box size was LARGER - 512
%          imsize = 518;
%          xvals = FULL_XY_TABLE.X{q};
%          xvals = xvals *(416/imsize);
%          xloom = xvals(600);
%          
%          yvals = FULL_XY_TABLE.Y{q};
%          yvals = yvals*(416/imsize);
%          yloom = yvals(600);
%          
%      elseif q >= 56 && q <89 || q>219
%          
%          imsize = 524;
%          xvals = FULL_XY_TABLE.X{q};
%          xvals = xvals *(416/imsize);
%          xloom = xvals(600);
%          
%          yvals = FULL_XY_TABLE.Y{q};
%          yvals = yvals*(416/imsize);
%          yloom = yvals(600);
%          
%      elseif q>=89% Newer experiments - box size = 416 
         imsize = 416;
         
         xvals = FULL_XY_TABLE.X{q};
         xloom = xvals(600);
         xshelt = xvals(rowbackshelt);
         yvals = FULL_XY_TABLE.Y{q};
         yloom = yvals(600);
         yshelt = yvals(rowbackshelt);
%      end
     
    imsize2 = 416; 
    plot(xvals(600), imsize2-yvals(600), 'Marker', '.', 'Color', col, 'MarkerSize', 25)
    hold on 
end 

ax = gca;
ax.XTick = []; 
ax.YTick = [];
axis square
axis([0 416 0 416])

%% BOXPLOTS - ALL TRIALS! 

% Col 17 = Distance Centre
% Col 18 = Distance Shelter
% Col 22 = Directedness

val = 18; % Column you want to assess. 

n = height(FULL_XY_TABLE);

dataWT = [];
dataHET = []; 

for i = 1:n
        if (FULL_XY_TABLE.GENO{i}) == 1 %&& string(FULL_XY_TABLE.Exp{i}) == "01_Loom" && (FULL_XY_TABLE.DAY{i}) == 1
            G = cell2mat(FULL_XY_TABLE{i, val});
            dataWT = vertcat(dataWT, G);
        elseif (FULL_XY_TABLE.GENO{i}) == 2 %&& string(FULL_XY_TABLE.Exp{i}) == "01_Loom" && (FULL_XY_TABLE.DAY{i}) == 1
            F = cell2mat(FULL_XY_TABLE{i, val});
            dataHET = vertcat(dataHET, F);
        end
end 

% PLOT 
n_wt = numel(dataWT);
n_het = numel(dataHET);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

geno = cell2mat(FULL_XY_TABLE.GENO); 
var = cell2mat(FULL_XY_TABLE{:, val});

figure
hold on 
scatter(x1, dataWT','SizeData', 175, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.05)
scatter(x2, dataHET', 'SizeData', 175,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05)
b = boxplot(var, geno, 'Colors', 'k');
set(b , 'LineWidth', 1.5)

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
ylim([0 6])

f = gcf;
f.Position = [704   288   275   491]; %[704   207   355   572]; 

% Stats
nanmean(dataWT)
nanmean(dataHET)
[p,h] = ranksum(dataWT, dataHET)
[h,p] = kstest2(dataWT, dataHET)



%% ana-animals - animal averages over all trials. FOR STATS & BOXPLOTS

all_animals = string(unique(FULL_XY_TABLE.Animal));
n_animals = numel(all_animals);

analysis_animals = zeros(n_animals, 10);

for j = 1:n_animals
    
    ani = all_animals{j};
    all_ani = find(string(FULL_XY_TABLE.Animal) == ani & cell2mat(FULL_XY_TABLE.DS_START) >10); % & cell2mat(FULL_XY_TABLE.DAY) == 5 % & string(FULL_XY_TABLE.Exp) == "01_Loom"
    
    if isempty(all_ani)
        all_ani2 = find(string(FULL_XY_TABLE.Animal) == ani);
        geno = FULL_XY_TABLE.GENO{all_ani2(1)};
    else 
        geno = FULL_XY_TABLE.GENO{all_ani(1)};
    end 
    
    % Assign animal number to  first column of table. 
    analysis_animals(j,1) = str2double(ani(3:end)); 
    
    analysis_animals(j,2) = geno; 
    
    cv = 3; 
    
    for k = [17, 18, 22, 24, 26, 13, 14] %[24, 25, 29, 13, 32, 19, 20, 21] %Ptchd1 [17, 18, 22, 13, 14]% Cul3 - [20, 21, 25, 27, 14, 15, 16, 17]%[17,18,22] % 24, 26
        
        if k == 24 %13 %27 %24
            all_vals = abs((FULL_XY_TABLE{all_ani, k}));
        else
            all_vals = cell2mat(FULL_XY_TABLE{all_ani, k});
        end 
        
        analysis_animals(j,cv) = nanmean(all_vals);
        cv = cv+1;
    end
    
end 

Animal = analysis_animals(:,1);
Geno = analysis_animals(:,2);
DC = analysis_animals(:,3);
DS = analysis_animals(:,4);
DRATIO = analysis_animals(:,5);
ANG = analysis_animals(:,6);
SPEED_AT = analysis_animals(:, 7);
% ACC_AT = analysis_animals(:,8);
MAXSP = analysis_animals(:,8);
T2M = analysis_animals(:,9);


ana_animals = table(Animal, Geno, DC, DS, DRATIO, MAXSP, T2M, ANG, SPEED_AT); %, ANG, SPEED_AT); ANG, SPEED_AT, ACC_AT,

%% ALL TRIALS

val = 4; % Column you want to assess. 

n = height(ana_animals);

dataWT = [];
dataHET = []; 

for i = 1:n
        if (ana_animals.Geno(i)) == 1 
            G = (ana_animals{i, val});
            dataWT = vertcat(dataWT, G);
        elseif (ana_animals.Geno(i)) == 2 
            F = (ana_animals{i, val});
            dataHET = vertcat(dataHET, F);
        end
end 

% PLOT 
n_wt = numel(dataWT);
n_het = numel(dataHET);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

geno = (ana_animals.Geno); 
var = table2array(ana_animals(:, val));

figure

scatter(x1, dataWT','SizeData', 175, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.05)
hold on
scatter(x2, dataHET', 'SizeData', 175,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05)
b = boxplot(var, geno, 'Colors', 'k');
set(b , 'LineWidth', 1.5)

xticks([1,2])
xticklabels({''})
ax = gca;
% ax.FontSize = 30;
box off
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 2;
xlim([0.5 2.5])

ylim([0 30])
% ylim([0 40])
% ylim([0 20])
ylim([0 0.4])

f = gcf;
f.Position = [704   288   275   491]; %[704   207   355   572]; 

nanmean(dataWT)
nanmean(dataHET)
[p,h] = ranksum(dataWT, dataHET)
% [p,h] = ttest2(dataWT, dataHET)



%% FIRST TRIALS 

val = 5; % Column you want to assess. 

n = height(ana_animals5);

dataWT = [];
dataHET = []; 

for i = 1:n
        if (ana_animals5.Geno(i)) == 1 
            G = (ana_animals5{i, val});
            dataWT = vertcat(dataWT, G);
        elseif (ana_animals5.Geno(i)) == 2 
            F = (ana_animals5{i, val});
            dataHET = vertcat(dataHET, F);
        end
end 

% PLOT 
n_wt = numel(dataWT);
n_het = numel(dataHET);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

geno = (ana_animals5.Geno); 
var = table2array(ana_animals5(:, val));

figure
b = boxplot(var, geno, 'Colors', 'k');
set(b , 'LineWidth', 1.5)
hold on 
scatter(x1, dataWT','SizeData', 175, 'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', 0.05)
scatter(x2, dataHET', 'SizeData', 175,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05)

xticks([1,2])
xticklabels({''})
ax = gca;
ax.FontSize = 30;
box off
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 2;
xlim([0.5 2.5])

ylim([0 12])
% ylim([0 40])
% ylim([0 40])

f = gcf;
f.Position = [704   288   275   491]; %[704   207   355   572]; 

nanmean(dataWT)
nanmean(dataHET)
[p,h] = ranksum(dataWT, dataHET)





%% ALL TRIALS  - scatter plot with correlations. 

n = height(FULL_XY_TABLE);

figure
for q = 1:n
    
    if  FULL_XY_TABLE.DS_START{q}>10 % FULL_XY_TABLE.ReturnToShelter{q} == 1 &
        genooo = FULL_XY_TABLE.GENO{q};
        
        if genooo == 1
            col = 'k';
        else
            col = 'r'; %[1 0.4471, 0.1255]; 
        end
        
        xval = cell2mat(FULL_XY_TABLE.DIST_RATIO2(q));
        yval = FULL_XY_TABLE.MAXSP{q};
        
        
        plot(xval, yval, 'Marker', '.', 'Color', col, 'MarkerSize', 15)
        hold on
    end
end

ax = gca;
box off

xlim([0 180])
xlim([0 0.3])
xlim([5 30])
xlim([0 32])
% ylim([20 110])

% maxsp or t2m 
ylim([0 180])
ylim([-1 11])
ylim([0 0.4])
ylim([0 32])

ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];

% DATA FOR CORRELATION

dataWTx = [];
dataWTy = [];
dataHETx = [];
dataHETy = [];

for q = 1:n
    if  FULL_XY_TABLE.DS_START{q}>10 % FULL_XY_TABLE.ReturnToShelter{q} == 1 &
        genooo = FULL_XY_TABLE.GENO{q};
        
        if genooo == 1
            xval = cell2mat(FULL_XY_TABLE.DIST_RATIO2(q));
            yval = FULL_XY_TABLE.MAXSP{q};
            dataWTx = vertcat(dataWTx, xval);
            dataWTy = vertcat(dataWTy, yval);
        else
            
            xval = cell2mat(FULL_XY_TABLE.DIST_RATIO2(q));
            yval = FULL_XY_TABLE.MAXSP{q};
            dataHETx = vertcat(dataHETx, xval);
            dataHETy = vertcat(dataHETy, yval);
        end
    end
end


[xData, yData] = prepareCurveData( dataWTx, dataWTy );

% Set up fittype and options.
ft = fittype( 'poly1' );
[fitresult{1}, gof(1)] = fit( xData, yData, ft );
h = plot(fitresult{1}, 'k');


[xData, yData] = prepareCurveData( dataHETx, dataHETy );
ft = fittype( 'poly1' );
[fitresult{2}, gof(2)] = fit( xData, yData, ft );
h = plot(fitresult{2}); %col = [1 0.4471, 0.1255];

legend off

f = gcf;
f.Position =  [52   404   557   262]; %[ 407   389   377   273];


% CORRELATION STATS 

[R, P, RL, RU] = corrcoef(dataWTx, dataWTy)
[rho, p] = corr(dataWTx,dataWTy, 'Type', 'Spearman')

[R, P, RL, RU] = corrcoef(dataHETx, dataHETy)
[rho, p] = corr(dataHETx,dataHETy, 'Type', 'Spearman')



%% CORRELATION - ANIMAL AVERAGES

n = height(ana_animals);

figure
for q = 1:n
    
    genooo = ana_animals.Geno(q);
    
    if genooo == 1
        col = 'k';
    else
        col = 'r'; %[1 0.4471, 0.1255];
    end 
    
    yval = ana_animals.T2M(q);
    xval = ana_animals.SPEED_AT(q);
    
    plot(xval, yval, 'Marker', '.', 'Color', col, 'MarkerSize', 25)
    hold on 
end 

ax = gca;
box off
xlim([90 150])

% T2m or maxsp
ylim([0 100])
ylim([0 5])
% ylim([40 80])

ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];
hold on 

% DATA FOR CORRELATION

dataWTx = [];
dataWTy = [];
dataHETx = [];
dataHETy = [];

for q = 1:n
        genooo = ana_animals.Geno(q);
        
        if genooo == 1
            xval = (ana_animals.SPEED_AT(q));
            yval = ana_animals.T2M(q);
            dataWTx = vertcat(dataWTx, xval);
            dataWTy = vertcat(dataWTy, yval);
        else
            xval = (ana_animals.SPEED_AT(q));
            yval = ana_animals.T2M(q);
            dataHETx = vertcat(dataHETx, xval);
            dataHETy = vertcat(dataHETy, yval);
        end
end

[xData, yData] = prepareCurveData( dataWTx, dataWTy );

% Set up fittype and options.
ft = fittype( 'poly1' );
[fitresult{1}, gof(1)] = fit( xData, yData, ft );
h = plot( fitresult{1}, 'k');


[xData, yData] = prepareCurveData( dataHETx, dataHETy );
ft = fittype( 'poly1' );
[fitresult{2}, gof(2)] = fit( xData, yData, ft );
h = plot( fitresult{2}, 'r');

legend off

f = gcf;
f.Position =  [52   404   557   262]; %[ 407   389   377   273];


% CORRELATION STATS 

[R, P, RL, RU] = corrcoef(dataWTx, dataWTy)
[rho, p] = corr(dataWTx,dataWTy, 'Type', 'Spearman')

[R, P, RL, RU] = corrcoef(dataHETx, dataHETy)
[rho, p] = corr(dataHETx,dataHETy, 'Type', 'Spearman')






%% Looking at the correlation between the speed of the mouse when the loom starts versus max speed / time to max. 


figure
for i = 1:height(all_xy_analysis)
    
    if all_xy_analysis.ReturnToShelter{i}==1
    g = string(all_xy_analysis.Geno{i});
    
    if g == "wt"
        col = 'k';
    elseif g == "het"
        col = 'r';
    end 
    
    xval = abs(all_xy_analysis.AngAtLoom(i));
    yval = all_xy_analysis.MaxSpEscape{i};
    
    if xval ~= 0 
    plot(xval, yval, 'Color', col, 'Marker', '.', 'MarkerSize', 15)
    end 
    hold on 
    end 
    
end 

ax = gca;
box off
xlim([0 180])
ylim([0 110])
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];
f = gcf;
f.Position =  [52   404   557   262]; %[ 407   389   377   273];

% % % % % % % % 

[xData, yData] = prepareCurveData( dataWTx, dataWTy );

% Set up fittype and options.
ft = fittype( 'poly1' );
[fitresult{1}, gof(1)] = fit( xData, yData, ft );
h = plot( fitresult{1}, 'k');


[xData, yData] = prepareCurveData( dataHETx, dataHETy );
ft = fittype( 'poly1' );
[fitresult{2}, gof(2)] = fit( xData, yData, ft );
h = plot( fitresult{2}, 'r');

legend off

% % % % % % . % % % % % % 

dataWTx = [];
dataWTy = [];
dataHETx = [];
dataHETy = [];

for q = 1:height(all_xy_analysis)
    
    if all_xy_analysis.ReturnToShelter{q} ==1
        
        g = string(all_xy_analysis.Geno{q});
        
        if g == "wt"
            xval = abs(all_xy_analysis.AngAtLoom(q));
            yval = all_xy_analysis.MaxSpEscape{q};
            dataWTx = vertcat(dataWTx, xval);
            dataWTy = vertcat(dataWTy, yval);
        elseif g == "het"
            xval = abs(all_xy_analysis.AngAtLoom(q));
            yval = all_xy_analysis.MaxSpEscape{q};
            dataHETx = vertcat(dataHETx, xval);
            dataHETy = vertcat(dataHETy, yval);
        end
    end 
end


% Least squares regression
b1 = dataWTx\dataWTy;

% Compute a linear regression: 
p = polyfit(dataWTx, dataWTy, 1)


% Correlation coefficients
% Pearsons correlation - linear and normal
% Spearmans - only assumption is general increase or decrease

%  the correlation coefficient only measures the strength of the relationship, not its magnitude.

% the R² value tells you the percentage of variation in Y that is explained by variation in X

% A = horzcat(dataWTx, dataWTy);
% [R, P] = corrcoef(A)
% [RHO,PVAL] = corr(dataWTx,dataWTy,'Type','Spearman')
% 
% WT 
[RHO,PVAL] = corr(dataWTx,dataWTy,'Type','Pearson')
% HET
[RHO,PVAL] = corr(dataHETx,dataHETy,'Type','Pearson')


% % % % %  Way of reporting these results:
% r(degrees of freedom) = the r stat, p = p value. 

% The r-square value in curve fitting is the squared version of the RHO
% value from the 'corr' function above. 


%% Angle of head/ body from shelter when the loom stimulus starts:

% Use all_xy_analysis
% Data found here: /Users/lauraburnett/Data_Analysis_Mac/Loom_Behaviour/Setd5/ALL4Cohorts
%load('220616_Setd5_4Cs_Angular_and_xy_analysis.mat', 'all_xy_analysis');
% For Ptchd1
% '/Users/lauraburnett/Data_Analysis_Mac/Loom_Behaviour/Ptchd1/Ptchd1_DLC'
% load 'DLC_ANALYSIS...'

%% ADD ANGLE DATA TO FULL_XY_TABLE

% all_ani = find(all_xy_analysis.Animal == "GN4472" | all_xy_analysis.Animal == "GN6558");
% all_xy_analysis(all_ani, :) = [];
% 
% all_ani = find(FULL_XY_TABLE.Animal == "GN4472" | FULL_XY_TABLE.Animal == "GN6558");
% FULL_XY_TABLE(all_ani, :) = [];

%  FULL_XY_TABLE.ANG = all_xy_analysis.AngAtLoom;
% FULL_XY_TABLE.RETURN = all_xy_analysis.ReturnToShelter;

 % folder: /Users/lauraburnett/Documents/Burnett_etal/DATA/SETD5/BEHAVIOUR/Position
% save('220617_FULL_XYLOOM_TABLE_withANG_C1-C4.mat', 'FULL_XY_TABLE');

FULL_XY_TABLE.ANG = dlc_analysis.AngAtLoom;
FULL_XY_TABLE.ANG360 = dlc_analysis.AngAtLoom360;
FULL_XY_TABLE.ANG_LOOM = dlc_analysis.AngLoomAtLoom;
FULL_XY_TABLE.ANG_LOOM360 = dlc_analysis.AngLoomAtLoom360;

FULL_XY_TABLE.ReturnToShelter = all_xy_analysis.ReturnToShelter;


%% HISTOGRAM - DISTRIBUTIONS - ALL TRIALS - ANGLES OF HEAD / BODY AT LOOM - WRT SHELTER

dataWT = []; 
dataHET = [];

for i = 1:height(FULL_XY_TABLE)
    
    if cell2mat(FULL_XY_TABLE.DC_START(i))>10 %& cell2mat(FULL_XY_TABLE.ReturnToShelter(i)) == 1 % & cell2mat(FULL_XY_TABLE.DAY(i)) == 5 
       
        if cell2mat(FULL_XY_TABLE.GENO(i)) == 1
%             xval = 1;
%             col = 'k';
            yval = cell2mat(FULL_XY_TABLE.SpeedATLoom(i));
            dataWT = [dataWT, yval];
        else
%             xval = 2;
%             col = 'm';
            yval = cell2mat(FULL_XY_TABLE.SpeedATLoom(i));
            dataHET = [dataHET, yval];
        end
        
    end
    %      plot(xval, yval, 'Marker', 'o', 'Color', col, 'MarkerSize', 10)
    %      hold on
end

% end 

% angls = 0:10:180;
% angls = 0:2:32;
% angls = -30:2:30;
angls = 0:0.02:0.4;
% max(dataWT)

% HISTOGRAM - VARIATION

figure
histogram(dataWT, angls, 'FaceColor', [0.3 0.3 0.3], 'Normalization', 'pdf')
hold on 
histogram(dataHET, angls, 'FaceColor', col,  'Normalization', 'pdf')
box off
ax = gca;
ax.TickDir = 'out'; 
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
f = gcf;
f.Position = [ 440   600   352   198];

yticks([0, 0.05, 0.1]); %, 0.04, 0.06])
xticks([0:5:30])

% STATS - Variance?
nanmean(dataWT)
nanmean(dataHET)
[h,p] = kstest2(dataWT, dataHET)


% rhow = ones(1, numel(dw));
% rhoh = ones(1, numel(dh))*1.5;
% 
% figure
% polarplot(dw, rhow, 'o', 'Color', [0.7 0.7 0.7]); hold on ; polarplot(dh, rhoh, 'o', 'Color', 'r')









