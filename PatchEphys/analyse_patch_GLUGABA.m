%% Analyse patch data according to their firing at baseline - GLU/ GABA.
% Glu = silent at baseline
% GABA = tonically firing. 

% Cell type details: 
% '/Users/lauraburnett/Data_Analysis_Mac/PATCH/CellType_Details.mat'

% Arrays: details_wt/ details_HET 
% First row = COHORT
% second row  = ANIMAL
% third row = CELL
% fourth row = GLU? ( 1 = GLU, 0 = GABA)

% DATA IN HERE: /Users/lauraburnett/Data_Analysis_Mac/PATCH/ANALYSIS/FiringFreq-Baseline

%% CELLS GROUPED BY 'TYPE' - GLU/GABA
% Plot the firing at baseline to determine if they are silent/ tonic. 
clear
close all

files = dir('*.abf');
nfiles = length(files);


close
fn = files(i).name;
[d,si,h]=abfload(fn);

figure; plot(d)
f = gcf;
f.Position = [105   521   1589   291];


%% BASELINE FIRNG BAR CHART with CELLS PLOTTED ON TOP

valsWT1 = data(:,1);
valsHET1 = data(:,2);
valsWT2 = data(:,3);
valsHET2 = data(:,4);

x = [1,2,3,4];
y = [nanmean(valsWT1), nanmean(valsHET1), nanmean(valsWT2), nanmean(valsHET2)];

col = [0.25 0.25 0.25];

figure
b = bar(x,y);
hold on
scatter(ones(numel(valsWT1),1), valsWT1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsHET1),1)*2, valsHET1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsWT2),1)*3, valsWT2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsHET2),1)*4, valsHET2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
b.FaceColor = [0.7 0.7 0.7];

box off
ylim([0 250])
xlim([0 5])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.FontSize = 14;

f = gcf;
f.Position = [440   391   306   407];

% STATS AND NUMBERS 

nWT1 = sum(~isnan(valsWT1))
nWT2 = sum(~isnan(valsWT2))
nHET1 = sum(~isnan(valsHET1))
nHET2 = sum(~isnan(valsHET2))

nanmean(valsWT1)
nanmean(valsWT2)
nanmean(valsHET1)
nanmean(valsHET2)

freq = data2(:,1);

[p,tbl,stats] = kruskalwallis(data)
tbl = multcompare(stats)








%% PLOT CURRENT INJ FIRING

close
 figure; plot(DTX_HET(:,15), 'LineWidth', 1.5)
 box off
 ylim([0 100])
 hold on; plot([0, 25], [0 100], 'k:')
 f = gcf;
 f.Position = [680   629   233   176];
 
 %%
 close
 figure; plot(valsWT2(:,18), 'Color', 'k', 'LineWidth', 1.5)
 box off
 ylim([0 100])
 hold on; plot([0, 25], [0 100], 'k:')
 f = gcf;
 f.Position = [680   629   233   176];
 
 
 
 %% 
 data(:, 22) = kmeans(data(:, 1:18), 10);
 
%  
%  evaluation = evalclusters(data(:, 1:18),"kmeans","CalinskiHarabasz","KList",1:30);
%  figure; plot(evaluation)
%  % Looking for value of x with highest value of y. 

 
%  figure
%  for i = 1:4
%     all1 = find(data(:,22)==i);
%     subplot(2,2,i); plot(1:1:23, data(all1, 1:23), 'k')
%  end 
 d2 = [data(:,21), data(:,24)];
d2 = array2table(d2);
[C, ia, ic] = unique(d2(:, [1,2]));
[C, ia, ic] = unique(data(:, [21,24]));

 %% % Plot based on cluster found by k means - colour and marker by geno and cell type. 
 
  figure
 for k = 1:42
     
     cc = data(k, 21);
     if cc == 1
         col = 'k';
     else 
         col = 'r';
     end
     
     gp = data(k,22); 
     
     subplot(2,5,gp);
     
     typee = data(k, 24);
     if typee == 1
         markk = '-';
     else
         markk = ':';
     end 
    
     plot(1:1:18, data(k, 1:18), strcat(col, markk), 'LineWidth', 1.2);
     hold on 
     ylim([0 120])
 end 
 
 %% % Line plot - firing to current injection 10pA steps 
 % 4 plots - based on geno and cell type
 % Each cell's response is in light colour underneath and then mean is on
 % top. 
 % Below is mean / SEM for both 
 
 % data
 % Col 21 = GENO
 % Col 22 = k means group 
 % Col 23 = Type - Glu / GABA / WT / HET
 % COl 24 = Glu/GABA
 
  figure
 for k = 1:42
     
     cc = data(k, 21);
     if cc == 1
         col = [0.8 0.8 0.8];
     else 
         col = [1 0.8 0.8];
     end
     
     typee = data(k, 23);
     if typee ==1
         subplot(2,2,1);
         title('GABA - HET')
         plot([0 20], [0 120], 'k:', 'LineWidth', 1.1)
     elseif typee ==2
         subplot(2,2,2);
         title('GLU - HET')
         plot([0 20], [0 120], 'k:', 'LineWidth', 1.1)
     elseif typee ==3
         subplot(2,2,3);
         title('GABA - WT')
         plot([0 20], [0 120], 'r:', 'LineWidth', 1.1)
     elseif typee ==4
         subplot(2,2,4);
         title('GLU - WT')
         plot([0 20], [0 120], 'r:', 'LineWidth', 1.1)
     end 
    
     plot(1:1:18, data(k, 1:18), 'Color', col, 'LineWidth', 1.2);
     hold on 
     ylim([0 120])
     xlim([0 18])
     xticks([6,11,16])
     box off
     ax = gca;
     ax.TickDir = 'out'; 
     ax.LineWidth = 1;
     ax.TickLength = [0.02 0.02];
 end 
 
 for ii = [1,2,3,4]
     
     allr = find(data(:, 23)==ii);
     d = data(allr, 1:20);
     md = nanmean(d);
     
     if ii ==1 
        subplot(2,2,1)
        hold on; plot(md, 'r', 'LineWidth', 1.5)
        ax = gca; ax.XAxis.Visible = 'off';
     elseif ii ==2
         subplot(2,2,2)
         hold on; plot(md, 'r', 'LineWidth', 1.5)
         ax = gca; ax.YAxis.Visible = 'off';
         ax.XAxis.Visible = 'off';
     elseif ii ==3
         subplot(2,2,3)
        hold on; plot(md, 'k', 'LineWidth', 1.5)
     elseif ii ==4
         subplot(2,2,4)
        hold on; plot(md, 'k', 'LineWidth', 1.5)
        ax = gca; ax.YAxis.Visible = 'off';
     end 
 end 
 
 f = gcf;
 f.Position = [ 676   372   344   299];
 
 
 
 
 
 
 %% MEAN /SEM plot - subplot - wt/het for each group
 
 figure
 
 x = 1:1:20; 
 
  for ii = [1,2,3,4]
     
     allr = find(data(:, 23)==ii);
     d = data(allr, 1:20);
     
     md = nanmean(d);
     n_gp = numel(d(:,1));
     
     semWT = nanstd(d)/sqrt(n_gp);
     y1 = md+semWT;
     y2 = md-semWT;
     
     if ii ==1 
         subplot(1,2,1) % HET - GABA
         
         plot(x, y1, 'w')
         hold on
         plot(x, y2, 'w')
         patch([x fliplr(x)], [y1 fliplr(y2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
         plot(md, 'r:', 'LineWidth', 1.3)
         box off
         
     elseif ii ==2
         subplot(1,2,2) % HET - GLU
      
         plot(x, y1, 'w')
         hold on
         plot(x, y2, 'w')
         patch([x fliplr(x)], [y1 fliplr(y2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
         plot(md, 'r', 'LineWidth', 1.3)
         
         box off
         ax = gca; ax.YAxis.Visible = 'off';

     elseif ii ==3
         subplot(1,2,1) % WT - GABA
       
         plot(x, y1, 'w')
         hold on
         plot(x, y2, 'w')
         patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
         plot(md, 'k:', 'LineWidth', 1.3)
         box off
         ylim([0 100])
        
         
     elseif ii ==4
         subplot(1,2,2) % WT - GLU
       
         plot(x, y1, 'w')
         hold on
         plot(x, y2, 'w')
         patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
         plot(md, 'k', 'LineWidth', 1.3)
         box off
         ylim([0 100])
         
     end 
  end 
 
 f = gcf;
 f.Position = [306   492   467   199];
  

% 
% xlabel('Current (pA)')
% ylabel('Firing Rate (Hz)')
% box off
% set(gca, 'FontSize', 20)
% ax = gca;
% ax.TickDir  = 'out';
% ax.LineWidth = 1.75;
% ax.TickLength = [0.03 0.03];
% xlim([0.5 24])
% ylim([0 100])
% 
%  
 
 
%% PROP OF CELLS GLU / GABA - stacked bar chart - proportions. 

allGLUW = find(details_WT(4, :)==1);
allGABAW = find(details_WT(4, :)==0);

allGLUH = find(details_HET(4, :)==1);
allGABAH = find(details_HET(4, :)==0);

n1 = numel(allGLUW);
n2 = numel(allGABAW);

% GABA CELLS
n3 = numel(allGLUH);
n4 = numel(allGABAH);


prop_TYPE = zeros(2,2);
% PROP GLU - WT
prop_TYPE(1,1) = n1/(n1+n2);
% PROP GABA - WT
prop_TYPE(1,2) = n2/(n1+n2);
% PROP GLU - HET
prop_TYPE(2,1) = n3/(n3+n4);
% PROP GABA - HET
prop_TYPE(2,2) = n4/(n3+n4);

b = bar(prop_TYPE, 'stacked');
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [1 1 1];
b(1).LineWidth = 1.2;
b(2).LineWidth = 1.2;
box off
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];
f = gcf;
f.Position =  [680   529   187   276];

%% PIE CHART


figure
expl = ones(2,1);
X = prop_TYPE(1,:);
ax1 = subplot(1,2,1);
pie(ax1,X, expl)
% title(ax1,'Setd5^+^/^+');
colormap(gray)

Y = prop_TYPE(2,:);
ax2 = subplot(1,2,2);
pie(ax2,Y, expl)
% title(ax2,'Setd5^-^/^-');
colormap(cmap)

cmap = [0 0 0;0.1  0 0; 0.2  0 0; 0.3  0 0; 0.4  0 0; 0.5 0 0; 0.6 0 0; 0.7 0 0; 0.8 0 0; 0.9 0 0;1,0, 0];





%% PLOT JUST MEAN LINE PLOT for each group. 

figure; plot(nanmean(vWT(allGLUW, :)), 'k', 'LineWidth', 1.5);
hold on; plot(nanmean(vWT(allGABAW, :)), 'b', 'LineWidth', 1.5)
hold on; plot(nanmean(vHET(allGLUH, :)), 'r', 'LineWidth', 1.5)
hold on; plot(nanmean(vHET(allGABAH, :)), 'm', 'LineWidth', 1.5)

f = gcf;
f.Position = [557   527   365   289];

box off
xlabel({''})
xticklabels({''})
ylabel({'Firing (Hz)'})

ax = gca;
ax.LineWidth = 1.5;
ax.TickDir = 'out';

%% MEAN /SEM plot

valsWGLU = vWT(allGLUW, :)';
valsWGAB = vWT(allGABAW, :)';

valsHGLU = vHET(allGLUH, :)';
valsHGAB = vHET(allGABAH, :)';

%  1 - GLU CELLS! 

% col = 'm';
v = 21;
inj_steps = curr_inj(3:3+v, 2);

nWT = numel(valsWGLU(1,:)); 
nHET = numel(valsHGLU(1,:));

mean_WT = nanmean(valsWGLU'); 
mean_HET = nanmean(valsHGLU'); 

x = (1:1:v);

semWT = nanstd(valsWGLU')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(valsHGLU')/sqrt(nHET); 
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

xticks([1:4:24])
xticklabels(string(inj_steps(1:4:24)))

xlabel('Current (pA)')
ylabel('Firing Rate (Hz)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.03 0.03];
xlim([0.5 24])
ylim([0 100])

f = gcf;
f.Position = [742   492   304   258];


% ADD GABA CELLS

nWT = numel(valsWGAB(1,:)); 
nHET = numel(valsHGAB(1,:));

mean_WT = nanmean(valsWGAB'); 
mean_HET = nanmean(valsHGAB'); 

x = (1:1:v);

semWT = nanstd(valsWGAB')/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = nanstd(valsHGAB')/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;


plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT(1:v)', 'b', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET(1:v)', 'Color', 'm', 'LineWidth', 1.3)




%% % % % % % % % % % % % % % % % % % % % % GOOD STATS - REPEATED MEASURES - WITHIN MODEL - WITH TUKEY MULTIPLE COMPARISONS  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% SORT DATA BASED ON CELL TYPE.

% Remove groups 1 and 3 - GABA cells - only want GLU
all13 = find(data(:, 23)==1 | data(:, 23)==3);
GABAdata = data(all13, :);

all24 = find(data(:, 23)==2 | data(:, 23)==4);
GLUdata = data(all24, :);

   
   %% REPEATED MEASURES ANOVA - WITHIN DESIGN
   
   datatbl = array2table(data, 'VariableNames', {'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9', 'y10', 'y11', 'y12', 'y13', 'y14', 'y15', 'y16', 'y17', 'y18', 'y19', 'y20', 'Geno','Group', 'Type'});
   wtbl = table('Size', [20 1], 'VariableType', {'double'}, 'VariableName', {'CurrentStep'});
  
   
   % Looking only at the effect of genotype on firing - grouping all cells together
    rm = fitrm(datatbl, 'y1-y20 ~ Geno', 'WithinDesign', wtbl);
   ranovatbl = ranova(rm, 'WithinMode', 'CurrentStep')
    
  % Multiple comparisons - across current steps for genotype - all cells together. The same as above but need to chage the 'Geno' column to categorical not double. 
   rm = fitrm(datatbl2, 'y1-y20 ~ Var24', 'WithinDesign', wtbl);
   tbl = multcompare(rm, 'Var24', 'By', 'CurrentStep')
   
   
   %% Looking at the effect of both genotype and cell type on firing. 
   rm = fitrm(datatbl, 'y1-y20 ~ Geno + Type', 'WithinDesign', wtbl);
   ranovatbl = ranova(rm, 'WithinMode', 'CurrentStep')
    
   
   %% RMANOVA - JUST FOR GLU / GABA data separately. 
   
   GABAdtbl = datatbl(all13, :);
   
   rm = fitrm(GABAdtbl, 'y1-y20 ~ Geno', 'WithinDesign', wtbl);
   ranovatbl = ranova(rm, 'WithinMode', 'CurrentStep')
    
   
   GLUdtbl = datatbl(all24, :);
   
   rm = fitrm(GLUdtbl, 'y1-y20 ~ Geno', 'WithinDesign', wtbl);
   ranovatbl = ranova(rm, 'WithinMode', 'CurrentStep')

   %% MULTIPLE COMPARISONS TESTS - TUKEY - Comparing Genotype at each current step. 

   
   % MULTIPLE COMPARSONS  - TUKEY - ONLY GLU DATA
  rm = fitrm(GLUdtbl2, 'y1-y20 ~ Var24', 'WithinDesign', wtbl);
  tbl = multcompare(rm, 'Var24', 'By', 'CurrentStep')
    
  
     % MULTIPLE COMPARSONS  - TUKEY - ONLY GLU DATA
  rm = fitrm(GABAdtbl2, 'y1-y20 ~ Var24', 'WithinDesign', wtbl);
  tbl = multcompare(rm, 'Var24', 'By', 'CurrentStep')
  
 
    %% STATS - TESTS for normality - sphericity etc. 
    
    %  1 - Test groups for normality with Shapiro Wilks test since samples
    %  are small (<50)  - otherwise would have used Kolmogorov Smirnov
    
    
    % 2 - Bartlett's Test / Levene's test for Homoscedasticity - same standard deviation. 
    
     [p, stats] = vartestn(yy', gp3)
    
    %%     % Split data 
%    
%     x = []; 
%     g = [];
%     for i = 1:31
%         x = [x GLUdata(i, 1:10)]; %concatenate 
%         val = GLUdata(i,23);
%         g = [g repmat(val, 1, 10)];
%     end 
%    
%    [p, tbl, stats] = kruskalwallis(x, g)  

    
 %% Plot firing-current relationship for DTX cells - wT versus HET - don't know type but important to show the responses of individual cells
 
 x = 1:1:21;
 
  figure
  
 % WT cells with DTX 
  subplot(2,3,4)
  hold on
  for i = 1:12
  plot(DTX_WT(1:21, i), 'Color', [0.8 0.8 1], 'LineWidth', 1.5)
  end 
  title('WT-DTX')
  plot([0 20], [0 120], 'k:', 'LineWidth', 1.1)
  xlim([0 20])
  ylim([0 120])
  
  subplot(2,3,1) % HET GABA
  hold on
  for i = 1:13
      if DTX_HET(22, i) == 0
          plot(DTX_HET(1:21, i), 'Color', [1 0.8 1], 'LineWidth', 1.5)
      end
  end
  plot([0 20], [0 120], 'k:', 'LineWidth', 1.1)
  xlim([0 20])
  ylim([0 120])
  title('HET-GABA - DTX')
 
  subplot(2,3,2) % HET GLU
  hold on
  for i = 1:13
      if DTX_HET(22, i) == 1
          plot(DTX_HET(1:21, i), 'Color', [1 0.8 1], 'LineWidth', 1.5)
      end
  end
  plot([0 20], [0 120], 'k:', 'LineWidth', 1.1)
  xlim([0 20])
  ylim([0 120])
  title('HET-GLU - DTX')
  
  subplot(2,3,3) % HET GLU
  hold on
  for i = 1:13
      if DTX_HET(22, i) ~= 1 & DTX_HET(22, i) ~= 0
          plot(DTX_HET(1:21, i), 'Color', [1 0.8 1], 'LineWidth', 1.5)
      end
  end
  plot([0 20], [0 120], 'k:', 'LineWidth', 1.1)
  xlim([0 20])
  ylim([0 120])
  title('HET- NA - DTX')
  
  
%  f = gcf;
%  f.Position = [567   249   227   437];
% %  
 
 
%% All 4 conditions - individual line traces for each cell
    
    figure
      
 % WT cells with DTX 
  subplot(2,2,1)
  hold on
  for i = 1:19
  plot(valsWT2(1:23, i), 'Color', 'k', 'LineWidth', 1.5)
  end 
  title('WT')
  plot([0 20], [0 120], 'r--', 'LineWidth', 2)
  xlim([0 20])
  ylim([0 120])
  
  subplot(2,2,2)
  hold on
  for i = 1:24
  plot(valsHET2(1:23, i), 'Color', 'r', 'LineWidth', 1.5)
  end 
   plot([0 20], [0 120], 'k--', 'LineWidth', 2)
    xlim([0 20])
  ylim([0 120])
    title('HET')
    
      
 % WT cells with DTX 
  subplot(2,2,3)
  hold on
  for i = 1:12
  plot(DTX_WT(:, i), 'Color', [0.7 0.7 1], 'LineWidth', 1.5)
  end 
  title('WT-DTX')
  plot([0 20], [0 120], 'k--', 'LineWidth', 2)
  xlim([0 20])
  ylim([0 120])
  
  subplot(2,2,4)
  hold on
  for i = 1:13
  plot(DTX_HET(:, i), 'Color', [0.8 0.6 1], 'LineWidth', 1.5)
  end 
   plot([0 20], [0 120], 'k--', 'LineWidth', 2)
    xlim([0 20])
  ylim([0 120])
    title('HET-DTX')
    

%%

































    
    
    

















%% SPIKE TABLE %%%%%%%%%%

% GLU
allWTG = find(spike_table.Geno == 1 & spike_table.GLU ==1); 
allHETG = find(spike_table.Geno == 2 & spike_table.GLU ==1); 

% GABA
allWTB =  find(spike_table.Geno == 1 & spike_table.GLU ==0);
allHETB =  find(spike_table.Geno == 2 & spike_table.GLU ==0);

%
st = 500;
nd = 1500;

vWT = spike_table.V(allWTG, st:nd);
vHET = spike_table.V(allHETG,  st:nd);

dvWT = spike_table.dvdt(allWTG,  st:nd);
dvHET = spike_table.dvdt(allHETG,  st:nd);

vWTB = spike_table.V(allWTB,  st:nd);
vHETB = spike_table.V(allHETB,  st:nd);
% 
dvWTB = spike_table.dvdt(allWTB,  st:nd);
dvHETB = spike_table.dvdt(allHETB,  st:nd);


%% PLOT ACTION POTENTIAL SHAPE

st = 1;
nd = 1001;
% x = (st:1:nd);
x = (1:1:nd-st+1);

nWT = numel(vWT(:,1)); 
nHET = numel(vHET(:,1));

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

nWT = numel(vWTB(:,1)); 
nHET = numel(vHETB(:,1));
% 
mean_WT = nanmean(vWTB); 
mean_HET = nanmean(vHETB); 
% 
% 
semWT = nanstd(vWTB)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
% 
semHET = nanstd(vHETB)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;
% 
% figure
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT', 'k-', 'LineWidth', 1.3)
plot(x, y3, 'w-')
hold on 
plot(x, y4, 'w-')
patch([x fliplr(x)], [y3 fliplr(y4)], 'r' , 'FaceAlpha', 0.2, 'EdgeColor', 'none') %[0.5 0.15 1]
plot(mean_HET', '-','Color', 'r', 'LineWidth', 1.3)


%% PLOT DVDT versus TIME


% WT - HET - overlaid
st = 1;
nd = 1001;
x = (st:1:nd);

nWT = numel(vWT(:,1)); 
nHET =  numel(vHET(:,1));

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

nWT = numel(vWTB(:,1)); 
nHET = numel(vHETB(:,1));

mean_WT = nanmean(dvWTB); 
mean_HET = nanmean(dvHETB); 


semWT = nanstd(dvWTB)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;

semHET = nanstd(dvHETB)/sqrt(nHET); 
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


%% PLOT PHASE PLANE


st = 250;
nd = 750;

x = nanmean(vWT(:, st:nd));
x2 = nanmean(vHET(:, st:nd));

nWT = numel(dvWT(:,1)); 
nHET = numel(dvHET(:,1));

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
x = nanmean(vWTB(:, :));
x2 = nanmean(vHETB(:, :));
% 
nWT = numel(dvWTB(:,1)); 
nHET = numel(dvHETB(:,1));
% 
mean_WT = nanmean(dvWT); 
mean_HET = nanmean(dvHET); 
% 
semWT = nanstd(dvWTB)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
% 
semHET = nanstd(dvHETB)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

% 
% 
% % WT
plot(x, y1, 'w-')
hold on
plot(x, y2, 'w-')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x, mean_WT', 'k-', 'LineWidth', 1.3)
% % HET
plot(x2, y3, 'w-')
hold on 
plot(x2, y4, 'w-')
patch([x2 fliplr(x2)], [y3 fliplr(y4)],  'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x2, mean_HET', '-','Color', 'r', 'LineWidth', 1.3)



%% DUCK PLOTS ACROSS SWEEPS - NEED FULL SPIKE_TABLE table - cannot use just rheobase one.. 

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

allWT = find(spike_table.Geno == 1 & spike_table.Geno1 ==1 & spike_table.Sweep ==sw); % & spike_table.SpikeN ==spn
allHET = find(spike_table.Geno == 2 & spike_table.Geno1 ==1 & spike_table.Sweep ==sw);

allWTDTX =  find(spike_table.Geno == 1 & spike_table.Geno1 ==0 & spike_table.Sweep ==sw);
allHETDTX =  find(spike_table.Geno == 2 & spike_table.Geno1 ==0 & spike_table.Sweep ==sw);

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

%% Analyse variables

allWT = find(spike_table.Geno == 1 & spike_table.Geno1 ==1 );  % sweep 4 = 10pA - sweep 23 = 200pA
allHET = find(spike_table.Geno == 2 & spike_table.Geno1 ==1); 

allWTDTX =  find(spike_table.Geno == 1 & spike_table.Geno1 ==0 );
allHETDTX =  find(spike_table.Geno == 2 & spike_table.Geno1 ==0 );

% k = [11, 15,16,17,18] 

k = 11;
    
    vWT = spike_table{allWT, k};
    vHET = spike_table{allHET, k};
    vWTD = spike_table{allWTDTX, k};
    vHETD = spike_table{allHETDTX, k};
    
    figure; plot([1,2,3,4], [nanmean(vWT), nanmean(vHET), nanmean(vWTD), nanmean(vHETD)], 'k.', 'MarkerSize', 20)
    
    
 %%
 
 %[11,15,16,17,18, 20]
 % 11 = MaxAmp
 % 15 = APThresh
 % 16 = Rise T
 % 17 = whp
 % 18 = Decay T
 % 20 = Min Amp
 % 21 = Max d2 rise
 
[cellid, cellgps] = findgroups(spike_table(:, [1,2,3,4,5]));
cellgps.idx = idx;
cellgps = sortrows(cellgps, 5);
    
values_per_cell = []; 

for k = 1:numel(cellgps(:,1))
    
    ddd = [];
    all_cell = find(cellid == k);
    
    mousee = spike_table.Ani(all_cell(1));
    celll = spike_table.Cell(all_cell(1));
    genoo = spike_table.Geno(all_cell(1));
    typee = spike_table.Geno1(all_cell(1));
    
    ddd(1) = mousee;
    ddd(2) = celll;
    ddd(3) = genoo;
    ddd(4) = typee;
    
    idn = 5;
    
    for p = [11,15,16,17,18, 20, 21]%[10, 14, 15, 16, 17, 19, 25] %23,25] % 22
       
        vals = spike_table{all_cell, p};
        ddd(idn) = nanmean(vals);
        idn = idn+1;
    end 

    values_per_cell = vertcat(values_per_cell, ddd);

end 
  
    
allWTG = find(values_per_cell(:,3)==1 & values_per_cell(:,4)==1);
allHETG = find(values_per_cell(:,3)==2 & values_per_cell(:,4)==1);

allWTB = find(values_per_cell(:,3)==1 & values_per_cell(:,4)==0);
allHETB = find(values_per_cell(:,3)==2 & values_per_cell(:,4)==0);

%%  Boxplot with the 4 different conditions (geno and cell type) 

 % 5 = MaxAmp
 % 6 = APThresh
 % 7 = Rise T
 % 8 = whp
 % 9 = Decay T
 % 10 = Min Amp
 % 11 = Max d2 rise

column = 11; % Column of 'values_per_cell'

wtvalsG = values_per_cell(allWTG, column);
hetvalsG = values_per_cell(allHETG, column);
wtvalsB = values_per_cell(allWTB, column);
hetvalsB = values_per_cell(allHETB, column);

n_wtG = numel(wtvalsG);
n_hetG = numel(hetvalsG);
n_wtB = numel(wtvalsB);
n_hetB = numel(hetvalsB);

x1 = ones(1, n_wtG);
x2 = ones(1, n_hetG)*3;
x3 = ones(1, n_wtB)*2;
x4 = ones(1, n_hetB)*4;

figure
scatter(x1, wtvalsG,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
hold on 
scatter(x2, hetvalsG, 'SizeData', 100,'MarkerEdgeColor', 'r', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(x3, wtvalsB,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
scatter(x4, hetvalsB, 'SizeData', 100,'MarkerEdgeColor', [1 0.6 0.6], 'MarkerFaceColor', [1 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvalsG', wtvalsB', hetvalsG', hetvalsB'], [ones(1,n_wtG), ones(1, n_wtB)*2, ones(1,n_hetG)*3, ones(1, n_hetB)*4], 'Colors', 'k');
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

% ylim([0 30])

f = gcf;
f.Position = [669   473   295   430]; % [1017 351  226  322]; %[1373  50  226 322]; %[ 705   451   216   327]; %[704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

%% Old boxplot

% figure
% scatter(x1, wtvalsG,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
% hold on 
% scatter(x2, hetvalsG, 'SizeData', 100,'MarkerEdgeColor', 'r', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
% scatter(x3, wtvalsB,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 1], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
% scatter(x4, hetvalsB, 'SizeData', 100,'MarkerEdgeColor', 'm', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
% b = boxplot([wtvalsG', hetvalsG', wtvalsB', hetvalsB'], [ones(1,n_wtG), ones(1,n_hetG)*2, ones(1, n_wtB)*3, ones(1, n_hetB)*4], 'Colors', 'k');
% set(b, 'linew', 1.25);
% 
% xticks([1,2])
% xticklabels({''})
% ax = gca;
% box off
% xlim([0.5 4.5])
% ax.XAxis.Visible = 'off'; 
% hold off
% ax.TickDir = 'out'; 
% ax.TickLength = [0.025 0.025];
% ax.LineWidth = 2;

%% STATS - Kruskal Wallis - non-parametric comparison between the 4 groups 
% [p, h] = ranksum(wtvalsG, hetvalsG)

% %%%%%% DATA %%%%%%%%
% wtvalsG = values_per_cell(allWTG, collumn);
% hetvalsG = values_per_cell(allHETG, collumn);
% wtvalsB = values_per_cell(allWTB, collumn);
% hetvalsB = values_per_cell(allHETB, collumn);
% 
% n_wtG = numel(wtvalsG);
% n_hetG = numel(hetvalsG);
% n_wtB = numel(wtvalsB);
% n_hetB = numel(hetvalsB);
% 
% x1 = ones(1, n_wtG);
% x2 = ones(1, n_hetG)*2;
% x3 = ones(1, n_wtB)*3;
% x4 = ones(1, n_hetB)*4;

% % % % % % % % % %  %% 

     y = vertcat(wtvalsG, wtvalsB, hetvalsG, hetvalsB);
     gp3 = vertcat(x1', x3', x2', x4'); % 4 TYPEs
     [p, tbl, stats] = kruskalwallis(y', gp3')  
     results = multcompare(stats, 'Dimension', [1], 'CType', 'dunn-sidak')

[p, h] = ranksum(wtvalsG, hetvalsG)


%% Intrinsic properties plot with mean lines and individual points

figure

nWT = numel(wtvalsG);
nHET = numel(hetvalsG);

mWT = nanmean(wtvalsG);
mHET = nanmean(hetvalsG);

scatter(ones(1,nWT), wtvalsG, 150, 'o','jitter', 'on', 'jitterAmount', 0.25, 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
hold on 
scatter(ones(1,nHET)*2, hetvalsG, 150, 'o', 'jitter', 'on', 'jitterAmount', 0.25,  'MarkerEdgeColor', [1 0.5 0.5], 'LineWidth', 1.5);
% Errorbar
errorbar(1, mWT, nanstd(wtvalsG)/sqrt(nWT), 'k', 'LineWidth', 1.5);
errorbar(2, mHET, nanstd(hetvalsG)/sqrt(nHET), 'r', 'LineWidth', 1.5);
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










%%  Intrinsic properties GLU ONLY 

allrows = find(cell2mat(patch_variables_table.CELL_TYPE) == 1);
p2 = patch_variables_table;
patch_variables_table = patch_variables_table(allrows, :);

%% BOX PLOT

allWT = find(patch_variables_table.Geno == 1);
allHET = find(patch_variables_table.Geno == 2);

% Column val
i = 8; 

% BOXPLOT WITH SCATTER POINTS
close

figure
wtvalsG = patch_variables_table{allWT, 4+i};
hetvalsG = patch_variables_table{allHET, 4+i};

var = patch_variables_table{:, 4+i};
grp = patch_variables_table{:,2};

b = boxplot(var, grp, 'Color', 'k');
set(b , 'LineWidth', 1.25)
hold on
scatter(ones(1,numel(allWT)), wtvalsG, 150, 'o','jitter', 'on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0.5 0.5 0.5]);
scatter(ones(1,numel(allHET))*2, hetvalsG, 150, 'ro', 'jitter', 'on', 'jitterAmount', 0.1);

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


ylim([0 17])
yticks([0:5:25])


[p,h] = ranksum(wtvalsG, hetvalsG)


%% STATS - KRUSKAL WALLIS - WT/HET- SPIKING ACROSS ALL CURRENT INPUTS - 4 GROUPS
   
%     % GLU CELLS 
%     n1 = numel(allGLUW);
%     
%     n2 = numel(allGABAW);
%     
%     % GABA CELLS
%     n3 = numel(allGLUH);
% 
%     n4 = numel(allGABAH);
%     
%     % COMBINE. 
%     y = vertcat(valsWGLU', valsWGAB', valsHGLU', valsHGAB');
%     
%     gp3 = vertcat(ones(1,n1)', (ones(1,n2)')*2, (ones(1,n3)')*3, ((ones(1,n4))')*4); % 4 TYPEs
%     [p, tbl, stats] = kruskalwallis(y', gp3')  
% %     results = multcompare(stats, 'Dimension', [1], 'CType', 'bonferroni')
%     % USE DUNN'S for multcompare
%     
%     g = vertcat(ones(17,1), ones(17,1)*2, ones(17,1)*3, ones(17,1)*4);
%     dunn(data2(:,1)', data2(:,2)')
    
    
    %% STATS - KW - diff method - only comparing between genotype not type. 
   
%     y = []; % GABA
%     y2 = []; % GLU
%     
%     gp1 = [];
%     gp2 = [];
%     
%     for ii = [1,2,3,4]
%         % ii = 1 - GABA HET
%         % ii = 2 - GLU HET
%         % ii = 3 - WT - GABA
%         % ii = 4 - WT - GLU
%         
%         allr = find(data(:, 23)==ii);
%         d = data(allr, 1:20);
%         n_gp = numel(d(:,1));
%         
%         if ii == 1 || ii == 3
%             y = vertcat(y, d);
%             gp1 = vertcat(gp1, ones(n_gp, 1)*ii);
%         elseif ii == 2 || ii ==4
%             y2 = vertcat(y2, d);
%             gp2 = vertcat(gp2, ones(n_gp, 1)*ii);
%         end
%         
%     end
    

%% Rheobase current - current to first spike.

rheo = []; 

for k = 1:numel(cellgps(:,1))
    
    ddd = [];
    all_cell = find(cellid == k);
    
    mousee = spike_table.Ani(all_cell(1));
    celll = spike_table.Cell(all_cell(1));
    genoo = spike_table.Geno(all_cell(1));
    typee = spike_table.Geno1(all_cell(1));
    
    ddd(1) = mousee;
    ddd(2) = celll;
    ddd(3) = genoo;
    ddd(4) = typee;
    
    % Find all the spikes for that cell
    all_spikes_for_cell = spike_table(all_cell, :);
    sweep_no = min(all_spikes_for_cell.Sweep);
    current_2_spike = (sweep_no-3)*10;
    ddd(5) = current_2_spike;

    rheo = vertcat(rheo, ddd);

end 
   
allWTG = find(rheo(:,3)==1 & rheo(:,4)==1);
allHETG = find(rheo(:,3)==2 & rheo(:,4)==1);

allWTB = find(rheo(:,3)==1 & rheo(:,4)==0);
allHETB = find(rheo(:,3)==2 & rheo(:,4)==0);

wtvalsG = rheo(allWTG, 5);
hetvalsG = rheo(allHETG, 5);
wtvalsB = rheo(allWTB, 5);
hetvalsB = rheo(allHETB, 5);

n_wtG = numel(wtvalsG);
n_hetG = numel(hetvalsG);
n_wtB = numel(wtvalsB);
n_hetB = numel(hetvalsB);

x1 = ones(1, n_wtG);
x2 = ones(1, n_hetG)*3;
x3 = ones(1, n_wtB)*2;
x4 = ones(1, n_hetB)*4;

figure
scatter(x1, wtvalsG,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
hold on 
scatter(x2, hetvalsG, 'SizeData', 100,'MarkerEdgeColor', 'r', 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(x3, wtvalsB,'SizeData', 100, 'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
scatter(x4, hetvalsB, 'SizeData', 100,'MarkerEdgeColor', [1 0.6 0.6], 'MarkerFaceColor', [1 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvalsG', wtvalsB', hetvalsG', hetvalsB'], [ones(1,n_wtG), ones(1, n_wtB)*2, ones(1,n_hetG)*3, ones(1, n_hetB)*4], 'Colors', 'k');
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

% ylim([0 100])

f = gcf;
f.Position = [669   473   295   430]; % [1017 351  226  322]; %[1373  50  226 322]; %[ 705   451   216   327]; %[704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

%% 
% [p,h] = ranksum(wtvalsG, hetvalsG)



% %% BOXPLOTS - Spike properties - rheobase spikes - all spikes 
% 
% % ALL SPIKES FROM ALL CELLS FROM RHEOBASE SWEEP
% allWT = find(spike_table.Geno == 1 & spike_table.DTX ==0); 
% allHET = find(spike_table.Geno == 2 & spike_table.DTX ==0); 
% 

%% ADD CELL TYPE TO SPIKE_TABLE WITH ALL SPIKES '220902_SPIKE_TABLE_NEW.mat'


for kk = 1:34
    
    ani = cellgps.Ani(kk);
    celll = cellgps.Cell(kk);
    genooo = cellgps.Geno(kk);
    cohortt = cellgps.Cohort(kk);
    
    Glu = cellgps.Geno1(kk); 
    
    allrows = find(spike_table.Ani == ani & spike_table.Cell == celll & spike_table.Geno == genooo & spike_table.Cohort == cohortt);
    
    spike_table.GLU(allrows) = Glu;

end 

% save('230503_Setd5_Patch_Spike_table_GLUGABA_ALLSWEEPS.mat', 'spike_table', 'cellgps');



%% Looking at average spike variables across all spikes in each sweep. 
% WT/ HET - GLU 

% X = current injection 
% y value = all the spike variable values for that genotype

column = 19;  % Max Amp

d1 = nan(1,21);
d2 = nan(1,21);
d3 = nan(1,21);
d4 = nan(1,21);

sem1 = nan(1,21);
sem2 = nan(1,21);
sem3 = nan(1,21);
sem4 = nan(1,21);


for ll = 4:24
    
    allGLUW = find(spike_table.Geno == 1 & spike_table.GLU == 1 & spike_table.Sweep == ll & spike_table.DTX == 0);
    allGLUH = find(spike_table.Geno == 2 & spike_table.GLU == 1 & spike_table.Sweep == ll & spike_table.DTX == 0);
    
    allGABAW = find(spike_table.Geno == 1 & spike_table.GLU == 0 & spike_table.Sweep == ll & spike_table.DTX == 0);
    allGABAH = find(spike_table.Geno == 2 & spike_table.GLU == 0 & spike_table.Sweep == ll & spike_table.DTX == 0);
    
    v1 = spike_table{allGLUW, column}';
    v2 = spike_table{allGABAW, column}';
    v3 = spike_table{allGLUH, column}';
    v4 = spike_table{allGABAH, column}';
    
    n1 = numel(v1);
    n2 = numel(v2);
    n3 = numel(v3);
    n4 = numel(v4);
    
    d1(ll-3) = nanmean(v1);
    d2(ll-3) = nanmean(v2);
    d3(ll-3) = nanmean(v3);
    d4(ll-3) = nanmean(v4);
    
    sem1(ll-3) = nanstd(v1')/sqrt(n1); 
    sem2(ll-3) = nanstd(v2')/sqrt(n2); 
    sem3(ll-3) = nanstd(v3')/sqrt(n3); 
    sem4(ll-3) = nanstd(v4')/sqrt(n4); 
    
end 

% max(spike_table.Sweep)

v = 21;
inj_steps = curr_inj(3:3+v, 2);

x = (1:1:v);

y1 = d1+sem1;
y2 = d1-sem1;
     
y3 = d2+sem2;
y4 = d2-sem2;

figure
plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(d1(1:v)', 'k', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(d2(1:v)', 'Color', col, 'LineWidth', 1.3)

% GABA

y1 = d3+sem3;
y2 = d3-sem3;
     
y3 = d4+sem4;
y4 = d4-sem4;

plot(x, y1(1:v), 'w')
hold on
plot(x, y2(1:v), 'w')
patch([x fliplr(x)], [y1(1:v) fliplr(y2(1:v))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(d3(1:v)', 'k--', 'LineWidth', 1.3)

plot(x, y3(1:v), 'w')
hold on 
plot(x, y4(1:v), 'w')
patch([x fliplr(x)], [y3(1:v) fliplr(y4(1:v))],  'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(d4(1:v)','r--', 'LineWidth', 1.3)

xticks([1:4:24])
xticklabels(string(inj_steps(1:4:24)))

xlabel('Current (pA)')
box off
set(gca, 'FontSize', 20)
ax = gca;
ax.TickDir  = 'out';
ax.LineWidth = 1.75;
ax.TickLength = [0.03 0.03];
xlim([-0.5 24])

% 
% f = gcf;
% f.Position = [742   492   304   258];


ylabel('Min Amp - AHP (mV)')
% ylim([0 100])


%% BOXPLOTS - all spikes in all sweeps. 

% ALL SPIKES FROM ALL CELLS FROM ALL SWEEPS
allGLUW = find(spike_table.Geno == 1 & spike_table.GLU == 1 & spike_table.DTX == 0);
allGLUH = find(spike_table.Geno == 2 & spike_table.GLU == 1  & spike_table.DTX == 0);

allGABAW = find(spike_table.Geno == 1 & spike_table.GLU == 0  & spike_table.DTX == 0);
allGABAH = find(spike_table.Geno == 2 & spike_table.GLU == 0 & spike_table.DTX == 0);

column = 19; 

wtvalsG = spike_table{allGLUW, column};
hetvalsG = spike_table{allGLUH, column};

wtvalsB = spike_table{allGABAW, column};
hetvalsB = spike_table{allGABAH, column};


%% PLOT - WT versus HET - GLU
n_wt = numel(wtvalsG);
n_het = numel(hetvalsG);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

figure
% scatter(x1, wtvals,'SizeData', 50, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
% scatter(x2, hetvals, 'SizeData', 50,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvalsG', hetvalsG'], [ones(1,n_wt), ones(1,n_het)*2], 'Colors', 'k');
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

f = gcf;
f.Position = [ 705   451   216   327]; %[704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

nanmean(wtvalsG)
nanmean(hetvalsG)
[p,h] = ranksum(wtvalsG, hetvalsG)
[h, p] = kstest2(wtvalsG, hetvalsG)



%% PLOT - GABA - WT versus HET 
n_wt = numel(wtvalsB);
n_het = numel(hetvalsB);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

figure
% scatter(x1, wtvals,'SizeData', 50, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
% scatter(x2, hetvals, 'SizeData', 50,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot([wtvalsB', hetvalsB'], [ones(1,n_wt), ones(1,n_het)*2], 'Colors', 'k');
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

f = gcf;
f.Position = [ 705   451   216   327]; %[704   507   187   272]; %[704   207   355   572]; 
f.Renderer = 'painters';

nanmean(wtvalsB)
nanmean(hetvalsB)
[p,h] = ranksum(wtvalsB, hetvalsB)
[h, p] = kstest2(wtvalsB, hetvalsB)


ylim([-80 -50])








