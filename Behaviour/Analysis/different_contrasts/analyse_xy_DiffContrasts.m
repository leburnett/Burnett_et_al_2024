%% xy_analyis Different Contrasts 

WT_C1 = []; 
WT_C2 = []; 
WT_C3 = []; 
WT_C4 = []; 

HET_C1 = []; 
HET_C2 = []; 
HET_C3 = []; 
HET_C4 = []; 


for i = 1:height(xy_analysis)
        if string(xy_analysis.Geno{i}) == "wt" && xy_analysis.Contrast{i} == 1
            G = cell2mat(xy_analysis.TimeToMaxSp(i));
            WT_C1 = vertcat(WT_C1, G);
        
        elseif string(xy_analysis.Geno{i}) == "wt" && xy_analysis.Contrast{i} == 2
            G = cell2mat(xy_analysis.TimeToMaxSp(i));
             WT_C2 = vertcat(WT_C2, G);
        
        elseif string(xy_analysis.Geno{i}) == "wt" && xy_analysis.Contrast{i} == 3
            G = cell2mat(xy_analysis.TimeToMaxSp(i));
            WT_C3 = vertcat(WT_C3, G);
    
        elseif string(xy_analysis.Geno{i}) == "wt" && xy_analysis.Contrast{i} == 0
            G = cell2mat(xy_analysis.TimeToMaxSp(i));
            WT_C4 = vertcat(WT_C4, G);
         
        elseif string(xy_analysis.Geno{i}) == "het" && xy_analysis.Contrast{i} == 1
            F = cell2mat(xy_analysis.TimeToMaxSp(i));
            HET_C1 = vertcat(HET_C1, F);
            
        elseif string(xy_analysis.Geno{i}) == "het" && xy_analysis.Contrast{i} == 2
            F = cell2mat(xy_analysis.TimeToMaxSp(i));
            HET_C2 = vertcat(HET_C2, F);
            
        elseif string(xy_analysis.Geno{i}) == "het" && xy_analysis.Contrast{i} == 3
            F = cell2mat(xy_analysis.TimeToMaxSp(i));
            HET_C3 = vertcat(HET_C3, F);
            
        elseif string(xy_analysis.Geno{i}) == "het" && xy_analysis.Contrast{i} == 0
            F = cell2mat(xy_analysis.TimeToMaxSp(i));
            HET_C4 = vertcat(HET_C4, F);
       
        end
end
  
%%  Multigroup boxplot with scatter plot points. 

%Remove contrast 2. 
% rows2 = find(cell2mat(xy_analysis.Contrast)==2);
% xy_analysis(rows2,:) = []; 

% Change the 'contrast' number for 3/4 .

% for j = 1:93
%     if cell2mat(xy_analysis.Contrast(j)) == 3
%         xy_analysis.Contrast{j} =2; 
%     elseif cell2mat(xy_analysis.Contrast(j)) == 4
%         xy_analysis.Contrast{j} =3; 
%     end 
% end 

allWT = find(string(xy_analysis.Geno) == "wt");
allHET = find(string(xy_analysis.Geno) == "het");

xy_analysisB = xy_analysis(allWT, :); 

xdata = cell2mat(xy_analysisB.TimeToMaxSp);
ydata = cell2mat(xy_analysisB.MaxSpEscape);
group = cell2mat(xy_analysisB.Contrast);

ymin = -0.5; 
ymax= 100;

scatterhist(xdata,ydata,'Group',group,'Kernel','on', 'Color', [0,0,0;0.5,.5,.5; .75, .75, .75], 'Direction', 'in', 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 25, 'Legend', 'off')
xlabel('Time to Max Sp - s')
ylabel('Max Speed - cm/s')
title('WT')
set(gca, 'FontSize', 14)
hold on 
plot([0 0], [ymin ymax], 'k:', 'LineWidth', 1.5)
plot([0.75 0.75], [ymin ymax], 'k:', 'LineWidth', 1)
plot([1.5 1.5], [ymin ymax], 'k:', 'LineWidth', 1)
plot([2.25 2.25], [ymin ymax], 'k:', 'LineWidth', 1)
plot([3 3], [ymin ymax], 'k:', 'LineWidth', 1)
plot([3.75 3.75], [ymin ymax], 'k:', 'LineWidth', 1)

% HET

xy_analysisC = xy_analysis(allHET, :); 

xdata2 = cell2mat(xy_analysisC.TimeToMaxSp);
ydata2 = cell2mat(xy_analysisC.MaxSpEscape);
group2 = cell2mat(xy_analysisC.Contrast);

ymin = -5; 
ymax= 70;

scatterhist(xdata2,ydata2,'Group',group2,'Kernel','on', 'Color', [1,0,0;1,.5,.5; 1, .8, .8], 'Direction', 'in', 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 25, 'Legend', 'off')
xlabel('Time to Max Sp - s')
ylabel('Max Speed - cm/s')
title('HET')
set(gca, 'FontSize', 14)
hold on 
plot([0 0], [ymin ymax], 'k:', 'LineWidth', 1.5)
plot([0.75 0.75], [ymin ymax], 'k:', 'LineWidth', 1)
plot([1.5 1.5], [ymin ymax], 'k:', 'LineWidth', 1)
plot([2.25 2.25], [ymin ymax], 'k:', 'LineWidth', 1)
plot([3 3], [ymin ymax], 'k:', 'LineWidth', 1)
plot([3.75 3.75], [ymin ymax], 'k:', 'LineWidth', 1)


% var = animal_xy_table.Freeze;
% group = string(animal_xy_table.Grp);
% figure
% b = boxplot(var, group, 'Colors', 'krbmk'); 
%   

%%

subplot(1,2,2)
histogram((HET_C1), 'BinWidth', 1, 'Normalization', 'pdf','DisplayStyle', 'stairs', 'LineWidth', 2.5, 'EdgeColor', [1 0 0])
hold on 
histogram((HET_C2), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2.5, 'EdgeColor', [1 0 0], 'EdgeAlpha', 0.4)
histogram((HET_C3), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2.5, 'EdgeColor', [1 0 0],  'EdgeAlpha', 0.2)
title('Time to Max Speed - HET')

subplot(1,2,1)
histogram((WT_C1), 'BinWidth', 1, 'Normalization', 'pdf','DisplayStyle', 'stairs', 'LineWidth', 2.5, 'EdgeColor', [0 0 0])
hold on 
histogram((WT_C2), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2.5, 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.4)
histogram((WT_C3), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2.5, 'EdgeColor', [0 0 0],  'EdgeAlpha', 0.2)
title('Time to Max Speed - WT')



%%
% PLOT 
n_wt_set = numel(WT_Setd5);
n_het_set = numel(HET_Setd5);
n_wt_cul = numel(WT_Cul3);
n_het_cul = numel(HET_Cul3);
n_wt_c57 = numel(WT_C57);

x1 = ones(1, n_wt_set)*1; 
x2 = ones(1, n_het_set)*2;
x3 = ones(1, n_wt_c57)*3; 
x4 = ones(1, n_het_cul)*4;
x5 = ones(1, n_wt_cul)*5; 

f = figure;
 f.Position = [25 300 350 450]; %l b w h
geno = string(animal_xy_table.Grp); 
var = ((animal_xy_table.TimeFrozen));
boxplot(var, geno, 'Colors', 'k')
hold on 
scatter(x1, WT_Setd5', 'filled','SizeData', 50, 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.8,'jitter', 'on', 'jitterAmount', 0.05)
scatter(x2, HET_Setd5', 'filled', 'SizeData', 50,'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.8,'jitter', 'on', 'jitterAmount', 0.05)
scatter(x3, WT_C57', 'filled','SizeData', 50, 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.8,'jitter', 'on', 'jitterAmount', 0.05)
scatter(x4, HET_Cul3', 'filled', 'SizeData', 50,'MarkerFaceColor', 'm', 'MarkerFaceAlpha', 0.8,'jitter', 'on', 'jitterAmount', 0.05)
scatter(x5, WT_Cul3', 'filled','SizeData', 50, 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.8,'jitter', 'on', 'jitterAmount', 0.05)
ylabel('Time - s')
title('Time Frozen')
set(gca, 'FontSize', 14)
xtickangle(45)
axis([0.5 5.5 0 1.3])




%% XYLOOM

n = height(xy_analysis);
details = ALL_XYLOOM_TABLE(:, [1,2,3,4,6,7,8,9]);

%Change contrast to more straightofrward values
for i = 1:n
    if cell2mat(xy_analysis.Contrast(i))==5
        xy_analysis.Contrast(i) = {2};
    elseif cell2mat(xy_analysis.Contrast(i))==4
        xy_analysis.Contrast(i) = {3};
    end 
end 

% Add contrast column to xy_loom
ALL_XYLOOM_TABLE.Contrast = xy_analysis.Contrast;


%% Mean +SEM - WT/HET


speed_WT = [];
speed_HET = []; 
for i = 1:n
    if string(ALL_XYLOOM_TABLE.Geno{i}) == "wt"  && ALL_XYLOOM_TABLE.Contrast{i} == 3% &&  ALL_XYLOOM_TABLE.Trial(i) < 10
        G = cell2mat(ALL_XYLOOM_TABLE{i,5});
        speed_WT = vertcat(speed_WT, G); 
    elseif string(ALL_XYLOOM_TABLE.Geno{i}) == "het"  && ALL_XYLOOM_TABLE.Contrast{i} == 3 %&&  ALL_XYLOOM_TABLE.Trial(i) < 10
        F = cell2mat(ALL_XYLOOM_TABLE{i,5});
        speed_HET = vertcat(speed_HET, F);
    end
end

% Plot ALL
plot(speed_HET', 'r')
hold on 
plot(speed_WT', 'k')

%% Plot mean + SEM

rowsb4 = 180;
lf = 46; 

nWT = numel(speed_WT(:,1)); 
nHET = numel(speed_HET(:,1));

mean_WT = mean(speed_WT); 
mean_HET = mean(speed_HET); 

x = (1:1:780);

semWT = std(speed_WT)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = std(speed_HET)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_WT', 'k', 'LineWidth', 1.3)

plot(x, y3, 'w')
hold on 
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(mean_HET', 'Color', col, 'LineWidth', 1.3)

plot([rowsb4 rowsb4], [0 80], 'k:', 'LineWidth', 1.5)
plot([rowsb4+lf rowsb4+lf], [0 80], 'k:', 'LineWidth', 1)
plot([rowsb4+lf*2 rowsb4+lf*2], [0 80], 'k:', 'LineWidth', 1)
plot([rowsb4+lf*3 rowsb4+lf*3], [0 80], 'k:', 'LineWidth', 1)
plot([rowsb4+lf*4 rowsb4+lf*4], [0 80], 'k:', 'LineWidth', 1)
plot([rowsb4+lf*5 rowsb4+lf*5], [0 80], 'k:', 'LineWidth', 1)
% title(strcat('Response Speed to Loom - Day ', string(k)'))
% title('Response Speed to Loom')
 xticks([60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840, 900, 960, 1020, 1080])
 xticklabels({'-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'})
 xlabel('Time - s')
 ylabel('Speed - cm/s')
 set(gca, 'FontSize', 24)
axis([0 780 0 70])
box off
% title('')
% xticks([])
% xticklabels({''})
% % yticks([])
% yticklabels({''})   


%% Add trial #

  all_animals = unique(ALL_XYLOOM_TABLE.Animal); 
  n_animals = numel(all_animals); 

for i = 1:n_animals
    ani = all_animals(i); 
    trial_number = 1; 
    
    for j = 1:n
        if string(ALL_XYLOOM_TABLE.Animal{j}) == ani 
            ALL_XYLOOM_TABLE.Trial(j) = trial_number; 
            trial_number = trial_number +1; 
        end   
    end
end 


%% Plots per animal per contrast
    ALL_XYLOOM_TABLE_2 = sortrows(ALL_XYLOOM_TABLE, 7);

  all_animals = unique(xy_analysis.Animal); 
  n_animals = numel(all_animals); 

for i = 1:n_animals
    ani = string(all_animals{i}); 

    
    speed_1 = [];
    speed_2 = [];
    speed_3 = []; 
    
    for i = 1:n
        if string(ALL_XYLOOM_TABLE_2.Animal{i}) == ani  && ALL_XYLOOM_TABLE_2.Contrast{i} == 1
            G = cell2mat(ALL_XYLOOM_TABLE_2{i,5});
            speed_1 = vertcat(speed_1, G);
        elseif string(ALL_XYLOOM_TABLE_2.Animal{i}) == ani  && ALL_XYLOOM_TABLE_2.Contrast{i} == 2
            F = cell2mat(ALL_XYLOOM_TABLE_2{i,5});
            speed_2 = vertcat(speed_2, F);
        elseif string(ALL_XYLOOM_TABLE_2.Animal{i}) == ani  && ALL_XYLOOM_TABLE_2.Contrast{i} == 3
            F = cell2mat(ALL_XYLOOM_TABLE_2{i,5});
            speed_3 = vertcat(speed_3, F);
        end
    end
    
    
    figure
    subplot(3,1,1)
    imagesc(speed_1)
    caxis([0 60])
    hold on
    plot([180 180], [0 150], 'w', 'LineWidth', 1.5)
    plot([225 225], [0 150], 'w:', 'LineWidth', 1.7)
    plot([270 270], [0 150], 'w:', 'LineWidth', 1.7)
    plot([315 315], [0 150], 'w:', 'LineWidth', 1.7)
    plot([360 360], [0 150], 'w:', 'LineWidth', 1.7)
    plot([405 405], [0 150], 'w', 'LineWidth', 1.2)
    xticks([60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840, 900, 960, 1020, 1080])
    xticklabels({'-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'})
    yticks([])
    ax = gca;
    ax.FontSize  =14;
    colormap(redblue)
    
     subplot(3,1,2)
    imagesc(speed_2)
    caxis([0 60])
    hold on
    plot([180 180], [0 150], 'w', 'LineWidth', 1.5)
    plot([225 225], [0 150], 'w:', 'LineWidth', 1.7)
    plot([270 270], [0 150], 'w:', 'LineWidth', 1.7)
    plot([315 315], [0 150], 'w:', 'LineWidth', 1.7)
    plot([360 360], [0 150], 'w:', 'LineWidth', 1.7)
    plot([405 405], [0 150], 'w', 'LineWidth', 1.2)
    xticks([60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840, 900, 960, 1020, 1080])
    xticklabels({'-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'})
    ax = gca;
    ax.FontSize  =14;
    colormap(redblue)
     yticks([])
    
     subplot(3,1,3)
    imagesc(speed_3)
    caxis([0 60])
    hold on
    plot([180 180], [0 150], 'w', 'LineWidth', 1.5)
    plot([225 225], [0 150], 'w:', 'LineWidth', 1.7)
    plot([270 270], [0 150], 'w:', 'LineWidth', 1.7)
    plot([315 315], [0 150], 'w:', 'LineWidth', 1.7)
    plot([360 360], [0 150], 'w:', 'LineWidth', 1.7)
    plot([405 405], [0 150], 'w', 'LineWidth', 1.2)
    xticks([60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840, 900, 960, 1020, 1080])
    xticklabels({'-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'})
    ax = gca;
    ax.FontSize  =14;
    colormap(redblue)
     yticks([])
    
    sgtitle(ani)
    

end 



