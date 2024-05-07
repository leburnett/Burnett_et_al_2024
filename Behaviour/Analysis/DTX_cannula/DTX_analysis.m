
% 30/03/23 - Analysis of DTX loom experiments 
% Setd5 animals with cannula
% Created by Burnett

%% 1 - REACTION TIME - SALINE VS DTX - for WT and HET animals 

all_animals = unique(xy_analysis.Animal);
n_animals = numel(all_animals);

valsWT1=[];
valsHET1 = [];
valsWT2=[];
valsHET2 = [];

figure
for i = 1:n_animals
    
    all_ani = find(string(xy_analysis.Animal) == all_animals{i});
    
    genoo = xy_analysis.Geno{all_ani(1)};
    if genoo == "wt"
        col = 'k';
    else
        col = 'r';
    end
    
    ani_before = find(string(xy_analysis.Animal) == all_animals{i} & cell2mat(xy_analysis.data) == 1);
    ani_after = find(string(xy_analysis.Animal) == all_animals{i} & cell2mat(xy_analysis.data) == 2);
    
    if ~isempty(ani_before) && ~isempty(ani_after)
        val_before = nanmean(cell2mat(xy_analysis.TimeToMaxSp(ani_before)));
        val_after = nanmean(cell2mat(xy_analysis.TimeToMaxSp(ani_after)));
        
        if genoo == "wt"
            valsWT1 = [valsWT1, val_before];
            valsWT2 = [valsWT2, val_after];
        else
            valsHET1 = [valsHET1, val_before];
            valsHET2 = [valsHET2, val_after];
        end 
        
        plot( [1+rand(1)/2,3+rand(1)/2], [val_before, val_after], '-', 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3) %'MarkerFaceColor', [1 0.8 0.8]
        hold on
        
%     elseif ~isempty(ani_before) && isempty(ani_after)
%         val_before = nanmean(cell2mat(xy_analysis.TimeToMaxSp(ani_before)));
%         plot(1+rand(1)/2, val_before, 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3)
%         hold on
%         if genoo == "wt"
%             valsWT1 = [valsWT1, val_before];
%         else
%             valsHET1 = [valsHET1, val_before];
%         end
        
%     elseif isempty(ani_before) && ~isempty(ani_after)
%                 val_before = nanmean(cell2mat(xy_analysis.TimeToMaxSp(ani_after)));
%                 plot(3+rand(1)/2, val_before, 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3)
%                 hold on
%         
%         if genoo == "wt"
%             valsWT2 = [valsWT2, val_after];
%         else
%             valsHET2 = [valsHET2, val_after];
%         end
    end

end 


%% Ammend genotype - xy_analysis 

het_animals = {'M21505', 'M21507', 'M30147', 'M30152', 'M22217', 'M22333', 'M30035', 'M30129', 'G30139', 'G30150'};

for i = 1:height(xy_analysis)
    ani = xy_analysis.Animal{i};
    if ismember(ani, het_animals)
        xy_analysis.Geno{i} = "het"; 
    end 
end 

%% BAR CHART - Reaction time (WT-B4, HET-B4, WT-DTX, HET-DTX) with points of individual animals on top. 

x = [1,2,3,4];
y = [nanmean(valsWT1), nanmean(valsWT2), nanmean(valsHET1), nanmean(valsHET2)];
y2 = [median(valsWT1), median(valsHET1), median(valsWT2), median(valsHET2)];

col = [0.25 0.25 0.25];

figure
b = bar(x,y);
hold on
scatter(ones(numel(valsWT1),1), valsWT1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05, 'LineWidth', 1.5)
scatter(ones(numel(valsHET1),1)*3, valsHET1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05, 'LineWidth', 1.5)
scatter(ones(numel(valsWT2),1)*2, valsWT2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05, 'LineWidth', 1.5)
scatter(ones(numel(valsHET2),1)*4, valsHET2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.05, 'LineWidth', 1.5)
b.FaceColor = [0.7 0.7 0.7];

for i = 1:numel(valsWT1)
    plot([1,2], [valsWT1(i), valsWT2(i)], '-', 'Marker', 'none', 'Color', col, 'LineWidth', 1.3)
end 

for i = 1:numel(valsHET1)
    plot([3,4], [valsHET1(i), valsHET2(i)], '-', 'Marker', 'none', 'Color', col, 'LineWidth', 1.3)
end 

box off
ylim([0 3])
xlim([0 5])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.FontSize = 14;

f = gcf;
f.Position = [440   391   306   407];


%% STATISTICS - REACTION TIME 

% 1 - Are the data normally distributed?
% Small samples therefore use shapiro wilk test 

[h, p] = swtest(valsWT1)
% 0.6581 % ACCEPT NULL

[h, p] = swtest(valsWT2)
% 0.5569 % ACCEPT NULL

[h, p] = swtest(valsHET1)
% 0.4875 % ACCEPT NULL

[h, p] = swtest(valsHET2)
% 0.0827 % ACCEPT NULL


% PAIRED! Non-parametric 
[p,h] = ranksum(valsWT1, valsWT2)
[p,h] = ranksum(valsHET1, valsHET2)
% PAIRED - parametric
[h, p] = ttest(valsWT1, valsWT2)
[h,p] = ttest(valsHET1, valsHET2)


% NOT-PAIRED - 
% non parametric version of not-paired. 
[h,p] = ranksum(valsWT1, valsHET1)
[h,p] = ranksum(valsWT1, valsHET2)
[h,p] = ranksum(valsHET1, valsWT2)

%Within genotypes
[h,p] = ranksum(valsWT1, valsWT2)
[h,p] = ranksum(valsHET1, valsHET2)


% Number of WT/ HET animals
allWT = find(string(xy_analysis.Geno) == "wt");
allHET = find(string(xy_analysis.Geno) == "het");

WTanimals = unique(xy_analysis.Animal(allWT)); % 9 animals 
HETanimals = unique(xy_analysis.Animal(allHET)); %  9 animals

% Number of animals with both BEFORE and AFTER. 
% 11 WT values - 5 PAIRED sample
% 9 HET values - 6 PAIRED SAMPLES. 

save('ReactionTime_DATA_DTX_Setd5.mat', 'valsWT1', 'valsWT2', 'valsHET1', 'valsHET2');

%% Repeated measures ANOVA analysis 

rxn1 = [(valsWT1),(valsHET1)];
rxn2 = [(valsWT2), (valsHET2)];
genooo = {'wt', 'wt', 'wt', 'wt', 'wt', 'het', 'het', 'het', 'het', 'het', 'het'};

t1 = array2table(rxn1', 'VariableNames', {'Rxn1'});
t2 = array2table(rxn2', 'VariableNames', {'Rxn2'});
gp = array2table(genooo', 'VariableNames', {'Geno'});

t3 = [t1, t2, gp];

rm = fitrm(t3, 'Rxn1-Rxn2~Geno') % 'WithinDesign', within, 'WithinModel', 'separatemeans' 
ranovatbl = ranova(rm)
tbl = multcompare(rm, 'Geno')



% OTHER OPTION

rxn = [(valsWT1),(valsHET1), (valsWT2), (valsHET2)];
genooo = {'wt', 'wt', 'wt', 'wt', 'wt', 'het', 'het', 'het', 'het', 'het', 'het', 'wt', 'wt', 'wt', 'wt', 'wt', 'het', 'het', 'het', 'het', 'het', 'het'};
DTX = {'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX'};

tA = array2table(rxn', 'VariableNames', {'Rxn'});
gp = array2table(genooo', 'VariableNames', {'Geno'});
gp2 = array2table(DTX', 'VariableNames', {'DTX'});

tB = [tA, gp, gp2];
rm = fitrm(tB, 'Rxn~DTX*Geno') 
ranovatbl = ranova(rm)
tbl = multcompare(rm, 'Geno', 'By', 'DTX')
tbl = multcompare(rm, 'DTX', 'By', 'Geno')

%% 2 - MAXIMUM SPEED -  SALINE VS DTX - for WT and HET animals 

all_animals = unique(xy_analysis.Animal);
n_animals = numel(all_animals);

valsWT1=[];
valsHET1 = [];
valsWT2=[];
valsHET2 = [];

figure
for i = 1:n_animals
    
    all_ani = find(string(xy_analysis.Animal) == all_animals{i});
    
    genoo = xy_analysis.Geno{all_ani(1)};
    if genoo == "wt"
        col = 'k';
    else
        col = 'r';
    end
    
    ani_before = find(string(xy_analysis.Animal) == all_animals{i} & cell2mat(xy_analysis.data) == 1);
    ani_after = find(string(xy_analysis.Animal) == all_animals{i} & cell2mat(xy_analysis.data) == 2);
    
    if ~isempty(ani_before) && ~isempty(ani_after)
        val_before = nanmean(cell2mat(xy_analysis.MaxSpEscape(ani_before)));
        val_after = nanmean(cell2mat(xy_analysis.MaxSpEscape(ani_after)));
        
        if genoo == "wt"
            valsWT1 = [valsWT1, val_before];
            valsWT2 = [valsWT2, val_after];
        else
            valsHET1 = [valsHET1, val_before];
            valsHET2 = [valsHET2, val_after];
        end 
        
        plot( [1+rand(1)/2,3+rand(1)/2], [val_before, val_after], '-', 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3) %'MarkerFaceColor', [1 0.8 0.8]
        hold on
        
%     elseif ~isempty(ani_before) && isempty(ani_after)
%         val_before = nanmean(cell2mat(xy_analysis.TimeToMaxSp(ani_before)));
%         plot(1+rand(1)/2, val_before, 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3)
%         hold on
%         if genoo == "wt"
%             valsWT1 = [valsWT1, val_before];
%         else
%             valsHET1 = [valsHET1, val_before];
%         end
        
%     elseif isempty(ani_before) && ~isempty(ani_after)
%                 val_before = nanmean(cell2mat(xy_analysis.TimeToMaxSp(ani_after)));
%                 plot(3+rand(1)/2, val_before, 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3)
%                 hold on
%         
%         if genoo == "wt"
%             valsWT2 = [valsWT2, val_after];
%         else
%             valsHET2 = [valsHET2, val_after];
%         end
    end

end 


%% BAR CHART - Maximum escape speed (WT-B4, HET-B4, WT-DTX, HET-DTX) with points of individual animals on top. 

x = [1,2,3,4];
y = [nanmean(valsWT1), nanmean(valsWT2), nanmean(valsHET1), nanmean(valsHET2)];
y2 = [median(valsWT1), median(valsHET1), median(valsWT2), median(valsHET2)];

col = [0.25 0.25 0.25];

figure
b = bar(x,y);
hold on
scatter(ones(numel(valsWT1),1), valsWT1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsHET1),1)*3, valsHET1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsWT2),1)*2, valsWT2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsHET2),1)*4, valsHET2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
b.FaceColor = [0.7 0.7 0.7];

for i = 1:numel(valsWT1)
    plot([1,2], [valsWT1(i), valsWT2(i)], '-', 'Marker', 'none', 'Color', col, 'LineWidth', 1.3)
end 

for i = 1:numel(valsHET1)
    plot([3,4], [valsHET1(i), valsHET2(i)], '-', 'Marker', 'none', 'Color', col, 'LineWidth', 1.3)
end 


box off
ylim([0 110])
xlim([0 5])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.FontSize = 14;

f = gcf;
f.Position = [440   391   306   407];


%% STATISTICS - Max Speed

% 1 - Are the data normally distributed?
% Small samples therefore use shapiro wilk test 

[h, p] = swtest(valsWT1)
% 0.1795 % ACCEPT NULL

[h, p] = swtest(valsWT2)
% 0.8481 % ACCEPT NULL

[h, p] = swtest(valsHET1)
% 0.7314 % ACCEPT NULL

[h, p] = swtest(valsHET2)
% 0.4584 % ACCEPT NULL


% PAIRED! Non-parametric 
[p,h] = ranksum(valsWT1, valsWT2)
[p,h] = ranksum(valsHET1, valsHET2)
% PAIRED - parametric
[h, p] = ttest(valsWT1, valsWT2)
[h,p] = ttest(valsHET1, valsHET2)


% NOT-PAIRED - 
% non parametric version of not-paired. 
[h,p] = ranksum(valsWT1, valsHET2)
[h,p] = ranksum(valsHET1, valsWT2)


% Number of animals with both BEFORE and AFTER. 
% 11 WT values - 5 PAIRED sample
% 9 HET values - 6 PAIRED SAMPLES. 

save('MaxSp_DATA_DTX_Setd5.mat', 'valsWT1', 'valsWT2', 'valsHET1', 'valsHET2');


%% ANOVA 

rxn1 = [(valsWT1),(valsHET1)];
rxn2 = [(valsWT2), (valsHET2)];
genooo = {'wt', 'wt', 'wt', 'wt', 'wt', 'het', 'het', 'het', 'het', 'het', 'het'};

t1 = array2table(rxn1', 'VariableNames', {'Rxn1'});
t2 = array2table(rxn2', 'VariableNames', {'Rxn2'});
gp = array2table(genooo', 'VariableNames', {'Geno'});

t3 = [t1, t2, gp];

rm = fitrm(t3, 'Rxn1-Rxn2~Geno') % 'WithinDesign', within, 'WithinModel', 'separatemeans' 
ranovatbl = ranova(rm)
tbl = multcompare(rm, 'Geno')

%

maxsp = [(valsWT1),(valsHET1), (valsWT2), (valsHET2)];
genooo = {'wt', 'wt', 'wt', 'wt', 'wt', 'het', 'het', 'het', 'het', 'het', 'het', 'wt', 'wt', 'wt', 'wt', 'wt', 'het', 'het', 'het', 'het', 'het', 'het'};
DTX = {'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX'};

tA = array2table(maxsp', 'VariableNames', {'MaxSp'});
gp = array2table(genooo', 'VariableNames', {'Geno'});
gp2 = array2table(DTX', 'VariableNames', {'DTX'});

tB = [tA, gp, gp2];
rm = fitrm(tB, 'MaxSp~DTX*Geno') 
ranovatbl = ranova(rm)
tbl = multcompare(rm, 'Geno', 'By', 'DTX')
tbl = multcompare(rm, 'DTX', 'By', 'Geno')
















%% 3 - HEATMAPS OF SPEED TRACES - using ALL_XYLOOM_TABLE
% Make ALL_XYLOOM tables for 'before / saline' responses and another for DTX

details = ALL_XYLOOM_B4(:, [1,2,3,4,6,7]);
details2 = ALL_XYLOOM_DTX(:, [1,2,3,4,6,7]);

ALL_XYLOOM_B4 = [];
ALL_XYLOOM_DTX = []; 

ALL_XYLOOM_B4 = vertcat(ALL_XYLOOM_B4, ALL_XYLOOM_TABLE(1:9, :));
ALL_XYLOOM_DTX = vertcat(ALL_XYLOOM_DTX, ALL_XYLOOM_TABLE(10:17, :));

save('230428_Setd5_ALL_XYLOOM_B4_DTX.mat', 'ALL_XYLOOM_B4', 'ALL_XYLOOM_DTX');


%% Plot HEATMAPS for WT/ HET - saline and DTX

% SORT BY T2M

ALL_XYLOOM_B4 = sortrows(ALL_XYLOOM_B4, 7);
ALL_XYLOOM_DTX = sortrows(ALL_XYLOOM_DTX, 7);

speed_WT_B4 = [];
speed_HET_B4 = []; 

speed_WT_DTX = [];
speed_HET_DTX = []; 

for i = 1:height(ALL_XYLOOM_B4)
    
    if string(ALL_XYLOOM_B4.Geno{i}) == "wt"
        G = cell2mat(ALL_XYLOOM_B4{i,5});
        speed_WT_B4 = vertcat(speed_WT_B4, G);
    elseif string(ALL_XYLOOM_B4.Geno{i}) == "het"
        F = cell2mat(ALL_XYLOOM_B4{i,5});
        speed_HET_B4 = vertcat(speed_HET_B4, F);
    end
end


for i = 1:height(ALL_XYLOOM_DTX)
    
    if string(ALL_XYLOOM_DTX.Geno{i}) == "wt"
        G = cell2mat(ALL_XYLOOM_DTX{i,5});
        speed_WT_DTX = vertcat(speed_WT_DTX, G);
    elseif string(ALL_XYLOOM_DTX.Geno{i}) == "het"
        F = cell2mat(ALL_XYLOOM_DTX{i,5});
        speed_HET_DTX = vertcat(speed_HET_DTX, F);
    end
end


figure;
% subplot(2,2,1)
imagesc(speed_WT_B4)
n_h = numel(speed_WT_B4(:,1));
hold on 
plot([rowsb4 rowsb4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf rowsb4+lf], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*2 rowsb4+lf*2], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*3 rowsb4+lf*3], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*4 rowsb4+lf*4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*5 rowsb4+lf*5], [0 n_h+1], 'w:', 'LineWidth', 1.5)
axis off
box off
% colormap(redblue)
colormap(gray)
caxis([0 80])
f = gcf; 
% f.Position = [100 600 750 n_h*8];
f.Position = [567   483   549  n_h*15]; 

% subplot(2,2,3)
figure
imagesc(speed_HET_B4)
n_h = numel(speed_HET_B4(:,1));
hold on 
plot([rowsb4 rowsb4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf rowsb4+lf], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*2 rowsb4+lf*2], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*3 rowsb4+lf*3], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*4 rowsb4+lf*4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*5 rowsb4+lf*5], [0 n_h+1], 'w:', 'LineWidth', 1.5)
axis off
box off
colormap(redblue)
caxis([0 80])
f = gcf; 
% f.Position = [100 600 750 n_h*8];
f.Position = [567   483   549  n_h*15]; 


figure
% subplot(2,2,2)
imagesc(speed_WT_DTX)
n_h = numel(speed_WT_DTX(:,1));
hold on 
plot([rowsb4 rowsb4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf rowsb4+lf], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*2 rowsb4+lf*2], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*3 rowsb4+lf*3], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*4 rowsb4+lf*4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*5 rowsb4+lf*5], [0 n_h+1], 'w:', 'LineWidth', 1.5)
axis off
box off
colormap(redblue)
caxis([0 80])
f = gcf; 
% f.Position = [100 600 750 n_h*8];
f.Position = [567   483   549  n_h*15]; 

figure
% subplot(2,2,4)
imagesc(speed_HET_DTX)
n_h = numel(speed_HET_DTX(:,1));
hold on 
plot([rowsb4 rowsb4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf rowsb4+lf], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*2 rowsb4+lf*2], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*3 rowsb4+lf*3], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*4 rowsb4+lf*4], [0 n_h+1], 'w:', 'LineWidth', 1.5)
plot([rowsb4+lf*5 rowsb4+lf*5], [0 n_h+1], 'w:', 'LineWidth', 1.5)
axis off
box off
colormap(redblue)
caxis([0 80])
f = gcf; 
% f.Position = [100 600 750 n_h*8];
f.Position = [567   483   549  n_h*15]; 




% If you want to add the colourbar, uncommment the lines below: 

% c = colorbar;
% c.Location = 'southoutside';
% c.Label.String = 'Speed (cm s^-1)';
% ax = gca;
% ax.FontSize = 18;




%% 4 - EXITS - saline vs dtx - loom trials 

all_animals = unique(exit_analysis.Animal);
n_animals = numel(all_animals);

valsWT1=[];
valsHET1 = [];
valsWT2=[];
valsHET2 = [];

figure
for i = 1:n_animals
    
    all_ani = find((exit_analysis.Animal) == all_animals(i));
    
    genoo = exit_analysis.geno(all_ani(1));
    if genoo == 1
        col = 'k';
    else
        col = 'r';
    end
    
    ani_before = find((exit_analysis.Animal) == all_animals(i) & cell2mat(exit_analysis.data) == 1);
    ani_after = find((exit_analysis.Animal) == all_animals(i) & cell2mat(exit_analysis.data) == 2);
    
    if ~isempty(ani_before) && ~isempty(ani_after)
        val_before = nanmean((exit_analysis.NumOut(ani_before)));
        val_after = nanmean((exit_analysis.NumOut(ani_after)));
        
        if genoo == 1
            valsWT1 = [valsWT1, val_before];
            valsWT2 = [valsWT2, val_after];
        else
            valsHET1 = [valsHET1, val_before];
            valsHET2 = [valsHET2, val_after];
        end 
        
        plot( [1+rand(1)/2,3+rand(1)/2], [val_before, val_after], '-', 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3) %'MarkerFaceColor', [1 0.8 0.8]
        hold on

    end

end 


%% BAR CHART - Exits during loom trials (WT-B4, HET-B4, WT-DTX, HET-DTX) with points of individual animals on top. 

x = [1,2,3,4];
y = [nanmean(valsWT1), nanmean(valsWT2), nanmean(valsHET1) nanmean(valsHET2)];
y2 = [median(valsWT1), median(valsHET1), median(valsWT2), median(valsHET2)];

col = [0.25 0.25 0.25];

figure
b = bar(x,y);
hold on
scatter(ones(numel(valsWT1),1), valsWT1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsHET1),1)*3, valsHET1,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsWT2),1)*2, valsWT2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
scatter(ones(numel(valsHET2),1)*4, valsHET2,'SizeData', 100, 'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.5)
b.FaceColor = [0.7 0.7 0.7];

for i = 1:numel(valsWT1)
    plot([1,2], [valsWT1(i), valsWT2(i)], '-', 'Marker', 'none', 'Color', col, 'LineWidth', 1.3)
end 

for i = 1:numel(valsHET1)
    plot([3,4], [valsHET1(i), valsHET2(i)], '-', 'Marker', 'none', 'Color', col, 'LineWidth', 1.3)
end 


box off
ylim([0 7])
xlim([0 5])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.FontSize = 14;

f = gcf;
f.Position = [440   391   306   407];


%% STATISTICS - EXITS 

% 1 - Are the data normally distributed?
% Small samples therefore use shapiro wilk test 

[h, p] = swtest(valsWT1)
%   0.2117 % ACCEPT NULL

[h, p] = swtest(valsWT2)
% 0.5141 % ACCEPT NULL

[h, p] = swtest(valsHET1)
% 0.2603 % ACCEPT NULL

[h, p] = swtest(valsHET2)
% 0.0127 % REJECT NULL


% PAIRED! Non-parametric 
[p,h] = ranksum(valsWT1, valsWT2)
[p,h] = ranksum(valsHET1, valsHET2)
% PAIRED - parametric
[h, p] = ttest(valsWT1, valsWT2)
[h,p] = ttest(valsHET1, valsHET2)


% NOT-PAIRED - 
% non parametric version of not-paired. 
[h,p] = ranksum(valsWT1, valsHET2)
[h,p] = ranksum(valsHET1, valsWT2)
[h,p] = ranksum(valsHET1, valsWT1)
[h,p] = ranksum(valsHET2, valsWT2)


save('EXITS_DATA_DTX_Setd5.mat', 'valsWT1', 'valsWT2', 'valsHET1', 'valsHET2');

%% ANOVA

rxn1 = [(valsWT1),(valsHET1)];
rxn2 = [(valsWT2), (valsHET2)];
genooo = {'wt', 'wt', 'wt', 'wt', 'wt', 'wt' 'het', 'het', 'het', 'het', 'het', 'het', 'het'};

t1 = array2table(rxn1', 'VariableNames', {'Rxn1'});
t2 = array2table(rxn2', 'VariableNames', {'Rxn2'});
gp = array2table(genooo', 'VariableNames', {'Geno'});

t3 = [t1, t2, gp];

rm = fitrm(t3, 'Rxn1-Rxn2~Geno') % 'WithinDesign', within, 'WithinModel', 'separatemeans' 
ranovatbl = ranova(rm)
tbl = multcompare(rm, 'Geno')

%

maxsp = [(valsWT1),(valsHET1), (valsWT2), (valsHET2)];
genooo = {'wt', 'wt', 'wt', 'wt', 'wt', 'wt','wt', 'het', 'het', 'het', 'het', 'het', 'het', 'het', 'het', 'wt', 'wt', 'wt', 'wt', 'wt', 'het', 'het', 'het', 'het', 'het', 'het'};
DTX = {'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'saline', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX', 'DTX'};

tA = array2table(maxsp', 'VariableNames', {'MaxSp'});
gp = array2table(genooo', 'VariableNames', {'Geno'});
gp2 = array2table(DTX', 'VariableNames', {'DTX'});

tB = [tA, gp, gp2];
rm = fitrm(tB, 'MaxSp~DTX*Geno') 
ranovatbl = ranova(rm)
tbl = multcompare(rm, 'Geno', 'By', 'DTX')
tbl = multcompare(rm, 'DTX', 'By', 'Geno')












%% 5 - EXITS - Stripey dot plots - during loom trials with saline vs with dtx. 

% Each row is the trial of an individual mouse.
% TOP = WT animals. 
% BOTTOM  = HET animals. 
% Uses 'exits' rather than 'exit_analysis'

% Set day and experiments. 
figure
data = 2;

n = height(exit_analysis);
all_animals = unique(exit_analysis.Animal);
n_animals = numel(all_animals);

allWT = find((exit_analysis.geno) == 1);
allHET = find((exit_analysis.geno) == 2);

all_animals_WT = unique(exit_analysis.Animal(allWT));
all_animals_HET = unique(exit_analysis.Animal(allHET));

nWT = numel(all_animals_WT);
nHET = numel(all_animals_HET);


for i = 1:2
    if i == 1 % WT  
        ax1 =  subplot(n_animals,1,1:nWT);
        y_val = 1;
%         rectangle('Position', [900 0 9100 18], 'FaceColor', [1 0.75 0.75, 0.3], 'EdgeColor', 'none')
        
        for j = 1:n_animals
            ani = all_animals(j);
            
            if ismember(ani, all_animals_WT)
                
                rectangle('Position', [0,y_val-0.25,15000,0.5], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none')
                    hold on
                    
                allrows = find(exits.Animal == ani & cell2mat(exits.data) == data);
                n_exits = numel(allrows);
                
                if n_exits >0
                    
                    
                    alltstamps = exits.FrameOut(allrows);
                    y_vals = ones(n_exits,1)*y_val;
                    hold on
                    
                    for k = 1:n_exits
                        if exits.LoomTrig(allrows(k))==1
                            mark = '.';
                            sz = 30;
                        elseif exits.LoomTrig(allrows(k))==0
                            mark = 'o';
                            sz = 9;
                        end
                        
                        plot(alltstamps(k), y_vals(k), mark, 'Color', 'k', 'MarkerSize', sz);
                    end
%                     y_val = y_val +1;
                end
                y_val = y_val +1;
            end
            
        end
        
        box off
        ax1.XAxis.Visible = 'off';
        ax1.YAxis.Visible = 'off';
%         ax1.YLim = ([0 nWT+1]);
%         ax1.YTick = [1:1:nWT];
        ax1.YLim = ([0 10]);
        ax1.YTick = [1:1:9];
        ax1.LineWidth = 1.2;
        ax1.FontSize = 18; 
        ax1.TickDir = 'out';
%         xlim([900 10000])
        
        
    elseif i ==2 % HET 
        ax2 =  subplot(n_animals,1,nWT+1:n_animals);
%         rectangle('Position', [900 0 9100 18], 'FaceColor', [1 0.75 0.75, 0.3], 'EdgeColor', 'none')
        y_val = 1;
        
        for j = 1:n_animals
            ani = all_animals(j);
            
            if ismember(ani, all_animals_HET)
                rectangle('Position', [0,y_val-0.25,15000,0.5], 'FaceColor', [1 0 0 0.3], 'EdgeColor', 'none') % CUl3 [1 0 1 0.1] % Ptchd1 [1 0.8 0.55 0.3]
                    hold on
                    
                allrows = find(exits.Animal == ani & cell2mat(exits.data) == data);
                n_exits = numel(allrows);
                
                if n_exits >0
                    
                    alltstamps = exits.FrameOut(allrows);
                    y_vals = ones(n_exits,1)*y_val;
                    hold on
                    
                    for k = 1:n_exits
                        if exits.LoomTrig(allrows(k))==1
                            mark = '.';
                            sz = 30;
                        elseif exits.LoomTrig(allrows(k))==0
                            mark = 'o';
                            sz = 9;
                        end
                        
                        plot(alltstamps(k), y_vals(k), mark, 'Color', [1 0 0], 'MarkerSize', sz);  % [1 0.6 0.25]
                    end
%                     y_val = y_val +1;
                end
                y_val = y_val +1;
            end
        end
        box off
        %     ax2.YLim = ([0 nHET+1]);
        %     ax2.YTick =(1:1:nHET);
        ax2.YLim = ([0 8]);
        ax2.YTick =(1:1:7);
        ax2.XTick = (15:60:240)*60; % frames
        ax2.XTickLabels = ({'0', '1', '2', '3', '4'});
        xlabel('Time (minutes)')
        ax2.LineWidth = 1.2;
        ax2.FontSize = 18;
        ax2.TickDir = 'out';
        ax2.YAxis.Visible = 'off';
        %     xlim([900 10000])
    end

end


f = gcf;
f.Position = [638   193   303   532]; %[680   605   677   493];


%% 
% exit_analysis

allWT = find((exit_analysis.geno) == 1);
allHET = find((exit_analysis.geno) == 2);

all_animals_WT = unique(exit_analysis.Animal(allWT));
all_animals_HET = unique(exit_analysis.Animal(allHET));

% xy_analysis 
allWT = find(string(xy_analysis.Geno) == "wt" & cell2mat(xy_analysis.data) == 2);
allHET = find(string(xy_analysis.Geno) == "het" & cell2mat(xy_analysis.data) == 2);

all_animals_WT = unique(xy_analysis.Animal(allWT));
all_animals_HET = unique(xy_analysis.Animal(allHET));

%%


























































%% EXITS - Extra code

%% STATS

[h,p] = ttest2(valsWT1, valsWT2)
[h,p] = ttest2(valsHET1, valsHET2)

[h,p] = ttest2(valsWT1, valsHET2)
[h,p] = ttest2(valsWT2, valsHET1)

% 

[h,p] = ranksum(valsWT1, valsWT2)
[h,p] = ranksum(valsHET1, valsHET2)

[h,p] = ranksum(valsWT1, valsHET2)
[h,p] = ranksum(valsWT2, valsHET1)


%% Exits made by animals during the LOOM experiments for days - first with saline - then second with DTX

% In exit_analysis - the column 'data' corresponds to before (1) or after
% (dtx). and data1 = the nM conc of DTX used. 

all_animals = unique(exit_analysis.Animal);
n_animals = numel(all_animals);

valsWT1=[];
valsHET1 = [];
valsWT2=[];
valsHET2 = [];

figure
for i = 1:n_animals
    
    all_ani = find((exit_analysis.Animal) == all_animals(i));
    
    genoo = exit_analysis.geno(all_ani(1));
    if genoo == 1
        col = 'k';
    else
        col = 'r';
    end
    
    ani_before = find((exit_analysis.Animal) == all_animals(i) & cell2mat(exit_analysis.data) == 1);
    ani_after = find((exit_analysis.Animal) == all_animals(i) & cell2mat(exit_analysis.data) == 2);
    
    if ~isempty(ani_before) && ~isempty(ani_after)
        val_before = nanmean((exit_analysis.NumOut(ani_before)));
        val_after = nanmean((exit_analysis.NumOut(ani_after)));
        
        if genoo == 1
            valsWT1 = [valsWT1, val_before];
            valsWT2 = [valsWT2, val_after];
        else
            valsHET1 = [valsHET1, val_before];
            valsHET2 = [valsHET2, val_after];
        end 
        
        plot( [1+rand(1)/2,3+rand(1)/2], [val_before, val_after], '-', 'Marker', 'o', 'Color', col, 'MarkerSize', 12, 'LineWidth', 1.3) %'MarkerFaceColor', [1 0.8 0.8]
        hold on

    end

end 


box off
xlim([0.25 4.25])
ylim([-0.5 7])
xticks([1.25 3.25])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.FontSize = 14;
xticklabels({'Saline', 'DTX'})
ylabel('Shelter exits')

f = gcf;
f.Position = [776   541   253   415]; 



%% BAR CHART (WT-B4, HET-B4, WT-DTX, HET-DTX) with points of individual animals on top. 

% valsWT1=[];
% valsHET1 = [];
% valsWT2=[];
% valsHET2 = [];

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
ylim([0 7])
xlim([0 5])
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.2;
ax.FontSize = 14;

f = gcf;
f.Position = [440   391   306   407];


%% STATISTICS

% valsWT1=[];
% valsHET1 = [];
% valsWT2=[];
% valsHET2 = [];

% 1 - Are the data normally distributed?
% Small samples therefore use shapiro wilk test 

[h, p] = swtest(valsWT1)
% 0.0108 % REJECT NULL

[h, p] = swtest(valsWT2)
% 0.2969 % ACCEPT NULL

[h, p] = swtest(valsHET1)
% 0.0237 % REJECT NULL

[h, p] = swtest(valsHET2)
% 0.1305 % ACCEPT NULL

% 2 - Not normal therefore Kruskal Wallis TEST - CANNOT DO THIS BECAUSE
% DATA ARE NOT INDEPENDENT.
% kruskalwallis(x) returns the p-value for the null hypothesis that the data in each column of the matrix x comes from the same distribution.

% x1 = [valsWT1, valsHET1]; % SALINE
% x2 = [valsWT2, valsHET2]; % DTX
% 
% % 4 groups
% x = [x1, x2]';
% g = [ones(1,6), ones(1,7)*2, ones(1,6)*3, ones(1,7)*4]';
% 
% % 2 groups 
% x = [x1', x2'];
% g = [ones(1,6), ones(1,7)*2]';
% g = [g,g];
% g2 = [ones(13,1)', (ones(13,1)*2)']';
% 
% gg = [g, g2];
% 
% % Compare saline vs dtx - by genotype
% [p,tbl,stats] = kruskalwallis(x, gg)
% multcompare(stats)


% PAIRED! Non-parametric 
[h,p] = ranksum(valsWT1, valsWT2)
[h,p] = ranksum(valsHET1, valsHET2)

% NOT-PAIRED - 
[h,p] = ttest2(valsWT1, valsHET2)
[h,p] = ttest2(valsWT2, valsHET1)
% non parametric version of not-paired. 
[h,p] = ranksum(valsWT1, valsHET2)
[h,p] = ranksum(valsHET1, valsWT2)


% Number of WT/ HET animals
allWT = find(exit_analysis.geno == 1);
allHET = find(exit_analysis.geno == 2);

WTanimals = unique(exit_analysis.Animal(allWT)); % 9 animals 
HETanimals = unique(exit_analysis.Animal(allHET)); %  9 animals

% Number of animals with both BEFORE and AFTER. 
% 6 WT values
% 7 HET values 
