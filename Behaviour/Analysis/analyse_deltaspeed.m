% Analyse Arrest
% 'DIP' in speed/ arrest behaviour in response to loom stimulus.
% Burnett - 13/04/22 

% Uses the data : 'Setd5_xy_return_XYLOOM_sortedl2m.mat' - 'xy_return'

%% Speed at versus Speed after. 
allWT = (find(string(xy_return.Geno)=="wt")); %& xy_return.Trial <= ntrials));
allHET = (find(string(xy_return.Geno)=="het")); 

rowsb4 = 175; 

for i = 1:n_ret
    
    % G = speed data
    G = cell2mat(ALL_XYLOOM{i,5});
    
    sp_at = nanmean(G(rowsb4-3:rowsb4+3));
    sp_immed = nanmean(G(rowsb4+15:rowsb4+30));% 45
    maxxx = max(G(rowsb4:end));
    
    xy_return.sp_at{i} = sp_at;
    xy_return.sp_immed{i} = sp_immed;
    xy_return.maxxx{i} = maxxx;

end 
      
% 30 frames = 500ms. 
% 10 frames = 166ms 

% nanmean(cell2mat(xy_return.sp_at(allWT)))
% nanmean(cell2mat(xy_return.sp_at(allHET)))
% nanmean(cell2mat(xy_return.sp_immed(allWT)))
% nanmean(cell2mat(xy_return.sp_immed(allHET)))
%     

%% Plot sp_at versus sp_imm - All trials - all animals - for a particular DAY

figure
for i = 1:height(xy_return)
    
    g = xy_return.Geno{i};
    t = xy_return.Trial(i);
    d = xy_return.Day(i);
    
%     if d==1
    
    if g == "wt"
    subplot(1,2,1)
    y1 = xy_return.sp_at{i};
    y2 = xy_return.sp_immed{i};
    plot([1 2], [y1 y2], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2)
    hold on;
    box off
    ylim([-5 70])
    xlim([0.5 2.5])
    xticks([1 2])    
    
    else
    subplot(1,2,2) 
    y1 = xy_return.sp_at{i};
    y2 = xy_return.sp_immed{i};
    plot([1 2], [y1 y2], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', [1 0.7 0.7], 'LineWidth', 1.2)
    hold on;
    box off  
    ylim([-5 70])
    xlim([0.5 2.5])
    xticks([1 2])
    
%     end 
    
    end 
 
end 
    

%% SETD5
% rowsb4 = 175; 
% 
% for i = 1:n_ret
%     
%     % G = speed data
%     G = cell2mat(ALL_XYLOOM{i,5});
%     
%     sp_at = nanmean(G(rowsb4-10:rowsb4+5));
%     sp_immed = nanmean(G(rowsb4+15:rowsb4+45));% 45
%     maxxx = max(G(rowsb4:end));
%     
%     xy_return.sp_at{i} = sp_at;
%     xy_return.sp_immed{i} = sp_immed;
%     xy_return.maxxx{i} = maxxx;
% 
% end 

%% CUL3

rowsb4 = 175; 

for i = 1:n_ret
    
    % G = speed data
    G = cell2mat(ALL_XYLOOM{i,5});
    
    sp_at = nanmean(G(rowsb4-10:rowsb4+3));
    sp_immed = nanmean(G(rowsb4+15:rowsb4+45));% 45
    maxxx = max(G(rowsb4:end));
    
    xy_return.sp_at{i} = sp_at;
    xy_return.sp_immed{i} = sp_immed;
    xy_return.maxxx{i} = maxxx;

end 

%% L2M = 'looms to maximum speed' 

lf = 46;

for i = 1:n_ret
    t2max = (xy_return.TimeToMaxSp{i}*60)+rowsb4; %+0.1667
    if t2max <= lf+rowsb4
        xy_return.L2M{i} = 1; 
    elseif t2max <= (lf*2)+rowsb4
        xy_return.L2M{i}= 2;
    elseif t2max <= (lf*3)+rowsb4
        xy_return.L2M{i} = 3;
    elseif t2max <= (lf*4)+rowsb4
        xy_return.L2M{i} = 4;
    elseif t2max <= (lf*5)+rowsb4
        xy_return.L2M{i} = 5;
    elseif t2max >(lf*5)+rowsb4 && xy_return.MaxSpEscape{i}>=25
        xy_return.L2M{i}= 6;
    elseif t2max >(lf*5)+rowsb4 && xy_return.MaxSpEscape{i}<25 % NO RESPONSE 
        xy_return.L2M{i} = 7;
    end 
end 



%% STATS - paried t-test

allWT1 = (find(string(xy_return.Geno)=="wt" & cell2mat(xy_return.L2M) == 1)); %& xy_return.Trial <= ntrials));
allHET1 = (find(string(xy_return.Geno)=="het" & cell2mat(xy_return.L2M) == 1)); %& xy_return.Trial <= ntrials));
 
allWT2 = (find(string(xy_return.Geno)=="wt" & cell2mat(xy_return.L2M) > 1)); %& xy_return.Trial <= ntrials));
allHET2 = (find(string(xy_return.Geno)=="het" & cell2mat(xy_return.L2M) > 1)); %& xy_return.Trial <= ntrials));


wt1 = cell2mat(xy_return.sp_at(allHET1));
wt2 = cell2mat(xy_return.sp_immed(allHET1));

% Observations are not independent... animal averages / regression? 
[p, h] = ttest2(wt1, wt2)

% Test for normality - Shapiro Wilk test - small sample number
[h, p]= swtest(wt2)



%% PLOT - Animal average of 'sp_at' and 'sp_imm' across trials where the animal responds within one loom, or after the first loom.  

all_animals = unique(xy_return.Animal);
n_animals = numel(all_animals);

valsWT1at = [];
valsWT1im = [];

valsWT2at = [];
valsWT2im = [];

valsHET1at = [];
valsHET1im = [];

valsHET2at = []; 
valsHET2im = []; 

 figure

for ii = 1:n_animals
    
    ani = all_animals{ii};
    
    % Trials respond within one loom 
    allANI1 = find(string(xy_return.Animal) == ani & cell2mat(xy_return.L2M)==1  & cell2mat(xy_return.maxxx)<100 & cell2mat(xy_return.sp_at)>5); 
    % Trials where the mouse responds after 1 loom 
    allANI2 = find(string(xy_return.Animal) == ani & cell2mat(xy_return.L2M)>1  & cell2mat(xy_return.maxxx)<100 & cell2mat(xy_return.sp_at)>5); 
    
    allani = vertcat(allANI1, allANI2);
    g = xy_return.Geno{allani(1)};
    
    if ~isempty(allANI1) % Respond within 1Loom
        
        y1 = nanmean(cell2mat(xy_return.sp_at(allANI1)));
        y2 = nanmean(cell2mat(xy_return.sp_immed(allANI1)));
        
        if g=="wt"
            % WT 1
            subplot(2,2,1)
            col = 'k';
            valsWT1at = [valsWT1at, y1];
            valsWT1im = [valsWT1im, y2];
            
        elseif g ~= "wt"
            subplot(2,2,2)
%             col = 'm';
            col = 'b'; %[1 0.45 0.13];
            valsHET1at = [valsHET1at, y1];
            valsHET1im = [valsHET1im, y2];
        end
        
        plot([1 2], [y1 y2], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 0.7)
        hold on;
        box off
        ylim([-5 60])
        xlim([0.5 2.5])
        xticks([1 2])
        ax = gca;
        ax.TickDir  = 'out';
        ax.LineWidth = 1.2;
        ax.TickLength = [0.03 0.03];
    end
    
    
    if ~isempty(allANI2) % Respond after 1 loom 
        
        y1 = nanmean(cell2mat(xy_return.sp_at(allANI2)));
        y2 = nanmean(cell2mat(xy_return.sp_immed(allANI2)));
        
        
        if g=="wt"
            % WT 1
            subplot(2,2,3)
            col = 'k';
%             col = [0.7 0.7 0.7];
            valsWT2at = [valsWT2at, y1];
            valsWT2im = [valsWT2im, y2];
            
        elseif g ~= "wt"
            subplot(2,2,4)
%             col = [1 0.7 0.7];
%             col = 'm';
            col = 'b'; %[1 0.45 0.13];
            valsHET2at = [valsHET2at, y1];
            valsHET2im = [valsHET2im, y2];
        end
        
        plot([1 2], [y1 y2], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 0.7)
        hold on;
        box off
        ylim([-5 60])
        xlim([0.5 2.5])
        xticks([1 2])
        ax = gca;
        ax.TickDir  = 'out';
        ax.LineWidth = 1.2;
        ax.TickLength = [0.03 0.03];
    end
    
    
end

f = gcf;
f.Position = [680   544   401   554]; %[680   374   471   724]; 


%%  STATS - ANIMALS

% valsWT1at = [];
% valsWT1im = [];
% 
% valsWT2at = [];
% valsWT2im = [];
% 
% valsHET1at = [];
% valsHET1im = [];
% 
% valsHET2at = []; 
% valsHET2im = []; 

[p, h] = ttest2(valsWT1at, valsWT1im)
[p, h] = ttest2(valsWT2at, valsWT2im)
[p, h] = ttest2(valsHET1at, valsHET1im)
[p, h] = ttest2(valsHET2at, valsHET2im)


[p, h] = ranksum([valsWT1at, valsWT2at], [valsHET1at, valsHET2at])
[p, h] = ranksum([valsWT1at], [valsHET1at])
[p, h] = ranksum([valsWT1im], [valsHET1im])
[p, h] = ranksum([valsWT2im], [valsHET2im])

%% Plot ALL TRIALS - not animal averages

figure
for i = 1:n_ret
    
    g = xy_return.Geno{i};
    t = xy_return.Trial(i);
    d = xy_return.Day(i);
    maxsp = xy_return.maxxx{i};
    spat = xy_return.sp_at{i};
    l2m = xy_return.L2M{i};
    
    if maxsp < 100 && spat>2 %&& spat < 25
%     if d==5
    
    if g == "wt"
        
        if l2m == 1
            col = 'k';
        else 
            col = [0.7 0.7 0.7];
        end 
        
        if l2m ==1
            subplot(2,2,1)
        else
            subplot(2,2,3)
        end
    y1 = xy_return.sp_at{i};
    y2 = xy_return.sp_immed{i};
%     y3 = xy_return.maxxx{i};
%     plot([1 2 3], [y1 y2 y3], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 1.2)
    plot([1 2], [y1 y2], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 0.7)
    hold on;
    box off
    ylim([-5 60])
    xlim([0.5 2.5])
    xticks([1 2])
        ax = gca;
        ax.TickDir  = 'out';
        ax.LineWidth = 1.2;
        ax.TickLength = [0.02 0.02];    
    
    else
        
        if l2m ==1
            subplot(2,2,2)
        else
            subplot(2,2,4)
        end
    
        if l2m == 1
            col = 'm';
        else 
            col = [1 0.7 1];
        end 
        
    y1 = xy_return.sp_at{i};
    y2 = xy_return.sp_immed{i};
%     y3 = xy_return.maxxx{i};
%     plot([1 2 3], [y1 y2 y3], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 1.2)
    plot([1 2], [y1 y2], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 0.7)
    hold on;
    box off  
    ylim([-5 60])
    xlim([0.5 2.5])
    xticks([1 2])
        ax = gca;
        ax.TickDir  = 'out';
        ax.LineWidth = 1.2;
        ax.TickLength = [0.02 0.02];
    
    end 
    
    end 
    
%     end 
 
end 

f = gcf;
f.Position = [680   374   471   724]; 



%% Speed at - speed im plots for trials sorted by L2M group 
% Like in Figure Suppl. 1k. 

n_ret = height(xy_return);

figure
for i = 1:n_ret
    
    g = xy_return.Geno{i};
    t = xy_return.Trial(i);
    d = xy_return.Day(i);
    maxsp = xy_return.maxxx{i};
    spat = xy_return.sp_at{i};
    l2m = xy_return.L2M{i};
    
    if maxsp < 100 && spat>5 %&& maxsp > 35
%     if d==5
    
    if g == "wt" && l2m<7
        
        if l2m == 1
            col = [0.6 0.6 0.6];
        else 
            col = [0.6 0.6 0.6];
        end 
        
        if l2m ==1
            subplot(2,6,1);
            ylim([-5 60])
        elseif l2m ==2 
            subplot(2,6,2);
        elseif l2m ==3
            subplot(2,6,3);
        elseif l2m ==4
            subplot(2,6,4);
        elseif l2m ==5
            subplot(2,6,5);
        elseif l2m > 5
            subplot(2,6,6);
        end
    y1 = xy_return.sp_at{i};
    y2 = xy_return.sp_immed{i};
%     y3 = xy_return.maxxx{i};
%     plot([1 2 3], [y1 y2 y3], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 1.2)
    plot([1 2], [y1 y2], '-', 'Color', col, 'LineWidth', 0.7, 'Marker', '.', 'MarkerSize', 12) % 'Marker', '.', 'MarkerSize', 30, 
    hold on;
    box off
    ylim([-5 60])
    xlim([0.5 2.5])
    xticks([1 2])
    xticklabels({''})
        ax = gca;
        ax.TickDir  = 'out';
        ax.LineWidth = 1.2;
        ax.TickLength = [0.03 0.03];    
    
    elseif g ~= "wt" && l2m<7
        
        if l2m ==1
            subplot(2,6,7);
            ylim([-5 60])
        elseif l2m ==2 
            subplot(2,6,8);
        elseif l2m ==3
            subplot(2,6,9);
        elseif l2m ==4
            subplot(2,6,10);
        elseif l2m ==5
            subplot(2,6,11);
        elseif l2m == 6
            subplot(2,6,12);
        end
    
        if l2m == 1
%             col = [1 0.4 1];
            col = [1 0.45 0.13];
        else 
%             col = [1 0.4 1];
            col = [1 0.45 0.13];
        end 
        
    y1 = xy_return.sp_at{i};
    y2 = xy_return.sp_immed{i};
%     y3 = xy_return.maxxx{i};
%     plot([1 2 3], [y1 y2 y3], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 1.2)
    plot([1 2], [y1 y2], '-', 'Color', col, 'LineWidth', 0.7, 'Marker', '.', 'MarkerSize', 12)
    hold on;
    box off  
    ylim([-5 60])
    xlim([0.5 2.5])
    xticks([1 2])
        ax = gca;
        ax.TickDir  = 'out';
        ax.LineWidth = 1.2;
        ax.TickLength = [0.03 0.03];
        xticklabels({''})
    
    end 
    
    end 
    
%     end 
 
end 

f = gcf;
f.Position = [259   436   720   352];%[432   111   692   549]; 


%% FIND DATA VALUES

L2M = 6;

wt_rows = find(string(xy_return.Geno)=="wt" & cell2mat(xy_return.maxxx)<100 & cell2mat(xy_return.sp_at)>5 & cell2mat(xy_return.L2M)==L2M);
wt_sp_at_values = xy_return{wt_rows, 32};
wt_sp_im_values = xy_return{wt_rows, 33};

het_rows = find(string(xy_return.Geno)=="het" & cell2mat(xy_return.maxxx)<100 & cell2mat(xy_return.sp_at)>5 & cell2mat(xy_return.L2M)==L2M);
het_sp_at_values = xy_return{het_rows, 32};
het_sp_im_values = xy_return{het_rows, 33};

%% STATS - per l2m group
[h,p] = ttest2(wt_sp_at_values, wt_sp_im_values)
n = numel(wt_sp_im_values)

%% Plot example traces to show how spat and spimm are found.
rowsb4 = 175;

row_valA = 152;
row_valB = 101;

close
figure; 
hold on 
rectangle('Position', [rowsb4-8 -130 10 250], 'FaceColor', [0.3 0.8 1 0.4], 'EdgeColor', 'none')
rectangle('Position', [rowsb4+20 -130 25 250], 'FaceColor', [0.3 0.3 1 0.4], 'EdgeColor', 'none')
plot(cell2mat(ALL_XYLOOM{row_valA,5}), 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2); 
plot((cell2mat(ALL_XYLOOM{row_valB,5})-30), 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2); 
plot([rowsb4-3 rowsb4-3], [-135 80], 'k:', 'LineWidth', 1)
ylim([-130 80])
% xticklabels({''})
yticklabels({''})
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
% title(strcat(string(row_valA), '-', string(row_valB)))

f = gcf;
f.Position = [430   522   682   283];
 
    