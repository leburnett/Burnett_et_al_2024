% Analyse Arrest
% 'DIP' in speed/ arrest behaviour in response to loom stimulus.
% Burnett - 13/04/22

% DATA HERE
save('/Users/lauraburnett/Documents/Burnett_etal/DATA/SETD5/BEHAVIOUR/ARREST/220414_all_xy_analysis_ARREST_Setd5.mat', 'all_xy_analysis');

%% Initialise variables
% all_xy_analysis = xy_analysis;

n = height(all_xy_analysis); 
% details = ALL_XYLOOM_TABLE(:, [1,2,3,4,6,7]);
rowsb4 = 175;  
% fps = frames per second
fps = 60; 
% lf = loom frames - how many frames each loom lasts for
lf = 46; 

% col = the colour to use for plotting HET variables. 
col = 'r'; 
% col = 'm'; % CUL3
% col = [255/255 114/255 32/255]; % PTCHD1

allWT = find(string(all_xy_analysis.Geno) == "wt");
allHET = find(string(all_xy_analysis.Geno) == "het");

all_animals_WT = unique(all_xy_analysis.Animal(allWT));
all_animals_HET = unique(all_xy_analysis.Animal(allHET));


%% VISUALISE INDIVIDUAL SPEED TRACES
i = i+1;
x = cell2mat(ALL_XYLOOM_TABLE{i,5}); 
vel = diff(x);
figure
plot(x)
hold on 
plot(vel)
plot([rowsb4 rowsb4], [-20 100], 'r:')

%% Add COLUMNS 
% IMMED SPEED CHANGE - Changes in speed at the time of loom versus in the period just after loom starts. 
% MAX SPEED CHANGE - maximum speed during escape / speed at loom. 

% 1 - Speed at loom = mean speed from 5 frames before loom presentation and 5 frames after loom presentation. 
% 2 - Speed immediately after = 300-800ms after loom starts. Frame 18 - 48 after. 0.5s time frame. The mean of the speed during these frames. 
% 3 - Max Speed = maximum absolute speed reached during escape to shelter. 

for i = 1:numel(ALL_XYLOOM_TABLE(:,1)) %numel(all_xy_analysis(:,1))
    
    xx = []; 
    
    x = cell2mat(ALL_XYLOOM_TABLE{i,5}); 
    ALL_XYLOOM_TABLE.Sp{i} = x;
    
    vel = diff(x);
    ALL_XYLOOM_TABLE.Vel{i} = vel;
    
    sp_at = nanmean(x(rowsb4-5:rowsb4+5));
    sp_immed = nanmean(x(rowsb4+18:rowsb4+48)); %nanmean
    
    all_xy_analysis.speedat{i} = sp_at;
    all_xy_analysis.speed_immed{i} = sp_immed;
    all_xy_analysis.DeltaImmed{i} = sp_immed/sp_at; 
%     all_xy_analysis.DeltaMax{i} = all_xy_analysis.MaxSpEscape{i}/sp_at;

    ALL_XYLOOM_TABLE.DeltaImmed{i} = sp_immed/sp_at; 

end 

%% Find out how many frames are < 2cm/s  in the period between loom and return to shelter - how long were these bouts for. 

for i = 1:numel(all_xy_analysis(:,1))
    
    xx = []; % Reset
    x = cell2mat(ALL_XYLOOM_TABLE{i,5});
    
    framestoshelter = all_xy_analysis.TimeToShelter{i}*60;
    
    speed_loom_to_shelter = x(rowsb4:rowsb4+framestoshelter);
    xx = find(speed_loom_to_shelter<2); %Find the frames - between loom start and mouse enters the shelter that the speed < 2cm/s
    
    % To check that the mouse reaches a speed of > 20cm/s - 'ESCAPE'
    % otherwise it's a very slow shuffle to the shelter. 
    xx2 = find(speed_loom_to_shelter>20);
    
    if ~isempty(xx)
        xx(2,1:end-1) = diff(xx(1,:));
        
        smallgap = find(xx(2,:)>1 & xx(2,:)<4);
        if ~isempty(smallgap)
            xx(2, smallgap) = 1;
            xx(2,1:end-1) = diff(xx(1,:));
        end
        
        runstarts = find(xx(2,:)>1);

        frames_frozen = diff(runstarts)+1; % length of bouts in frames where mouse <2cm/s
        all_xy_analysis.FramesFrozenBouts{i} = frames_frozen;
        
        %Remove periods of <2cm/s which last less than 15 frames
        frames_frozen2 = frames_frozen(frames_frozen>15);
        if ~isempty(xx)
            all_xy_analysis.freezeimmed{i} = 1;
        else
            all_xy_analysis.freezeimmed{i} = 0;
        end
        
    else
        all_xy_analysis.FramesFrozenBouts{i} = 0;
    end
   
    % Colum of 'Y/N' if mouse runs with > 20cm/s back to shelter. 
    if ~isempty(xx2)
        all_xy_analysis.morethan20{i} = 1;
    else 
        all_xy_analysis.morethan20{i} = 0;
    end 
    
end 


%% SCATTER - Delta Speed - Immediate versus Max 

figure
for i = 1:n
    if (all_xy_analysis.speedat{i})>2 %& all_xy_analysis.Trial(i)<=3
    if all_xy_analysis.Geno{i} == "wt"
        col = 'k';
    else
        col = 'r';
    end 
    
    x1 = ((all_xy_analysis.DeltaImmed{i})); % Immed
    y1 = (all_xy_analysis.DeltaMax{i}); % Max
    plot(y1, x1, 'o', 'Color', col)
    hold on 
    end 
end 
plot([0 40], [1 1], 'k:')
plot([1 1], [0 40], 'k:')
plot([0 40], [0 40], 'k:')
box off
xlabel('Max Escape Speed / Speed at Loom')
ylabel('Mean Immediate Speed / Speed at Loom')
ylim([0 25])
xlim([0 30])
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;
ax.FontSize = 16;


%% Delta Speed - Immediate versus Max Absolute Speed

figure
for i = 1:n
%     if (all_xy_analysis.speedat{i})>2 %& all_xy_analysis.Trial(i)<=3
    if all_xy_analysis.Geno{i} == "wt"
        col = 'k';
    else
        col = 'r';
    end 
    
    x1 = log(cell2mat(all_xy_analysis.DeltaImmed(i)));
    y1 = cell2mat(all_xy_analysis.MaxSpEscape(i));
    plot(y1, x1, 'o', 'Color', col)
    hold on 
%     end 
end 
plot([0 150], [0 0], 'k:')
% plot([1 1], [0 70], 'k:')
box off
xlabel('Max Escape Speed (cm s^-1)')
ylabel('Speed Immediate/ Speed at Loom')
% ylim([0 20])
xlim([0 75])
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;
ax.FontSize = 16;


%% LOG VERSION - ALL TRIALS - POOLED ACROSS ANIMALS
figure

for i = 1:n
    if (all_xy_analysis.speedat{i})>2 & (all_xy_analysis.MaxSpEscape{i})>20  & (all_xy_analysis.ReturnToShelter{i})==1   %& all_xy_analysis.Trial(i)<=5
    if all_xy_analysis.Geno{i} == "wt"
        col = 'k';
    else
        col = 'r';
    end 
    
    x1 = log((all_xy_analysis.DeltaImmed{i}));
    y1 = (all_xy_analysis.MaxSpEscape{i});
    plot(y1,x1, 'o', 'Color', col)
    hold on 
    end 
end 
plot([0 150], [0 0], 'k:')

box off
xlabel('Max Speed (cm s^-1)')
ylabel('Log(Speed Immediately After/Speed At Loom)')
 ylim([-4 4])
 xlim([0 100])
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;
ax.FontSize = 16;

% Create textbox
% annotation(gcf,'textbox',[0.191476190476191 0.829383886255924 0.135054421768707 0.0497630331753559],'String',{'< 5 Trials'},'FitBoxToText','off');

%%
for k = 1:height(xy_return)
%     xy_return.LogDeltaImmed(k) = log(cell2mat(xy_return.DeltaImmed(k)));
%     xy_return.MS(k) = cell2mat(xy_return.MaxSpEscape(k));
xy_return.G(k) = xy_return.Geno{k};
end 


%% MARGINAL HISTOGRAM AND SCATTER PLOT 

all10 = find(cell2mat(all_xy_analysis.ReturnToShelter) ==1 & cell2mat(all_xy_analysis.speedat)>2 & cell2mat(all_xy_analysis.MaxSpEscape)>20); %& cell2mat(all_xy_analysis.speedat)>2

xvals = log(cell2mat(all_xy_analysis.DeltaImmed(all10))); 
yvals = cell2mat(all_xy_analysis.MaxSpEscape(all10));
gp = string(all_xy_analysis.Geno(all10));

% [1 0.2 1; 0.3 0.3 0.3]
%1 0.4471 0.1255
% cul3 - [1 0.2 1; 0.3 0.3 0.3]
% [0.3 0.3 0.3; 1 0.3 0.3]

scatterhist(xvals, yvals, 'Group', gp, 'Kernel', 'overlay', 'Color', [1 0.3 0.3; 0.4 0.4 0.4], 'Style', 'bar', 'LineStyle', '-', 'Location', 'SouthWest', 'Direction', 'out', 'MarkerSize', 7)
scatterhistogram(xy_return, 'LogDeltaImmed', 'MS', 'GroupVariable', 'G', 'HistogramDisplayStyle', 'smooth', 'LineStyle', '-', 'LineWidth', 1.5, 'MarkerStyle', 'o', 'LegendVisible', 'off', 'Color', [1 0.3 0.3; 0.4 0.4 0.4])
box off
hold on 
plot([0 0], [0 120], 'k:', 'LineWidth', 1.5)
% xlabel('Max Speed (cm s^-1)')
% ylabel('Log(Speed Immediately After/Speed At Loom)')
xlim([-4 4])
ylim([0 120])
ax = gca;
ax.TickDir = 'out';
% ax.LineWidth = 1.5;
% ax.FontSize = 16;

f = gcf;
f.Position = [535   300   660   651]; 


%% MARGINAL HISTOGRAM AND SCATTER PLOT - Only first 10 trials per animals

% all10 = find(all_xy_analysis.Trial <=10);
all10 = find(cell2mat(all_xy_analysis.ReturnToShelter) ==1 & cell2mat(all_xy_analysis.speedat)>2);

xvals = log(cell2mat(all_xy_analysis.DeltaImmed(all10))); 
yvals = cell2mat(all_xy_analysis.MaxSpEscape(all10));
gp = string(all_xy_analysis.Geno(all10));

scatterhist(yvals, xvals, 'Group', gp, 'Kernel', 'overlay', 'Color', [1 0.2 0.2; 0.3 0.3 0.3], 'Style', 'bar', 'LineStyle', '-', 'Location', 'SouthWest', 'Direction', 'out');

% s = scatterhistogram(yvals, xvals, 'GroupData', gp,  'MarkerStyle', 'o', 'Color', [1 0.1 0.1; 0.3 0.3 0.3], 'HistogramDisplayStyle', 'bar', 'MarkerFilled', 'off', 'ScatterPlotLocation', 'SouthWest');
% s.BinWidths = [15; 0.5];

box off
hold on 
plot([0 150], [0 0], 'k:', 'LineWidth', 1.5)
xlabel('Max Speed (cm s^-1)')
ylabel('Log(Speed Immediately After/Speed At Loom)')
ylim([-4 4])
xlim([0 110])
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.6;
ax.FontSize = 16;
ax.TickLength = [0.02 0.02];

f = gcf;
f.Position = [535   300   660   651]; 


%% STATS 
WTrows = find(cell2mat(all_xy_analysis.ReturnToShelter) ==1 & cell2mat(all_xy_analysis.speedat)>2 & string(all_xy_analysis.Geno)=="wt");
HETrows = find(cell2mat(all_xy_analysis.ReturnToShelter) ==1 & cell2mat(all_xy_analysis.speedat)>2 & string(all_xy_analysis.Geno)=="het");

WTvals = log(cell2mat(all_xy_analysis.DeltaImmed(WTrows)));
HETvals = log(cell2mat(all_xy_analysis.DeltaImmed(HETrows)));

[h,p] = kstest2(WTvals, HETvals)
[p, h] = ranksum(WTvals, HETvals)

 %% LOG MEDIAN PER ANIMAL

 ani_vals = zeros(28, 3);
 
figure
% WT
for i = 1:14
    ani = all_animals_WT{i};
    all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial<=10 & cell2mat(all_xy_analysis.ReturnToShelter) ==1); %& all_xy_analysis.Trial<5
    col = [0.4 0.4 0.4];
    x1 = log(nanmedian(cell2mat(all_xy_analysis.DeltaImmed(all_ani))));
    y1 = nanmedian(cell2mat(all_xy_analysis.MaxSpEscape(all_ani)));
    plot(y1, x1, 'o', 'Color', col)
    hold on 
    
    ani_vals(i,1) = x1;
    ani_vals(i,2) = y1;
    ani_vals(i,3) = 1;
    
end 

% HET
for i = 1:14
    ani = all_animals_HET{i};
    all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial<=10 & cell2mat(all_xy_analysis.ReturnToShelter) ==1);
    col = 'r';
    x1 = log(nanmedian(cell2mat(all_xy_analysis.DeltaImmed(all_ani))));
    y1 = nanmedian(cell2mat(all_xy_analysis.MaxSpEscape(all_ani)));
    plot(y1, x1, 'o', 'Color', col)
    hold on 
    
    ani_vals(i+14,1) = x1;
    ani_vals(i+14,2) = y1;
    ani_vals(i+14,3) = 2;
end 

% plot([0 30], [1 1], 'k:')
% plot([1 1], [0 70], 'k:')
box off
xlabel('Median Max Escape Speed')
ylabel('Log (Median Immed Sp/ Speed at Loom)')
plot([0 150], [0 0], 'k:')
 ylim([-2.5 2.5])
 xlim([0 120])
 title('<= 10 Trials')
%  title('1st trial')
 ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;
ax.FontSize = 16;
f = gcf;
f.Position = [440   428   419   370]; 

%%

xvals = ani_vals(:,1);
yvals = ani_vals(:,2);
gp = ani_vals(:,3);

scatterhist(yvals, xvals, 'Group', gp, 'Kernel', 'overlay', 'Color', [0.3 0.3 0.3; 1 0.2 0.2], 'Style', 'bar', 'LineStyle', '-', 'Location', 'SouthWest', 'Direction', 'out');
hold on 
box off
xlabel('Median Max Escape Speed')
ylabel('Log (Median Immed Sp/ Speed at Loom)')
plot([0 150], [0 0], 'k:')
 ylim([-2.5 2.5])
 xlim([0 100])
 title('<= 10 Trials')
%  title('1st trial')
 ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;
ax.FontSize = 16;
f = gcf;
f.Position = [440   428   419   370]; 


%% AVERAGE PER ANIMAL
% 
figure
% WT
for i = 1:14
    ani = all_animals_WT{i};
    all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial<5);
    col = 'k';
    x1 = nanmean(cell2mat(all_xy_analysis{all_ani,24}));
    y1 = nanmean(cell2mat(all_xy_analysis{all_ani, 25}));
    plot(y1, x1, 'o', 'Color', col)
    hold on 
end 

% HET
for i = 1:14
    ani = all_animals_HET{i};
    all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial<5);
    col = 'r';
    x1 = nanmean(cell2mat(all_xy_analysis{all_ani,24}));
    y1 = nanmean(cell2mat(all_xy_analysis{all_ani, 25}));
    plot(y1, x1, 'o', 'Color', col)
    hold on 
end 
% 
% plot([0 30], [1 1], 'b')
% plot([1 1], [0 70], 'b')
% box off
% xlabel('Speed Change - Max')
% ylabel('Speed Change - Immediate')
% ylim([0 20])
% xlim([0 30])
% 

%%

% %% MEDIAN
% 
% figure
% % WT
% for i = 1:14
%     ani = all_animals_WT{i};
%     all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial< 5); %& all_xy_analysis.Trial<5
%     col = [0.4 0.4 0.4];
%     x1 = nanmedian(cell2mat(all_xy_analysis{all_ani,37}));
%     y1 = nanmedian(cell2mat(all_xy_analysis{all_ani, 38}));
%     plot(y1, x1, 'o', 'Color', col)
%     hold on 
% end 
% 
% % HET
% for i = 1:14
%     ani = all_animals_HET{i};
%     all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial< 5);
%     col = 'r';
%     x1 = nanmedian(cell2mat(all_xy_analysis{all_ani,37}));
%     y1 = nanmedian(cell2mat(all_xy_analysis{all_ani, 38}));
%     plot(y1, x1, 'o', 'Color', col)
%     hold on 
% end 
% 
% plot([0 30], [1 1], 'k:')
% plot([1 1], [0 70], 'k:')
% box off
% xlabel('Speed Change - Max')
% ylabel('Speed Change - Immediate')
% ylim([0 10])
% xlim([0 20])


% %% MAX 
% 
% figure
% % WT
% for i = 1:14
%     ani = all_animals_WT{i};
%     all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial< 5); %& all_xy_analysis.Trial<5
%     col = 'k';
%     x1 = max(cell2mat(all_xy_analysis{all_ani,37}));
%     y1 = max(cell2mat(all_xy_analysis{all_ani, 38}));
%     plot(y1, x1, 'o', 'Color', col)
%     hold on 
% end 
% 
% % HET
% for i = 1:14
%     ani = all_animals_HET{i};
%     all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial< 5);
%     col = 'r';
%     x1 = max(cell2mat(all_xy_analysis{all_ani,37}));
%     y1 = max(cell2mat(all_xy_analysis{all_ani, 38}));
%     plot(y1, x1, 'o', 'Color', col)
%     hold on 
% end 
% 
% plot([0 30], [1 1], 'b')
% plot([1 1], [0 70], 'b')
% box off
% xlabel('Speed Change - Max')
% ylabel('Speed Change - Immediate')
% ylim([0 10])
% xlim([0 20])
% 
% 
% %% MIN
% 
% figure
% % WT
% for i = 1:14
%     ani = all_animals_WT{i};
%     all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial< 5); %& all_xy_analysis.Trial<5
%     col = 'k';
%     x1 = min(cell2mat(all_xy_analysis{all_ani,37}));
%     y1 = min(cell2mat(all_xy_analysis{all_ani, 38}));
%     plot(y1, x1, 'o', 'Color', col)
%     hold on 
% end 
% 
% % HET
% for i = 1:14
%     ani = all_animals_HET{i};
%     all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & all_xy_analysis.Trial< 5);
%     col = 'r';
%     x1 = min(cell2mat(all_xy_analysis{all_ani,37}));
%     y1 = min(cell2mat(all_xy_analysis{all_ani, 38}));
%     plot(y1, x1, 'o', 'Color', col)
%     hold on 
% end 
% 
% plot([0 30], [1 1], 'b')
% plot([1 1], [0 70], 'b')
% box off
% xlabel('Speed Change - Max')
% ylabel('Speed Change - Immediate')
% ylim([0 10])
% xlim([0 20])



%% Boxplot + DOTS - log(delat_immed) - animal averages. 

nWT = numel(all_animals_WT);
nHET = numel(all_animals_HET);

WTVALS = zeros(nWT,1);
HETVALS = zeros(nHET, 1);

% WT
for i = 1:nWT
    ani = all_animals_WT{i};
    all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2 & cell2mat(all_xy_analysis.ReturnToShelter) ==1); %& all_xy_analysis.Trial<5
    x1 = log(nanmean(cell2mat(all_xy_analysis.DeltaImmed(all_ani))));
%     x1 = (nanmean(cell2mat(all_xy_analysis.MaxSpEscape(all_ani))));
    WTVALS(i) = x1;
end 

% HET
for i = 1:nHET
    ani = all_animals_HET{i};
    all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.speedat) >2  & cell2mat(all_xy_analysis.ReturnToShelter) ==1); %& all_xy_analysis.Trial<=10
    x1 = log(nanmean(cell2mat(all_xy_analysis.DeltaImmed(all_ani))));
%     x1 = (nanmean(cell2mat(all_xy_analysis.MaxSpEscape(all_ani))));
    HETVALS(i) = x1;
end 

x1 = ones(nWT,1);
x2 = ones(nHET,1)*2;

VALS = vertcat(WTVALS, HETVALS);
GRP = vertcat(x1,x2);

% Figure
figure
scatter(x1, WTVALS,'SizeData', 150, 'MarkerEdgeColor', [0.2 0.2 0.2], 'jitter', 'on', 'jitterAmount', 0.05)
hold on
scatter(x2, HETVALS', 'SizeData', 150,'MarkerEdgeColor', [0.2 0.2 1] , 'jitter', 'on', 'jitterAmount', 0.05) % [1 0.2 1] % [1 0.7 0.5]
b = boxplot(VALS, GRP, 'Colors', 'k');
set(b, 'LineWidth', 1.5)

xticks([1,2])
xticklabels({''})
ax = gca;
ax.FontSize = 25;
box off
xlim([0.5 2.5])
ylim([20 85])
% axis([0.5 2.5 -2.5 2.5])
% ax.XAxis.Visible = 'off';
% ylabel('Number of Exits')
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.LineWidth = 1.5;
f = gcf;
f.Position = [680   844   250   400];

%% STATS
[p,h] = ranksum(WTVALS, HETVALS)

% Setd5Emx1 % % % % % %

% LOG DELTA IMMED
% p =  0.5338

% MaxSpEscape 
% p = 0.4452

% Reation Time
% p = 0.3869


% CUL3 % % % % % %

% LOG DELTA IMMED
% p =  0.0039 

% MaxSpEscape 
% p = 6.2894X10-4

% Reation Time
% p = 0.0011

% PTCHD1 % % % % % %

% LOG DELTA IMMED
% p =  0.0039 



%% POOLED TRIALS

allWT = find(string(all_xy_analysis.Geno) == "wt" & cell2mat(all_xy_analysis.speedat) >2 & cell2mat(all_xy_analysis.ReturnToShelter) ==1  & cell2mat(all_xy_analysis.MaxSpEscape) >20); 
allHET = find(string(all_xy_analysis.Geno) == "het" & cell2mat(all_xy_analysis.speedat) >2 & cell2mat(all_xy_analysis.ReturnToShelter) ==1  & cell2mat(all_xy_analysis.MaxSpEscape) >20); 

WTVALS = ((cell2mat(all_xy_analysis.TimeToMaxSp(allWT))));
HETVALS = ((cell2mat(all_xy_analysis.TimeToMaxSp(allHET))));

[p,h] = ranksum(WTVALS, HETVALS)

%%



%% STATISTICAL TESTS 

% Number of WT/ HET trials in xy_return - double check above! Should be the same :) 
allWT = find(string(all_xy_analysis.Geno) == "wt");
allHET = find(string(all_xy_analysis.Geno) == "het");

dataWT = log(cell2mat(all_xy_analysis.DeltaImmed(allWT)));
dataHET = log(cell2mat(all_xy_analysis.DeltaImmed(allHET)));

% 1 - plot data
% figure; hist(data);

% 2 - test normality with Ks test. - more than 30 samples.
[h,p] =kstest(dataWT)
[h,p] =kstest(dataHET)

% 2 - test normality with K-S test. - less than 30 samples.
[H, pValue] = swtest(dataWT)
[H, pValue] = swtest(dataHET)

% 3 - Levene test for homogenity of variance. 
nWT = numel(dataWT);
nHET = numel(dataHET);

col1 = vertcat(dataWT, dataHET);
col2 = vertcat(ones(nWT,1), ones(nHET,1)*2);

X = horzcat(col1, col2);

%Levene's test for variance:
Levenetest(X)

% Multiple different tests comparing the mean:
% Mann whitney U test:
mwwtest(dataWT', dataHET')
% Wilcoxon RankSum
ranksum(dataWT, dataHET)
% Welch t-test with unequal variance
[p, h] = ttest2(dataWT, dataHET, 'Vartype', 'unequal')


%% Finding animal avreage of log delta immed
all_animals = unique(all_xy_analysis.Animal);

speed_WT = [];
speed_HET = []; 

for j = 1:numel(all_animals)
    ani = all_animals{j};
    speed_ani = [];
    
    for i = 1:n
        if  string(all_xy_analysis.Animal(i)) == ani
            G = log(cell2mat(all_xy_analysis.DeltaImmed(i)));
            speed_ani = vertcat(speed_ani, G);
            genoo = string(all_xy_analysis.Geno{i});
        end
    end
    
    if genoo == "wt"
        speed_WT = vertcat(speed_WT, nanmean(speed_ani));
    else
        speed_HET = vertcat(speed_HET, nanmean(speed_ani));
    end
    
end


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

% Wilcoxon RankSum
nanmean(speed_WT)
nanmean(speed_HET)
ranksum(speed_WT, speed_HET)








