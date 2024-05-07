%% Adaptation analysis

% Main script used to plot the figures and analyse the adaptation of the loom response of the mice from Burnett et
% al. 

% Generated by Burnett 26/01/22

% This script requires the following files/arrays/tables:

% ALL_XYLOOM_TABLE
% all_xy_analysis

% ALL_XYLOOM
% xy_return

%%  This script looks at adaptation in the ESCAPE response - both the VIGOUR of the response and the SPEED of initiating a response:

% 1 - across trials
% 2 - across days

% FOR ADAPTATION OF EXIT BEHAVIOUR - i.e. with BANANA - see 'analyse_exits.m'

%% ADAPTATION ACROSS DAYS 

%% ACROSS DAYS - DOT and ERRORBAR LINE 
% Figure 2d

% Full code: 
data_per_day = zeros(2,5); 
sem_per_day = zeros(2,5); 
n_per_day  =zeros(2,5); 
% stats_per_day  =zeros(1,5); 
% 
allWT = find(string(xy_return.Geno) =="wt");
allHET = find(string(xy_return.Geno) =="het");

val = 13; %Column of xy_return you want to asses: Return to shelter. 

for j = 1:5
    allWT = find(string(all_xy_analysis.Geno) =="wt" & (all_xy_analysis.Day) == j); %  & cell2mat(all_xy_analysis.ReturnToShelter) == 1 & cell2mat(all_xy_analysis.speedat)>2); 
    allHET = find(string(all_xy_analysis.Geno) =="het" & (all_xy_analysis.Day) == j); % & cell2mat(all_xy_analysis.ReturnToShelter) == 1 & cell2mat(all_xy_analysis.speedat)>2); 
    data_WT = (cell2mat(all_xy_analysis{allWT, val})); 
    data_HET = (cell2mat(all_xy_analysis{allHET, val})); 
    data_per_day(1,j) = nanmean(data_WT); 
    data_per_day(2,j) = nanmean(data_HET); 
    n_per_day(1,j) = numel(data_WT); 
    n_per_day(2,j) = numel(data_HET); 
    sem_per_day(1,j) = std(data_WT)/sqrt(numel(data_WT));
    sem_per_day(2,j) = std(data_HET)/sqrt(numel(data_HET));
%     if numel(data_WT)>1 && numel(data_HET)>1
%     [p, h] = ranksum(data_WT, data_HET);
%     stats_per_day(1,j) = p; 
%     end 
end 

% [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] 

% PTCHD1
% figure
% errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
% hold on 
% errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] , 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
% errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
% errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0.4 0.2],'Marker', 'none', 'LineWidth', 1.75)

% CUL3
figure
errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0.4 0.4], 'MarkerFaceColor', [1 0.7 0.7] , 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0.4 0.4],'Marker', 'none', 'LineWidth', 1.75)

[fitresult, gof] = createFits(1:1:5, data_per_day(2,:), col);
[fitresult, gof] = createFits(1:1:5, data_per_day(1,:), 'k');
legend off

xticks([1,2,3,4,5])
box off
ax = gca;
ax.FontSize = 22;
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.TickLength = [0.025 0.025];
xlim([0.5 5.5])
% ylim([30 80])
% ylim([0 2])
% ylim([-2 2.5])
% ylabel('Escape Probability (%)')
% ylabel('Maximum Speed (cm s^-1)')
% ylabel('Looms to Escape')
% ylabel('Reaction Time (s)')
% ylabel('Log(Speed Change)')
% xticklabels({})
xlabel('Day')
f = gcf;
% f.Position = [487   540   204   310];
% f.Position = [ 607   361   235   308]; % small
f.Position = [ 928   491   233   209]; % square


%% P7 - DOTBOX PLOT - Number of Trials - WT/HET
% XY_RETURN DOES NOT CONTAIN ALL THE TRIALS - ONLY WHERE RESPOND WITHIN THE 5 LOOMS. 
% Need xy_analysis.
% Figure 1e

n = height(all_xy_analysis);

all_animals = string(unique(all_xy_analysis.Animal));

val = 13; % Which column to assess in all_xy_analysis : T2M

speed_WT = [];
speed_HET = []; 

for j = 1:numel(all_animals)
    ani = all_animals(j);
    speed_ani = [];
    
    for i = 1:n
        if all_xy_analysis.Animal(i) == ani
            G = cell2mat(all_xy_analysis{i, val});
            speed_ani = vertcat(speed_ani, G);
            genoo = string(all_xy_analysis.Geno{i});
        end
    end
    
    if genoo == "wt"
        speed_WT = vertcat(speed_WT, numel(speed_ani));
    else
        speed_HET = vertcat(speed_HET, numel(speed_ani));
    end
    
end

 %%%%% Combine arrays %%%% 
speed_ALL = horzcat(speed_WT, speed_HET);

% PLOT 
figure
n_wt = numel(speed_WT(:,1));
n_het = numel(speed_HET);

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

scatter(x1, speed_WT(:,1),'SizeData', 200, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, speed_HET', 'SizeData', 200,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot(speed_ALL, 'Colors', 'k');
set(b, 'linew', 1.25);

xticks([1,2])
xticklabels({''})
ax = gca;
ax.FontSize = 30;
box off
axis([0.5 2.5 0 80])
ax.XAxis.Visible = 'off'; 
hold off
ax.TickDir = 'out'; 
ylabel('Number of Trials')
ax.LineWidth = 1.8;

f = gcf;
f.Position = [704   207   355   572]; 

[p, h] = ranksum(speed_WT, speed_HET)
nanmean(speed_WT)
nanmean(speed_HET)
% 0.0333



%% P9 - Dot for Mean + SEM error vertical bar - with LINE FIT - DAY  VERSUS VARIABLE 
% Figure Suppl.3 F & M

allWT = find(string(xy_return.Geno) =="wt");
allHET = find(string(xy_return.Geno) =="het");

% WT 
n_trialsWT = 5; 
dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM

for i = 1:5
    vals = find(xy_return.Day == i & string(xy_return.Geno) == "wt");
    data = cell2mat(xy_return.MaxSpEscape(vals)); 
    dWT(i) = mean(data);
    sWT(i) = std(data)/sqrt(numel(vals));
end

% HET 
n_trialsHET = 5;
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

for i = 1:5
    vals = find(xy_return.Day == i & string(xy_return.Geno) == "het");
    data = cell2mat(xy_return.MaxSpEscape(vals));
    dHET(i) = mean(data);
    sHET(i) = std(data)/sqrt(numel(vals));
end


% PLOT 

figure
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.8 0.05], 'MarkerFaceColor', [1 0.8 0.05], 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75)
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.05],'Marker', 'none', 'LineWidth', 1.75)

% CREATE FIT LINES - LINEAR FIT - y = ax+b % % % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % % 
% 
% [fitresult, gof] = createFits(dWT, dHET);
% legend off

axis([0.5 5.5 30 90])
box off
% yticks([1,2,3,4,5])
xticks([1,2,3,4,5])
xlabel('Day')
ylabel('Max. Speed (cm s^-1)')
ax = gca;
ax.FontSize = 18; 
hold off
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.025 0.025];
f =gcf;
f.Position =  [2143 217 193 268]; %[2143 217 287  268]; %[2143 212 340 273];
% NARROW = [2143 217 193 268];

% Save the figure; 
% saveas(f, 'Dot_ERROR_LINESOFFIT_LINEAR_T2M_Trial.svg')


%% P2 - INSET PLOT - 'errorbar' SLOPE OF LINE + 95% CI

% Find slope of WT
sWT = coeffvalues(fitresult{1,1});
ciWT = confint(fitresult{1,1});

sHET = coeffvalues(fitresult{2,1});
ciHET = confint(fitresult{2,1});

% For errorbar need data point and error. 

figure
errorbar(1,sWT(1),(ciWT(2,1)-sWT(1)), 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
hold on 
errorbar(1,sHET(1),(ciHET(2,1)-sHET(1)),'Color', [1 0.4 0.4],  'LineWidth', 1);
rectangle('Position', [0.85,sWT(1),0.3,0.015],'FaceColor', 'k', 'EdgeColor', 'none')
rectangle('Position', [0.85,sHET(1),0.3,0.015],'FaceColor', 'r', 'EdgeColor', 'none')
box off
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ylim([-0.15 0.15])
% yticks([ -0.1,-0.05, 0,0.05, 0.1,0.15, 0.2])
% ylim([-0.1 0.15])
ax.XAxis.Visible = 'off';
xlim([0.75 1.25])
ax.TickLength = [0.025 0.025];
f =gcf;
f.Position = [  767   634   130   196];

%Save the fig.
% saveas(f, 'MaxSp_Trial_Error_Of_Slope.svg')




%% STATS - for finding animal averages. 

all_animals = unique(xy_return.Animal);
n_animals = numel(all_animals);

analysis_animals = zeros(n_animals, 7);

for j = 1:n_animals
    ani = all_animals{j};
    all_ani = find(string(xy_return.Animal) == ani & xy_return.Day == 5);
    if isempty(all_ani)
        all_ani2 = find(string(xy_return.Animal) == ani);
        geno = xy_return.Geno{all_ani2(1)};
    else 
        geno = xy_return.Geno{all_ani(1)};
    end 
    
    % Assign animal number to  first column of table. 
    analysis_animals(j,1) = str2double(ani(3:end)); 
    
    if geno =="wt"
        g = 1;
    elseif geno =="het"
        g = 2;
    end 
    analysis_animals(j,2) = g; 
    cv = 3; 
    
    for k = [7,8,9,12,13]
        all_vals = cell2mat(xy_return{all_ani, k});
%         if k == 8 || k ==9 || k == 18 || k ==17
%             exit_analysis_animals(j, cv) = sum(all_vals);
%         else
            analysis_animals(j,cv) = nanmean(all_vals);
%         end
        cv = cv+1;
    end
    
end 

Animal = analysis_animals(:,1);
Geno = analysis_animals(:,2);
MeanB4 = analysis_animals(:,3);
DShelt = analysis_animals(:,4);
Return = analysis_animals(:,5);
MaxSp = analysis_animals(:,6);
T2M = analysis_animals(:,7);

ana_animals5 = table(Animal, Geno, MeanB4, DShelt, Return, MaxSp, T2M);


%% For pooling trials across animals: 

% % % % % % % % % % % % % Tests: 

% % % Reaction time - T2M

% 1 - WT versus HET - Day 1
allWT = find(string(xy_return.Geno) =="wt" & xy_return.Day == 1);
allHET = find(string(xy_return.Geno) =="het" & xy_return.Day == 1);

WTVALS = cell2mat(xy_return.TimeToMaxSp(allWT));
HETVALS = cell2mat(xy_return.TimeToMaxSp(allHET));

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)

% 1 - WT versus HET - Day 5
allWT = find(string(xy_return.Geno) =="wt" & xy_return.Day == 5);
allHET = find(string(xy_return.Geno) =="het" & xy_return.Day == 5);

WTVALS = cell2mat(xy_return.TimeToMaxSp(allWT));
HETVALS = cell2mat(xy_return.TimeToMaxSp(allHET));

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)


% Look at trend in data across days - Day 1 to Day 5 - WT 

%
allWT = find(string(xy_return.Geno) == "wt");
X = fitlm(xy_return.Day(allWT),  cell2mat(xy_return.TimeToMaxSp(allWT)))

% Simple linear regression was used to test if the day studied significantly predicted the maximum distance travelled.
% The fitted regression model was: Maximum Distance for WT = 18.242 + -1.0793*(hours studied).
% The overall regression was statistically significant (R2 = 0.0605, F(2, 66) = 4.25, p < 0.0431).

% figure; plot(X)
% figure; plot(X2)

% Look at trend in data across days - Day 1 to Day 5 - WT 
allHET = find(string(xy_return.Geno) == "het");
X = fitlm(xy_return.Trial(allHET),  cell2mat(xy_return.TimeToMaxSp(allHET)))

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % Maximum Escape Speed

% 1 - WT versus HET - Day 1
allWT = find(string(xy_return.Geno) =="wt" & xy_return.Day == 1);
allHET = find(string(xy_return.Geno) =="het" & xy_return.Day == 1);

WTVALS = cell2mat(xy_return.MaxSpEscape(allWT));
HETVALS = cell2mat(xy_return.MaxSpEscape(allHET));

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)


% 1 - WT versus HET - Day 5
allWT = find(string(xy_return.Geno) =="wt" & xy_return.Day == 5);
allHET = find(string(xy_return.Geno) =="het" & xy_return.Day == 5);

WTVALS = cell2mat(xy_return.MaxSpEscape(allWT));
HETVALS = cell2mat(xy_return.MaxSpEscape(allHET));

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)


% Look at trend in data across days - Day 1 to Day 5 - WT 

% 
allWT = find(string(xy_return.Geno) == "wt");
X = fitlm(xy_return.Day(allWT),  cell2mat(xy_return.MaxSpEscape(allWT)))

% Simple linear regression was used to test if the day studied significantly predicted the maximum distance travelled.
% The fitted regression model was: Maximum Distance for WT = 18.242 + -1.0793*(hours studied).
% The overall regression was statistically significant (R2 = 0.0605, F(2, 66) = 4.25, p < 0.0431).

% figure; plot(X)
% figure; plot(X2)

% Look at trend in data across days - Day 1 to Day 5 - WT 
allHET = find(string(xy_return.Geno) == "het");
X = fitlm(xy_return.Day(allHET),  cell2mat(xy_return.MaxSpEscape(allHET)))

%% Across presentations
% Look at trend in data across presentations

allHET = find(string(xy_return.Geno) == "het");
X = fitlm(xy_return.Day(allHET),  cell2mat(xy_return.TimeToMaxSp(allHET)))

allWT = find(string(xy_return.Geno) == "wt");
X = fitlm(xy_return.Trial(allWT),  cell2mat(xy_return.TimeToMaxSp(allWT)))

% Max Speed

allWT = find(string(xy_return.Geno) == "wt");
X = fitlm(xy_return.Trial(allWT),  cell2mat(xy_return.MaxSpEscape(allWT)))

allHET = find(string(xy_return.Geno) == "het");
X = fitlm(xy_return.Trial(allHET),  cell2mat(xy_return.MaxSpEscape(allHET)))



%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% %  ONLY THE FIRST TRIAL 

% 1 - MAxSP
allWT = find(string(xy_return.Geno) =="wt" & xy_return.Trial == 1);
allHET = find(string(xy_return.Geno) =="het" & xy_return.Trial == 1);

WTVALS = cell2mat(xy_return.MaxSpEscape(allWT));
HETVALS = cell2mat(xy_return.MaxSpEscape(allHET));

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)

% 1 - Reaction Time
allWT = find(string(xy_return.Geno) =="wt" & xy_return.Trial == 1);
allHET = find(string(xy_return.Geno) =="het" & xy_return.Trial == 1);

WTVALS = cell2mat(xy_return.TimeToMaxSp(allWT));
HETVALS = cell2mat(xy_return.TimeToMaxSp(allHET));

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)


%% BANANA EXP
allWT = find(exit_analysis.geno == 1);
allHET = find(exit_analysis.geno == 2);

% Number of Bouts:
WTVALS = exit_analysis.NumBOUTS(allWT);
HETVALS = exit_analysis.NumBOUTS(allHET);

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)

% Probability of getting Banana
WTVALS = exit_analysis.GotBanana(allWT);
HETVALS = exit_analysis.GotBanana(allHET);

nanmean(WTVALS)
nanmean(HETVALS)
[p,h] = ranksum(WTVALS, HETVALS)

% Trend of RT/ MAx SP over repetitions
for j = 1:height(ALL_XYLOOM_TABLE)
    val = ALL_XYLOOM_TABLE.T2M{j};
    v2 = (val-180)/60;
    ALL_XYLOOM_TABLE.Time2Max{j} = v2;
end 

xx = find(cell2mat(ALL_XYLOOM_TABLE.Time2Max) <0);
ALL_XYLOOM_TABLE(xx,:) = [];

allWT = find(string(ALL_XYLOOM_TABLE.Geno) == "wt");
X = fitlm(ALL_XYLOOM_TABLE.Trial(allWT),  cell2mat(ALL_XYLOOM_TABLE.Time2Max(allWT)))

allHET = find(string(ALL_XYLOOM_TABLE.Geno) == "het");
X = fitlm(ALL_XYLOOM_TABLE.Trial(allHET),  cell2mat(ALL_XYLOOM_TABLE.Time2Max(allHET)))

% Max Speed

allWT = find(string(ALL_XYLOOM_TABLE.Geno) == "wt");
X = fitlm(ALL_XYLOOM_TABLE.Trial(allWT),  cell2mat(ALL_XYLOOM_TABLE.MaxSp(allWT)))

allHET = find(string(ALL_XYLOOM_TABLE.Geno) == "het");
X = fitlm(ALL_XYLOOM_TABLE.Trial(allHET),  cell2mat(ALL_XYLOOM_TABLE.MaxSp(allHET)))






%%  NUMBER OF LOOMS TRIGGERED ACROSS DAYS - Average PER ANIMAL
% Figure 2c

% Full code: 
all_animals = unique(all_xy_analysis.Animal);

% Number of WT/ HET trials in total 
allWT = find(string(all_xy_analysis.Geno) == "wt");
allHET = find(string(all_xy_analysis.Geno) == "het");

% Number of WT/ HET animals in total. 
all_animals_WT = unique(all_xy_analysis.Animal(allWT));
all_animals_HET = unique(all_xy_analysis.Animal(allHET));

data_per_day = zeros(2,5); 
sem_per_day = zeros(2,5); 
n_per_day  =zeros(2,5); 
stats_per_day  =zeros(1,5); 

val = 12; %Column of xy_return you want to asses: Return to shelter. 

for j = 1:5
    
    data_WT = [];
    data_HET = [];
    for k = 1:n_animals
        ani = all_animals(k);
        allani = find(string(all_xy_analysis.Animal)==ani & (all_xy_analysis.Day)==j); %  & cell2mat(all_xy_analysis.ReturnToShelter) == 1 & cell2mat(all_xy_analysis.speedat)>2);
        
        if contains(ani, all_animals_WT)
            data_WT =  [data_WT, numel(allani)];
        else
            data_HET =  [data_HET, numel(allani)];
        end
    end 
    
    data_per_day(1,j) = nanmean(data_WT); 
    data_per_day(2,j) = nanmean(data_HET); 
    n_per_day(1,j) = numel(data_WT); 
    n_per_day(2,j) = numel(data_HET); 
    sem_per_day(1,j) = std(data_WT)/sqrt(numel(data_WT));
    sem_per_day(2,j) = std(data_HET)/sqrt(numel(data_HET));
%     if numel(data_WT)>1 && numel(data_HET)>1
%     [p, h] = ranksum(data_WT, data_HET);
%     stats_per_day(1,j) = p; 
%     end 
end 

% [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] 

% % PTCHD1
% figure
% errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
% hold on 
% errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] , 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
% errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
% errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0.4 0.2],'Marker', 'none', 'LineWidth', 1.75)

% CUL3
figure
errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 12, 'LineWidth', 1.75)
hold on 
errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0.7 1], 'MarkerFaceColor', [1 0.7 1] , 'MarkerEdgeColor','none', 'MarkerSize', 12, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:5, data_per_day(1,:), sem_per_day(1,:), 'o', 'CapSize', 0, 'Color', [0 0 0], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:5, data_per_day(2,:), sem_per_day(2,:), 'o', 'CapSize', 0, 'Color', [1 0 1],'Marker', 'none', 'LineWidth', 1.75)

xticks([1,2,3,4,5])
box off
ax = gca;
ax.FontSize = 22;
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.TickLength = [0.025 0.025];
xlim([0.5 5.5])
ylim([0 6])
% ylim([0 2])
% ylim([-2 2.5])
% ylabel('Escape Probability (%)')
% ylabel('Maximum Speed (cm s^-1)')
% ylabel('Looms to Escape')
% ylabel('Reaction Time (s)')
% ylabel('Log(Speed Change)')
% xticklabels({})
% xlabel('Day')
f = gcf;
% f.Position = [487   540   204   310];
% f.Position = [ 607   361   235   308]; % small
f.Position = [ 928   491   233   209]; % square
% f.Position = [889   274   282   302] % for fig 2



%% FOR STATISTICS - Number of looms triggered over days - animal average. - CORRELATIONS!!! 

    data_WT2 = [];
    data_HET2 = [];
    
for j = 1:5
    
    for k = 1:n_animals
        ani = all_animals(k);
        allani = find(string(all_xy_analysis.Animal) ==ani & (all_xy_analysis.Day)==j); %  & cell2mat(all_xy_analysis.ReturnToShelter) == 1 & cell2mat(all_xy_analysis.speedat)>2);
        
        if contains(ani, all_animals_WT)
            data_WT2 =  vertcat(data_WT2, [j, numel(allani)]);
        else
            data_HET2 =  vertcat(data_HET2, [j, numel(allani)]);
        end
    end 
    
end 

%%  ADD FIT LINES TO THE NUMBER OF TRIALS ACROSS DAYS - ERRORBAR PLOT

    data_WT2 = [];
    data_HET2 = [];
    
for j = 1:5
    
    for k = 1:2
        if k ==1
        allani = find(string(all_xy_analysis.Geno) =="wt" & (all_xy_analysis.Day)==j); %  & cell2mat(all_xy_analysis.ReturnToShelter) == 1 & cell2mat(all_xy_analysis.speedat)>2);
            data_WT2 =  vertcat(data_WT2, [j, numel(allani)/9]);
        elseif k ==2 
            allani = find(string(all_xy_analysis.Geno) =="het" & (all_xy_analysis.Day)==j); %  & cell2mat(all_xy_analysis.ReturnToShelter) == 1 & cell2mat(all_xy_analysis.speedat)>2);
            data_HET2 =  vertcat(data_HET2, [j, numel(allani)/9]);
        end
    end 
    
end 

% WT LINE
X = data_WT2(:, 1);
Y = data_WT2(:, 2);
[xData, yData] = prepareCurveData(X, Y);

% Set up fittype and options.
ft = fittype( 'poly1' );
[fitresult{1}, gof(1)] = fit( xData, yData, ft );
h = plot( fitresult{1}, 'k');


% HET LINE
X = data_HET2(:, 1);
Y = data_HET2(:, 2);
[xData, yData] = prepareCurveData(X, Y);
ft = fittype( 'poly1' );
[fitresult{2}, gof(2)] = fit( xData, yData, ft );
h = plot(fitresult{2}, 'm');

legend off

% [rho, p] = corr(X,Y, 'Type', 'Pearson')
[R, P, RL, RU] = corrcoef(X, Y)
[rho, p] = corr(X,Y, 'Type', 'Spearman')


%% Correlation of variables over days - xy_return

allWT = find(string(xy_return.Geno) == "wt");
allHET = find(string(xy_return.Geno) == "het");

% Versus LOOMS TO ESCAPE 
% % WT
% X = cell2mat(xy_return.Looms2Max(allWT));
% Y = cell2mat(xy_return.MaxSpEscape(allWT));
% 
% % HET
% X = cell2mat(xy_return.Looms2Max(allHET));
% Y = cell2mat(xy_return.MaxSpEscape(allHET));


% VERSUS DAY 
% WT
X = (xy_return.Day(allWT));
Y = cell2mat(xy_return.TimeToMaxSp(allWT));

% HET
X = (xy_return.Day(allHET));
Y = cell2mat(xy_return.TimeToMaxSp(allHET));

% STATS 
% [rho, p] = corr(X,Y, 'Type', 'Pearson')
n = numel(Y)
[R, P, RL, RU] = corrcoef(X, Y)
[rho, p] = corr(X,Y, 'Type', 'Spearman')



%% Correlation of variables over trials - BANANA CHIP EXPERIMENTS - ALL_XYLOOM_TABLE

% allWT = find(string(ALL_XYLOOM_TABLE.Geno) == "wt");
% allHET = find(string(ALL_XYLOOM_TABLE.Geno) == "het");
% 
% % WT
% X = ALL_XYLOOM_TABLE.Trial(allWT);
% Y = cell2mat(ALL_XYLOOM_TABLE.T2M(allWT));
% 
% % HET
% X = ALL_XYLOOM_TABLE.Trial(allHET);
% Y = cell2mat(ALL_XYLOOM_TABLE.T2M(allHET));

% OR FROM XY_ANALYSIS - Cul3 

%   all_animals = unique(xy_analysis.Animal); 
%   n_animals = numel(all_animals); 
% 
% for i = 1:n_animals
%     ani = all_animals(i); 
%     trial_number = 1; 
%     
%     for j = 1:n
%         if string(xy_analysis.Animal{j}) == ani 
%             xy_analysis.Trial(j) = trial_number; 
%             trial_number = trial_number +1; 
%         end   
%     end
% end 

allWT = find(string(xy_analysis.Geno) == "wt");
allHET = find(string(xy_analysis.Geno) == "het");

% WT
X = xy_analysis.Trial(allWT);
Y = cell2mat(xy_analysis.TimeToMaxSp(allWT));

% HET
X = xy_analysis.Trial(allHET);
Y = cell2mat(xy_analysis.TimeToMaxSp(allHET));

% STATS 
% [rho, p] = corr(X,Y, 'Type', 'Pearson')
n = numel(Y)
[R, P, RL, RU] = corrcoef(X, Y)
[rho, p] = corr(X,Y, 'Type', 'Spearman')





