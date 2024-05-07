% Extra code that generates plots that didn't make it into the paper

% Originally found in 'analyse_adaptation.m'

%% P1 - Dot for Mean + SEM error vertical bar - with LINE FIT - TRIAL VERSUS LOOMS 2 ESCAPE 
col = 'm';

allWT = find(string(xy_return.Geno) =="wt");
allHET = find(string(xy_return.Geno) =="het");

% WT 
n_trialsWT = max(xy_return.Trial(allWT));
dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM

for i = 1:n_trialsWT
    vals = find(xy_return.Trial == i & string(xy_return.Geno) == "wt");
    data = cell2mat(xy_return.Looms2Max(vals)); 
    dWT(i) = mean(data);
    sWT(i) = std(data)/sqrt(numel(vals));
end

% HET 
n_trialsHET = max(xy_return.Trial(allHET));
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

for i = 1:n_trialsHET
    vals = find(xy_return.Trial == i & string(xy_return.Geno) == "het");
    data = cell2mat(xy_return.Looms2Max(vals));
    dHET(i) = mean(data);
    sHET(i) = std(data)/sqrt(numel(vals));
end

dHET = dHET(1:n_trialsWT);
sHET = sHET(1:n_trialsWT);
n_trialsHET = n_trialsWT; 

% PLOT 

%Pthcd1 = . [1 0.8 0.05], setd5 = [1 0 0], cul3 = [1, 0 1]

figure
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0], 'MarkerFaceColor', col , 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.7 0.5],'Marker', 'none', 'LineWidth', 1.75)

% CREATE FIT LINES - LINEAR FIT - y = ax+b % % % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % % 

[fitresult, gof] = createFits(dWT, dHET, col);
legend off

axis([0 n_trialsWT+1 0 6])
box off
% yticks([1,2,3,4,5])
xlabel('Trial')
ylabel('Looms to Escape')
ax = gca;
ax.FontSize = 18; 
hold off
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.025 0.025];
f =gcf;
f.Position = [2143 212 340 273];

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
errorbar(1,sWT(1),(ciWT(2,1)-sWT(1)), 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
hold on 
errorbar(1,sHET(1),(ciHET(2,1)-sHET(1)),'Color', [1 0.4 0.4],  'LineWidth', 1.5);
rectangle('Position', [0.85,sWT(1),0.3,0.01],'FaceColor', 'k', 'EdgeColor', 'none')
rectangle('Position', [0.85,sHET(1),0.3,0.01],'FaceColor', col, 'EdgeColor', 'none')
box off
ax = gca;
ax.FontSize = 18; 
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ylim([-4 0.7])
% yticks([ -0.1,-0.05, 0,0.05, 0.1,0.15, 0.2])
% ylim([-0.1 0.15])
ax.XAxis.Visible = 'off';
xlim([0.75 1.25])
ax.TickLength = [0.025 0.025];
f =gcf;
f.Position = [  767   634   130   196];

%Save the fig.
% saveas(f, 'MaxSp_Trial_Error_Of_Slope.svg')


%% P6 - ERRORBAR - VALUES ACROSS DAYS - ONLY TRIALS WHERE RETURN TO SHELTER - All trials not Animal Averages
 
% Full code: 
data_per_day = zeros(2,5); 
sem_per_day = zeros(2,5); 
n_per_day  =zeros(2,5); 
stats_per_day  =zeros(1,5); 
% 
allWT = find(string(xy_return.Geno) =="wt");
allHET = find(string(xy_return.Geno) =="het");

val = 13; %Column of xy_return you want to asses: Return to shelter. 

for j = 1:5
    allWT = find(string(xy_return.Geno) =="wt" & (xy_return.Day) == j); 
    allHET = find(string(xy_return.Geno) =="het" & (xy_return.Day) == j); 
    data_WT = cell2mat(xy_return{allWT, val}); 
    data_HET = cell2mat(xy_return{allHET, val}); 
%        data_WT = table2array(xy_return(allWT, val)); 
%     data_HET = table2array(xy_return(allHET, val)); 
    data_per_day(1,j) = nanmean(data_WT); 
    data_per_day(2,j) = nanmean(data_HET); 
    n_per_day(1,j) = numel(data_WT); 
    n_per_day(2,j) = numel(data_HET); 
    sem_per_day(1,j) = std(data_WT)/sqrt(numel(data_WT));
    sem_per_day(2,j) = std(data_HET)/sqrt(numel(data_HET));
    if numel(data_WT)>1 && numel(data_HET)>1
    [p, h] = ranksum(data_WT, data_HET);
    stats_per_day(1,j) = p; 
    end 
end 


% figure
% hold off
% errorbar((data_per_day(1,:)), (sem_per_day(1,:)), 'k', 'LineWidth', 1.4, 'Marker', '.', 'MarkerSize', 30, 'MarkerFaceColor', 'k')
% hold on 
% errorbar((data_per_day(2,:)), (sem_per_day(2,:)), 'Color', col, 'LineWidth', 1.4, 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w')

figure
hold off
errorbar((data_per_day(1,:)), (sem_per_day(1,:)), 'k', 'LineWidth', 1.5) %'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'none', 'CapSize', 0
hold on 
errorbar((data_per_day(2,:)), (sem_per_day(2,:)), 'Color', col, 'LineWidth', 1.5) %'Marker', 'o', 'MarkerSize', 10,  'MarkerFaceColor', 'none', 'CapSize', 0

xticks([1,2,3,4,5])
box off
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickDir = 'out';
xlim([0.5 5.5])
% ylim([30 90])
% ylim([0 2.5])
ylim([-0.05 2.2])
% ylabel('Escape Probability (%)')
% ylabel('Speed (cm s^-1)')
% ylabel('Looms to Escape')
ylabel('Time (s)')
% ylabel('ylabel')
% xticklabels({})
xlabel('Day')
f = gcf;
f.Position = [680   844   319   254];



%% P8 - Cumulative Histograms - Looms2Escape across days 

% WT
allD1 = find(string(all_xy_analysis.Geno) =="wt" & all_xy_analysis.Day == 1);
allD2 = find(string(all_xy_analysis.Geno) =="wt" & all_xy_analysis.Day == 2);
allD3 = find(string(all_xy_analysis.Geno) =="wt" & all_xy_analysis.Day == 3);
allD4 = find(string(all_xy_analysis.Geno) =="wt" & all_xy_analysis.Day == 4);
allD5 = find(string(all_xy_analysis.Geno) =="wt" & all_xy_analysis.Day == 5);

figure
histogram(cell2mat(all_xy_analysis.Looms2Max(allD1)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0 0 0],  'LineWidth', 2)
hold on 
% histogram(cell2mat(all_xy_analysis.Looms2Max(allD2)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [1 0.25 0],  'LineWidth', 2)
histogram(cell2mat(all_xy_analysis.Looms2Max(allD3)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.5 0.5 0.5],  'LineWidth', 2)
% histogram(cell2mat(all_xy_analysis.Looms2Max(allD4)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [1 0.75 0],  'LineWidth', 2)
histogram(cell2mat(all_xy_analysis.Looms2Max(allD5)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.8 0.8 0.8],  'LineWidth', 2)
format_figure('Looms to Escape', 'Cumul. Probability',  [])
legend('Day1', 'Day3', 'Day5')
xticks([1,2,3,4,5,6])
f = gcf;
f.Position = [795   660   287   298];

% HET
allD1 = find(string(all_xy_analysis.Geno) =="het" & all_xy_analysis.Day == 1);
allD2 = find(string(all_xy_analysis.Geno) =="het" & all_xy_analysis.Day == 2);
allD3 = find(string(all_xy_analysis.Geno) =="het" & all_xy_analysis.Day == 3);
allD4 = find(string(all_xy_analysis.Geno) =="het" & all_xy_analysis.Day == 4);
allD5 = find(string(all_xy_analysis.Geno) =="het" & all_xy_analysis.Day == 5);

figure
histogram(cell2mat(all_xy_analysis.Looms2Max(allD1)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [1 0 0],  'LineWidth', 2)
hold on 
% histogram(cell2mat(all_xy_analysis.Looms2Max(allD2)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [1 0.25 0],  'LineWidth', 1.2)
histogram(cell2mat(all_xy_analysis.Looms2Max(allD3)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [1 0.5 0],  'LineWidth', 2)
% histogram(cell2mat(all_xy_analysis.Looms2Max(allD4)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [1 0.75 0],  'LineWidth', 1.2)
histogram(cell2mat(all_xy_analysis.Looms2Max(allD5)),0.5:1:6.5, 'Normalization','cdf', 'DisplayStyle', 'stairs', 'EdgeColor', [1 0.9 0],  'LineWidth', 2)
format_figure('Looms to Escape', 'Cumul. Probability',  [])
legend('Day1', 'Day3', 'Day5')
xticks([1,2,3,4,5,6])
f = gcf;
f.Position = [795   660   287   298];



%% P1 - Dot for Mean + SEM error vertical bar - with LINE FIT - TRIAL VERSUS LOOMS 2 ESCAPE 
col = 'r';

allWT = find(string(ALL_XYLOOM_TABLE.Geno) =="wt");
allHET = find(string(ALL_XYLOOM_TABLE.Geno) =="het");

% WT 
n_trialsWT = max(ALL_XYLOOM_TABLE.Trial(allWT));
dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM

for i = 1:n_trialsWT
    vals = find(ALL_XYLOOM_TABLE.Trial == i & string(ALL_XYLOOM_TABLE.Geno) == "wt");
    data = log(cell2mat(ALL_XYLOOM_TABLE.DeltaImmed(vals))); 
    dWT(i) = mean(data);
    sWT(i) = std(data)/sqrt(numel(vals));
end

% HET 
n_trialsHET = max(ALL_XYLOOM_TABLE.Trial(allHET));
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

for i = 1:n_trialsHET
    vals = find(ALL_XYLOOM_TABLE.Trial == i & string(ALL_XYLOOM_TABLE.Geno) == "het");
    data = log(cell2mat(ALL_XYLOOM_TABLE.DeltaImmed(vals)));
    dHET(i) = mean(data);
    sHET(i) = std(data)/sqrt(numel(vals));
end

% dHET = dHET(1:n_trialsWT);
% sHET = sHET(1:n_trialsWT);
% n_trialsHET = n_trialsWT; 

% PLOT 
dHET = dHET(1:15);
sHET = sHET(1:15);
n_trialsHET = 15; 

dWT = dWT(1:2);
sWT = sWT(1:2);
n_trialsWT = 2; 

%Pthcd1 = . [1 0.8 0.05], setd5 = [1 0 0], cul3 = [1, 0 1]

figure ; hold on
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.7 0.7], 'MarkerFaceColor', [1 0.7 0.7], 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.4],'Marker', 'none', 'LineWidth', 1.75)

% CREATE FIT LINES - LINEAR FIT - y = ax+b % % % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % % 

[fitresult, gof] = createFits(dWT, dHET, col);
legend off

axis([0 15.5 -2.5 2])
box off
% yticks([1,2,3,4,5])
% xlabel('Trial')
% ylabel('Looms to Escape')
ax = gca;
% ax.FontSize = 18; 
hold off
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];
f =gcf;
f.Position = [993   645   385   198]; %[2143 212 340 273];


%%

[rho, p] = corr([1:1:15]',dHET', 'Type', 'Pearson')
[rho, p] = corr([1:1:2]',dWT', 'Type', 'Pearson')

[R, P, RL, RU] = corrcoef(1:1:15, dHET(1:15))

[rho, p] = corr([1:1:15]',dHET', 'Type', 'Spearman')

% save('220904_Setd5_ALL_XYLOOM_REWARD.mat', 'ALL_XYLOOM_TABLE');


allWT = find(ALL_XYLOOM_TABLE.Trial <12 & string(ALL_XYLOOM_TABLE.Geno) == "wt"); 




%% CHECKING REWARD TRIALS 



  all_animals = unique(ALL_XYLOOM_TABLE.Animal); 
  n_animals = numel(all_animals); 

for i = 1:n_animals
    ani = all_animals(i); 
    trial_number = 1; 
    
    for j = 1:height(ALL_XYLOOM_TABLE)
        if string(ALL_XYLOOM_TABLE.Animal{j}) == ani 
            ALL_XYLOOM_TABLE.Trial(j) = trial_number; 
            trial_number = trial_number +1; 
        end   
    end
end 


details = ALL_XYLOOM_TABLE(:, [1,2,3,4,6,7,8]);


% Trend of RT/ MAx SP over repetitions
for j = 1:height(ALL_XYLOOM_TABLE)
    val = ALL_XYLOOM_TABLE.T2M{j};
    v2 = (val-180)/60;
    ALL_XYLOOM_TABLE.T2M{j} = v2;
end 


%%

all_animals = unique(ALL_XYLOOM_TABLE.Animal);

figure
hold on 

for k = 1:numel(all_animals)

    ani= all_animals{k};
    allANI = find(string(ALL_XYLOOM_TABLE.Animal) == ani);

    n_trialsWT = max(ALL_XYLOOM_TABLE.Trial(allANI));
    dWT = zeros(1, n_trialsWT); %data
    sWT = zeros(1, n_trialsWT); % SEM
    
    for i = 1:n_trialsWT
        vals = find(ALL_XYLOOM_TABLE.Trial(allANI) == i);
        data = (cell2mat(ALL_XYLOOM_TABLE.T2M(allANI(vals))));
        dWT(i) = mean(data);
        sWT(i) = std(data)/sqrt(numel(vals));
    end
    
    geno = string(ALL_XYLOOM_TABLE.Geno{allANI(1)});
    
    if geno == "wt"
        errorbar(1:1:n_trialsWT, dWT, sWT, 'o-', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
        errorbar(1:1:n_trialsWT, dWT, sWT, 'o-', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)

    elseif geno == "het"
        errorbar(1:1:n_trialsWT, dWT, sWT, 'o-', 'CapSize', 0, 'Color', [1 0.7 0.7], 'MarkerFaceColor', [1 0.7 0.7], 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
        errorbar(1:1:n_trialsWT, dWT, sWT, 'o-', 'CapSize', 0, 'Color', [1 0.4 0.4],'Marker', 'none', 'LineWidth', 1.75)
    end 
    
end 


axis([0 11.5 0 12])
box off

%%

