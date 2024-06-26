% Analyse reward trials with banana chip
% Generated by Burnett

% Requires 'exit_analysis' from 'Setd5_RewardTrials_N16_Exit_Analysis.mat'

%% Number of looms triggered
% Figure 2h

col = 'r'; 

% How many animals: 
all_animals = (unique(exit_analysis.Animal));

% Which column to assess:
col_idx = 17; % 'NumBOUTS' 

speed_WT = [];
speed_HET = []; 

for j = 1:numel(all_animals)
    
    % Which animal to assess:
    ani = all_animals(j);
    
    % All trials of that animal
    all_ani = find(exit_analysis.Animal == ani); 
    
    % All the values from those trials
    ani_vals = exit_analysis{all_ani, col_idx};
    
    % The average of these values
    av_ani_val = nanmean(ani_vals);
    
    % Genotype of the animal
    genoo = (exit_analysis.geno(all_ani(1)));
    
    if genoo == 1
        speed_WT = vertcat(speed_WT, av_ani_val);
    else
        speed_HET = vertcat(speed_HET, av_ani_val);
    end
    
end

 %%%%% Combine arrays %%%% 
speed_ALL = horzcat(speed_WT, speed_HET);

% PLOT 
figure
n_wt = numel(speed_WT(:,1));
n_het = numel(speed_HET(:,1));

x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

scatter(x1, speed_WT,'SizeData', 200, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
hold on 
scatter(x2, speed_HET, 'SizeData', 200,'MarkerEdgeColor', col, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 1.2)
b = boxplot(speed_ALL, 'Colors', 'k');
set(b, 'linew', 1.5);

% STATISTICS

[p, h] = ranksum(speed_WT, speed_HET)

meanWT = nanmean(speed_WT)
meanHET = nanmean(speed_HET)



%% STRIPE / DOT PLOTS - Plot of Banana Trial for each mouse. Dots for exits from shelter. - * - * - * - * - * - * - * - * - *
% Figure 2i 

% Each row is the trial of an individual mouse.
% TOP = WT animals. 
% BOTTOM  = HET animals. 
% Uses 'exits' rather than 'exit_analysis'

n = height(exit_analysis);
all_animals = unique(exit_analysis.Animal);
n_animals = numel(all_animals);

allWT = find((exit_analysis.geno) == 1);
allHET = find((exit_analysis.geno) == 2);

all_animals_WT = unique(exit_analysis.Animal(allWT));
all_animals_HET = unique(exit_analysis.Animal(allHET));

nWT = numel(all_animals_WT);
nHET = numel(all_animals_HET);

figure
for i = 1:2
    if i == 1 % WT  
        ax1 =  subplot(n_animals,1,1:nWT);
        y_val = 1;
        
        for j = 1:n_animals
            ani = all_animals(j);
            
            if ismember(ani, all_animals_WT)
                rectangle('Position', [0,y_val-0.25,15000,0.5], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none')
                hold on
                allrows = find(exits.Animal == ani);
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
                    
                end
                y_val = y_val +1;
            end
            
        end
        
        box off
        ax1.XAxis.Visible = 'off';
        ax1.YLim = ([0 nWT+1]);
        ax1.YTick = [1:1:nWT];
        ax1.LineWidth = 1.2;
        ax1.FontSize = 18; 
        ax1.TickDir = 'out';
        ax1.YAxis.Visible = 'off';
        
        
        
    elseif i ==2 % HET 
        ax2 =  subplot(n_animals,1,nWT+1:n_animals);
 
        y_val = 1;
        
        for j = 1:n_animals
            ani = all_animals(j);
            
            if ismember(ani, all_animals_HET)
                rectangle('Position', [0,y_val-0.25,15000,0.5], 'FaceColor', [1 0.4471, 0.1255, 0.1], 'EdgeColor', 'none') % [1 0 1 0.1]
                hold on
                allrows = find(exits.Animal == ani);
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
                        
                        plot(alltstamps(k), y_vals(k), mark, 'Color', col, 'MarkerSize', sz);
                    end
                    
                end
                y_val = y_val +1;
            end
        end
    box off
    ax2.YLim = ([0 nHET+1]);
    ax2.YTick =(1:1:nHET);
    ax2.XTick = (0:60:240)*60; % frames
    ax2.XTickLabels = ({'0', '1', '2',  '3',  '4'});
    xlabel('Time (minutes)')
    ax2.LineWidth = 1.2;
    ax2.FontSize = 18; 
    ax2.TickDir = 'out';
    ax2.YAxis.Visible = 'off';
    end

end

f = gcf;
f.Position = [680   605   677   493]; 



%% TRAJECTORIES during banana experimen. 

% From exit analysis/exits for MULTILOOM - therefore only one row per animal! 
% Run through 'exits' table for multiloom.
% Find the  frame when the mouse passes the threshold 

% Plots the trajectories for the exits where the mouse passes the
% threshold.
% Plots the position of the mouse from when it leaves the shelter til when
% it re-enters the shelter. 
% Plots dot at position when loom happened. 
% Traj before loom = thick and dark
% Traj after = thin adn light grey. 

% WTS

figure
hold on 
for i = 1:height(exits)
    
    %Find fixed variables for this exit.
    ani = exits.Animal(i);
    geno = exits.Geno(i);
    date = exits.Date(i);
    
    if geno == 1
        % Find which row in exit_analysis corresponds to this animal
        row_ana = find(exit_analysis.Animal == ani);
        %     num_exits=numel(find(exits.Animal == ani));
        
        % Specify the distrance from the shelter that triggers the loom.
%         if contains(string(date), "191")
%             thrD = 15;
%             IM_SIZE = 518;
%         elseif contains(string(date), "200")
%             thrD = 15;
%             IM_SIZE = 524;
%         elseif contains(string(date), "220")
%             thrD = 15;
%             IM_SIZE = 416;
%         end
        thrD = 15; 
        IM_SIZE = 416;

        % Set the X, Y, and dist from the shelter for that mouse fromexit_analysis
        xVALS = cell2mat(exit_analysis.X(row_ana));
        xVALS = xVALS/(IM_SIZE/416);
        yVALS = cell2mat(exit_analysis.Y(row_ana));
        yVALS = yVALS/(IM_SIZE/416);
        DVALS = cell2mat(exit_analysis.DShelt(row_ana));
        
        fout = exits.FrameOut(i);
        fin = exits.FrameIn(i);
        
        % Find the frame when the mouse passes the loom threshold
        dloom = DVALS(fout:fin);
        overT = find(dloom>thrD);
        if ~isempty(overT)
            
            lr = overT(1);
            xxWT = smooth(xVALS(fout:fin));
            yyWT = smooth(yVALS(fout:fin));
            nframes= exits.FramesExit(i);
            
            for q = 1:nframes
                x = xxWT(q);
                y = 416 - yyWT(q);
                x2 = xxWT(q+1);
                y2 = 416 - yyWT(q+1);
                if q<lr
                    plot([x, x2],[y,y2],'k', 'LineWidth', 1.4)
                elseif q>=lr
                    plot([x, x2],[y,y2],'Color', [0.7 0.7 0.7], 'LineWidth', 0.75)
                end
                hold on
                plot(xxWT(lr), 416-yyWT(lr), 'Color', [1 0.4471, 0.1255], 'Marker', '.', 'MarkerSize', 20)
            end
        end
        
    end
end

% box off
axis([0 416 0 416])
axis square
ax = gca;
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';
box on
xticks([])
yticks([])

% % % % % % % % % % 

% HETS

figure
hold on 
for i = 1:height(exits)
    
    %Find fixed variables for this exit.
    ani = exits.Animal(i);
    geno = exits.Geno(i);
    date = exits.Date(i);
    
    if geno == 2
        % Find which row in exit_analysis corresponds to this animal
        row_ana = find(exit_analysis.Animal == ani);
        %     num_exits=numel(find(exits.Animal == ani));
        
        % Specify the distrance from the shelter that triggers the loom.
%         if contains(string(date), "191")
%             thrD = 15;
%             IM_SIZE = 518;
%         elseif contains(string(date), "2002")
%             thrD = 15;
%             IM_SIZE = 524;
%         elseif contains(string(date), "220") ||contains(string(date), "2004") 
            thrD = 15;
            IM_SIZE = 416;
%         end
        
        % Set the X, Y, and dist from the shelter for that mouse fromexit_analysis
        xVALS = cell2mat(exit_analysis.X(row_ana));
        xVALS = xVALS/(IM_SIZE/416);
        yVALS = cell2mat(exit_analysis.Y(row_ana));
        yVALS = yVALS/(IM_SIZE/416);
        DVALS = cell2mat(exit_analysis.DShelt(row_ana));
        
        fout = exits.FrameOut(i);
        fin = exits.FrameIn(i);
        
        % Find the frame when the mouse passes the loom threshold
        dloom = DVALS(fout:fin);
        overT = find(dloom>thrD);
        if ~isempty(overT)
            
            lr = overT(1);
            xxWT = smooth(xVALS(fout:fin));
            yyWT = smooth(yVALS(fout:fin));
            nframes= exits.FramesExit(i);
            
            for q = 1:nframes
                x = xxWT(q);
                y = 416 - yyWT(q);
                x2 = xxWT(q+1);
                y2 = 416 - yyWT(q+1);
                if q<lr
                    plot([x, x2],[y,y2], 'Color', [1 0.4471, 0.1255], 'LineWidth', 1.4)
                elseif q>=lr
                    plot([x, x2],[y,y2],'Color', [0.7 0.7 0.7], 'LineWidth', 0.75)
                end
                hold on
                plot(xxWT(lr), 416-yyWT(lr), 'k.', 'MarkerSize', 20)
            end
        end
        
    end
end

% box off
axis([0 416 0 416])
axis square
ax = gca;
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';
box on
xticks([])
yticks([])


%% P1B - Dot for Mean + SEM error vertical bar - with LINE FIT - TRIAL VERSUS REACTION TIME (T2M) 
% Figure 2j

allWT = find(string(all_xy_analysis.Geno) =="wt");
allHET = find(string(all_xy_analysis.Geno) =="het");

% WT 
n_trialsWT = max(all_xy_analysis.Trial(allWT));
% n_trialsWT = 12; 
dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM

for i = 1:n_trialsWT
    vals = find(all_xy_analysis.Trial == i & string(all_xy_analysis.Geno) == "wt"  & cell2mat(all_xy_analysis.ReturnToShelter)==1); % & cell2mat(all_xy_analysis.speedat)>2 & cell2mat(all_xy_analysis.ReturnToShelter)==1);
    data = cell2mat(all_xy_analysis.MaxSpEscape(vals)); 
    dWT(i) = nanmean(data);
    sWT(i) = nanstd(data)/sqrt(numel(vals));
end

% HET 
n_trialsHET = max(all_xy_analysis.Trial(allHET));
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

for i = 1:n_trialsHET
    vals = find(all_xy_analysis.Trial == i & string(all_xy_analysis.Geno) == "het" & cell2mat(all_xy_analysis.ReturnToShelter)==1); % & cell2mat(all_xy_analysis.speedat)>2); %& cell2mat(all_xy_analysis.ReturnToShelter)==1);
    data = cell2mat(all_xy_analysis.MaxSpEscape(vals));
    dHET(i) = nanmean(data);
    sHET(i) = nanstd(data)/sqrt(numel(vals));
end

% dHET = dHET(1:n_trialsWT);
% sHET = sHET(1:n_trialsWT);
% n_trialsHET = n_trialsWT; 

% PLOT 

%Pthcd1 = . [1 0.8 0.05], setd5 = [1 0 0], cul3 = [1, 0 1]
% Ptchd1 - [1 0.7 0.5]

figure
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 1], 'MarkerFaceColor', [1 0.7 1] , 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 1],'Marker', 'none', 'LineWidth', 1.75)

% CREATE FIT LINES - LINEAR FIT - y = ax+b % % % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % % 

[fitresult, gof] = createFits(dWT, dHET, col);
legend off

axis([0 n_trialsWT+1 0 2.25])
box off
% yticks([1,2,3,4,5])
% xlabel('Trial')
xlabel('')
% ylabel('Reaction Time (s)')
ax = gca;
ax.FontSize = 20; 
hold off
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.025 0.025];
f =gcf;
% f.Position = [2143 212 340 273];
f.Position = [2395 174 301 180]; %small
% f.Position = [ 819   552   588   180]; %long

%
[R, P, RL, RU] = corrcoef(1:1:n_trialsWT, dWT)
[R, P, RL, RU] = corrcoef(1:1:5, dHET(1:5))

[rho, p] = corr([1:1:n_trialsWT]', dWT', 'Type', 'Spearman')
[rho, p] = corr([1:1:5]', dHET(1:5)', 'Type', 'Spearman')

%% STATS - ACROSS TRIALS!!

dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

pDATA = zeros(1, n_trialsWT);

for i = 1:n_trialsWT
    valsWT = find(all_xy_analysis.Trial == i & string(all_xy_analysis.Geno) == "wt"); % & cell2mat(all_xy_analysis.speedat)>2 & cell2mat(all_xy_analysis.ReturnToShelter)==1);
    dataWT = cell2mat(all_xy_analysis.MaxSpEscape(valsWT)); 
    dWT(i) = mean(dataWT);
    sWT(i) = std(dataWT)/sqrt(numel(valsWT));

    valsHET = find(all_xy_analysis.Trial == i & string(all_xy_analysis.Geno) == "het"); % & cell2mat(all_xy_analysis.speedat)>2); %& cell2mat(all_xy_analysis.ReturnToShelter)==1);
    dataHET = cell2mat(all_xy_analysis.MaxSpEscape(valsHET));
    dHET(i) = mean(dataHET);
    sHET(i) = std(dataHET)/sqrt(numel(valsHET));
    
    [p, h] = ranksum(dataWT, dataHET);
    pDATA(i) = p;
    
end


%% Dot for Mean + SEM error vertical bar - MAX SPEED across trials -REWARD TRIALS
% Figure 2j

allWT = find(string(ALL_XYLOOM_TABLE.Geno) =="wt");
allHET = find(string(ALL_XYLOOM_TABLE.Geno) =="het");

% WT 
% n_trialsWT = max(all_xy_analysis.Trial(allWT));
n_trialsWT = 12;
dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM

for i = 1:n_trialsWT
    vals = find(ALL_XYLOOM_TABLE.Trial == i & string(ALL_XYLOOM_TABLE.Geno) == "wt"); % & cell2mat(all_xy_analysis.speedat)>2 
    data = cell2mat(ALL_XYLOOM_TABLE.MaxSp(vals)); 
    dWT(i) = mean(data);
    sWT(i) = std(data)/sqrt(numel(vals));
end

% HET 
n_trialsHET = max(ALL_XYLOOM_TABLE.Trial(allHET));
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

for i = 1:n_trialsHET
    vals = find(ALL_XYLOOM_TABLE.Trial == i & string(ALL_XYLOOM_TABLE.Geno) == "het");
    data = cell2mat(ALL_XYLOOM_TABLE.MaxSp(vals));
    dHET(i) = mean(data);
    sHET(i) = std(data)/sqrt(numel(vals));
end

dHET = dHET(1:n_trialsWT);
sHET = sHET(1:n_trialsWT);
n_trialsHET = n_trialsWT; 

% PLOT 

%Pthcd1 = . [1 0.8 0.05], setd5 = [1 0 0], cul3 = [1, 0 1]
% Ptchd1 - [1 0.7 0.5] -  [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] 

figure
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] , 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.2],'Marker', 'none', 'LineWidth', 1.75)

% CREATE FIT LINES - LINEAR FIT - y = ax+b % % % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % % 

% [fitresult, gof] = createFits(dWT, dHET, col);
% legend off

axis([0 n_trialsWT+1 20 100])
box off
% yticks([1,2,3,4,5])
% xlabel('Trial')
xlabel('')
% ylabel('Max. Escape Speed (cms^-1)')
ylabel('')
ax = gca;
ax.FontSize = 20; 
hold off
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.025 0.025];
f =gcf;
% f.Position = [2143 212 340 273];
f.Position = [2395 174 301 180]; %small
% f.Position = [ 819   552   588   180]; %long

%% Dot for Mean + SEM error vertical bar - REACTION TIME across trials - REWARD TRIALS
% Figure 2j

allWT = find(string(ALL_XYLOOM_TABLE.Geno) =="wt");
allHET = find(string(ALL_XYLOOM_TABLE.Geno) =="het");

% WT 
% n_trialsWT = max(all_xy_analysis.Trial(allWT));
n_trialsWT = 12;
dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM

for i = 1:n_trialsWT
    vals = find(ALL_XYLOOM_TABLE.Trial == i & string(ALL_XYLOOM_TABLE.Geno) == "wt"); % & cell2mat(all_xy_analysis.speedat)>2 
    data = cell2mat(ALL_XYLOOM_TABLE.T2M(vals)); 
    dWT(i) = mean(data);
    sWT(i) = std(data)/sqrt(numel(vals));
end

% HET 
n_trialsHET = max(ALL_XYLOOM_TABLE.Trial(allHET));
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

for i = 1:n_trialsHET
    vals = find(ALL_XYLOOM_TABLE.Trial == i & string(ALL_XYLOOM_TABLE.Geno) == "het");
    data = cell2mat(ALL_XYLOOM_TABLE.T2M(vals));
    dHET(i) = mean(data);
    sHET(i) = std(data)/sqrt(numel(vals));
end

dHET = dHET(1:n_trialsWT);
sHET = sHET(1:n_trialsWT);
n_trialsHET = n_trialsWT; 

% PLOT 

%Pthcd1 = . [1 0.8 0.05], setd5 = [1 0 0], cul3 = [1, 0 1]
% Ptchd1 - [1 0.7 0.5] -  [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] 

figure
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.2], 'MarkerFaceColor', [1 0.7 0.5] , 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75) %[1 0.8 0.8]
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.2],'Marker', 'none', 'LineWidth', 1.75)

% CREATE FIT LINES - LINEAR FIT - y = ax+b % % % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % % 

% [fitresult, gof] = createFits(dWT, dHET, col);
% legend off

axis([0 n_trialsWT+1 20 100])
box off
xlabel('')
ylabel('')
ax = gca;
ax.FontSize = 20; 
hold off
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.025 0.025];
f =gcf;

