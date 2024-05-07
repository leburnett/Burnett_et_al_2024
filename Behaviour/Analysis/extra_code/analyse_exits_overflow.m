% analyse_exits_overflow

%% Add trial number to 'exit_analysis' and then make plots looking at adaptation across trials

  all_animals = unique(exit_analysis.Animal); 
  n_animals = numel(all_animals); 

for i = 1:n_animals
    ani = all_animals(i); 
    trial_number = 1; 
    
    for j = 1:height(exit_analysis)
        if exit_analysis.Animal(j) == ani 
            exit_analysis.Trial(j) = trial_number; 
            trial_number = trial_number +1; 
        end   
    end
end 

%% P10  - MEAN and SEM - Error line - across TRIALS


allWT = find((exit_analysis.geno) ==1);
allHET = find((exit_analysis.geno) ==2);

% WT 
n_trialsWT = max(exit_analysis.Trial(allWT)); 
dWT = zeros(1, n_trialsWT); %data
sWT = zeros(1, n_trialsWT); % SEM

for i = 1:n_trialsWT
    vals = find(exit_analysis.Trial == i & (exit_analysis.geno) == 1);
    data = (exit_analysis.NumOut(vals)); 
    dWT(i) = mean(data);
    sWT(i) = std(data)/sqrt(numel(vals));
end

% HET 
n_trialsHET = max(exit_analysis.Trial(allHET)); 
dHET = zeros(1, n_trialsHET); %data
sHET = zeros(1, n_trialsHET); % SEM

for i = 1:n_trialsHET
    vals = find(exit_analysis.Trial == i & (exit_analysis.geno) == 2);
    data = (exit_analysis.NumOut(vals));
    dHET(i) = mean(data);
    sHET(i) = std(data)/sqrt(numel(vals));
end


% PLOT 

figure
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor', 'none',  'MarkerSize', 10, 'LineWidth', 1.75)
hold on 
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0 0], 'MarkerFaceColor', [1 0.8 0.8], 'MarkerEdgeColor','none', 'MarkerSize', 10, 'LineWidth', 1.75)
errorbar(1:1:n_trialsWT, dWT, sWT, 'o', 'CapSize', 0, 'Color', [0.4 0.4 0.4], 'Marker', 'none', 'LineWidth', 1.75)
errorbar(1:1:n_trialsHET, dHET, sHET, 'o', 'CapSize', 0, 'Color', [1 0.4 0.4],'Marker', 'none', 'LineWidth', 1.75)

% CREATE FIT LINES - LINEAR FIT - y = ax+b % % % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % %  % % % % % % % % 
% % 
% [fitresult, gof] = createFits(dWT, dHET);
% legend off

axis([0 16 0 5])
box off
% yticks([1,2,3,4,5])
xticks([1:1:15])
xlabel('Experiment')
ylabel('Shelter Exits')
ax = gca;
ax.FontSize = 18; 
hold off
ax.TickDir = 'out'; 
ax.LineWidth = 1.2;
ax.TickLength = [0.02 0.02];
f =gcf;
f.Position = [2143 217 287  268]; %[2143 212 340 273];
% NARROW = [2143 217 193 268];


%% 4 - PLOT - Distance from the shelter at all times during trial -  WT/ HET siblings. 

n = 15000; %8000;
n1 = 30;
n2 = 32;
dist_loom_trig = 7.8; % 9.2;%9.2 7.8 - this will change for the different cohorts. 

d = exit_analysis.DShelt{n1};
figure; plot(smooth(d), 'k', 'LineWidth', 1.2)
axis([0 n -10 30])
hold on
plot([0 xl], [0 0], 'k', 'LineWidth', 2)
plot([0 n], [dist_loom_trig dist_loom_trig], 'k:', 'LineWidth', 1)

d = exit_analysis.DShelt{n2};
plot(smooth(d), 'r', 'LineWidth', 1.2)
box off
ax = gca;
ax.TickDir = 'out';
ax.FontSize = 18;
xlabel('Time')
xticks([])
ylabel('Distance from Shelter (cm)')
ax.LineWidth = 1.2;
f = gcf;
f.Position = [680   821   528   277];


%% 2 - MEAN SEM PLOT - Similar to boxplot but with dot and errorbar vetical line - all trials!

% Set 'v' as the column number of the variable you are interested inanalysing. 
v = 8;

valsWT = [];
valsHET = [];
for i = 1:n
    if exit_analysis.geno(i) == 1 
        val = exit_analysis{i, v};
        valsWT = [valsWT, val];
    elseif exit_analysis.geno(i) == 2 
        val = exit_analysis{i, v};
        valsHET = [valsHET, val];
    end
end

n_wt = numel(valsWT);
n_het = numel(valsHET);
x1 = ones(1, n_wt);
x2 = ones(1, n_het)*2;

% PLOT
x = [1,2];
y = [nanmean(valsWT), nanmean(valsHET)];

err(1) =  nanstd(valsWT)/sqrt(numel(valsWT));
err(2) =  nanstd(valsHET)/sqrt(numel(valsHET));

colvals = {[0 0 0]; [1 0 0]};

figure
hold on
for i = 1:2
    scatter(x(i), y(i), 250, 'MarkerFaceColor', colvals{i}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')
    hold on
    errorbar(x(i), y(i), err(i),'Color',  colvals{i}, 'CapSize', 0, 'LineWidth', 1.1)
end
% 
% hold on
% scatter(x1, valsWT','SizeData', 120, 'MarkerEdgeColor', [0.6 0.6 0.6], 'jitter', 'on', 'jitterAmount', 0.2, 'MarkerEdgeAlpha', 0.5)
% scatter(x2, valsHET', 'SizeData', 120,'MarkerEdgeColor', [1 0 0], 'jitter', 'on', 'jitterAmount', 0.2,  'MarkerEdgeAlpha', 0.5)
ax = gca;
ax.FontSize = 25;
box off
axis([0.5 2.5 -0.1 1.1])
% ax.XAxis.Visible = 'on';
% xticks([1, 2])
% xticklabels(["Setd5^+/+", "Setd5^+/-"])
ax.XAxis.Visible = 'off';
xtickangle(45)
ax.TickDir = 'out';
ax.LineWidth = 2;
ylabel('Probability of Getting Reward')

f = gcf;
f.Position = [680   844   250   400];



%% TRAJECTORIES during Loom experiments

% From exit analysis/exits for DAY1 Loom Exps
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
    exp = exits.Exp(i);
    
    if exp == "01_Loom"
        
        if geno == 1
            % Find which row in exit_analysis corresponds to this animal
            row_ana = find(exit_analysis.Animal == ani & exit_analysis.Exp == "01_Loom");
            loom_r = cell2mat(exit_analysis.LoomRows(row_ana(1)));
            loom_r = loom_r(1);
            
            %     num_exits=numel(find(exits.Animal == ani));
            
            % Specify the distrance from the shelter that triggers the loom.
            if contains(string(date), "191")
                thrD = 15;
                IM_SIZE = 518;
                c = 1;
            elseif contains(string(date), "200")
                thrD = 15;
                IM_SIZE = 524;
                c = 2;
            elseif contains(string(date), "210") 
                thrD = 15;
                IM_SIZE = 416;
                c = 3;
            end
            
            % Set the X, Y, and dist from the shelter for that mouse fromexit_analysis
            xVALS = cell2mat(exit_analysis.X(row_ana));
            xVALS = xVALS/(IM_SIZE/416);
            yVALS = cell2mat(exit_analysis.Y(row_ana));
            yVALS = yVALS/(IM_SIZE/416);
            DVALS = cell2mat(exit_analysis.DShelt(row_ana));
            
            fout = exits.FrameOut(i);
            fin = exits.FrameIn(i);
            
            if loom_r > fout && loom_r <fin
                % Find the frame when the mouse passes the loom threshold
%                 dloom = DVALS(fout:fin);
%                 overT = find(dloom>thrD);
%                 if ~isempty(overT)
                    
                    lr = loom_r-fout; %overT(1);
                    xxWT = smooth(xVALS(fout:fin));
                    yyWT = smooth(yVALS(fout:fin));
                    nframes= exits.FramesExit(i);
                    
                    for q = 1:nframes
                        x = xxWT(q);
                        y = 416 - yyWT(q);
                        x2 = xxWT(q+1);
                        y2 = 416 - yyWT(q+1);
                        if q>lr
                            plot([x, x2],[y,y2],'k', 'LineWidth', 1)
                        elseif q<=lr
%                             plot([x, x2],[y,y2],'Color', [0.7 0.7 0.7], 'LineWidth', 0.75)
                        end
                        hold on
                        plot(xxWT(lr), 416-yyWT(lr), 'r.', 'MarkerSize', 20)
                    end
%                 end
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
    exp = exits.Exp(i);
    
    if exp == "03_Loom" 
        
        if geno == 2
            % Find which row in exit_analysis corresponds to this animal
            row_ana = find(exit_analysis.Animal == ani & exit_analysis.Exp == "03_Loom");
            loom_r = cell2mat(exit_analysis.LoomRows(row_ana(1)));
            loom_r = loom_r(1);
            
            %     num_exits=numel(find(exits.Animal == ani));
            
            % Specify the distrance from the shelter that triggers the loom.
            if contains(string(date), "191")
                thrD = 15;
                IM_SIZE = 518;
                c = 1;
            elseif contains(string(date), "200")
                thrD = 15;
                IM_SIZE = 524;
                c = 2;
            elseif contains(string(date), "210")
                thrD = 15;
                IM_SIZE = 416;
                c = 3;
            end
            
            % Set the X, Y, and dist from the shelter for that mouse fromexit_analysis
            xVALS = cell2mat(exit_analysis.X(row_ana));
            xVALS = xVALS/(IM_SIZE/416);
            yVALS = cell2mat(exit_analysis.Y(row_ana));
            yVALS = yVALS/(IM_SIZE/416);
            DVALS = cell2mat(exit_analysis.DShelt(row_ana));
            
            fout = exits.FrameOut(i);
            fin = exits.FrameIn(i);
            
            if loom_r > fout && loom_r <fin
                % Find the frame when the mouse passes the loom threshold
%                 dloom = DVALS(fout:fin);
%                 overT = find(dloom>thrD);
%                 if ~isempty(overT)
                    
                    lr = loom_r-fout; %overT(1);
                    xxWT = smooth(xVALS(fout:fin));
                    yyWT = smooth(yVALS(fout:fin));
                    nframes= exits.FramesExit(i);
                    
                    for q = 1:nframes
                        x = xxWT(q);
                        y = 416 - yyWT(q);
                        x2 = xxWT(q+1);
                        y2 = 416 - yyWT(q+1);
                        if q>lr
                            plot([x, x2],[y,y2],'r', 'LineWidth', 1)
                        elseif q<=lr
%                             plot([x, x2],[y,y2],'Color', [0.7 0.7 0.7], 'LineWidth', 0.75)
                        end
                        hold on
                        plot(xxWT(lr), 416-yyWT(lr), 'k.', 'MarkerSize', 20)
                    end
%                 end
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

[p, h] = ttest(valsWT, valsHET)