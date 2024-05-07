%%  Analyse Laser Power Trials from Optogenetics experiments
% To be used after 'analyse_OPTO.m'

% Created 16/12/21 - Burnett

%% Initialise variables

% Convert xy_optoA back into xy_opto for plots:
xy_opto = xy_optoA;
xy_opto_info = xy_opto_infoA;
xy_analysis = xy_analysisA; 

n = numel(xy_opto(:,1));
col = 'r';
frames_b4 = 595;
total_lighton_frames = 60;  % light was on for 1s. - 60 frames.

all_animals = unique(string(xy_opto_info.Animal));
n_animals = numel(all_animals);



%% P1 - Individual plot for each animal - LINE PLOT - mean speed response of each laser power. 

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(string(xy_opto_info.Animal) == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.EC)); 
    n_EC = numel(all_EC);
    figure
    total_frames_light_on = 60;
    interpulse_frames = total_frames_light_on/10;
    
    rectangle('Position', [frames_b4 0 total_frames_light_on 65], 'FaceColor', [0.67 0.84 0.9, 0.2], 'EdgeColor', [0.67 0.84 0.9, 0.2]) %xy wh
    hold on
    rectangle('Position', [frames_b4 65 total_frames_light_on 5], 'FaceColor', [0.67 0.84 0.9], 'EdgeColor', [0.67 0.84 0.9]) %xy wh
    
%     added = 0;
%     lw = 0.5; 
%     for k = 1:10
%         plot([frames_b4+added frames_b4+added], [55 65],'LineStyle', '-', 'Color', [0.37 0.54 0.6], 'LineWidth', lw)
%         added = added+interpulse_frames;
%         lw = lw+0.1; 
%     end
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Animal) == ani);
        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        
    end 
    title(ani)
    axis([500 800 0 70])
    box off
    ylabel('Speed - cm/s')
    xticks([frames_b4-60,frames_b4, frames_b4+60, frames_b4+120, frames_b4+180])
        xticklabels({'-1', '0', '1', '2', '3'})
        ax = gca;
        ax.FontSize = 20;
        ax.TickDir = 'out';
        ax.LineWidth = 1;
end 


%% P2 - Plot for each GENOTYPE - LINE PLOT - mean speed response of each laser power.

% 1 - WT TRIALS % % % % % % % % % % % % % % % % % % % % % % % % % 

figure
total_frames_light_on = 60;
interpulse_frames = total_frames_light_on/10;

rectangle('Position', [frames_b4 0 total_frames_light_on 55], 'FaceColor', [0.67 0.84 0.9, 0.2], 'EdgeColor', [0.67 0.84 0.9, 0.2]) %xy wh
hold on
rectangle('Position', [frames_b4 55 total_frames_light_on 5], 'FaceColor', [0.67 0.84 0.9], 'EdgeColor', [0.67 0.84 0.9]) %xy wh
added = 0;

for k = 1:10
    plot([frames_b4+added frames_b4+added], [55 60],'LineStyle', '-', 'Color', [0.37 0.54 0.6], 'LineWidth', 1)
    added = added+interpulse_frames;
end

all_ani = find(string(xy_opto_info.Geno) == "wt");
n_trials_ani = numel(all_ani);

all_EC = unique(cell2mat(xy_opto_info.EC));
n_EC = numel(all_EC);

for j = [1,2,3 5] % 8, 10, 11, 12] % [1,2,3,4,5,8,10,11,12,13]  %,8,11,13,14]
    EC_str = all_EC(j,1);
    all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "wt");
    col = 1 -[j/n_EC j/n_EC j/n_EC];
    
    if ~isempty(all_ani_EC) && numel(all_ani_EC)>1
        av_resp = mean(xy_opto(all_ani_EC, :));
        %                   av_resp(611:620) = linspace(av_resp(611), av_resp(621), 10);
        %                     av_resp(556:565) = linspace(av_resp(556), av_resp(566), 10);
        plot(av_resp, 'Color', col, 'LineWidth', 1.4);
        hold on
    elseif ~isempty(all_ani_EC) && numel(all_ani_EC)==1
        av_resp = (xy_opto(all_ani_EC, :));
        plot(av_resp, 'Color', col, 'LineWidth', 1.4);
        hold on
    end
    %             xticklabels({''})
end
axis([500 800 0 60])
box off
ylabel('Speed - cm/s')

xticks([frames_b4-60,frames_b4, frames_b4+60, frames_b4+120, frames_b4+180])
xticklabels({'-1', '0', '1', '2', '3'})
ax = gca;
ax.FontSize = 20;
ax.TickDir = 'out';
ax.LineWidth = 1;

% HET TRIALS % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

figure
total_frames_light_on = 60;
interpulse_frames = total_frames_light_on/10;

rectangle('Position', [frames_b4 0 total_frames_light_on 65], 'FaceColor', [0.67 0.84 0.9, 0.2], 'EdgeColor', [0.67 0.84 0.9, 0.2]) %xy wh
hold on
rectangle('Position', [frames_b4 55 total_frames_light_on 5], 'FaceColor', [0.67 0.84 0.9], 'EdgeColor', [0.67 0.84 0.9]) %xy wh

added = 0;
for k = 1:10
    plot([frames_b4+added frames_b4+added], [55 60],'LineStyle', '-', 'Color', [0.37 0.54 0.6], 'LineWidth', 1)
    added = added+interpulse_frames;
end

all_ani = find(string(xy_opto_info.Geno) == "het");
n_trials_ani = numel(all_ani);

all_EC = unique(cell2mat(xy_opto_info.EC));
n_EC = numel(all_EC);
for j = [1,3,4,5,8,10,11,12]
    EC_str = all_EC(j,1);
    all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "het");
    
    col = 1 -[j/n_EC j/n_EC j/n_EC];
    
    if ~isempty(all_ani_EC)&& numel(all_ani_EC)>1
        av_resp = mean(xy_opto(all_ani_EC, :));
        plot(av_resp, 'Color', col, 'LineWidth', 1.4);
        hold on
    elseif ~isempty(all_ani_EC) && numel(all_ani_EC)==1
        av_resp = (xy_opto(all_ani_EC, :));
        plot(av_resp, 'Color', col, 'LineWidth', 1.4);
        hold on
    end
end

axis([500 800 0 60])
box off
ylabel('Speed - cm/s')
xticks([frames_b4-60,frames_b4, frames_b4+60, frames_b4+120, frames_b4+180])
xticklabels({'-1', '0', '1', '2', '3'})

ax = gca;
ax.FontSize = 20;
ax.TickDir = 'out';
ax.LineWidth = 1;



%% Individual frequencies. - GOOD PLOTS!!!! 

laser_val = 14; 
col = 'r';

figure
total_frames_light_on = 60; 
interpulse_frames = total_frames_light_on/10; 
rectangle('Position', [frames_b4 0 total_frames_light_on 65], 'FaceColor', [0.67 0.84 0.9, 0.2], 'EdgeColor', [0.67 0.84 0.9, 0.2]) %xy wh
hold on
rectangle('Position', [frames_b4 55 total_frames_light_on 5], 'FaceColor', [0.67 0.84 0.9], 'EdgeColor', [0.67 0.84 0.9]) %xy wh 

added = 0; 
    for k = 1:10
        plot([frames_b4+added frames_b4+added], [55 65],'LineStyle', '-', 'Color', [0.37 0.54 0.6], 'LineWidth', 1)
        added = added+interpulse_frames; 
    end 
    

for i = 1:2
    if i ==1
        all_ani = find(string(xy_opto_info.Geno) == "wt");
        n_trials_ani = numel(all_ani);
        
        all_EC = unique(cell2mat(xy_opto_info.EC));
        n_EC = numel(all_EC);
        
        for j = laser_val %1:n_EC
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "wt");
            
            speed_WT = (xy_opto(all_ani_EC, :));

            nWT = numel(speed_WT(:,1));
            mean_WT = mean(speed_WT);
            
            x = (1:1:1200);
            
            semWT = std(speed_WT)/sqrt(nWT);
            y1 = mean_WT+semWT;
            y2 = mean_WT-semWT;

            plot(x, y1, 'w')
            hold on
            plot(x, y2, 'w')
            patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
            plot(mean_WT', 'k', 'LineWidth', 1.3)

        end
%         axis([500 800 0 60])
        ylim([0 60])
%         xlim([500 950])
        xlim([500 800])
        box off
        ylabel('Speed - cm/s')
        xticks([frames_b4-60,frames_b4, frames_b4+60, frames_b4+120, frames_b4+180, frames_b4+240, frames_b4+300])
        xticklabels({'-1', '0', '1', '2', '3', '4', '5'})
        ax = gca;
        ax.FontSize = 20;
        ax.TickDir = 'out'; 
        ax.LineWidth = 1; 
        
    elseif i ==2
        all_ani = find(string(xy_opto_info.Geno) == "het");
        n_trials_ani = numel(all_ani);
        
        all_EC = unique(cell2mat(xy_opto_info.EC));
        n_EC = numel(all_EC);
        
        for j = laser_val%1:n_EC
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "het");
            speed_HET = (xy_opto(all_ani_EC, :));
            
            nHET = numel(speed_HET(:,1));
            mean_HET = mean(speed_HET);
            
            x = (1:1:1200);

            semHET = std(speed_HET)/sqrt(nHET);
            y3 = mean_HET+semHET;
            y4 = mean_HET-semHET;

            plot(x, y3, 'w')
            hold on
            plot(x, y4, 'w')
            patch([x fliplr(x)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
            plot(mean_HET', 'Color', col, 'LineWidth', 1.3)
        end
    end
end

%%




%% Psychometric curve plots. 


% Branco paper - each point has an errorbar (s.d. not sem)
% escape probability = logistic fit. 
% reaction time and escape speed = linear fit. 

% Use xy_analysis for different variables. 

% errorbar(x,y,err, '.')

n_xy_ana= height(xy_analysis);

allWT = find(string(xy_analysis.Geno) == "wt");
allHET = find(string(xy_analysis.Geno) == "het");
ECs = unique(cell2mat(xy_analysis.EC));
n_EC = numel(ECs);
column = 18; 

%

xvals = ECs; 
yvalsWT = zeros(n_EC,1);
errWT = zeros(n_EC,1);
yvalsHET = zeros(n_EC,1);
errHET = zeros(n_EC,1);

%

for i = 1:n_xy_ana
    
    for j = 1:n_EC
        EC_val = ECs(j);
        allwt = find(cell2mat(xy_analysis.EC(allWT)) == EC_val);
        allhet = find(cell2mat(xy_analysis.EC(allHET)) == EC_val);
        
        valsWT = cell2mat(table2array(xy_analysis(allwt,column)));
        valsHET = cell2mat(table2array(xy_analysis(allhet,column)));
        
        yvalsWT(j) = nanmean(valsWT);
        yvalsHET(j) = nanmean(valsHET);
        
        errWT(j) = nanstd(valsWT)/sqrt(numel(valsWT));
        errHET(j) = nanstd(valsHET)/sqrt(numel(valsHET));
    end 

end 

%%
figure
errorbar(xvals, yvalsWT, errWT, 'k.', 'MarkerSize', 20)
hold on
errorbar(xvals, yvalsHET, errHET, 'r.', 'MarkerSize', 20)
box off













































%% Trouble Shooting

% T1 - Removing animal.   
% all_ani = find(xy_opto_info.Animal == "MJ1415");
% 
% xy_opto(all_ani, :) = []; 
% xy_opto_info(all_ani, :) = [];
% xy_analysis(all_ani, :) = []; 








