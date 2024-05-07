%% Analyse OPTOGENETICS experiments in PCF arena. 
% Made for Setd5/VGlut2 OPTO experiments from May 2021 - for Burnett et al.

% Created by Burnett - Nov 2021

%% 1 - Initiailise variables
% Open xy_opto and xy_opto_info

% Number of trials. 
n = numel(xy_opto_info(:,1)); 
col = 'r'; 
frames_b4 = 595;
total_lighton_frames = 60;  % light was on for 1s. - 60 frames. 

all_animals = unique(string(xy_opto_info.Animal)); 
n_animals = numel(all_animals); 

%% 2 - Add trial number - cumulative 

for i = 1:n_animals
    ani = all_animals(i); 
    trial_number = 1; 
    
    for j = 1:n
        if string(xy_opto_info.Animal{j}) == ani 
            xy_opto_info.Trial(j) = trial_number; 
            trial_number = trial_number +1; 
        end   
    end
end 

%% 3 - Add 'DAY' - IRRELEVANT OF DATE. 
for i = 1:n 
    if string(xy_opto_info.Date(i)) == "210521" || string(xy_opto_info.Date(i)) == "210830" || string(xy_opto_info.Date(i)) == "210906" || string(xy_opto_info.Date(i)) == "210913" 
        xy_opto_info.Day(i) = 1;
    elseif string(xy_opto_info.Date(i)) == "210522" || string(xy_opto_info.Date(i)) == "210831" || string(xy_opto_info.Date(i)) == "210907"
        xy_opto_info.Day(i) = 2;
    elseif string(xy_opto_info.Date(i)) == "210523" || string(xy_opto_info.Date(i)) == "210901" || string(xy_opto_info.Date(i)) == "210908"
        xy_opto_info.Day(i) = 3;
    elseif string(xy_opto_info.Date(i)) == "210524" || string(xy_opto_info.Date(i)) == "210916" 
        xy_opto_info.Day(i) = 4;
    elseif string(xy_opto_info.Date(i)) == "210525"
        xy_opto_info.Day(i) = 5;
    end 
end 

%% 4 - Add 'TrialPerDay' 

for i = 1:n_animals
    ani = all_animals(i); 
    
    for k = 1:5 % Run through days
        day = k;
        trial_number = 1; 
            
        for j = 1:n
            if xy_opto_info.Day(j)== day && string(xy_opto_info.Animal{j}) == ani
                    xy_opto_info.TrialPerDay(j) = trial_number;
                    trial_number = trial_number +1;
            end
        end
        
    end
end 

save('XY_OPTO_ARRAY_INFO.mat', 'xy_opto_info');


%% Combining multiple files

% Analyse PULSE WIDTH

b1 = xy_analysisC;
b2 = xy_opto_infoC;
b3 = xy_optoC;

xy_analysisC = vertcat(a1,b1, xy_analysisC);
xy_opto_infoC = vertcat(a2,b2, xy_opto_infoC);
xy_optoC = vertcat(a3,b3,xy_optoC);


%% Sort xy_opto based on EC.

% LASER POWER TRIALS
% xy_optoA = xy_opto(1:198, :);
% xy_opto_infoA = xy_opto_info(1:198, :); 
% xy_analysisA = xy_analysis(1:198, :); 

% xy_optoA = xy_opto(1:204, :);
% xy_opto_infoA = xy_opto_info(1:204, :); 

% xy_optoA = xy_opto(1:154, :);
% xy_opto_infoA = xy_opto_info(1:154, :); 
%  xy_analysisA = xy_analysis(1:154, :);  

% xy_optoA = xy_opto([1:287,767:end], :);
% xy_opto_infoA = xy_opto_info([1:287, 767:end], :); 

% xy_optoA = xy_opto(154:end, :);
% xy_opto_infoA = xy_opto_info(154:end, :); 
% xy_analysisA = xy_analysis(154:end, :);

% xy_optoA = xy_opto([1:164, 219:336, 385:457], :);
% xy_opto_infoA = xy_opto_info([1:164, 219:336, 385:457], :); 
% xy_analysisA = xy_analysis([1:164, 219:336, 385:457], :);

xy_optoA = xy_opto([1:195], :);
xy_opto_infoA = xy_opto_info([1:195], :); 
xy_analysisA = xy_analysis([1:195], :);

% SORT
xy_optoA(:,1201) = cell2mat(xy_opto_infoA.EC);
xy_opto_infoA = sortrows(xy_opto_infoA, 9); % 9 for EC
xy_optoA = sortrows(xy_optoA, 1201); 
% Remove last column 
xy_optoA = xy_optoA(:, 1:1200); 

allWT = find(string(xy_opto_infoA.Geno) == "wt");
allHET = find(string(xy_opto_infoA.Geno) == "het");

figure
subplot(1,2,1)
ax1 = imagesc(xy_optoA(allWT, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('WT')
subplot(1,2,2)
ax2 =imagesc(xy_optoA(allHET, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('HET')

% figure
% subplot(1,2,1)
% ax1 = plot(mean((xy_optoA(allWT, :))));
% subplot(1,2,2)
% ax2 =plot(mean((xy_optoA(allHET, :))));

% xy_opto(557:843, :) = [];
% xy_opto_info(557:843, :) = []; 
% xy_analysis(557:843, :) = []; 
% 
% xy_opto = vertcat(xy_opto, xy_optoA);
% xy_opto_info = vertcat(xy_opto_info, xy_opto_infoA);
% xy_analysis = vertcat(xy_analysis, xy_analysisA);
% 
% save('LaserPower_ALLCs.mat', 'xy_opto', 'xy_opto_info', 'xy_analysis'); 

% % % % % % % %

% PULSE WIDTH TRIALS
% xy_optoB = xy_opto(199:315, :);
% xy_opto_infoB = xy_opto_info(199:315, :); 
% xy_analysisB = xy_analysis(199:315, :);

% xy_optoB = xy_opto(205:446, :);
% xy_opto_infoB = xy_opto_info(205:446, :); 

% xy_optoB = xy_opto(155:338, :);
% xy_opto_infoB = xy_opto_info(155:338, :); 

% xy_optoB = xy_opto(289:533, :);
% xy_opto_infoB = xy_opto_info(289:533, :); 

% xy_optoB = xy_opto([165:218, 337:384, 458:537], :);
% xy_opto_infoB = xy_opto_info([165:218, 337:384, 458:537], :); 
% xy_analysisB = xy_analysis([165:218, 337:384, 458:537], :);

% xy_optoB = xy_opto(155:338, :);
% xy_opto_infoB = xy_opto_info(155:338, :); 
% xy_analysisB = xy_analysis(155:338, :); 

xy_optoB = xy_opto(196:417, :);
xy_opto_infoB = xy_opto_info(196:417, :); 
xy_analysisB = xy_analysis(196:417, :); 

% SORT
xy_optoB(:,1201) = cell2mat(xy_opto_infoB.T_pulse);
xy_opto_infoB = sortrows(xy_opto_infoB, 7); % 9 for EC
xy_optoB = sortrows(xy_optoB, 1201); 
% Remove last column 
xy_optoB = xy_optoB(:, 1:1200); 

allWT = find(string(xy_opto_infoB.Geno) == "wt");
allHET = find(string(xy_opto_infoB.Geno) == "het");

figure
subplot(1,2,1)
ax1 = imagesc(xy_optoB(allWT, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('WT')
subplot(1,2,2)
ax2 =imagesc(xy_optoB(allHET, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('HET')

% xy_optoB = xy_opto;
% xy_opto_infoB = xy_opto_info; 
% xy_analysisB = xy_analysis;


% % % % % % % % 

% FREQ TRIALS
% xy_optoC = xy_opto(316:end, :);
% xy_opto_infoC = xy_opto_info(316:end, :); 
% xy_analysisC = xy_analysis(316:end, :); 

% xy_optoC = xy_opto(447:end, :);
% xy_opto_infoC = xy_opto_info(447:end, :); 

% xy_optoC = xy_opto(339:end, :);
% xy_opto_infoC = xy_opto_info(339:end, :); 

% xy_optoC = xy_opto(534:766, :);
% xy_opto_infoC = xy_opto_info(534:766, :); 

% xy_optoC = xy_opto(339:end, :);
% xy_opto_infoC = xy_opto_info(339:end, :); 
% xy_analysisC = xy_analysis(339:end, :); 

% xy_optoC = xy_opto(538:end, :);
% xy_opto_infoC = xy_opto_info(538:end, :); 
% xy_analysisC = xy_analysis(538:end, :); 

xy_optoC = xy_opto(418:end, :);
xy_opto_infoC = xy_opto_info(418:end, :); 
xy_analysisC = xy_analysis(418:end, :); 


% SORT
xy_optoC(:,1201) = cell2mat(xy_opto_infoC.FreqPulse);
xy_opto_infoC = sortrows(xy_opto_infoC, 7); % 9 for EC
xy_optoC = sortrows(xy_optoC, 1201); 
% Remove last column 
xy_optoC = xy_optoC(:, 1:1200); 

allWT = find(string(xy_opto_infoC.Geno) == "wt");
allHET = find(string(xy_opto_infoC.Geno) == "het");

figure
subplot(1,2,1)
ax1 = imagesc(xy_optoC(allWT, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('WT')
subplot(1,2,2)
ax2 =imagesc(xy_optoC(allHET, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('HET')


save('LaserPower_C2_N5.mat', 'xy_analysisA', 'xy_optoA', 'xy_opto_infoA'); 
save('PulseWidth_C2_N5.mat', 'xy_analysisB', 'xy_optoB', 'xy_opto_infoB'); 
save('Frequency_C2_N5.mat', 'xy_analysisC', 'xy_optoC', 'xy_opto_infoC'); 

%%

c1 = xy_optoA; 
c2 = xy_opto_infoA; 
c3 = xy_analysisA; 

xy_opto = vertcat(a1,b1,c1, xy_optoA); 
xy_opto_info = vertcat(a2,b2,c2, xy_opto_infoA);
xy_analysis = vertcat(a3,b3,c3,xy_analysisA); 

%%



% Add EC column to xy_opto; 
xy_opto(:,1201) = cell2mat(xy_opto_info.T_pulse);

% Sort based on EC. 
xy_opto_info = sortrows(xy_opto_info, 7); % 9 for EC
xy_opto = sortrows(xy_opto, 1201); 
% Remove last column 
xy_opto = xy_opto(:, 1:1200); 

figure
subplot(1,2,1)
ax1 = imagesc(xy_opto(allWT, :));
subplot(1,2,2)
ax2 =imagesc(xy_opto(allHET, :));
% linkaxes([ax1, ax2], 'x');

allWT = find(string(xy_opto_info.Geno) == "wt");
allHET = find(string(xy_opto_info.Geno) == "het");

%%




%% Make heatmap of all trials - ordered by EC for each animal. 

% FOR EACH ANIMAL MAKE SEPARATE PLOT. 

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_info.Animal == ani); 
    max_EC = find(max(xy_opto_info.EC(all_ani)))
    n_trials_ani = numel(all_ani);
    
    subplot(1,4,i); imagesc(xy_opto(all_ani, :)); title(ani); caxis([0 70]); 
    hold on 
    plot([frames_b4 frames_b4], [0 n_trials_ani+10], 'w:', 'LineWidth', 1.2)
    plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 n_trials_ani+10], 'w:', 'LineWidth', 1.2)
    xticklabels({''})
end 


%% Line plot for each animal - where trials of same EC are averaged. 

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_info.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.T_pulse)); 
    n_EC = numel(all_EC); 
    subplot(n_animals,1, i)
    
    for j = n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.T_pulse) == EC_str & xy_opto_info.Animal == ani);
        
        if EC_str == 0 || EC_str == 1
            v2 = 0.85; 
            col = [v2 v2 v2]; 
        elseif EC_str > 0
%             v2 = 1-(j/n_EC);
            v2 = 1-(EC_str/20); 
            col = [v2 v2 v2]; 
        end 

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 60], 'k:', 'LineWidth', 1.2)
%         plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 60], 'k:', 'LineWidth', 1.2)
        xticklabels({''})
    end 
    title(ani)
    axis([400 900 0 50])
    box off
%     ylabel('Speed - cm/s')
%     xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
%     xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
%     ax = gca;
%     ax.FontSize = 14; 
end 


%%
WT_RESP = [];
HET_RESP =[];

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_info.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.EC)); 
    n_EC = numel(all_EC); 
    EC_Str = 3.25; 
    geno = xy_opto_info.Geno{all_ani(1)};
    
    for j = n_EC
%         EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_Str & xy_opto_info.Animal == ani);

        if ~isempty(all_ani_EC) 
            if numel(all_ani_EC)>1
            av_resp = mean(xy_opto(all_ani_EC, :));
            else
                av_resp = xy_opto(all_ani_EC, :);
            end 
        end 
        
        if string(geno) == "het"
        HET_RESP = vertcat(HET_RESP, av_resp);
        elseif string(geno) == "wt"
        WT_RESP = vertcat(WT_RESP, av_resp);
        end 
       
    end 
end 
        figure
        plot(mean(WT_RESP), 'k',  'LineWidth', 1.2)
        hold on 
        plot(mean(HET_RESP), 'r',  'LineWidth', 1.2)
        plot([frames_b4 frames_b4], [0 60], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+65 frames_b4+65], [0 60], 'k:', 'LineWidth', 1.2)
        xticklabels({''})
        axis([500 800 0 40])
        box off
        ax = gca; 
        ax.FontSize = 16; 
        title(string(EC_Str))
        

%% Individual plots for each EC for each animal to combine post to make step plots. 


for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_info.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.T_pulse)); 
    n_EC = numel(all_EC); 

    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.T_pulse) == EC_str & xy_opto_info.Animal == ani);
        
        figure
        
        if EC_str == 0 || EC_str == 1
            v2 = 0.85; 
            col = [v2 v2 v2]; 
        elseif EC_str > 0
%             v2 = 1-(j/n_EC);
            v2 = 1-(EC_str/20); 
            col = [v2 v2 v2]; 
        end 

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            box off
            axis off
            axis([575 725 -10 75])
            hold on
            plot([600 660], [-5 -5], 'k', 'LineWidth', 3);
            title(EC_str)
        end 
%         plot([frames_b4 frames_b4], [0 100], 'k:', 'LineWidth', 1.2)
%         plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 100], 'k:', 'LineWidth', 1.2)
%         xticklabels({''})
    end 
%     title(ani)
%     axis([500 800 0 75])
%     box off
%     ylabel('Speed - cm/s')
%     xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
%     xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
%     ax = gca;
%     ax.FontSize = 14; 
end 


%% Plot with each EC in a subplot with the average trace for each animal plotted ontop of each other. 

all_EC = unique(cell2mat(xy_opto_info.T_pulse));
n_EC = numel(all_EC);
 figure       
for j = 2:n_EC-3  % For each EC value. 
    EC_str = all_EC(j,1);
    subplot(1,13,j-1)
    
    % Colour 
%     if EC_str == 0 || EC_str == 1
%         v2 = 0.85;
% %         col = [v2 v2 v2];
%     elseif EC_str > 0
%         %             v2 = 1-(j/n_EC);
%         v2 = 1-(EC_str/5);
% %         col = [v2 v2 v2];
%     end
    
    for i = 1:n_animals
        ani = string(all_animals{i});
        all_ani_EC = find(cell2mat(xy_opto_info.T_pulse) == EC_str & xy_opto_info.Animal == ani);
        
        if ani == "MJ0580" || ani == "MJ0582" % If HET make RED not black.
            col = [1 0 0];
        else 
            col = [0 0 0];
        end
        
        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            box off
            axis off
            axis([575 850 -10 75])
            hold on
            plot([600 660], [-5 -5], 'k', 'LineWidth', 3);
            title(EC_str)
        end
        
    end
end




%% Combine animals - split into genotype - make heatmaps and mean+SEM for different EC. 

    figure
    all_geno = find(string(xy_opto_info.Geno) == "wt"); 
    n_trials_geno = numel(all_geno);
    ax1 = subplot(1,8,1:3); imagesc(xy_opto(all_geno, :)); title("WT"); caxis([0 70]); 
    hold on 
    plot([frames_b4 frames_b4], [0 n_trials_geno+10], 'w:', 'LineWidth', 1.2)
    plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 n_trials_geno+10], 'w:', 'LineWidth', 1.2)
    xticklabels({''})
%     colormap(redblue)
    caxis([0 80])
    
    subplot(1,8,4)
    imagesc(cell2mat(xy_opto_info{all_geno, 9}))
    axis off 
    box off
    colormap(gray)

    
    all_geno = find(string(xy_opto_info.Geno) == "het"); 
    n_trials_geno = numel(all_geno);
    ax2 = subplot(1,8,5:7); imagesc(xy_opto(all_geno, :)); title("HET"); caxis([0 70]); 
    hold on 
    plot([frames_b4 frames_b4], [0 n_trials_geno+10], 'w:', 'LineWidth', 1.2)
    plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 n_trials_geno+10], 'w:', 'LineWidth', 1.2)
    xticklabels({''})
%     colormap(redblue)
    caxis([0 80])

    subplot(1,8,8)
    imagesc(cell2mat(xy_opto_info{all_geno, 9}))
    axis off
    box off
    linkaxes([ax1,ax2],'x')
  


  
%% For step plots - geno based. 

    all_EC = unique(cell2mat(xy_opto_info.EC)); 
    n_EC = numel(all_EC); 

    for j = 2:n_EC-3
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "het");
        n_trials = numel(all_ani_EC);
        
        figure

        if ~isempty(all_ani_EC)
            
            for jj = 1:n_trials
                plot(xy_opto(all_ani_EC(jj), 575:725), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.2)
                hold on 
            end 
            plot(mean(xy_opto(all_ani_EC, 575:725)), 'Color', 'r', 'LineWidth', 1.5)
            plot([25 85], [-5 -5], 'k', 'LineWidth', 3);
            box off
            axis off
            axis([0 200 -10 100])
            title(EC_str)
        end 

    end 

%% WT/HET mean plotted for step plot

  all_EC = unique(cell2mat(xy_opto_info.EC)); 

    for j =[2,4,6,8,9,10,11]
        
        EC_str = all_EC(j,1);
        figure
        
        all_WT = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "wt");
        all_HET = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "het");
        
        plot(mean(xy_opto(all_WT, 575:725)), 'Color', 'k', 'LineWidth', 1.5)
        hold on
        plot(mean(xy_opto(all_HET, 575:725)), 'Color', 'r', 'LineWidth', 1.5)
        plot([25 85], [-5 -5], 'k', 'LineWidth', 3);
        box off
        axis off
        axis([0 200 -10 70])
        title(EC_str)
    end



%%

    all_EC = unique(cell2mat(xy_opto_info.EC)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = [1:6,8:n_EC]
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Geno) == "wt");
        
        if EC_str == 0 || EC_str == 1
            v2 = 0.85; 
            col = [v2 v2 v2]; 
        elseif EC_str > 0
%             v2 = 1-(j/n_EC);
            v2 = 1-(EC_str/5)+0.15; 
            col = [v2 v2 v2]; 
        end 

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 100], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 100], 'k:', 'LineWidth', 1.2)
        xticklabels({''})
    end 
    title("WT")
    axis([500 800 0 75])
    box off
    ylabel('Speed - cm/s')
    xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
    xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
    ax = gca;
    ax.FontSize = 14; 


    %% TIME PULSE
    
    all_EC = unique(cell2mat(xy_opto_info.T_pulse)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC
        
        subplot(1,7,j)
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.T_pulse) == EC_str & string(xy_opto_info.Geno) == "het");
        
        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', 'r', 'LineWidth', 1.2);
            hold on 
        end 
        
        plot([frames_b4 frames_b4], [0 40], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 40], 'k:', 'LineWidth', 1.2)
        xticklabels({''})
        axis([500 800 0 40])
        box off
        title(EC_str)
    end 
   

%% xy analysis



% Want to find :

% 1 - probability of escape

% 2 - reaction time

% 3 - escape speed


% First, make subsets of xy_analysis with only trials per each animal

  all_animals = unique(string((xy_opto_info.Animal)));
  n_animals = numel(all_animals); 


for i = 1:n_animals
    ani = all_animals(i); 
    all_ani = find(string(xy_analysis.Animal) == ani);
    xy_analysis2 = xy_analysis(all_ani, :); 

 all_EC = unique(cell2mat(xy_analysis2.EC)); 
 n_EC = numel(all_EC); 

 for j = 1:n_EC
     EC_val = all_EC(j,1);
     EC_rows = find(cell2mat(xy_analysis2.EC) == EC_val);
     val = mean(cell2mat(xy_analysis2{EC_rows, 28})); 
     all_EC(j,2) = val;
 end 
 figure
 if i == 1 | i == 3
 plot(all_EC(:,1), all_EC(:,2), 'k.-')
 else 
      plot(all_EC(:,1), all_EC(:,2), 'r.-')
 end 
 title(strcat('Change Speed - ', ani))
 box off
%  axis([0 5 -0.1 1.1])
     
end 



%%

  all_animals = unique(string(xy_analysis.Animal)); 
  n_animals = numel(all_animals); 

  % 14 - change speed. 
  % 15 . = LSI_during. 
  % 21 - Freeze
  % 25 - return to shelter. 
  % 27 - speed to shelter. 
  % 28 - maxspescape. 
  % 29 - time to maxsp. 
  
  column = 28;  
  
  for i = 1:n_animals
      ani = all_animals(i);
%       all_ani = find(string(xy_analysis.Animal) == ani & cell2mat(xy_analysis.ReturnToShelter) == 1);
       all_ani = find(string(xy_analysis.Animal) == ani & cell2mat(xy_analysis.Freeze) == 1);
      xy_analysis2 = xy_analysis(all_ani, :);
      geno = string(xy_analysis2.Geno(1));
      
      all_EC = unique(cell2mat(xy_analysis2.EC));
      n_EC = numel(all_EC);
      
      xvals = []; 
      yvals = []; 
      
      for j = 1:n_EC
          EC_val = all_EC(j,1);
          EC_rows = find(cell2mat(xy_analysis2.EC) == EC_val);
          n_ECrows = numel(EC_rows);
          
          dt = cell2mat(xy_analysis2{EC_rows, column});
          val = mean(dt);
          sem_val = std(dt)/sqrt(numel(dt));
          all_EC(j,2) = val;
          all_EC(j,3) = sem_val;
          all_EC(j,4) = numel(EC_rows);
          
%           xvals = [xvals, ones(1,n_ECrows)*EC_val];
%           yvals = [yvals, dt']; 
      end
      
%       [Qpre, p, sm, varcov] = fit_logistic(xvals,yvals);
      
      subplot(2,6,i)
      if geno == "wt"
          errorbar(all_EC(:,2)', all_EC(:,3)', 'k', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 28)
      elseif geno == "het"
          errorbar(all_EC(:,2)', all_EC(:,3)', 'r', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 28)
      end
      
%       % Add logistic curve.
%       plot(xvals,yvals,'o')
%        plot(xvals,Qpre)    % best fitting logistic
%       % Add text values of 'p'
%       text(1,0.1,string(p))
       
       
%       title(strcat('Norm Max Speed -  ', ani))
      title((ani))

      box off
      hold on
      xticks(1:1:numel(all_EC(:,1)))
      xticklabels(string(all_EC(:,1)))
      %  axis([0 5 -0.1 1.1])
      axis([0 numel(all_EC(:,1))+1 0 100])
%       ylabel('Speed - cm/s')
      ax = gca;
      ax.FontSize = 14;
      xlabel('Arbitrary Laser Power')
      
  end
  
%%


 
  column = 28;  
  
for i = 1:n_animals
      ani = all_animals(i);
%       all_ani = find(string(xy_analysis.Animal) == ani ); & cell2mat(xy_analysis.MaxSpEscape) >20
       all_ani = find(string(xy_analysis.Animal) == ani);
      xy_analysis2 = xy_analysis(all_ani, :);
      geno = string(xy_analysis2.Geno(1));
      
      all_EC = unique(cell2mat(xy_analysis2.EC));
      n_EC = numel(all_EC);
      
      xvals = []; 
      yvals = []; 
      
      for j = 1:n_EC
          EC_val = all_EC(j,1);
          EC_rows = find(cell2mat(xy_analysis2.EC) == EC_val);
          n_ECrows = numel(EC_rows);
          
          dt = cell2mat(xy_analysis2{EC_rows, column});
          val = nanmean(dt);
          sem_val = nanstd(dt)/sqrt(numel(dt));
          all_EC(j,2) = val;
          all_EC(j,3) = sem_val;
          all_EC(j,4) = numel(EC_rows);
          
          xvals = [xvals, ones(1,n_ECrows)*EC_val];
          yvals = [yvals, dt']; 
      end
      
      [Qpre, p, sm, varcov] = fit_logistic(xvals,yvals);
      
      figure
      subplot(1,2,1)
      if geno == "wt"
          errorbar(all_EC(:,2)', all_EC(:,3)', 'k', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 28)
      elseif geno == "het"
          errorbar(all_EC(:,2)', all_EC(:,3)', 'r', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 28)
      end
      box off
      hold on
      xticks(1:1:numel(all_EC(:,1)))
      xticklabels(string(all_EC(:,1)))
      %  axis([0 5 -0.1 1.1])
      axis([0 numel(all_EC(:,1))+1 0 100])
%       ylabel('Speed - cm/s')
      ax = gca;
      ax.FontSize = 18;
      xlabel('Laser Power')

      subplot(1,2,2)
      % Add logistic curve.
      plot(xvals,yvals,'o')
       plot(xvals,Qpre)    % best fitting logistic
      % Add text values of 'p'
      text(0.5,24.5,string(p))
       ax = gca;
      ax.FontSize = 18;
      xlabel('Laser Power')
      axis([0 4.5 0 100])
      
      sgtitle((ani))
end 
      

  
% for i = 1:726
%     if isempty(xy_analysis.TimeHalfHeight{i})
%         xy_analysis.TimeHalfHeight{i} = NaN;
%     end 
% end 
    


%% Errorbar

 col = 'r';  %Ptchd1

 % 1, 1.75, 1.9, 2, 2.25, 2.5, 2.75, 3. 
 EC_vals = [1, 1.75, 1.9, 2, 2.25, 2.5, 2.75, 3];

% 

for i = 1:n_animals
    ani = all_animals(i);
    
    data_per_day = zeros(1,8);
    sem_per_day = zeros(1,8);
    n_per_day  =zeros(1,8);
    
    for j = 1:8
        EC = EC_vals(j);
        
        allWT = find(string(xy_analysis.Animal) == ani & cell2mat(xy_analysis.EC) == EC);
        
        data_WT = cell2mat(xy_analysis.SpeedToShelter(allWT));
        
        data_per_day(1,j) = nanmean(data_WT);
        
        n_per_day(1,j) = numel(data_WT);
        
        sem_per_day(1,j) = std(data_WT)/sqrt(numel(data_WT));
        
    end
    
    if i == 1 || i ==3 
        col = 'k';
    else 
        col = 'r';
    end 
    
    % % Errorbar
    figure
    errorbar(data_per_day(1,:), sem_per_day(1,:), col, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 28)
    xticks([1,2,3,4,5,6,7, 8])
    xticklabels({EC_vals});
    xlabel('EC')
    ylabel('Speed - cm/s')
%     ylabel('Probability')
%     ylabel('Time - s')
    set(gca, 'FontSize', 15)
    axis([0.5 8.5 0 30])
    box off
    title(strcat('Speed to Shelter - ', ani))
end







%%   PLOT POSITION OF MOUSE WHEN LIGHT STARTED. 

oob = []; 

figure
for i = 1:569 
    
    xval = XY_OPTO((2*i)-1, 600); 
    yval = 416 - XY_OPTO(2*i, 600); 
    
    if xval >=175 && yval <=230
        oob = [oob, i];
        col = 'r.';
    else 
        col = 'k.';
    end 
    
    plot(xval, yval, col)
    hold on 
    
end 


xy_analysis(oob, :) = []; 
xy_opto(oob, :) = []; 
xy_opto_info(oob, :) = []; 


%% Combine 

b1 = xy_analysis; 
b2 = xy_opto; 
b3 = xy_opto_info; 


xy_analysis = vertcat(a1, b1, xy_analysis);
xy_opto = vertcat(a2, b2, xy_opto);
xy_opto_info = vertcat(a3, b3, xy_opto_info); 

save('3C_xy_opto.mat', 'xy_analysis', 'xy_opto', 'xy_opto_info');

day1 = find(xy_opto_info.Day == 1);
day2 = find(xy_opto_info.Day == 2);
day3 = find(xy_opto_info.Day == 3);




%%

xy_optoA = xy_opto(day1, :);
xy_opto_infoA = xy_opto_info(day1, :); 

% SORT
xy_optoA(:,1201) = cell2mat(xy_opto_infoA.EC);
xy_opto_infoA = sortrows(xy_opto_infoA, 9); % 9 for EC
xy_optoA = sortrows(xy_optoA, 1201); 
% Remove last column 
xy_optoA = xy_optoA(:, 1:1200); 

allWT = find(string(xy_opto_infoA.Geno) == "wt");
allHET = find(string(xy_opto_infoA.Geno) == "het");

figure
subplot(1,2,1)
ax1 = imagesc(xy_optoA(allWT, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('WT')
subplot(1,2,2)
ax2 =imagesc(xy_optoA(allHET, :));
hold on 
plot([600 600], [0 1000], 'w:', 'LineWidth', 2)
xlim([450 850])
caxis([0 150])
title('HET')




n = numel(xy_optoA(:,1)); 
col = 'r'; 
frames_b4 = 600;
total_lighton_frames = 60;  % light was on for 1s. - 60 frames. 

% EC - per animal

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_infoA.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_infoA.EC)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & xy_opto_infoA.Animal == ani);
       
        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_optoA(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
        xticklabels({''})
    end 
    title(ani)
    axis([500 800 0 70])
    box off
    ylabel('Speed - cm/s')
    xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
    xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
    ax = gca;
    ax.FontSize = 16; 
end 


% EC - per geno

for i = 1:2
    
    if i  == 1
        geno = 'wt'; 
    else 
        geno = 'het';
    end 
    
    all_ani = find(string(xy_opto_infoA.Geno) == geno); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_infoA.EC)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC %[1,2,3,4,5,6,7,10,12,16,17]
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & string(xy_opto_infoA.Geno) == geno);
       
        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_optoA(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
        xticklabels({''})
%         tttt = input(string(EC_str)); 
    end 
%     title(ani)
    axis([500 800 0 60])
    box off
    ylabel('Speed - cm/s')
    xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
    xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
    ax = gca;
    ax.FontSize = 16; 
end 


% vals = {}; 
% 
% for i = 1:2
%     
%     if i  == 1
%         geno = 'wt'; 
%     else 
%         geno = 'het';
%     end 
%     
%     all_ani = find(string(xy_opto_infoA.Geno) == geno); 
%     n_trials_ani = numel(all_ani);
%     
%     all_EC = unique(cell2mat(xy_opto_infoA.EC)); 
%     n_EC = numel(all_EC); 
%     
%     for j = 1:n_EC
%         EC_str = all_EC(j,1); 
%         all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & string(xy_opto_infoA.Geno) == geno);
%         vals(j, 1) = {all_ani_EC};     %numel(all_ani_EC); 
%     end 
% %   
% end 


%%

% Pulse Width - per geno

for i = 1:2
    
    if i  == 1
        geno = 'wt'; 
    else 
        geno = 'het';
    end 
    
    all_ani = find(string(xy_opto_info.Geno) == geno); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.T_pulse)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.T_pulse) == EC_str & string(xy_opto_info.Geno) == geno);
       
        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
        xticklabels({''})
        tttt = input(string(EC_str)); 
    end 
%     title(ani)
    axis([500 800 0 70])
    box off
    ylabel('Speed - cm/s')
    xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
    xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
    ax = gca;
    ax.FontSize = 16; 
end 


for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_info.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.T_pulse)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.T_pulse) == EC_str & xy_opto_info.Animal == ani);
       
        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_opto(all_ani_EC, :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
        xticklabels({''})
    end 
    title(ani)
    axis([500 800 0 70])
    box off
    ylabel('Speed - cm/s')
    xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
    xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
    ax = gca;
    ax.FontSize = 16; 
end 


%13
% 20 ,25



%% OPLO



for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_info.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.EC)); 
    n_EC = numel(all_EC); 
    subplot(n_animals,1, i)
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & xy_opto_info.Animal == ani);
        
       if j ==1 
           col = [0.6 0.6 0.6];
       elseif j ==2
           col = [0 0 0];
       end 

        if ~isempty(all_ani_EC)
            av_resp = (xy_opto(all_ani_EC(10), :));
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 60], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+30 frames_b4+30], [0 60], 'r:', 'LineWidth', 1.2)
        plot([frames_b4+90 frames_b4+90], [0 60], 'k', 'LineWidth', 0.75)
        plot([frames_b4+30+46 frames_b4+30+46], [0 60], 'r:', 'LineWidth', 1.2)
        plot([frames_b4+30+46*2 frames_b4+30+46*2], [0 60], 'r:', 'LineWidth', 1.2)
        plot([frames_b4+30+46*3 frames_b4+30+46*3], [0 60], 'r:', 'LineWidth', 1.2)
        plot([frames_b4+30+46*4 frames_b4+30+46*4], [0 60], 'r:', 'LineWidth', 1.2)
        plot([frames_b4+30+46*5 frames_b4+30+46*5], [0 60], 'r', 'LineWidth', 1)

        xticklabels({''})
    end 
    title(ani)
    axis([400 900 0 70])
    box off

end 






