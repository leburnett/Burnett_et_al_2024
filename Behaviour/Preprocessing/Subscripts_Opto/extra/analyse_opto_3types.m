
xy_optoA = xy_opto; 
xy_opto_infoA = xy_opto_info;
xy_analysisA = xy_analysis;

% xy_optoC = xy_opto; 
% xy_opto_infoC = xy_opto_info;
% xy_analysisC = xy_analysis;

% xy_optoA = vertcat(xy_opto, xy_optoA); 
% xy_opto_infoA = vertcat(xy_opto_info, xy_opto_infoA);
% xy_analysisA = vertcat(xy_analysis, xy_analysisA);


%%

n = height(xy_opto_info);
all_animals = unique(xy_opto_info.Animal);
n_animals = numel(all_animals);

allWT = find(string(xy_opto_info.Geno) == "wt");
allHET = find(string(xy_opto_info.Geno) == "het");

all_animals_WT = unique(xy_opto_info.Animal(allWT));
all_animals_HET = unique(xy_opto_info.Animal(allHET));



%% HEAMTAPS - Sorted by time to max speed for ech laser power - per animal

n = numel(xy_optoA(:,1));
col = 'r';
frames_b4 = 595;
total_lighton_frames = 60;  % light was on for 1s. - 60 frames.

all_animals = unique(string(xy_opto_infoA.Animal));
n_animals = numel(all_animals);

EC_VALS = [1.9, 2, 2.25, 2.5, 2.75, 3];
  
% EC

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_infoA.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
%     all_EC = unique(cell2mat(xy_opto_infoA.EC(all_ani))); 
%     n_EC = numel(all_EC); 
    
    figure
    
    for j = 1:6
        
        % Make a plot where each laser power is a new subplot
        subplot(6, 1, j)
        
        % %Define the laser power to assess and find all the rows
        % corresponding to that animal at that laser power. 
%         EC_str = all_EC(j,1); 
        EC_str = EC_VALS(j);
        all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & xy_opto_infoA.Animal == ani);

        if ~isempty(all_ani_EC)
            %Extract the speed traces for these trials. 
            av_resp = (xy_optoA(all_ani_EC, :));
            imagesc(av_resp)
            colormap(redblue)
            caxis([0 70])
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 90], 'w:', 'LineWidth', 1.5)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'w:', 'LineWidth', 1.5)
        xticks([])
        yticks([])
        ylabel(EC_str)
        xlim([500 800])
    end 
    sgtitle(ani)
%     axis([500 800 0 90])
    box off
%     ylabel('Speed - cm/s')
%     xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
%     xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
    f = gcf;
    f.Position = [440 73 299 705]; 

end 







%% BY GENOTYPE!!! % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% 
% % Rows removed. 
% save('210924_LaserPower.mat', 'xy_opto', 'xy_analysis', 'xy_opto_info');
% savefig(gcf, '210924_WT_LaserPower_595.fig');
% savefig(gcf, '210924_HET_LaserPower_595.fig');
% 


%% MEAN SEM - Line plot - ONLY MAX EC for each animal - genotype. 

xy_optoA = xy_opto;
xy_opto_infoA = xy_opto_info;
xy_analysisA = xy_analysis; 

%

laser_val = 6; 

figure
rectangle('Position', [frames_b4 52 60 5], 'FaceColor', [0.67 0.84 0.9], 'EdgeColor', [0.67 0.84 0.9]) %xy wh 
hold on 

for i = 1:2
    if i ==1
        all_ani = find(string(xy_opto_infoA.Geno) == "wt");
        n_trials_ani = numel(all_ani);
        
        all_EC = unique(cell2mat(xy_opto_infoA.EC));
        n_EC = numel(all_EC);
        
        for j = laser_val %1:n_EC
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & string(xy_opto_infoA.Geno) == "wt");
            
            speed_WT = (xy_optoA(all_ani_EC, :));

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
        axis([500 800 0 60])
        box off
        ylabel('Speed - cm/s')
        xticks([frames_b4-60,frames_b4, frames_b4+60, frames_b4+120, frames_b4+180])
        xticklabels({'-1', '0', '1', '2', '3'})
        ax = gca;
        ax.FontSize = 20;
        ax.TickDir = 'out'; 
        
    elseif i ==2
        all_ani = find(string(xy_opto_infoA.Geno) == "het");
        n_trials_ani = numel(all_ani);
        
        all_EC = unique(cell2mat(xy_opto_infoA.EC));
        n_EC = numel(all_EC);
        
        for j = laser_val%1:n_EC
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & string(xy_opto_infoA.Geno) == "het");
            speed_HET = (xy_optoA(all_ani_EC, :));
            
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

  plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
  plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)


  
  
  
  
  
  
  
  
  
  
  
  
  
%% Heatmap plots - ordered by EC - WT/ HET. 

% - Need to make this with normalised laser power too! 

    allWT = find(string(xy_opto_infoA.Geno) =="wt");
    allHET = find(string(xy_opto_infoA.Geno) =="het");
    
    WT_array = xy_optoA(allWT, :);
    HET_array = xy_optoA(allHET, :);
    
    imagesc(WT_array(:, 500:850))
    figure;imagesc(HET_array(:, 500:850))
    
 % ORDERED BY EC. 
 
 xy_opto2 = xy_optoA; 
 xy_opto_info2 = xy_opto_infoA; 
 
 xy_opto2(:,1201) = cell2mat(xy_opto_infoA.EC); 
 xy_opto2 = sortrows(xy_opto2, 1201); 
 xy_opto_info2 = sortrows(xy_opto_info2, 9); 
 
    allWT = find(string(xy_opto_info2.Geno) =="wt");
    allHET = find(string(xy_opto_info2.Geno) =="het");
    
    WT_array = xy_opto2(allWT, :);
    HET_array = xy_opto2(allHET, :);
    
    % PLOT % % % % % % % % % % % % %
    subplot(1,2,1); imagesc(WT_array(:, 550:850))
    hold on; plot([50 50], [0 450], 'w:', 'LineWidth', 1.5); plot([150 150], [0 450], 'w:', 'LineWidth', 1.5); yticklabels({''}); xticklabels('')
    subplot(1,2,2);imagesc(HET_array(:, 550:850))
    hold on; plot([50 50], [0 450], 'w:', 'LineWidth', 1.5); plot([150 150], [0 450], 'w:', 'LineWidth', 1.5);  yticklabels({''}); xticklabels('')

%% 





%% HEATMAPS - only the highest EC values - sorted by T2M 

xy_opto2 = xy_optoA;
xy_opto_info2 = xy_opto_infoA;

% xy_opto2(:,1201) = cell2mat(xy_opto_infoA.T2M);
% 
% xy_opto2 = sortrows(xy_opto2, 1201);
% xy_opto_info2 = sortrows(xy_opto_info2, 20);

allWT = find(string(xy_opto_info2.Geno) =="wt" & cell2mat(xy_opto_info2.EC) == 3);
allHET = find(string(xy_opto_info2.Geno) =="het" & cell2mat(xy_opto_info2.EC) == 3);

WT_array = xy_opto2(allWT, :);
HET_array = xy_opto2(allHET, :);

%% How many animals from each:

numel(unique(xy_opto_info2.Animal(allWT)))
numel(unique(xy_opto_info2.Animal(allHET)))

% WT = N = 4, n = 60
% HET = N = 4, n = 60

% PLOT % % % % % % % % % % % % %


% Timing:

% 600-660 = LASER ON.  On for 1s! 
% 3s before = 600-180 = 420
% 7s after = 600 + (60*7) = 1020

% Draw lines for laser at 180 and 240. 

% windw = [550:850]; 
windw = [420:1020]; 

figure
subplot(1,2,1); imagesc(WT_array(:, windw))
% hold on; plot([50 50], [0 450], 'w:', 'LineWidth', 1.5); plot([150 150], [0 450], 'w:', 'LineWidth', 1.5); yticklabels({''}); xticklabels('')
hold on; plot([180 180], [0 450], 'w:', 'LineWidth', 1.5); plot([240 240], [0 450], 'w:', 'LineWidth', 1.5); yticklabels({''}); xticklabels('')

subplot(1,2,2);imagesc(HET_array(1:58, windw))
% hold on; plot([50 50], [0 450], 'w:', 'LineWidth', 1.5); plot([150 150], [0 450], 'w:', 'LineWidth', 1.5);  yticklabels({''}); xticklabels('')
hold on; plot([180 180], [0 450], 'w:', 'LineWidth', 1.5); plot([240 240], [0 450], 'w:', 'LineWidth', 1.5);  yticklabels({''}); xticklabels('')


% WT_array([10, 13, 21 64], :) = []; 
% HET_array([1:8], :) = []; 


%% Find frame at which speed > 20 cm/s 

rows = []; 
for i = 1:numel(WT_array(:, 1))
    
    array = WT_array(i, :); 
    val_20 = min(find(array(609:900)>20));
    
    if isempty(val_20)
        rows= [rows, i];
    else
        WT_array(i, 1202) = val_20;
    end
    
end 

% Remove rows where the mouse ~= > 20 cm/s
WT_array(rows, :) = []; 

rows = []; 
for i = 1:numel(HET_array(:, 1))
    
    array = HET_array(i, :); 
    
    val_20 = min(find(array(609:900)>20));
    
    if isempty(val_20)
        rows= [rows, i];
    else
        HET_array(i, 1202) = val_20;
    end
    
end 

HET_array(rows, :) = []; 


WT_array2 = sortrows(WT_array, 1202);
HET_array2 = sortrows(HET_array, 1202);

% figure
% subplot(1,2,1); imagesc(WT_array2(:, windw))
% % hold on; plot([50 50], [0 450], 'w:', 'LineWidth', 1.5); plot([150 150], [0 450], 'w:', 'LineWidth', 1.5); yticklabels({''}); xticklabels('')
% hold on; plot([180 180], [0 450], 'w:', 'LineWidth', 1.5); plot([240 240], [0 450], 'w:', 'LineWidth', 1.5); yticklabels({''}); xticklabels('')
% 
% subplot(1,2,2);imagesc(HET_array2(:, windw))
% % hold on; plot([50 50], [0 450], 'w:', 'LineWidth', 1.5); plot([150 150], [0 450], 'w:', 'LineWidth', 1.5);  yticklabels({''}); xticklabels('')
% hold on; plot([180 180], [0 450], 'w:', 'LineWidth', 1.5); plot([240 240], [0 450], 'w:', 'LineWidth', 1.5);  yticklabels({''}); xticklabels('')
% 

%%

figure;
% subplot(2,1,1)
imagesc(WT_array2(:, windw))
n_h = numel(WT_array2(:, 1));
hold on; plot([180 180], [0 450], 'w:', 'LineWidth', 2); plot([240 240], [0 450], 'w', 'LineWidth', 2); yticklabels({''}); xticklabels('')
axis off
box off
colormap(redblue)
caxis([0 70])
f = gcf; 
f.Position = [100 600 750 n_h*8];

figure;
% subplot(2,1,1)
imagesc(HET_array2(:, windw))
n_h = numel(HET_array2(:,1));
hold on; plot([180 180], [0 450], 'w:', 'LineWidth', 2); plot([240 240], [0 450], 'w', 'LineWidth', 2);  yticklabels({''}); xticklabels('')
axis off
box off
colormap(redblue)
caxis([0 70])
f = gcf; 
f.Position = [100 600 750 n_h*8];


%%






%% Create table with average speed for each Laser Power for each animal. 

av_animal_speed = {}; 

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(string(xy_opto_info.Animal) == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_info.EC)); 
    n_EC = numel(all_EC); 

    geno_str = xy_opto_info.Geno{all_ani(1)}; 
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_info.EC) == EC_str & string(xy_opto_info.Animal) == ani); 
        
        if ~isempty(all_ani_EC)
            n_ani = numel(all_ani_EC);
            if n_ani ==1 
                av_resp = xy_opto(all_ani_EC, :);
            elseif n_ani >1 
                av_resp = mean(xy_opto(all_ani_EC, :));
            end 
        end 
        av_animal_speed{j+n_EC*(i-1),1} = ani; 
        av_animal_speed{j+n_EC*(i-1),2} = EC_str;
        av_animal_speed{j+n_EC*(i-1),3} = geno_str;
        av_animal_speed{j+n_EC*(i-1),4} = av_resp; 
        av_resp = []; 
    end 
end 

Animal = av_animal_speed(:,1);
EC = av_animal_speed(:,2);
Geno = av_animal_speed(:,3);
Speed = av_animal_speed(:,4);

opto_animal = table(Animal, EC, Geno, Speed);

%%


for i = 1:2
    if i ==1
        all_ani = find(string(opto_animal.Geno) == "wt");
        all_EC = unique(cell2mat(opto_animal.EC));
        n_EC = numel(all_EC);
        figure
        
        for j = 1:n_EC%[2, 3, 4, 8, 11, 12, 15, 16, 17, 19] % 1:n_EC
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(opto_animal.EC) == EC_str & string(opto_animal.Geno) == "wt");
            
            col = 1 -[j/n_EC j/n_EC j/n_EC];
            
            if ~isempty(all_ani_EC)&& numel(all_ani_EC)>1
                av_resp = nanmean(cell2mat(opto_animal.Speed(all_ani_EC)));
                plot(av_resp, 'Color', col, 'LineWidth', 1.2);
                hold on
            elseif ~isempty(all_ani_EC)&& numel(all_ani_EC)==1
                av_resp = (cell2mat(opto_animal.Speed(all_ani_EC)));
                plot(av_resp, 'Color', col, 'LineWidth', 1.2);
                hold on
            end
            
            plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
            plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
            xticklabels({''})
        end
        title('WT')
        axis([500 800 0 45])
        box off
        ylabel('Speed - cm/s')
        xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
        xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
        ax = gca;
        ax.FontSize = 14;
           
        
        
        
    elseif i ==2
        all_ani = find(string(opto_animal.Geno) == "het");
        all_EC = unique(cell2mat(opto_animal.EC));
        n_EC = numel(all_EC);
        figure

        for j = 1:n_EC %[2, 3, 4, 8, 9,12, 15, 16, 17, 19]
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(opto_animal.EC) == EC_str & string(opto_animal.Geno) == "het");
            col = 1 -[j/n_EC j/n_EC j/n_EC];
            
            if ~isempty(all_ani_EC)&& numel(all_ani_EC)>1
                av_resp = nanmean(cell2mat(opto_animal.Speed(all_ani_EC)));
                plot(av_resp, 'Color', col, 'LineWidth', 1.2);
                hold on
            elseif ~isempty(all_ani_EC)&& numel(all_ani_EC)==1
                av_resp = (cell2mat(opto_animal.Speed(all_ani_EC)));
                plot(av_resp, 'Color', col, 'LineWidth', 1.2);
                hold on
            end
            
            plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
            plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
            xticklabels({''})
            
        end
        title('HET')
        axis([500 800 0 45])
        box off
        ylabel('Speed - cm/s')
        xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
        xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
        ax = gca;
        ax.FontSize = 14;
        
    end
end

title('')
xlabel('')
ylabel('')
xticklabels({''})
yticklabels({''})

p = 1; 












    
    %%
    
%     
%        ani = string(all_animals{i}); 
%     all_ani = find(xy_opto_infoA.Animal == ani); 
%     n_trials_ani = numel(all_ani);
%     
%     all_EC = unique(cell2mat(xy_opto_infoA.EC)); 
%     n_EC = numel(all_EC); 
%  
%     val = 0; 
%     plotnum = 1; 
%     
%     figure
%     for j = 6:n_EC
%         subplot(6,1,plotnum)
%         EC_str = all_EC(j,1);
%         all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & xy_opto_infoA.Animal == ani);
%         n_trials = numel(all_ani_EC);
%         
%         if EC_str == 0 || EC_str == 1
%             v2 = 0.85;
%             col = [v2 v2 v2];
%         elseif EC_str > 0
%             v2 = (1.3-(j/n_EC))+(0.02*j);
%             %             v2 = 1-(EC_str/20);
%             col = [v2 v2 v2];
%         end
%         
%         %         col = 1 -[j/n_EC j/n_EC j/n_EC];
%         
%         for jj = 1:n_trials
%             plot(xy_optoA(all_ani_EC(jj), :), 'Color', col, 'LineWidth', 1.2);
%             hold on
%         end 
%         
% %         if ~isempty(all_ani_EC)
% %             av_resp = mean(xy_optoA(all_ani_EC, :));
% %             plot(av_resp, 'Color', col, 'LineWidth', 1.2);
% %             hold on
% %         end
%         
%         xticklabels({''})
%         title(ani)
%         axis([400 900 0 110])
%         box off
%         ylabel('Speed - cm/s')
%         xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
%         xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
%         plot([frames_b4 frames_b4], [0 120], 'k:', 'LineWidth', 1.2)
%         text(800, 80-val, string(EC_str));
%         plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 120], 'k:', 'LineWidth', 1)
%         inp = input('');
%         val = val+4;
%         plotnum = plotnum+1; 
%         
%     end
    
    
%%


%    
%     ani = string(all_animals{i}); 
%     all_ani = find(xy_opto_infoA.Animal == ani); 
%     n_trials_ani = numel(all_ani);
%     
%     all_EC = unique(cell2mat(xy_opto_infoA.EC)); 
%     n_EC = numel(all_EC); 
%  
%     val = 0; 
%     plotnum = 1; 
%     
%     figure
%     for j = 1:n_EC
%         subplot(14,1,plotnum)
%         EC_str = all_EC(j,1);
%         all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & xy_opto_infoA.Animal == ani);
%         n_trials = numel(all_ani_EC);
%         
%         col = 'r';
% %         if j == 4
% %             v2 = 0.9;
% %             col = [v2 v2 v2];
% %         elseif j ~= 4
% %             v2 = (1.25-(j/14))+(0.02*j);
% %             %             v2 = 1-(EC_str/20);
% %             col = [v2 v2 v2];
% %         end
% %         
%         %         col = 1 -[j/n_EC j/n_EC j/n_EC];
%         
%         for jj = 1:n_trials
%             plot(xy_optoA(all_ani_EC(jj), :), 'Color', col, 'LineWidth', 1.2);
%             hold on
%         end 
%         
% %         if ~isempty(all_ani_EC)
% %             av_resp = mean(xy_optoA(all_ani_EC, :));
% %             plot(av_resp, 'Color', col, 'LineWidth', 1.2);
% %             hold on
% %         end
%         
%         xticklabels({''})
%         title(ani)
%        
%         box off
%         ylabel('Speed - cm/s')
%         xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
%         xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
%         plot([frames_b4 frames_b4], [0 120], 'k:', 'LineWidth', 1.2)
%         text(800, 80-val, string(EC_str));
%         plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 120], 'k:', 'LineWidth', 1)
%         inp = input('');
%         val = val+4;
%         plotnum = plotnum+1; 
%         axis([500 800 0 110])
%         
%     end
    
    

%% Individual plots for each light power - WT/ HET overlaid. 

    
    all_EC = unique(cell2mat(xy_opto_infoA.EC)); 
    n_EC = numel(all_EC); 
    
    for j = 1:n_EC
        
      figure
      EC_str = all_EC(j,1);
      
      for jj = 1:2
          ani = string(all_animals{jj});
          all_ani = find(xy_opto_infoA.Animal == ani);
          n_trials_ani = numel(all_ani);
          
          all_ani_EC = find(cell2mat(xy_opto_infoA.EC) == EC_str & xy_opto_infoA.Animal == ani);
          n_trials = numel(all_ani_EC);
          
          if jj ==1
              col = [0.6 0.6 0.6];
          elseif jj == 2
              col = [1 0 0];
          end
          
          for jjj = 1:n_trials
              plot(xy_optoA(all_ani_EC(jjj), :), 'Color', col, 'LineWidth', 1);
              hold on
          end
          axis([400 850 0 90])
          plot([frames_b4 frames_b4], [0 120], 'k:', 'LineWidth', 1.2)
          title(EC_str)
          box off
          
      end
      
    end 
      
      
      
%         
figure
plot(mean(xy_optoC(16:20, :)), 'LineWidth', 1)
hold on 
plot(mean(xy_optoC(21:25, :)), 'LineWidth', 1)
plot(mean(xy_optoC(26:30, :)), 'LineWidth', 1)
plot(mean(xy_optoC(89:98, :)), 'LineWidth', 1)
% 



%% Individual plots for each pulse width - WT/ HET overlaid. 

    
    all_EC = unique(cell2mat(xy_opto_infoB.T_pulse)); 
    n_EC = numel(all_EC); 
    
    for j = 1:n_EC
        
      figure
      EC_str = all_EC(j,1);
      
      for jj = 1:2
          ani = string(all_animals{jj});
          all_ani = find(xy_opto_infoB.Animal == ani);
          n_trials_ani = numel(all_ani);
          
          all_ani_EC = find(cell2mat(xy_opto_infoB.T_pulse) == EC_str & xy_opto_infoB.Animal == ani);
          n_trials = numel(all_ani_EC);
          
          if jj ==1
              col = [0.6 0.6 0.6];
          elseif jj == 2
              col = [1 0 0];
          end
          
          for jjj = 1:n_trials
              plot(xy_optoB(all_ani_EC(jjj), :), 'Color', col, 'LineWidth', 1);
              hold on
          end
          axis([400 850 0 90])
          plot([frames_b4 frames_b4], [0 120], 'k:', 'LineWidth', 1.2)
          title(EC_str)
          box off
          
      end
      
    end 
      
    
    
    
%% Individual plots for each frequency - WT/ HET overlaid. 

    
    all_EC = unique(cell2mat(xy_opto_infoC.FreqPulse)); 
    n_EC = numel(all_EC); 
    
    for j = 1:n_EC
        
      figure
      EC_str = all_EC(j,1);
      
      for jj = 1:2
          ani = string(all_animals{jj});
          all_ani = find(string(xy_opto_infoC.Animal) == ani);
          n_trials_ani = numel(all_ani);
          
          all_ani_EC = find(cell2mat(xy_opto_infoC.FreqPulse) == EC_str & string(xy_opto_infoC.Animal) == ani);
          n_trials = numel(all_ani_EC);
          
          if jj ==1
              col = [0.6 0.6 0.6];
          elseif jj == 2
              col = [1 0 0];
          end
          
          for jjj = 1:n_trials
              plot(xy_optoC(all_ani_EC(jjj), :), 'Color', col, 'LineWidth', 1);
              hold on
          end
          axis([400 850 0 90])
          plot([frames_b4 frames_b4], [0 120], 'k:', 'LineWidth', 1.2)
          title(EC_str)
          box off
          
      end
      
    end

    %%

% for i = 447:569
%     if xy_opto_info.Animal{i} == "MJ1347"
%         xy_opto_info.Geno{i} = "wt";
%         xy_analysis.Geno{i} = "wt"; 
%     elseif xy_opto_info.Animal{i} == "MJ1348"
%         xy_opto_info.Geno{i} = "het";
%         xy_analysis.Geno{i} = "het"; 
%     end 
% end 
% 
% save('XY_OPTO_ARRAY_INFO.mat', 'xy_opto_info');
% save('XY_analysis.mat', 'xy_analysis');



%%



xy_optoA = xy_opto;
xy_opto_infoA = xy_opto_info;
xy_analysisA = xy_analysis; 

for i = 1:2
    if i ==1
        all_ani = find(string(xy_opto_infoA.Geno) == "wt");
        n_trials_ani = numel(all_ani);
        
        all_EC = unique(cell2mat(xy_opto_infoA.T_pulse));
        n_EC = numel(all_EC);
        figure
        
        for j = 1:n_EC %[1,2,3,6,5,8,10,11,12] %1:n_EC
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(xy_opto_infoA.T_pulse) == EC_str & string(xy_opto_infoA.Geno) == "wt");
            
            col = 1 -[j/n_EC j/n_EC j/n_EC];
            rectangle('Position', [frames_b4 45 60 5], 'FaceColor', [0.67 0.84 0.9], 'EdgeColor', [0.67 0.84 0.9]) %xy wh 

            if ~isempty(all_ani_EC)
                av_resp = mean(xy_optoA(all_ani_EC, :));
                if j == 6
                av_resp(550:570) = linspace(av_resp(550), av_resp(570), 21);
                end 
                plot(av_resp, 'Color', col, 'LineWidth', 1.4);
                hold on
            end
            plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1)
            plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
            xticklabels({''})
        end
        axis([500 800 0 50])
        box off
        ylabel('Speed - cm/s')
%         xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
%         xticks([frames_b4-90, frames_b4-60, frames_b4-30, frames_b4, frames_b4+30, frames_b4+60, frames_b4+90, frames_b4+120, frames_b4+150, frames_b4+180])
%         xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'}) 
         xticks([frames_b4-60,frames_b4, frames_b4+60, frames_b4+120, frames_b4+180])
        xticklabels({'-1', '0', '1', '2', '3'})
        ax = gca;
        ax.FontSize = 20;
        ax.TickDir = 'out';
        
%         x.YAxis.Visible = 'off';
%         x.XAxis.Visible = 'off';
        
    elseif i ==2
        all_ani = find(string(xy_opto_infoA.Geno) == "het");
        n_trials_ani = numel(all_ani);
        
        all_EC = unique(cell2mat(xy_opto_infoA.T_pulse));
        n_EC = numel(all_EC);
        figure
        rectangle('Position', [frames_b4 45 60 5], 'FaceColor', [0.67 0.84 0.9], 'EdgeColor', [0.67 0.84 0.9]) %xy wh 
       hold on 
        for j = 1:n_EC %[1,3,4,5,8,10,11,13] %1:n_EC
            EC_str = all_EC(j,1);
            all_ani_EC = find(cell2mat(xy_opto_infoA.T_pulse) == EC_str & string(xy_opto_infoA.Geno) == "het");
            
            col = 1 -[j/n_EC j/n_EC j/n_EC];
            
            if ~isempty(all_ani_EC)
                av_resp = mean(xy_optoA(all_ani_EC, :));
                plot(av_resp, 'Color', col, 'LineWidth', 1.2);
                hold on
            end
%             xticklabels({''})
        end
        plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
        axis([500 800 0 50])
        box off
        ylabel('Speed - cm/s')
%         xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
%         xticks([frames_b4-90, frames_b4-60, frames_b4-30, frames_b4, frames_b4+30, frames_b4+60, frames_b4+90, frames_b4+120, frames_b4+150, frames_b4+180])    
%         xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
       xticks([frames_b4-60,frames_b4, frames_b4+60, frames_b4+120, frames_b4+180])
        xticklabels({'-1', '0', '1', '2', '3'})

        ax = gca;
        ax.FontSize = 20;
        ax.TickDir = 'out'; 
        
%         ax.YAxis.Visible = 'off';
%         ax.XAxis.Visible = 'off';
        
    end
end




%%



%%


% Specify which trials are LASER POWER TRIALS - xy_optoA
xy_optoA = xy_opto([1:195], :);
xy_opto_infoA = xy_opto_info([1:195], :); 
xy_analysisA = xy_analysis([1:195], :);

% SORT these trials by laser power. 
xy_optoA(:,1201) = cell2mat(xy_opto_infoA.EC);
xy_opto_infoA = sortrows(xy_opto_infoA, 9); % 9 for EC
xy_optoA = sortrows(xy_optoA, 1201); 

% Remove last column once the table/array is sorted. 
xy_optoA = xy_optoA(:, 1:1200); 

% Find all WT/ HET trials. 
allWT = find(string(xy_opto_infoA.Geno) == "wt");
allHET = find(string(xy_opto_infoA.Geno) == "het");

% Create a heatmap of WT/HET trials sorted by laser power. 
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


% % % % % % % %

% PULSE WIDTH TRIALS - xy_optoB
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

% % % % % % % % 

% FREQ TRIALS - xy_optoC
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





%% P1 - HEATMAP of all trials - ordered by EC for each animal. 

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_opto_info.Animal == ani); 
    max_EC = find(max(xy_opto_info.EC(all_ani)));
    n_trials_ani = numel(all_ani);
    
    subplot(1,4,i); imagesc(xy_opto(all_ani, :)); title(ani); caxis([0 70]); 
    hold on 
    plot([frames_b4 frames_b4], [0 n_trials_ani+10], 'w:', 'LineWidth', 1.2)
    plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 n_trials_ani+10], 'w:', 'LineWidth', 1.2)
    xticklabels({''})
end 


%% P2 - Line plot for each animal - where trials of same EC are averaged. 

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
  


  
%% For step plots - geno based - mean trace and all trials in the background. 

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

%% Mean / SEM line plots for each Laser Power

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
   


%%

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



%%

% Frequency - per geno

for i = 1:2
    
    if i  == 1
        geno = 'wt'; 
    else 
        geno = 'het';
    end 
    
    all_ani = find(string(xy_opto_infoC.Geno) == geno); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_infoC.FreqPulse)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_infoC.T_pulse) == EC_str & string(xy_opto_infoC.Geno) == geno);
       
        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_optoC(all_ani_EC, :));
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
    all_ani = find(xy_opto_infoC.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_opto_infoC.T_pulse)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_opto_infoC.T_pulse) == EC_str & xy_opto_infoC.Animal == ani);
       
        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            av_resp = mean(xy_optoC(all_ani_EC, :));
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



