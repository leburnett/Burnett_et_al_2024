function analyse_traj_comp()

% Analyse the directedness of the trajectory - minimum euclidean distance between
% position when loom started and shelter versus actual distance mouse
% travelled on escape back to shelter. 

global exp_name 

% End result will be an array
% Col 1 = row when loom starts.
% Col 2 = row when enters shelter.
% Col 3 = minimum distance from pos of start of loom to shelter.
% Col 4 = actual distance travelled
% Col 5 = distance travelled/ minimum distance.

% dthresh constant. 
% boxC = [430 440]; 
boxC = [430 78]; % y = 518-y
dbox_thresh = 90*(32/518); 

if exist(strcat('loom_row_all_', exp_name, '.mat'))   
    
    % Load loom rows
    load(strcat('loom_row_all_', exp_name, '.mat'), 'loom_row_all')
    n_looms = numel(loom_row_all(1,:));
    
    %Make array
    trajcomp = zeros(n_looms, 5);
    
    % Load response array
    load(strcat('response_array_', exp_name, '.mat'), 'response_array')
    
    % Load xy_array
    load(strcat('XY_array_', exp_name, '.mat'), 'xy_array')
    
    
    for i = 1:n_looms
        
        row1 = loom_row_all(1,i); %row when loom starts;
        row2 = response_array(i,3); %row when mouse enters the shelter.
        
        trajcomp(i,1) = row1;
        
        trajcomp(i,2) = row2;
        
        x_loom = xy_array(row1, 3);
        y_loom = 518 - xy_array(row1, 4);
        
        trajcomp(i,3) = (pdist([x_loom, y_loom; boxC(1), boxC(2)])*(32/518))- dbox_thresh;% minimum distance from shelter edge when loom started
        
        trajcomp(i,4) =  sum(xy_array(row1:row2,5))/60; % add all the distances between frames from the xy_array between these two rows to find the actual distance travelled.
        
        trajcomp(i,5) = trajcomp(i,4) / trajcomp(i,3);
        
    end
    
    save(strcat('traj_comp_', exp_name, '.mat'), 'trajcomp')
    
    %SAVE GLOBALLY AS WELL - NOT JUST IN FOLDER.
    save(fullfile(strcat('C:\Data_analysis\DATA\Setd5\ALL_DAY_SUMMARIES\TRAJ\traj_comp_', exp_name, '.mat')), 'trajcomp')
    
else 
end 
    clearvars
end 


%% Visualise this graphically. 

% figure 
% x = xy_array(row1,3);
% y = xy_array(row1,4);
% plot(x,518-y,'r.', 'MarkerSize', 15)
% axis([0 520 0 520])
% hold on     
% x2 = xy_array(row2,3);
% y2 = xy_array(row2,4);
% plot(x2,518-y2,'k.', 'MarkerSize', 15)
% 
% % Plot the actual trajectory. 
% for q = row1:row2 
%      x = xy_array(q,3);
%      y = 518 - xy_array(q,4); 
%      x2 = xy_array(q+1, 3);
%      y2 = 518 - xy_array(q+1, 4);
%      plot([x, x2],[y,y2],'k')
%      hold on 
% end 
% 
% x = xy_array(row1,3);
% y = xy_array(row1,4);
% plot([x, boxC(1)], [518-y, boxC(2)], 'k:')
% 
% hold on 
% viscircles(boxC, radius, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 0.5)
% title('Computing trajectory directedness')
% savefig(gcf, 'Traj_Directedness.fig')
% close

