function dynamic_traj(stim_start, stim_stop) %% Plot the trajectory of the mouse and the dot.
 %% Plot the trajectory of the mouse and the dot.
 % Animated graphic of the mouse and the dot position. s
 % Created by Burnett 10/10/19. 
 
 % stim_start can go from 1. 
 % stim_stop can finish at num_stim 
 
 global exp_name
 
 load(strcat('ALLDATA_', exp_name, '.mat'), 'ALLDATA');

DATA_size = size(ALLDATA);
DATA_size = DATA_size(1);

load(strcat('STIMROWS_', exp_name, '.mat'), 'STIMROWS');

num_stim = size(STIMROWS);
num_stim = num_stim(1);
 
 if isempty(stim_start) || ~exist('stim_start')
     stim_start = 1;
 end 
 
 if isempty(stim_stop) || ~exist('stim_stop')
     stim_stop = num_stim;
 end 
 
alpha_zero = 0.004385936788792; %Previously -0.006637070683990;
basler_image_size = 512; %previously 925 

start_row = (STIMROWS{stim_start,2}+1);
end_row = STIMROWS{stim_stop,2};


%% Draw Plot 

ax = axes;
axis([0 550 0 550]);
hold on 

% dot_x2 = 0;
% dot_y2 = 0; 
mouse_x = 0;
mouse_y = 0;

h = animatedline('Color', 'k', 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 2);
% h2 = animatedline('LineStyle', 'none', 'Marker', '.', 'Color', 'r');
% p = plot(ax, dot_x2, dot_y2, 'r.', 'MarkerSize', 25);
q = plot(ax, mouse_x, mouse_y, 'ko', 'MarkerSize', 25);

%% Plot for the entire video. 

 for k = 5000:end_row%start_row:end_row
     
    %Plot dot 
%     dot_x = ALLDATA.dot_x{k};
%     dot_y = ALLDATA.dot_y{k};
%        
%     %Transform dot coords from python. 
%     [alpha1, rad1] = cart2pol (dot_x, dot_y);
%     alpha1 = alpha1 + alpha_zero ;
%     rad1 = rad1*basler_image_size/1000;
%     [dot_x2, dot_y2] = pol2cart(alpha1, rad1);
%    
%     p.XData = dot_x2;
%     p.YData = dot_y2;
%     addpoints(h2, dot_x2, dot_y2);
%     drawnow limitrate
%     pause(0.01)
   
     %plot mouse position
    mouse_x = ALLDATA.Centroid_x(k);
    mouse_y = ALLDATA.Centroid_y(k);
    
    q.XData = mouse_x;
    q.YData = mouse_y;
    addpoints(h, mouse_x, mouse_y);
    drawnow limitrate
    pause(0.01)
 end 


%% Plot for one part of the experiment 

% % global StimNum
% 
% % StimNum = log_data.StimNum; 
% StimNum = 4;
% load(strcat('STIM', string(StimNum), '_DATA.mat'));
% 
% DATA_size = size(stimDATA);
% DATA_size= DATA_size(1);
% 
% %% Draw Plot 
% 
% ax = axes;
% axis([0 550 0 550]);
% hold on 
% h = animatedline('Color', 'k', 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 2);
% h2 = animatedline('LineStyle', 'none', 'Marker', '.', 'Color', 'r');
% p = plot(ax, dot_x2, dot_y2, 'r.', 'MarkerSize', 25);
% q = plot(ax, mouse_x, mouse_y, 'ko', 'MarkerSize', 25);
% 
%  for k = 1:DATA_size
%      
%     %Plot dot 
%     dot_x = stimDATA{k,3};
%     dot_y = stimDATA{k,4};
%        
%     %Transform dot coords from python. 
%     [alpha1, rad1] = cart2pol (dot_x, dot_y);
%     alpha1 = alpha1 + alpha_zero ;
%     rad1 = rad1*basler_image_size/1000;
%     [dot_x2, dot_y2] = pol2cart(alpha1, rad1);
%    
%     p.XData = dot_x2;
%     p.YData = dot_y2;
%     addpoints(h2, dot_x2, dot_y2);
%     drawnow limitrate
%     pause(0.01)
%    
%      %plot mouse position
%     mouse_x = stimDATA{k,11};
%     mouse_y = stimDATA{k,12};
%     
%     q.XData = mouse_x;
%     q.YData = mouse_y;
%     addpoints(h, mouse_x, mouse_y);
%     drawnow limitrate
%     pause(0.01)
%  end 




%%
% 
%  end_row_stim1 = log_data.DATArow(1);
% A = find(rows(:,1)>end_row_stim1);
% end_stim1 = min(A)-1;
% 
% end_row_stim2 = log_data.DATArow(2);
% A = find(rows(:,1)>end_row_stim2);
% end_stim2 = min(A)-1;
% 
% end_row_stim3 = log_data.DATArow(3);
% A = find(rows(:,1)>end_row_stim3);
% end_stim3 = min(A)-1;
% 
% end_row_stim4 = log_data.DATArow(4);
% A = find(rows(:,1)>end_row_stim4);
% end_stim4 = min(A)-1;
% 
% end_row_stim5 = log_data.DATArow(5);
% A = find(rows(:,1)>end_row_stim5);
% end_stim5 = min(A)-1;
 
 
% % for i = 1: 5849
% %     if mouse_coords(i,1)==500.5
% %         mouse_coords(i,1)= NaN; 
% %     end 
% %     if mouse_coords(i,2)==500.5
% %         mouse_coords(i,2)=NaN;
% %     end 
% % end
 
%% Old Scripts


% h = animatedline('Color', 'k', 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 2);
% % h2 = animatedline('LineStyle', 'none', 'Marker', '.', 'Color', 'r');
% 
% for i = 1:1000
% %     x = mouse_coords(i,1);
% %     y = mouse_coords(i,2);
%     x = ORIENT_DATA.Var4(i);
%     y = ORIENT_DATA.Var5(i);
%     addpoints(h, x, y);
%     drawnow
% %     x2 = dot_coords(i,1);
% %     y2 = dot_coords(i,2);
% %     addpoints(h2, dot_coords(i,1), dot_coords(i,2));
% %     drawnow
% end 
% 
% mouse_x = num2cell(ORIENT_DATA.Var4(1:3740));
% mouse_y = num2cell(ORIENT_DATA.Var5(1:3740));


%% 
% 
% DATA_size = size(ORIENT_DATA);
% DATA_size= DATA_size(1);

% dot_x = cell2mat(ORIENT_DATA(:,2));
% dot_y = cell2mat(ORIENT_DATA(:,3));
% mouse_x = cell2mat(ORIENT_DATA(:,4));
% mouse_y = cell2mat(ORIENT_DATA(:,5));
% rows = cell2mat(ORIENT_DATA(:,6));

% for i = 1:DATA_size
%     [alpha1, rad1] = cart2pol (dot_x(i), dot_y(i));
%     alpha1 = alpha1 + alpha_zero ;
%     rad1 = rad1*basler_image_size/1000;
%     [dot_x2(i,1), dot_y2(i,1)] = pol2cart(alpha1, rad1);
% end 
% 
% plot(mouse_x, mouse_y, 'k')
% hold on 
% plot(dot_x2, dot_y2, 'r')
% hold off

% 
% %Green Screen - no dot. 
% plot(mouse_x(1:end_stim1), mouse_y(1:end_stim1), 'k')
% hold on 
% plot(dot_x2(1:end_stim1), dot_y2(1:end_stim1), 'r')
% hold off
% 
% 
% %Moving dot 1
% plot(mouse_x(end_stim1+1:end_stim2), mouse_y(end_stim1+1:end_stim2), 'k')
% hold on 
% plot(dot_x2(end_stim1+1:end_stim2), dot_y2(end_stim1+1:end_stim2), 'r')
% hold off
% 
% % Moving dot 2
% plot(mouse_x(end_stim2+1:end_stim3), mouse_y(end_stim2+1:end_stim3), 'k')
% hold on 
% plot(dot_x2(end_stim2+1:end_stim3), dot_y2(end_stim2+1:end_stim3), 'r')
% hold off
% 
% %Moving dot 3
% plot(mouse_x(end_stim3+1:end_stim4), mouse_y(end_stim3+1:end_stim4), 'k')
% hold on 
% plot(dot_x2(end_stim3+1:end_stim4), dot_y2(end_stim3+1:end_stim4), 'r')
% hold off
% 
% %Moving dot 4
% plot(mouse_x(end_stim4+1:end_stim5), mouse_y(end_stim4+1:end_stim5), 'k')
% hold on 
% plot(dot_x2(end_stim4+1:end_stim5), dot_y2(end_stim4+1:end_stim5), 'r')
% hold off
% 
% %Green 
% plot(mouse_x(end_stim5+1:end), mouse_y(end_stim5+1:end), 'k')
% hold on 
% plot(dot_x2(end_stim5+1:end), dot_y2(end_stim5+1:end), 'r')
% hold off
% 
% 
% 
% x = ORIENT_DATA(i).Var4;
% y = ORIENT_DATA(i).Var5;
% 
% % dot_coords(:,1) = dot_x2;
% % dot_coords(:,2) = dot_y2;
% % find(dot_coords(:,1)>417 & dot_coords(:,1)<421 & dot_coords(:,2)>250 & dot_coords(:,2)<259);
% % 
% % mouse_coords(:,1) = mouse_x;
% % mouse_coords(:,2) = mouse_y;
% % 
