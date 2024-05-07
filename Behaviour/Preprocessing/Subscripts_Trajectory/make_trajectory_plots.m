function make_trajectory_plots()
% Function that makes animated trajectory plots. Other 'static' plots are
% made in 'make_EXP_SUM_TABLE.m'

%% Can add 'gif' function into these scripts to generate gifs of the mouse trajectory. 

% Different versions for dot stimuli and moving stripes.  
% 1- Animated plot. Mouse centroid and dot.
 
 global exp_name 

load(strcat('EXP_SUM_TAB', exp_name, '.mat'), 'EXP_SUMMARY_TABLE');
load(strcat('ALLDATA_', exp_name, '.mat'), 'ALLDATA');

DATA_size = size(ALLDATA);
DATA_size = DATA_size(1);

load(strcat('STIMROWS_', exp_name, '.mat'), 'STIMROWS');

num_stim = size(STIMROWS);
num_stim = num_stim(1);

alpha_zero = 0.004385936788792; %Previously -0.006637070683990;
basler_image_size = 512; %previously 925 

%%  1- Plot animated trajectory with ellipse and centroid
% Include either the dot or add stripes image for moving stripes. 


%%%%%%%%%%%%% FOR DOT STIMULI.
%%% This animated plot draws the mouse centroid position as black small dots. 
%%% The dot is drawn as a red dot. 
% 
% ax = axes;
% axis([0 550 0 550]);
% hold on 
% dot_x2 = 0;
% dot_y2 = 0; 
% mouse_x = 0;
% mouse_y = 0;
% 
% h = animatedline('Color', 'k', 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 2);
% h2 = animatedline('LineStyle', 'none', 'Marker', '.', 'Color', 'r');
% p = plot(ax, dot_x2, dot_y2, 'r.', 'MarkerSize', 25);
% q = plot(ax, mouse_x, mouse_y, 'ko', 'MarkerSize', 25);
% 
%  for k = 1:DATA_size
% 
%     %Plot dot
%     dot_x = ALLDATA.dot_x{k};
%     dot_y = ALLDATA.dot_y{k};
% 
%     %Transform dot coords from python. 
%     [alpha1, rad1] = cart2pol (dot_x, dot_y);
%     alpha1 = alpha1 + alpha_zero ;
%     rad1 = rad1*basler_image_size/1000;
%     [dot_x2, dot_y2] = pol2cart(alpha1, rad1);
% 
%   
%     p.XData = dot_x2;
%     p.YData = dot_y2;
%     addpoints(h2, dot_x2, dot_y2);
%     drawnow limitrate
%     pause(0.01)
%     
%      %plot mouse position
% 
%     mouse_x = ALLDATA.Centroid_x(k);
%     mouse_y = ALLDATA.Centroid_y(k);
% 
%     q.XData = mouse_x;
%     q.YData = mouse_y;
%     addpoints(h, mouse_x, mouse_y);
%     drawnow limitrate
%     pause(0.01)
% 
%  end 
% 
% %%%%%%%%% DOT STIMULI WITH ELLIPSE AND CENTROID
% 
% % Animated plot of an ellipse of the mouse body and a red dot. 
% centre = [260 254];
% radius = 232;
% 
%  for k = 1:DATA_size
%      
%     imshow(ones(550));
%     hold on 
% 
%     %Plot dot 
%     dot_x = ALLDATA.dot_x{k};
%     dot_y = ALLDATA.dot_y{k};
%     
%     %Transform dot coords from python. 
% 
%     [alpha1, rad1] = cart2pol (dot_x, dot_y);
%     alpha1 = alpha1 + alpha_zero ;
%     rad1 = rad1*basler_image_size/1000;
%     [dot_x2, dot_y2] = pol2cart(alpha1, rad1);
% 
%      %plot mouse position
% 
%     Xc = ALLDATA.Centroid_x(k);
%     Yc = ALLDATA.Centroid_y(k);   
%     t = linspace(0,2*pi,50);
%     a = ALLDATA.MajorAxis(k);
%     b = ALLDATA.MinorAxis(k);
%     phi = ALLDATA.phi(k);
%     x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%     y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%     
% % % % Plot ellipse
% 
%     plot(dot_x2, dot_y2, 'r.', 'MarkerSize', 5)
%     plot(x,y,'k','Linewidth',0.5)
%     viscircles(centre, radius, 'Color', 'k', 'LineWidth', 0.5);
%     hold off
%     pause(0.01)
% 
%  end
 
%%%%%%%%%%%%%%%% Moving Stripes Stimuli----- 

% Want to work with colours but at the moment. The Blue dot moves to
% different positions for still, clock and counter. 
 
ax = axes;
axis([0 550 0 550]);

hold on 
dot_x2 = 0;
dot_y2 = 0; 
mouse_x = 0;
mouse_y = 0;

h = animatedline('Color', 'k', 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 2);
h2 = animatedline('LineStyle', 'none', 'Marker', '.', 'Color', 'r');
p = plot(ax, dot_x2, dot_y2, 'b.', 'MarkerSize', 30);
q = plot(ax, mouse_x, mouse_y, 'ko', 'MarkerSize', 25);

 for k = 1:DATA_size

    %Plot Stripes as dot. 
   stim_name = ALLDATA.Stimulus{k};
      if stim_name == "Acclim" 
%          colour = 'w';
         dot_x = 25;
         dot_y = 25;
     elseif stim_name == "Still"
%          colour = 'y';
         dot_x = 25;
         dot_y = 100; 
     elseif stim_name == "Clock"
%          colour = 'g';
          dot_x = 25;
         dot_y = 300;
     elseif stim_name == "Counter"
%          colour = 'b'; 
          dot_x = 20;
         dot_y = 450;
     else 
%          colour = 'w';
          dot_x = 25;
         dot_y = 25;
     end 
     
    p.XData = dot_x;
    p.YData = dot_y;
    addpoints(h2, dot_x, dot_y);
    drawnow limitrate
    pause(0.01)

    %plot mouse position
    mouse_x = ALLDATA.Centroid_x(k);
    mouse_y = ALLDATA.Centroid_y(k);

    q.XData = mouse_x;
    q.YData = mouse_y;
    addpoints(h, mouse_x, mouse_y);
    drawnow limitrate
    pause(0.01)

 end 

 


%% Plot end plot - full trajectory for the entire video - just centroid

 %% MAKE PLOT - Trajectory -  1 second bins. 
% 
% Ts = ceil(DATA_size/40); %Basler acquires at 40 fps - averaging 40 frames. 
% Ts_array = (1:40:DATA_size);
% bin_values_position = [];
% 
% for j = 1:Ts
%     if j == 1
%         start_row =1;
%     elseif j>1
%         start_row = Ts_array(j-1);
%     end 
%     if j == Ts
%         end_row = DATA_size;
%     elseif j <Ts
%         end_row = Ts_array(j);
%     end 
%     bin_values_position(j,1) = mean(ALLDATA{start_row:end_row, 12});
%     bin_values_position(j,2) = mean(ALLDATA{start_row:end_row, 13});
%      
% end 


%% Colourful binned full trajectory. 
% Uncomment 'drawnow' to make it animated. 
% 
% figure 
% hold on 
% cmap = hsv(Ts);
% 
% for j = 1: Ts-1
%     x1 = bin_values_position(j,1);
%     x2 = bin_values_position(j+1,1); 
%     y1 = bin_values_position(j,2); 
%     y2 =  bin_values_position(j+1,2);
%  plot([x1 x2],[y1 y2], 'color', cmap(j,:));
% %  drawnow 
% end 
% savefig(gcf, strcat('TRAJPLOT_COL_BIN_FULL_', exp_name));
% 
% 



%% Colourful non-binned full trajectory. 
% Uncomment 'drawnow' to make it animated. 
% 
% figure 
% hold on 
% numb = size(ALLDATA);
% numb = numb(1);
% cmap = hsv(numb);
% 
% for j = 1: numb
%     x1 = ALLDATA{j,12};
%     x2 = ALLDATA{j+1,12}; 
%     y1 = ALLDATA{j,13}; 
%     y2 =  ALLDATA{j+1,13};
%  plot([x1 x2],[y1 y2], 'color', cmap(j,:))
% %  drawnow 
% end 



%% Plot end plot for individual stimuli 
% A -  Loop to make and save the plots. 

% 
% %  if isempty(stim_start) || ~exist('stim_start')
% %      stim_start = 1;
% %  end 
% % 
% %  if isempty(stim_stop) || ~exist('stim_stop')
% %      stim_stop = num_stim;
% %  end 
% 
% num_stim = size(log_data);
% num_stim = num_stim(1);
% 
% for i = 1:num_stim-1 
% 
%  stim_start = i;
%  stim_stop = i+1;
%  
% start_row = (STIMROWS{stim_start,2}+1);
% end_row = STIMROWS{stim_stop,2};
% 
% 
% figure 
% hold on 
% numb = size(ALLDATA);
% numb = numb(1);
% cmap = hsv(numb);
% 
% for j = start_row: end_row
%     x1 = ALLDATA{j,12};
%     x2 = ALLDATA{j+1,12}; 
%     y1 = ALLDATA{j,13}; 
%     y2 =  ALLDATA{j+1,13};
%  plot([x1 x2],[y1 y2], 'color', cmap(j,:))
% %  drawnow 
% end 
% stim_name = log_data.ExpName{i};
% savefig(gcf, strcat('TRAJPLOT_', string(i), string(stim_name)))
% close
% 
% end 


%% To visualise the trajectory of the mouse for each stimulus. 
% num_stim = size(log_data);
% num_stim = num_stim(1);
% centre = [260 254];
% radius = 232;
% 
% %Enter stimulus number here. 
% stimulus_number =1; 
% 
% for i = stimulus_number 
% 
%  stim_start = i;
%  stim_stop = i+1;
%  stim_name = log_data.ExpName{i};
%  
%  start_row = (STIMROWS{stim_start,2}+1);
%  end_row = STIMROWS{stim_stop,2};
% 
% figure 
% hold on 
% viscircles(centre, radius, 'Color', 'k', 'LineWidth', 0.3);
% numb = size(ALLDATA);
% numb = numb(1);
% cmap = hsv(numb);
% 
% for j = start_row: end_row
%     x1 = ALLDATA{j,12};
%     x2 = ALLDATA{j+1,12}; 
%     y1 = ALLDATA{j,13}; 
%     y2 =  ALLDATA{j+1,13};
%     plot([x1 x2],[y1 y2], 'color', cmap((j),:))
%     title(stim_name)
%    drawnow 
% end 
% end 




%% PLOTTING TRIGGER THRESHOLD IMAGE
% 
% im = 'img_13188421340226.jpg';
% im = imread(im); 
% im= im(4:521, 23:540);
% imshow(im);
%  boxC = [430, 440];
%  radius= 90; 
%  hold on 
% viscircles(boxC, radius)
% viscircles(boxC, 200, 'Color','b')
% viscircles(boxC, 300, 'Color','g')
% viscircles(boxC, 250, 'Color','m')




end 

