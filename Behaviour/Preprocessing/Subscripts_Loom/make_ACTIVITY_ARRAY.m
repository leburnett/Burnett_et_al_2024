function make_ACTIVITY_ARRAY()

% Created by Burnett 25/03/20

global exp_name output_folder %track_req 

    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    num_bas = numel(Info(:,1));
    load(strcat('XY_array_', exp_name, '.mat'), 'xy_array');
    
%% Activity Array - 

    % Col 1 - Total # rows for entire experiment. 
    % COl 2 - Total time for experiment (s)
    % Col 3 - Total distance covered - (100ms time bins)  
    % Col 4 - Average Speed - (100ms time bins)  - greater than 2 cm/s and
    % less than 40cm/s. 
    % Col 5 - Maximum speed across entire experiment.(100ms time bins)
    % Col 6 - Percentage of time that speed < 2cm/s 
    % Col 7 - Percentage of time mouse spends OUTSIDE of the box.
    % Col 8 - Percentage of time mouse spends in the CENTRE of the arena.
    % Col 9 - Percentage of time mouse spends at the EDGE of the arena.
    
    % Values for FIRST 10 seconds of recording 
    % Col 10 - Total Dist 
    % Col 11 - Av Speed 
    % Col 12 - Max Speed 
    % Col 13 - Percentage of time that speed < 2cm/s 
    % Col 14 - Percentage of time mouse spends OUTSIDE of the shelter.
    % Col 15 - Percentage of time mouse spends in the CENTRE of the arena.
    % Col 16 - Percentage of time mouse spends at the EDGE of the arena. 
    
    bin_window = 6; %average over 100ms.  - for moving mean 
    activity_array = zeros(1,16);
    
    %Col - 1 
    activity_array(1,1) = num_bas; 
    activity_array(1,2) = num_bas/60; %time in seconds. 
    
    xy_100ms = movmean(xy_array(:,5), bin_window); %xy_array col 5 is already a movmean average. 
    
    activity_array(1,3) = sum(xy_100ms)/100; %Distance in m over the entire recording. 
    
    xy_mid = xy_100ms(xy_100ms> 2 & xy_100ms<40); 
    activity_array(1,4) = mean(xy_mid); %average speed not slow and not really fast. 
    
    activity_array(1,5) = max(xy_100ms);
    
    % Percentage of time for which speed is <2cm/s.
    P = find(xy_100ms<2); 
      if isempty(P)
             activity_array(1,6) = 0;
      else
            activity_array(1,6) = (numel(P)/num_bas)*100;
      end 
      
     % The percentage of time for which the mouse is outside of the box. 
      A = find(xy_array(:,8)>0);
      activity_array(1,7)= (numel(A)/num_bas)*100; % Percent of time.  
      
      % The percentage of of total time for which the mouse is in the CENTRE of the arena
       B = find(xy_array(:,7)<7); % Within 7cm radius of centre of arena. 
      activity_array(1,8)= (numel(B)/num_bas)*100; % Percent of time.  
      
      % The percentage of of total time for which the mouse is at the EDGE of the arena 
      C = find(xy_array(:,7)>14 & xy_array(:,8)>0); %Further than 14cm from the centre of the arena. 
      activity_array(1,9)= (numel(C)/num_bas)*100; % Percent of time. 

      % % % Looking at only the First 600 frames of the video - 10s - the
      % 'acclim' time before the looms. 
      
     xy_100ms_600 = movmean(xy_array(1:600,5), bin_window); 
     activity_array(1,10) = sum(xy_100ms_600); 
     activity_array(1,11) = mean(xy_100ms_600);
     activity_array(1,12) = max(xy_100ms_600); 
     
      % The percentage of of total time for which speed <2
      Q = find(xy_100ms_600<2); 
      if isempty(Q)
             activity_array(1,13) = 0;
      else
            activity_array(1,13) = (numel(Q)/600)*100;
      end
      
     % The percentage of of total time for which the mouse is outside the
     % shelter. 
      A2 = find(xy_array(1:600,8)>0);
      activity_array(1,14)=(numel(A2)/600)*100; %Time in seconds. 
      
      % The percentage of of total time for which the mouse is in the CENTRE of the arena
       B2 = find(xy_array(1:600,7)<7); % Within 7cm radius of centre of arena. 
      activity_array(1,15)= (numel(B2)/600)*100; %Time in seconds. 
      
      % The percentage of of total time for which the mouse is at the EDGE of the arena 
      C2 = find(xy_array(1:600,7)>14 & xy_array(1:600,8)>0); %Further than 14cm from the centre of the arena. 
      activity_array(1,16)= (numel(C2)/600)*100; %Time in seconds.
      
      
      %Save once locally within the experiment folder. 
      save(strcat('activity_array_', exp_name, '.mat'), 'activity_array');
     
    % Save in 'Exp Project' folder. 
    activity_folder = strcat(output_folder, 'SUMMARIES\ACTIVITY\ACTIVITY_ARRAYS\');
    
    if ~exist(activity_folder,'dir')
        mkdir(activity_folder)
    end

    save(fullfile(strcat(activity_folder,'ACTIVITY_ARRAY_', exp_name,'.mat')), 'activity_array');  
%         save(fullfile(strcat('C:\Data_analysis\DATA\2003_SERT\SUMMARIES\ACTIVITY_ARRAYS', 'activity_array_', exp_name, '.mat')), 'activity_array');
        
    clearvars
        
end 
