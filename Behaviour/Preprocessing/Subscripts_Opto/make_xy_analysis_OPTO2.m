function make_xy_analysis_OPTO2()

 %Created by Burnett 015/10/20 

% Run through the folders, if they had a loom, make an - 'xy_loom' array
% but with LONGER beforehand. 

% 1 second light stimulation - 10ms pulses - 10Hz - different light intensities. 

global exp_name dir3_path inpath

if exist(strcat(dir3_path, '/', 'XY_analysis.mat'), 'file')
    load(strcat(dir3_path, '/', 'XY_analysis.mat'),'xy_analysis')
end 

if exist(strcat(dir3_path, '/', 'PeriLoom_TrialSpeed.mat'), 'file')
    load(strcat(dir3_path, '/', 'PeriLoom_TrialSpeed.mat'),'xy_loom', 'peri_acc')
end 

% Make empty cell array in which to enter data. 
if ~exist('xy_analysis', 'var')
    xy_analysis = table(); 
    current_size = 0;
    xy_loom = []; 
    peri_acc = [];
else 
     current_size = length(xy_analysis.Date); 
end 

length_of_stimulus = 1; % Time in seconds that the light was on. 
rows_light_on = length_of_stimulus*60; 

if exist(strcat('INFO_', exp_name,'.mat'), 'file') 
    
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    load(strcat('LOGDATA_', exp_name, '.mat'), 'log_data');
    
    date = exp_name(1:6);
    animal = exp_name(8:13);
    exp = exp_name(15:end);
    
    NumPulses = log_data.NumPulses(2);
    T_pulse = log_data.T_pulse(2);
    FreqPulse = log_data.FreqPulse(2);
    EC = log_data.EC(2);
    
     opto_rows = find(cell2mat(Info(:,2))==1);
     num_looms = numel(opto_rows);
     
%      if num_looms ~= 6
     % Add the first loom! 
%      opto_rows(2:num_looms+1,1) = opto_rows(1:num_looms,1);
% %      opto_rows(1,1) = opto_rows(2,1)-1803; 
%      num_looms = numel(opto_rows);
% %      end 
    
     num_images = length(Info);
      

  for i = 1:num_looms
         if isempty(xy_analysis)
             p = 1; 
         else 
             p =  current_size + i ; %which row to enter data into in the table.
         end 
         
         opto = char(strcat('0',string(i)));
        
         row_light_on = opto_rows(i,1); 
         row = row_light_on - 599; % Start 10s (599 frames) before the loom stimulus starts. Row light on  == 600 
         row_end = row_light_on + 600;  % end 10s after loom starts. 
         
         if row<0 
             row = 1;
         end 
         
         if row_end > num_images
            row_end = num_images;
         end 
         
         %Make array of speed 
         xy_array = [];   
         for k = 2:num_images
             A = pdist([Info{k-1,3}, Info{k-1,4}; Info{k,3}, Info{k,4}]);
             xy_array(k,1) = A*(32/416)*60; %Cm/s
         end 
%          xy_array = filloutliers(xy_array, 'linear');
         xy_array = movmean(xy_array, 5); 
         
         %Make array of acceleration. 
         ACC = diff(xy_array(:,1));
         
         max_val = max(xy_array);
         height_rect = max_val + (max_val/10);
    
       
         %%  Make 'peri-stimulus' Speed array 
         
          XY_loom = [];
          XY_loom = xy_array(row:row_end,1); 
          
          M = movmean(XY_loom, 5); %col 5 = speed! Movmean over 83ms.
    
          if numel(M) < 1200 && row ~=1
              diff_val = 1200- numel(M); 
              M(numel(M)+1:1200)= zeros(1,diff_val);
          elseif numel(M) < 1200 && row ==1
              diff_val = 1200 - numel(M); 
              added_zeros = zeros(1,diff_val)'; 
              M = vertcat(added_zeros, M);
          end
          
          % Add XY_LOOM to xy_loom 
          xy_loom = vertcat(xy_loom, M'); 
          
          %% Make 'peri-stimulus' ACC array 
          PERI_ACC = [];
          PERI_ACC = ACC(row:row_end,1); 
          
          M2 = movmean(PERI_ACC, 5); %col 5 = speed! Movmean over 83ms.
    
          if numel(M2) < 1200 && row ~=1
              diff_val = 1200- numel(M2); 
              M2(numel(M2)+1:1200)= zeros(1,diff_val);
          elseif numel(M2) < 1200 && row ==1
              diff_val = 1200 - numel(M2); 
              added_zeros = zeros(1,diff_val)'; 
              M2 = vertcat(added_zeros, M2);
          end
          
          % Add XY_LOOM to xy_loom 
          peri_acc = vertcat(peri_acc, M2'); 
          
          
          %% Make xy_analysis table. 
          
            xy_analysis.Date{p} = date;
            xy_analysis.Animal{p} = animal;
            xy_analysis.Exp{p} = exp;
            xy_analysis.Loom{p} = opto;
            xy_analysis.LoomRow{p} = opto_rows(i,1); 

          
          %% Baseline speed
          
          % Speed in the 10s before
          mean_10s = mean(M(1:599)); 
          xy_analysis.Mean_10s{p} = mean_10s;
          
          % Speed in the 5s before
          mean_5s = mean(M(300:599));
          xy_analysis.Mean_5s{p} = mean_5s;
          
          %% Speed during and after stimulation
          
          % Speed in the 5s DURING stimulation
          mean_during = mean(M(600:660));
          xy_analysis.Mean_During{p} = mean_during;
          
          % Speed in the 5s after stimulation
          mean_after5 = mean(M(661:961));
          xy_analysis.Mean_After5{p} = mean_after5;
          
           % Speed in the 10s after stimulation
          mean_after10 = mean(M(661:end));
          xy_analysis.Mean_After10{p} = mean_after10;
              

          %% Maxsp
          
          % MaxSpSpeed in the 5s before
          maxsp_5s_pre = max(M(300:599));
          xy_analysis.MaxSp_5s_PRE{p} = maxsp_5s_pre;
         
          % MaxSpeed in the 5s DURING stimulation
          maxsp_5s_during = max(M(600:660));
          xy_analysis.MaxSp_5s_LIGHT{p} = maxsp_5s_during;
          
          % MaxSpeed in the 5s after stimulation
          maxsp_5s_post = max(M(661:961));
          xy_analysis.MaxSp_5s_POST{p} = maxsp_5s_post;
          
  
           %% LSI
          
          % Maxspeed during / baseline 5s
          LSI_during = maxsp_5s_during/mean_5s; 
          
          % Did the mouse change speed wrt opto stim? Was MaxSpeed during 4
          % x higher than average before.     
          if LSI_during > 4
              xy_analysis.ChangeSpeed{p} = 1;
          else 
              xy_analysis.ChangeSpeed{p} = 0;
          end 
          
          xy_analysis.LSI_during{p} = LSI_during;

          % Maxspeed during / speed 5s post
          LSI_after = maxsp_5s_during/mean_after5; 
          xy_analysis.LSI_after{p} = LSI_after;
          
          % speed 5s pre/ speed 5s post
          LSI_speed_change = mean_5s/mean_after5; 
          xy_analysis.LSI_speed_change{p} = LSI_speed_change;
          
          
         
          %% Acc - DURING LOOM
        
          % Did the mouse change acceleration DURING OPTO STIM?
          rows_acc = find((M2(600:600+rows_light_on))>10); % Positive change in acceleration.    
          
          if isempty(rows_acc)
              xy_analysis.ChangeAcc{p} = 0;
          else 
              xy_analysis.ChangeAcc{p} = 1;
          end 
          
            
          % Max Acc in the 5s before
          maxAcc_5s_pre = max(M2(300:599));
          xy_analysis.MaxAcc_5s_PRE{p} = maxAcc_5s_pre;
         
          % Max Acc in the 5s DURING+ After stimulation
          maxacc_5s_during = max(M2(600:900));
          xy_analysis.MaxAcc_5s_LIGHT{p} = maxacc_5s_during;
        
  
          
          %% Freezing
           
          % Did the mouse freeze during opto presentation?
         
           rows_frozen = find(M(600:900)<1.5); 
           fr_row_diff = diff(rows_frozen(:,1))';
           seq = ones(1,30); % for at least 1/2 s
           Index = strfind(fr_row_diff, seq);
           
           if isempty(Index)
               xy_analysis.Freeze{p} = 0;
               xy_analysis.T2Fr{p} = 0;
               xy_analysis.TFr{p} = 0;
           else 
               xy_analysis.Freeze{p} = 1;
               row_freeze = Index(1);
               %t2fr
               T2Fr = row_freeze/60; 
               xy_analysis.T2Fr{p} = T2Fr;  
               
               rows_more5 = find(fr_row_diff(row_freeze:end)>5);
               if ~isempty(rows_more5)
                    end_fr_row = rows_more5(1);
               else
                   end_fr_row =300; % the end row.  
               end 
               time_frozen = (end_fr_row-row_freeze)/60;
               xy_analysis.TFr{p} = time_frozen;
               
           end 
               

       
          
          %% Distance and Speed
          shelter_size = 65;
          shelterC = [330, 325];
          sz_ratio = 32/416; 
          
          % Distance from the EDGE of the shelter.
          dbox_thresh = shelter_size*sz_ratio; % 90pixels from centre ~6cm.
          for i = 1:num_images
              D = pdist([Info{i,3}, Info{i,4}; shelterC(1),shelterC(2)]);
              Q2(i) = (D*sz_ratio)-dbox_thresh; %distance from centre of the box. - cm.
          end
          Q2 = Q2(row:row_end);
          Q = movmean(Q2, 5);
%           Q = Q2(row:row_end);

          % Distance from shelter when loom happened. 
          dist_at_loom_start = Q(600);
          xy_analysis.DShelt_Start{p} = dist_at_loom_start;
%           data{i,7} = dist_at_loom_start;
          
          in_shelter = find(Q(600:end)<0); 
         
          if isempty(in_shelter)
              xy_analysis.ReturnToShelter{p} = 0;
              xy_analysis.TimeToShelter{p} = 0;
              xy_analysis.SpeedToShelter{p} = 0;
%               data{i,8} = 0; 
%               data{i,9} = 0;
%               data{i,10} = 0;
              row_shelter = 1200;
               
          else
              xy_analysis.ReturnToShelter{p} = 1;

              rows_to_shelter = in_shelter(1); % +600 for actual row.
              row_shelter = rows_to_shelter+600; 
          
               % Time to shelter
               t_2_shelter = rows_to_shelter/60; 
%                data{i,9} = t_2_shelter;
               xy_analysis.TimeToShelter{p} = t_2_shelter;

               
               % Speed loom to shelter. 
               speed_2_shelter = dist_at_loom_start/t_2_shelter;
%                data{i,10} =speed_2_shelter;
               xy_analysis.SpeedToShelter{p} = speed_2_shelter;              
          end 
          
         
          
          %% Maxsp
          
          % Maxspeed - loom to shelter
          max_sp_during_escape = max(M(600:row_shelter));
          xy_analysis.MaxSpEscape{p} = max_sp_during_escape; 
%           data{i,11} =max_sp_during_escape;
          
          % Time to MaxSp
          t_max = find(M == max_sp_during_escape); 
          
          if numel(t_max)>1
              t_max = t_max(1);
          end 
          
          t_max = (t_max - 600)/60; %time from loom til max speed in s.
          xy_analysis.TimeToMaxSp{p} = t_max;
%           data{i,12} =t_max;
          
          % Maxspeed/ baseline 10s
          norm_maxsp_10 = max_sp_during_escape/mean_10s; 
          xy_analysis.MaxSpNorm10{p} = norm_maxsp_10;
%           data{i,13} =norm_maxsp_10;
          
          % Maxspeed/ baseline 5s
          norm_maxsp_5 = max_sp_during_escape/mean_5s; 
          xy_analysis.MaxSpNorm5{p} = norm_maxsp_5;
%           data{i,14} =norm_maxsp_5;
          
          %% Acc - DURING LOOM
          
          % Did the mouse change acceleration DURING loom presentation?
          rows_acc = find((ACC(600:825))>10); % Positive change in acceleration. 
          rows_dec = find((ACC(600:825))<-10);   
          
          if isempty(rows_acc)
              xy_analysis.ChangeAcc{p} = 0;
              xy_analysis.TimeToAccChange{p} = 0;   
          else 
              xy_analysis.ChangeAcc{p} = 1;
              %When did it first start to acc?
              first_start_acc = rows_acc(1)-1;
              t_start_acc = first_start_acc/60;
              xy_analysis.TimeToAccChange{p} = t_start_acc;
          end 
          
           if isempty(rows_dec)
              xy_analysis.DecAcc{p} = 0;
              xy_analysis.TimeToDecAcc{p} = 0;   
          else 
              xy_analysis.DecAcc{p} = 1;
              %When did it first start to acc?
              first_start_dec = rows_dec(1)-1;
              t_start_dec = first_start_dec/60;
              xy_analysis.TimeToDecAcc{p} = t_start_dec;
          end 
          
           % How fast did it accelerate?
           max_acc = max(ACC(600:825));
           xy_analysis.MaxAcc{p} = max_acc;
           
           % Time to Max Acceleration 
           row_max = find(ACC==max_acc);
           time2maxACC = (row_max - 600)/60; 
           xy_analysis.TimeToMaxAcc{p} = time2maxACC; 
           

          
          %% Freezing
           
          % Did the mouse freeze during loom presentation?
          if row_shelter<825 
              rows_frozen = find(M(600:row_shelter)<1.5);
          else 
              rows_frozen = find(M(600:825)<1.5); 
          end 
              
          if isempty(rows_frozen)
              xy_analysis.Freeze{p} = 0;
              xy_analysis.TimeToFreeze{p} = 0;
              xy_analysis.TimeFrozen{p} = 0;
%               data{i,21} = 0; 
%               data{i,22} = 0;
%               data{i,23} = 0;
          else
              xy_analysis.Freeze{p} = 1;
%               data{i,21} = 1; 
          % When did it start to freeze? 
                rowfr = rows_frozen(1)-1; 
                t_start_freeze = rowfr/60; 
                xy_analysis.TimeToFreeze{p} = t_start_freeze;
%                 data{i,22} = t_start_freeze; 

          % How long did it freeze for?
                row_when_freezes = rowfr + 600; 
                num_rows_frozen = find(M(row_when_freezes:end)<1.5);
                n1 = numel(num_rows_frozen);
                num_rows_frozen(1,2) = 1;
                    for q = 2:n1
                        num_rows_frozen(q,2) = num_rows_frozen(q,1)-num_rows_frozen(q-1,1); 
                    end 
                when_stops_freezing = find(num_rows_frozen(:,2)>15); %CHANGED THIS TO 15 - was WHEN ~=1. (1/4s)
                
                if isempty(when_stops_freezing) %frozen for rest of recording. 
                    rows_frozen2 = 1200 - row_when_freezes; 
                    time_frozen = rows_frozen2/60; 
                else 
%                     rowstopped = when_stops_freezing(1)+row_when_freezes; 
                    time_frozen = when_stops_freezing(1)/60; 
                end 
               % How long mouse is frozen for. 
               xy_analysis.TimeFrozen{p} = time_frozen;
%                 data{i,23} = time_frozen;
          
          end 
          
          %%
           

           
           xy_analysis.NumPulses{p} = NumPulses;
           xy_analysis.T_pulse{p} = T_pulse;
           xy_analysis.FreqPulse{p} = FreqPulse;
           xy_analysis.EC{p} = EC;

%            het_animals = ["MJ0580", "MJ0582", "MJ1186", "MJ1347", "MJ1415", "MJ1417", "MJ1579", "MJ1584", "MJ1585", "MJ1586", "MJ1964", "MJ1970", "MJ2142"];
            het_animals  =["MJ2371", "MJ2552", "MJ2547", "MJ2549"];

           if contains(xy_analysis.Animal{p}, het_animals)
               xy_analysis.Geno{p} = "het";
           else
               xy_analysis.Geno{p} = "wt";
           end

  
  end
  
   save(strcat(dir3_path, '/', 'XY_analysis.mat'), 'xy_analysis') 
   save(strcat(dir3_path,'/', 'PeriLoom_TrialSpeed.mat'), 'xy_loom', 'peri_acc') 

end 

end 