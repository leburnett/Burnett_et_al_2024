function make_xy_analysis()

%Created by Burnett 09/09/20 

% Run through the folders, if they had a loom, make an - 'xy_loom' array
% but with LONGER beforehand. 
global exp_name dir3_path %xy_analysis

if exist(strcat(dir3_path,'\', '201013_Setd5_C1_XY_analysis.mat'), 'file')
    load(strcat(dir3_path,'\', '201013_Setd5_C1_XY_analysis.mat'),'xy_analysis')
end 

if exist(strcat(dir3_path,'\', '201013_Setd5_C1_PeriLoom_TrialSpeed.mat'), 'file')
    load(strcat(dir3_path,'\', '201013_Setd5_C1_PeriLoom_TrialSpeed.mat'),'xy_loom')
end 

% Make empty cell array in which to enter data. 
if ~exist('xy_analysis', 'var')
    xy_analysis = table(); 
    current_size = 0;
     % Make longer and smoothed 'xy_loom' array as well
     xy_loom = []; 
else 
     current_size = length(xy_analysis.Date); 
end 
% data = {}; 

if exist(strcat('ALL_LOOM_ROWS_', exp_name,'.mat'), 'file') 

    load((strcat('ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS'); 
    load(strcat('XY_array_', exp_name, '.mat'), 'xy_array');
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    
    date = exp_name(1:6);
    animal = exp_name(8:13);
    exp = exp_name(15:end);
    
     bin_window = 5; %average over 83 ms. 
     loom_rows = ALL_LOOM_ROWS(1,:)'; 
     num_looms = numel(loom_rows);

  for i = 1:num_looms
         if isempty(xy_analysis)
             p = 1; 
         else 
             p =  current_size + i ; %which row to enter data into in the table.
         end 
         
         loom = char(strcat('0',string(i)));
        
         row = loom_rows(i,1) - 599; % Start 10s (599 frames) before the loom stimulus starts 
         row3 = loom_rows(i,1) + 600;  % end 10s after loom starts. 
         
         if row<0 
             row = 1;
         end 
         
         if row3 > length(xy_array)
            row3 = length(xy_array);
         end 
         
         % Speed array 
          XY_loom = [];
          XY_loom = xy_array(row:row3,:); 
          
          M = movmean(XY_loom(:,5), bin_window); %col 5 = speed! Movmean over 83ms.
    
          if numel(M) < 1200 && row ~=1
              diff = 1200- numel(M); 
              M(numel(M)+1:1200)= zeros(1,diff);
          elseif numel(M) < 1200 && row ==1
              diff = 1200 - numel(M); 
              added_zeros = zeros(1,diff)'; 
              M = vertcat(added_zeros, M);
          end
          
          M2 = smooth(M, 6);
          
          % Add XY_LOOM to xy_loom 
          xy_loom = vertcat(xy_loom, M2); 
          
          Q = XY_loom(:,8); %Distance from edge of shelter. 
          Q = smooth(Q, 6); 
          
          if numel(Q) < 1200 && row ~=1
              diff = 1200- numel(Q); 
              Q(numel(Q)+1:1200)= zeros(1,diff);
          elseif numel(Q) < 1200 && row ==1
              diff = 1200 - numel(Q); 
              added_zeros = zeros(1,diff)'; 
              Q = vertcat(added_zeros, Q);
          end 
          
          % Acceleration
          ACC = movmean(XY_loom(:,6), bin_window); % Acceleration. 
          ACC = smooth(ACC, bin_window);
          
          if numel(ACC) < 1200 && row ~=1
              diff = 1200- numel(ACC); 
              ACC(numel(ACC)+1:1200)= zeros(1,diff);
          elseif numel(ACC) < 1200 && row ==1
              diff = 1200 - numel(ACC); 
              added_zeros = zeros(1,diff)'; 
              ACC = vertcat(added_zeros, ACC);
          end 
          
            xy_analysis.Date{p} = date;
            xy_analysis.Animal{p} = animal;
            xy_analysis.Exp{p} = exp;
            xy_analysis.Loom{p} = loom;
            xy_analysis.LoomRow{p} = loom_rows(i,1); 

%           data{i,1} = {date};
%           data{i,2} = {animal};
%           data{i,3} = {exp};
%           data{i,4} = {loom};
          
          %% Baseline speed
          
          % Speed in the 10s before
          mean_10s = mean(M(1:600)); 
          xy_analysis.Mean_10s{p} = mean_10s;
%           data{i,5} = mean_10s; 
          
          % Speed in the 5s before
          mean_5s = mean(M(300:600));
          xy_analysis.Mean_5s{p} = mean_5s;
%           data{i,6} = mean_5s; 
          
          %% Distance and Speed
          
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
           
           % % % % % 
           
            % How fast did it deccelerate?
           max_dec = min(ACC(600:825));
           xy_analysis.MaxDec{p} = max_dec;
           
           % Time to Max Decceleration 
           row_min = find(ACC==max_dec);
           time2maxDEC = (row_min - 600)/60; 
           xy_analysis.TimeToMaxDec{p} = time2maxDEC; 

           %% ACC - DURING ESCAPE - til back in shelter
           
%            if xy_analysis.ReturnToShelter{p} == 1 %data{i,8} == 1 %returned to shelter. 
%                
%             % Did the mouse change acceleration end of loom til return to shelter?
%             rows_acc_escape = find(ACC(600:row_shelter)>5); 
%           
%             if isempty(rows_acc_escape)
%                 xy_analysis.AccPostLoom{p} = 0;
%                 xy_analysis.TimeAccPostLoom{p} = 0;
% %                 data{i,18} = 0;
% %                 data{i,19} = 0;
%             else 
%                 xy_analysis.AccPostLoom{p} = 1;
% %                 data{i,18} = 1;
%                 %When did it first start to acc?
%                 first_start_acc_escape = rows_acc_escape(1)-1;
%                 t_start_acc_escape = first_start_acc_escape/60;
%                 xy_analysis.TimeAccPostLoom{p} = t_start_acc_escape;
% %                 data{i,19} = t_start_acc_escape; 
%             end 
%           
%            % How fast did it accelerate?
%            max_acc_escape = max(ACC(600:row_shelter));
%            xy_analysis.MaxAccPostLoom{p} = max_acc_escape;
% %            data{i,20} = max_acc_escape; 
%            
%            end 
          
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
          
          
          % Col 15 -  CONTRAST
                
%                 contrast = [Info{loom_rows(i,1), 10}];
%                 
% %                 if contrast == [0 0 0] %High contrast
% %                     contrast_type = 1;
% %                 elseif contrast == [0 35 70] %Low contrast
% %                     contrast_type = 2;
% %                 elseif contrast == [0 100 0] %GREEN only
% %                     contrast_type = 3;
% %                 elseif contrast == [0 0 100] %UV only
% %                     contrast_type = 4;
% %                 else
% %                     contrast_type = 5;
% %                 end
%           
%                 % OLD - diff contrasts
%                         if contrast == [0 0 0]
%                             contrast_type = 1;
%                         elseif contrast == [0 20 40]
%                             contrast_type = 2;
%                         elseif contrast == [0 35 70]
%                             contrast_type = 3;
%                         elseif contrast == [0 40 80]
%                             contrast_type = 4;
%                         else
%                             contrast_type = 5;
%                         end
%                 
%                 xy_analysis.Contrast{p} = contrast_type;
          
          
%           %% Geno - setd5
          if xy_analysis.Animal{p} == "MJ0245" || xy_analysis.Animal{p} == "GN4244" || xy_analysis.Animal{p} == "GN4473" || xy_analysis.Animal{p} == "GN6560"  || xy_analysis.Animal{p} == "GN6562"  || xy_analysis.Animal{p} == "GN6610"|| xy_analysis.Animal{p} == "GN7269"|| xy_analysis.Animal{p} == "GN7476"|| xy_analysis.Animal{p} == "GN7614"|| xy_analysis.Animal{p} == "GN7790" 
              xy_analysis.Geno{p} = "het";
          else 
              xy_analysis.Geno{p} = "wt";
          end 

        %% Geno - Cul3 
%           if xy_analysis.Animal{p} == "GN2375" || xy_analysis.Animal{p} == "GN2377" || xy_analysis.Animal{p} == "GN2382" || xy_analysis.Animal{p} == "GN2628"  || xy_analysis.Animal{p} == "GN2637"  || xy_analysis.Animal{p} == "GN2754"|| xy_analysis.Animal{p} == "GN2902"|| xy_analysis.Animal{p} == "GN2901"|| xy_analysis.Animal{p} == "GN2900"|| xy_analysis.Animal{p} == "GN2832" || xy_analysis.Animal{p} == "GN2829"
%               xy_analysis.Geno{p} = "het";
%           else 
%               xy_analysis.Geno{p} = "wt";
%           end 
          

 
        save(strcat(dir3_path,'\', '201013_Setd5_C1_XY_analysis.mat'), 'xy_analysis') 
        save(strcat(dir3_path,'\', '201013_Setd5_C1_PeriLoom_TrialSpeed.mat'), 'xy_loom') 
          % For how long did the mouse freeze once it returned to the
          % shelter? 
   
%           xy_loom = vertcat(xy_loom, M); %speed array 
  end
% xy_analysis = [xy_analysis; data]; 
% clearvars -except xy_analysis
  

end 

end 