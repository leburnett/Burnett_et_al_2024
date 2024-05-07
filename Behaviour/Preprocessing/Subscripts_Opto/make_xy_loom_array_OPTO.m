function make_xy_loom_array_OPTO()

% Modified from make_xy_loom_array/ make_xy_analysis  - for OPTO. 
% 25/05/21

global exp_name dir3_path 

if exist(strcat(dir3_path, '/', 'XY_OPTO_ARRAY_INFO.mat'), 'file')
    load(strcat(dir3_path, '/', 'XY_OPTO_ARRAY_INFO.mat'),'xy_opto_info')
end 

if exist(strcat(dir3_path, '/', 'XY_OPTO_ARRAY.mat'), 'file')
    load(strcat(dir3_path, '/', 'XY_OPTO_ARRAY.mat'),'xy_opto', 'XY_OPTO')
end 

% Make empty cell array in which to enter data. 
if ~exist('xy_opto_info', 'var')
    xy_opto_info = table(); 
    current_size = 0;
    xy_opto = []; 
    XY_OPTO= []; 
else 
     current_size = length(xy_opto_info.Date); 
end 

%%

    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    load(strcat('LOGDATA_', exp_name, '.mat'), 'log_data'); 
%     load(strcat('xy_array.mat'), 'xy_array');

    dt = cell2mat(Info(:, [3,4]));
    xy_array = zeros(numel(dt(:,1)), 1);
    
     for k = 2:numel(dt(:,1))
        A = pdist([dt(k-1,1), dt(k-1,2); dt((k),1), dt((k),2)]);
        xy_array(k,1) = A*(32/416)*60; % (cm/s)
     end
     
%      speed_threshold = 150; %cm/s 
     xy_array(:,1) = movmean(xy_array(:,1), 5);

    date = exp_name(1:6);
    animal = exp_name(8:13);
    exp = exp_name(15:end);
    
    NumPulses = log_data.NumPulses(2);
    T_pulse = log_data.T_pulse(2);
    FreqPulse = log_data.FreqPulse(2);
    EC = log_data.EC(2);
   
%     % MANUALLY INPUT - MUST CHANGE IS LIGHT STIMULUS IS DIFFERENT.
%     length_of_stimulus = 1; % Time in seconds that the light was on. 
%     rows_light_on = length_of_stimulus*60;  
        
     opto_rows = find(cell2mat(Info(:,2))==1);
     num_looms = numel(opto_rows);
     
     num_images = length(Info);

     XY_opto = [];
     XY_opto2 = []; 
     
  for i = 1:num_looms
         
         if isempty(xy_opto_info)
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
         
         %% Make xy_loom_info array 
         
         xy_opto_info.Date{p} = date;
         xy_opto_info.Animal{p} = animal;
         xy_opto_info.Exp{p} = exp;
         xy_opto_info.Loom{p} = opto;
         xy_opto_info.LoomRow{p} = opto_rows(i,1);
         xy_opto_info.NumPulses{p} = NumPulses;
         xy_opto_info.T_pulse{p} = T_pulse;
         xy_opto_info.FreqPulse{p} = FreqPulse;
         xy_opto_info.EC{p} = EC;
         
%          het_animals = ["MJ0580", "MJ0582", "MJ1186", "MJ1347", "MJ1415", "MJ1417", "MJ1579", "MJ1584", "MJ1585", "MJ1586", "MJ1964", "MJ1970", "MJ2142"];
         het_animals = ["MJ2371", "MJ2552", "MJ2547", "MJ2549"];
         
         if contains(xy_opto_info.Animal{p}, het_animals)
             xy_opto_info.Geno{p} = 2; %"het";
         else
             xy_opto_info.Geno{p} = 1; % "wt";
         end
         
         %% Speed array
         
          XY_opto = xy_array(row:row_end,1); 
          
          M = movmean(XY_opto, 5); %col 5 = speed! Movmean over 83ms.
    
          if numel(M) < 1200 && row ~=1
              diff_val = 1200- numel(M); 
              M(numel(M)+1:1200)= zeros(1,diff_val);
          elseif numel(M) < 1200 && row ==1
              diff_val = 1200 - numel(M); 
              added_zeros = zeros(1,diff_val)'; 
              M = vertcat(added_zeros, M);
          end
          
          % Add XY_LOOM to xy_loom 
          xy_opto = vertcat(xy_opto, M');
          
          %% Position array 
          
          XY_opto2 = cell2mat(Info(row:row_end,3:4))'; % col 3 and 4 = adjusted x and y. 
         
          if numel(XY_opto2(1,:)) < 1200 && row ~=1
              diff = 1200- numel(XY_opto2(1,:)); 
              num_orig = numel(XY_opto2(1,:));
              XY_opto2(1,num_orig+1:1200)= zeros(1,diff); 
              XY_opto2(2,num_orig+1:1200)= zeros(1,diff);
          elseif numel(XY_opto2(1,:)) < 1200 && row ==1
              diff = 1200- numel(XY_opto2(1,:)); 
              added_zeros2 = zeros(2,diff);
              XY_opto2= [added_zeros2, XY_opto2]; 
          end 
          
          XY_OPTO = vertcat(XY_OPTO, XY_opto2); % x,y position array. 
  end
 
  
    save(strcat(dir3_path, '/', 'XY_OPTO_ARRAY.mat'),'xy_opto', 'XY_OPTO')
    save(strcat(dir3_path, '/', 'XY_OPTO_ARRAY_INFO.mat'),'xy_opto_info')

    
end        

