function make_DLC_LOOM_array()
% MAKE 'CROPPPED' DLC files. 

% Run through folders. 
% Find the ones with ALL_LOOM_ROWS
% Use 'exp_name' to 'find' the DLC files in pre-specified DLC folders. 

% 'CUT' the DLC arrays - 5s (300) before the loom row and 10s (600) after the loom row. 

% Save in 'individual' loom files '01_Loom_01' etc. 

global exp_name 

if exist(strcat('ALL_LOOM_ROWS_', exp_name,'.mat'), 'file') 

    pos_folder = 'C:\Data_analysis\DATA\2007_Cul3_C2\SUMMARIES\DLC\POS\';
    speed_folder = 'C:\Data_analysis\DATA\2007_Cul3_C2\SUMMARIES\DLC\SPEED\';
    
    pos_saving_folder = 'C:\Data_analysis\DATA\2007_Cul3_C2\SUMMARIES\DLC\LOOMS\POS\';
    speed_saving_folder = 'C:\Data_analysis\DATA\2007_Cul3_C2\SUMMARIES\DLC\LOOMS\SPEED\';
    
    
    %Run through all the files 
load((strcat('ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS'); 
     loom_rows = ALL_LOOM_ROWS(1,:)'; 
     num_looms = numel(loom_rows);
   
     % Load the DLC files. 
     if exist(fullfile(strcat(pos_folder,'DLC_LOOM_POS_', exp_name,'DLC.mat')), 'file') 
        load(fullfile(strcat(pos_folder,'DLC_LOOM_POS_', exp_name,'DLC.mat')), 'dlc_tbl'); 
        load(fullfile(strcat(speed_folder,'DLC_LOOM_SPEED_', exp_name,'DLC.mat')), 'dlc_speed');
        

for j = 1:num_looms
           
%         dlc_pos_array = zeros(900, 24); 
%         dlc_speed_array = zeros(900,8);
        
         row = loom_rows(j,1) - 299; % Start 5s (180 frames) before the loom stimulus starts 
         row3 = loom_rows(j,1) + 600;  % end 10s after loom starts. 
         
         if row<0 
             row = 1;
         end 
         
         if row3 > height(dlc_speed)
            row3 = height(dlc_speed);
         end 
         
         dlc_pos_array = dlc_tbl(row:row3, :); 
         dlc_speed_array = dlc_speed(row:row3, :); 
         
          save(fullfile(strcat(pos_saving_folder,'DLC_LOOM_POS_', exp_name, '_0', string(j), '.mat')), 'dlc_pos_array'); 
         save(fullfile(strcat(speed_saving_folder,'DLC_LOOM_SPEED_', exp_name,'_0', string(j),'.mat')), 'dlc_speed_array');
         
         
end 
          
 
     else 
     end 

end 

end 

        

  %          %pos
%           if numel(dlc_pos_array) < 900 && row ~=1
%               diff = 900 - height(dlc_pos_array); 
%               dlc_pos_array(height(dlc_pos_array)+1:1081)= zeros(1,diff);
%           elseif numel(M) < 780 && row ==1
%               diff = 780 - numel(M); 
%               added_zeros = zeros(1,diff); 
%               M = [added_zeros, M];
%           end  