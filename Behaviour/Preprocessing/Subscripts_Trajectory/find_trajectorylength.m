function find_trajectorylength()

global exp_name %animal_num exp_date exp_name_short

%% Load ALLDATA

load(strcat('ALLDATA_', exp_name, '.mat'), 'ALLDATA');
load(strcat('LOGDATA_', exp_name, '.mat'), 'log_data');
load(strcat('STIMROWS_', exp_name, '.mat'), 'STIMROWS');

% DATA_size is the number of jpg images. Same as num_bas
DATA_size = size(ALLDATA); 
DATA_size = DATA_size(1);

% num_stim is the number of stimuli presented in this experiment. 
num_stim = size(STIMROWS);
num_stim = num_stim(1);
     


%% DISTANCE COVERED FOR EACH FRAME 
% Make table: ALL_DIST. Contains distance information between frames for ALL FRAMES. 
% In DATA remember that is contains 'no stimulus' rows before and after.  

% This table contains:
% Col 1: the JPG names. 
% Col 2/3: centroid location from TRACKING_DATA
% Col 4: Distance between the centroid locations in one row to the row
% after. 
% Col 5: Cumulative distance 
  
     ALL_DIST = table();
     d_mouse = 0;
     
     %Add JPG names to column 1. 
     ALL_DIST(:,1)= ALLDATA(:,1);

     
     %% Find distance between each row and add to Col4. 
     for i = 1:DATA_size-1
 %     d_mouse(i,1) = (pdist([ALLDATA{i,11}, ALLDATA{i,12}; ALLDATA{(i+1),11}, ALLDATA{(i+1),12}]));
       ALL_DIST{i,2} = ALLDATA{i,12}; %Centroid x
       ALL_DIST{i,3} = ALLDATA{i,13}; %Centroid y 
             A = pdist([ALLDATA{i,12}, ALLDATA{i,13}; ALLDATA{(i+1),12}, ALLDATA{(i+1),13}]); 
             d_mouse = d_mouse+A;
       ALL_DIST{i,4} = d_mouse;  % Distance 
     end 
     
     %% Add cumulative distance to ALL_DIST Col5. 
      cum_d = 0; 
      for i = 1:DATA_size-1
          row_dist = ALL_DIST{i,4};
          cum_d = cum_d + row_dist;
          ALL_DIST{i,5} = cum_d; % cumulative distance 
      end   
      
      ALL_DIST.Properties.VariableNames = {'JPG', 'Centroid_x', 'Centroid_y', 'Distance', 'Cum_D'};
    
save(strcat('ALL_DIST_', exp_name, '.mat'), 'ALL_DIST');   

%% DISTANCE PER STIMULUS 
% Generate table STIMDIST

% This tables contains distance information where each row corresponds to a
% different stimulus presentation. 

% Col 1 = Stimulus #
% Col 2 = Sum of the distance between the centre of the mouse position
% calculated from 'Track_mouse_position' for the number of frames
% corresponding to that stimulus.
% Col 3 = Number of rows for that stimulus 
% Col 4 = Start row
% Col 5 = End row
% Col 6 = normalised distance for that stimulus. (Distance / # rows). 


    STIMDIST = table(); 
     
%% Finding the distance according to 'stimuli' 
    % K is the number of the stimuli 
 
    for k = 1:num+2 %the number of stimuli 
      
       STIMDIST{k,1} = strcat('STIM',string(k));
       stim_name = STIMROWS{k,3};
%        STIMDIST{k,2} = string(stim_name);
       
       mouse_d = 0; 
       A = 0;
       
        if k ==1 
         start_row = 1;
         end_row = STIMROWS{k,2};
        elseif k > 1
            start_row = STIMROWS{k-1,2}+1;
            end_row = STIMROWS{k,2};
        end 
        
     for i = start_row:end_row-1
            A = pdist([ALLDATA{i,12}, ALLDATA{i,13}; ALLDATA{(i+1),12}, ALLDATA{(i+1),13}]); 
            mouse_d = mouse_d+A;
     end 
    
        STIMDIST{k,2} = mouse_d; % distance 
        STIMDIST{k,3} = end_row-start_row+1; % number of rows
        STIMDIST{k,4} = start_row;
        STIMDIST{k,5} = end_row;
          
    end 
    
    % 4th column is the normalised distance. distance/ number of rows. 
    for j = 1:num+2
        STIMDIST{j,6} = STIMDIST{j,2}/STIMDIST{j,3};
    end 
  
    STIMDIST.Properties.VariableNames = {'StimNum', 'distance', 'nRows', 'startRow', 'endRow', 'NormDist'}; 
    
    TOTAL_DIST = sum((STIMDIST{:,2})); %double - sum all distances together. 
  

%% Generate table 'EXP_SUMMARY_TABLE'
% First ROW is TOTAL values for that 'experiment'
% Col 1: Date 
% Col 2: Animal #
% Col 3: ExpName
% Col 4: Stim #
% Col 5: Stim type
% Col 6: Start ROW in DATA/ALLDATA
% Col 7: End ROW in DATA/ALLDATA
% Cpl 8: nRows
% Col 9: Distance
% Col 10: Norm Dist

EXP_SUMMARY_TABLE = table();

animal_num = '19RB01';
exp_date = '190930';
exp_name_short = '02_Dots';

% Fill in values for the first row corresponding to the WHOLE VIDEO/EXPERIMENT. 

EXP_SUMMARY_TABLE{1,1}= {exp_date};
EXP_SUMMARY_TABLE{1,2}= {animal_num};
EXP_SUMMARY_TABLE{1,3}= {exp_name_short};
EXP_SUMMARY_TABLE{1,4}= {'ALL'};
EXP_SUMMARY_TABLE{1,5}= {'Total'};
EXP_SUMMARY_TABLE{1,6}= 1;
EXP_SUMMARY_TABLE{1,7}= DATA_size;
EXP_SUMMARY_TABLE{1,8}= DATA_size; %nRows
EXP_SUMMARY_TABLE{1,9}= TOTAL_DIST; 
EXP_SUMMARY_TABLE{1,10}= TOTAL_DIST/EXP_SUMMARY_TABLE{1,8}; 
EXP_SUMMARY_TABLE{1,11}= TOTAL_DIST/(EXP_SUMMARY_TABLE{1,8}*25); 

for i =2:num+1 
    EXP_SUMMARY_TABLE{i,1}= {exp_date};
    EXP_SUMMARY_TABLE{i,2}= {animal_num};
    EXP_SUMMARY_TABLE{i,3}= {exp_name_short};
    EXP_SUMMARY_TABLE{i,4}= {strcat('STIM', string(i-1))}; 
    EXP_SUMMARY_TABLE{i,5}= log_data{i-1,1};
    EXP_SUMMARY_TABLE{i,6}= STIMDIST{i,4};
    EXP_SUMMARY_TABLE{i,7}= STIMDIST{i,5};
    EXP_SUMMARY_TABLE{i,8}= STIMDIST{i,3}; %nRows
    EXP_SUMMARY_TABLE{i,9}= STIMDIST{i,2}; %distance 
    EXP_SUMMARY_TABLE{i,10}= STIMDIST{i,6}; 
    EXP_SUMMARY_TABLE{i,11}= EXP_SUMMARY_TABLE{i,9}/(EXP_SUMMARY_TABLE{i,8}*25); %speed (dist/nRows*25ms)
end 


EXP_SUMMARY_TABLE.Properties.VariableNames = {'Date', 'Animal_Num', 'Exp_Name','StimulusNumber', 'StimulusType', 'StartRow', 'EndRow', 'nRows', 'Distance', 'NormDist', 'Speed'};

save(strcat('EXP_SUM_TAB_', exp_name, '.mat'), 'EXP_SUMMARY_TABLE');

%% Make ANIMAL_SUMMARY_TABLE with all the EXP_SUM in one giant table. 

%% Plot the distances             
        
%%% Improve plot with distances. Add axis titles etc. 

x = (1:1:num+1);
y = STIMDIST{:,4}; %normalised distance. 
bar(x,y);
 








end 












%% 
% From ALLDATA calculate: 

% the length of the mouse trajectory

% the speed of the mouse - will need to bin into certain sections. 

% Find these values before and after the presentation of the stimulus. 
% In the 'grey screen' time
% When stripes are stationary 
% When stripes move clock/anti 
% When mouse moves with / against stripes. 

 %%  To find the distance the mouse covered in a particular period of time.
 
 %Each ROW of DATA and ALLDATA corresponds to 25ms. (1000ms/40) - 40fps. 

%% Extra code: 

% if stim_start == 0
%     start_row = 1;
% elseif stim_start > 0
%     start_row = STIMROWS(stim_start,2);
% end 
%   
% if stim_stop == DATA_size+1
%     end_row = DATA_size;
% else
%     end_row = STIMROWS(stim_stop,2);
% end

    % Dot position
    % dot_x = ALLDATA{k,3};
    % dot_y = ALLDATA{k,4};
       
    %Transform dot coords from python. 
%     [alpha1, rad1] = cart2pol (dot_x, dot_y);
%     alpha1 = alpha1 + alpha_zero ;
%     rad1 = rad1*basler_image_size/1000;
%     [dot_x2, dot_y2] = pol2cart(alpha1, rad1);
 
   % Mouse position
   % mouse_x = ALLDATA{k,11};
   % mouse_y = ALLDATA{k,12};
   
   
   
   %%% 
%        
%         
%        stim_start = k-1;
%        stim_stop = k;
%        
%         if stim_start == 0
%             start_row = 1;
%         elseif stim_start > 0
%             start_row = STIMROWS(stim_start,2);
%         end 
%   
%         if stim_stop == num+1
%             end_row = DATA_size;
%         else
%             end_row = STIMROWS(stim_stop,2);
%         end