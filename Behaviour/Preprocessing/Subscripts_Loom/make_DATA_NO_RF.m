function make_DATA_NO_RF()
% Created by Burnett 11/03/22

% FOR MAKING DATA WHEN PROBLEM WITH RF.

% -- Must run 'make_log_data.m' before running this script to generate 'LOGDATA' and 'MATDATA'. 
% -- Generates INFO and STIMROWS as outputs. 
% -- Process data from log files to generate table 'INFO' which contains the jpg image name and the corresponding mouse position and stimulus properties for that frame. 
% -- Aligns rows based on detection from arduino and when RF were commanded in the vbl_array from MATLAB and linearly interpolates the vbl_array rows between the 'anchored' red frames. 
 

% IF - running through 'process_dome_data.m' 
global inpath exp_name animal_num exp_date 


load(strcat('LOGDATA_', exp_name, '.mat'), 'log_data');
% num = height(log_data);
total_stim_time = log_data.t_cumul(end);

%% STAGE 1: FIND RFS FROM ARDUINO 


% 
% 
% 
% %% RED FRAMES RECORDED BY PHOTODIODE - - - - - - - Find out how many RFs the photodiode recorded from the arduino file. 
% % inpath = 'C:\Data_analysis\DATA\Setd5\191216\GN3959\01_Loom';
ardfile = dir(fullfile(inpath,'*.txt')); %gets txt files in folder
ard_name = ardfile(1).name;
ard = dlmread(ard_name, ',', 60,0); %Ignore the first 50 lines. 
% 
% %% FILTER ARDUINO FILE - DORIC SHOULD ALWAYS START BEFORE MATLAB. 
% 
% % Make any value <1000 in the 4th column == 0
% ard(ard(:,4)<800,4) = 0;
% ard(ard(:,4)>=800,4) = 100;
% num_ard = numel(ard(:,1));
% 
% %Generate sequence of zeroes and a 100 and isolate where this is found in the fourth column of the 'ard' file. 
% seq1 = zeros(30,1);
% seq1(31, 1) = 100;
% 
% col4_trans = ard(:,4)'; %Transpose column 4
% seq1 = seq1'; %Transpose seq
% Index = strfind(col4_trans, seq1);
% Index_RFs = Index';
% Index_RFs = Index_RFs+30; %Array of ‘ard’ row values of all RFs. 
%  
% % 'Index_RFs' should start with the first RF of the Grey Screeen stimulus. 
% % The first GAP should be at the END of the Grey stimulus. 
% 
% % Find when MATLAB presentation of stimulus started and stopped. 
% start_RF = Index_RFs(1); %row value of first RF in ‘ard’
% end_RF = Index_RFs(end); % row value of last Rf in ‘ard’
% 
% %Find the number of Rfs from the arduino file. 
% num_RF = numel(Index_RFs);
% 
% for j = 1:num_RF
%     row = Index_RFs(j); %ard row value of RF
%     ard(row+1:row+21,4) = 0; % Set 20 rows after the ‘RF row’ found in Index_Orig to zero. 
%     Index_RFs(j,2)= ard(row,1); % add time from ard to Index_RFs
% end 

%The TIME of the first ROW - i.e. the presentation of the FIRST RF, i.e.
%Index_RFs(1,2) will be the time in us from when the arduino file first
%starts recording til when the first RF was recorded. 

%% Plot the RFs. 
% 
% figure
% for j = 1:num_RF
%     row = Index_RFs(j);
%     x = [ard(row, 1), ard(row,1)];
%     y = [0, ard(row, 4)];
%     plot(x,y,'k');
%     hold on 
% end

% Add title and save figure.
% title_new = regexprep(ard_name, '_', '-');
% title_new = title_new(1:end-4);
% title(title_new)
% savefig(gcf, strcat('RF_',title_new, '.fig'));

%% Find the IRFI. 
%  
% for k = 1:num_RF-1
%   Index_RFs(k,3) = Index_RFs(k+1,2)-Index_RFs(k,2);
% end 

% Col 3 Index_RFs - is now the IRFI in us. Should be ~83000 == 83ms. 
 % this is the time differece from the RF of this row til the NEXT RF after
 % it! 
 
%%  Plot IRFI 
% plot(Index_RFs(:,3))
% title(title_new)
% savefig(gcf, strcat('IRFI_',title_new, '.fig'));
% 
% if height(log_data)==1   % For Acclim experiments!! 
%     
%     n_stim = 1; 
%     
%     % Here there should not be "(Index_RFs(:,3)>250000)" . 
%     stim_ard_rows = zeros(1,1);
%     % Col 1 = 'ard' row values for the first RF in each stimulus. 
%     % Col 2 = 'ard' row values for the last RF in each stimulus.
% 
%     stim_ard_rows(1,1)= start_RF; %(1,1) = row value in 'ard' of FIRST RF
%     stim_ard_rows(1,2)= end_RF; % THIS IS NOT THE END OF THE STIMULUS BUT THE LAST RF. - could be a few omre rows after last RF before stim ends. 
% 
%     row1 = stim_ard_rows(1,1);
%     row2 = stim_ard_rows(1,2);
%     stim_ard_rows(1,3)=numel(find((ard(row1:row2,3))==1)); % Number of images between these points.
% 
% else   % If NOT ACCLIM and THERE IS >1 STIM.  
%     
%     % Find when each new stimulus started from finding the row values where the big gaps in RFs are.
%     
%     % end_of_stim_rows = find(Index_RFs(:,3)>300000); %array - row value in Index_RFs for LAST RF BEFORE BREAK (0.5s).
%     start_of_stim_rows = find(Index_RFs(:,3)>350000)+1; %The value 250000 depends on HOW LONG THE RF GAP IS!! (0.25s)
%     % start_of_stim_rows - Col1 = ROW value of the RF GAP in IndexRFs. 
%     % For the WT_01 Exps : 120 - 660- 780- 1320- 1440 -1980
%     
%     n_stim = numel(start_of_stim_rows(:,1));
%     % For the 'old' loom stimulus wth Green-Loom-Green-Loom-Green
%     % start_of_stim_rows should have 4 rows since there are 4 gaps between 5
%     % stim. n_stim should = 4. See comments below for more details.
%    % Whereas with the MULTILOOM EXPERIMENTS - with 12 stimuli - there
%    % should be 11 gaps. 
%     
%    % Col 3 - IRFI between these gaps and the RF before. 
%     for j = 1:n_stim
%         start_of_stim_rows(j,3) = Index_RFs(start_of_stim_rows(j,1)-1,3);
%     end
%     
%     % When there are weird blocks of lots of red at the beg /end of stimulus
% %     large_RF_gaps = find(start_of_stim_rows(:,3)> 500000);
% %     start_of_stim_rows(large_RF_gaps, :)=[]; %REMOVE THESE ROWS!
%     
% %     if exp_name_short(end-4:end)=="_Loom"
% %         if n_stim>4
% %             fprintf("More than 4 RF blocks\n");
% %             % Also need to change End_RF!!
% %             if Index_RFs(end,3)>100000
% %                 end_RF = Index_RFs(end-1,1);
% %             end
% %         end
% %     end
%     
%     % Refine the number of stimuli now that have removed rows from
%     % 'start_of_stim_rows'
% %     n_stim = numel(start_of_stim_rows(:,1));
%     
%     %%%%%%%%%% KEY NOTES FOR TROUBLESHOOTING
%     % For the Setd5 Loom experiments 'start_of_stim_rows' should have 4 values!
%     % "find(Index_RFs(:,3)>250000)" will find the LAST RF before the RF GAP. These will be:
%     % - at the end of Grey1
%     % - at the end of loom1
%     % - at the end of Grey Mid
%     % - at the end of Loom2.
%     
%     % Adding +1 finds the ard row value for the FIRST RF of:
%     % - the LOOM1 STIM (~180)
%     % - the Grey Mid Stim (~900)
%     % - the LOOM2 STIM (~1620)
%     % - the Grey end stim. (~2340)
%     
%     % if 0.25s - gap of RFs ~ 250,000 us. SHOULD MAKE THIS LONGER!!! If drops
%     % 2/3 RFs would also be this long!
%     
%     % ( +1 ) to this value since the beginning of the next stimulus is the first RF
%     % after this break. The RF corresponding to this row is the last RF of the
%     % stimulus before the RF BREAK.
%     
%     % Grey stim for 15s = 900 frames. 900/5 = 180 RFs. Therefore the first
%     % Index_Rfs row after the first RF gap should ~= 180.
%     
%     for i = 1:n_stim
%         %     index_end_row = end_of_stim_rows(i);
%         %     end_of_stim_rows(i,2)= Index_RFs(index_end_row,1);
%         index_start_row = start_of_stim_rows(i);
%         start_of_stim_rows(i,2)= Index_RFs(index_start_row,1); %row in ard file of the first RF of that stimulus
%     end
%     
%     stim_ard_rows = zeros(n_stim+1,1);
%     % Col 1 = 'ard' row values for the first RF in each stimulus.
%     % Col 2 = 'ard' row values for the last RF in each stimulus.
%     
%     % For Setd5 loom experiments - should have 5 ROWS.
%     % Row 1 = Grey1, Row 2 = Loom1, Row 3 = Grey Mid, Row 4 = Loom2, Row 5 =
%     % Grey End.
%     
%     stim_ard_rows(1,1)= start_RF; %(1,1) = row value in 'ard' of FIRST RF
%     stim_ard_rows(n_stim+1,2)= end_RF; % THIS IS NOT THE END OF THE STIMULUS BUT THE LAST RF. - could be a few omre rows after last RF before stim ends.
%     
%     for j = 1:n_stim
%         stim_ard_rows(j+1,1) = start_of_stim_rows(j,2); %Col 1 = row in ard file that corresponds to first RF of stimulus
%         stim_ard_rows(j,2)= (start_of_stim_rows(j,2)-1); % Col 2 = row in ard file that corresponds to LAST RF of that stimulus
%     end
%     
%     
%     for j = 1 :n_stim+1
%         row1 = stim_ard_rows(j,1);
%         row2 = stim_ard_rows(j,2);
%         stim_ard_rows(j,3)=numel(find([ard(row1:row2,3)]==1)); % Number of images between these points.
%     end
% 
% end 
%Col 3 - should be 900, 3600, 3600, 3600, 900! 
% Col 3 - for multiloom of 20s stimuli = 1200. 

% for j = 1:n_stim+1
%     stim_ard_rows(j,3) = stim_ard_rows(j,2)- stim_ard_rows(j,1);
% end 

%% ROWS FROM ARDUINO WHERE RFS WERE RECORDED.

%Remove 'noisy signal' from RF value in arduino file by changing all values
%25 rows after the 'start row' of the RF into '0'. 
% for i = 1:num_RFs
%     row = Index_RFs(i);
%     ard(row+1:row+31,4)= 0;
% end 
%     
%Find the number of rows between the red frames. 
% for i = 1:num_RF-1
%     Index_RFs(i,2) = Index_RFs(i+1,1)-Index_RFs(i,1); 
% end

%Find the de facto irfi by finding the mean time between the recorded rfs, with a gap of less than 95ms. 
% rows_lessthan95 = find(Index_RFs(:,3)<95000);
% nrows = numel(rows_lessthan95);
% 
% irfi_sum = 0;
% for j = 1:nrows
%     ind = rows_lessthan95(j,1);
%     val = Index_RFs(ind,3); 
%     irfi_sum = val + irfi_sum; 
% end 
%      
% irfi = round(irfi_sum/nrows); % Will force RFs to be added this distance between the last recorded red frame - to fill in spaces of dropped RFs. 
 
%% DROPPPED FRAMES. 
% Find the rows in Index_Orig where column 3 is >95000
% rows_morethan95 = find(Index_RFs(:,3)>95000); %This gives the row index of the row in index_orig where the gap is greater than 100 rows. 
% n_more95 = numel(rows_morethan95);
% 
% for i = 1: n_more95
%     rows_morethan95(i,2) = Index_RFs(rows_morethan95(i),3);
% end 
% 
% for k = 1:n_more95
%     ind2 = rows_morethan95(k); %Row at which start of RF is before gap.
%     nrows_diff = rows_morethan95(k,2);
%     ard_row = Index_RFs(ind2,1);
%     num_RFs_to_add = round(nrows_diff/irfi); 
%     if num_RFs_to_add == 2
%         ard(ard_row+60,4)=100;
%     elseif num_RFs_to_add == 3
%         ard(ard_row+60,4)=100;
%         ard(ard_row+(60*2),4)=100;
%     elseif num_RFs_to_add == 4
%         ard(ard_row+60,4)=100;
%         ard(ard_row+(60*2),4)=100;
%         ard(ard_row+(60*3),4)=100;
%     elseif num_RFs_to_add == 5
%         ard(ard_row+60,4)=100;
%         ard(ard_row+(60*2),4)=100;
%         ard(ard_row+(60*3),4)=100;
%         ard(ard_row+(60*4),4)=100;
%     end 
% end 
% 
% %% ROW VALUES OF ARDUINO FILE WHERE RFS HAPPENED.           ** "Index_StartRF" **
% col4_trans_02 = ard(:,4)'; %Transpose column 4
% Index_RFs_2 = strfind(col4_trans_02, seq1);
% Index_RFs_2 = Index_RFs_2';
% Index_RFs_2 = Index_RFs_2+30; 
% 
% %% TOTAL NUMBER OF RFS IDENTIFIED FROM PROCESSED ARDUINO FILE.               ** "nRFs_Ard" **
% nRFs_Ard = numel(Index_RFs_2);

%For TROUBLESHOOTING: Find the # of rows difference between when the RFs took place. 

% for i = 1:nRFs_Ard-1
%     Index_StartRF(i,2) = Index_StartRF(i+1,1) -Index_StartRF(i,1);
% end 

%Plot graph of difference in rows between RF rows. 
% figure 
% plot(Index_StartRF(:,2))

% Find when MATLAB presentation of stimulus started and stopped. 
% start_RF = Index_StartRF(1,1);
% end_RF = Index_StartRF(end,1);

%% STAGE 2: COMBINE VBL_ARRAYS














%% Generate 'ALL_DATA_VBL' cell array with all the information from the '.mat'/vbl_array files for each stimulus. 
% Load MAT file. (my_matfiles)
matfiles = dir(fullfile(inpath,'MATDATA*')); %gets all log files in folder
load(matfiles(1).name);  %loads my_matfiles - all the names of mat files
num = length(my_matfiles); %num = 5

fps_matlab = 60;
tot_all_vbl = total_stim_time*fps_matlab; 
cum_T = 0;

% Define 'DATA_ALL_VBL'
DATA_ALL_VBL = cell(round(tot_all_vbl), 13); 

for i = 1:num %run through the mat files. %%%%% CHANGE -1 BACK WHEN DON'T HAVE TRACKING FILE 
    
    stimulus_type = string(log_data.stim(i));  
      
    mat= load(my_matfiles(i).name);
    mat= mat.vbl_array;
    end_of_stim = find(mat(:,1)==0); %After stimulus is finished vbl_array is filled with '0'
    end_of_stim = end_of_stim(1);
    
%     if end_of_stim>tot_all_vbl
%         end_of_stim = tot_all_vbl;
%     end 
    mat = mat(1:end_of_stim-1, :);
 
      if stimulus_type=="CL_green_screen" || stimulus_type=="CL_black_screen" || stimulus_type=="CL_graded_screen" || stimulus_type=="CL_green_screen_1S"

        for j = 1:end_of_stim-1
            DATA_ALL_VBL(j+cum_T, 1)= {mat(j,1)};
            DATA_ALL_VBL(j+cum_T, 2)= {mat(j,2)};
            DATA_ALL_VBL(j+cum_T, 3)= log_data.stim(i);
            DATA_ALL_VBL(j+cum_T, 4)= {50}; % dot x 
            DATA_ALL_VBL(j+cum_T, 5)= {50}; % dot y 
            DATA_ALL_VBL(j+cum_T, 6)= {mat(j, 3)}; % mouse x
            DATA_ALL_VBL(j+cum_T, 7)= {mat(j, 4)}; % mouse y 
            DATA_ALL_VBL(j+cum_T, 8)= {0};
            DATA_ALL_VBL(j+cum_T, 9)= {[log_data.BkR(i), log_data.BkG(i), log_data.BkB(i)]}; % BKG RGB
            DATA_ALL_VBL(j+cum_T, 10)= {0}; % DOT RGB. 
            DATA_ALL_VBL(j+cum_T, 11)= {animal_num};
            DATA_ALL_VBL(j+cum_T, 12)= {exp_date};
            DATA_ALL_VBL(j+cum_T, 13)= {log_data.ExpName(i)};
        end       
        
      elseif  stimulus_type=="CL_looming_1S_TRIGGER_simplified" || stimulus_type=="CL_looming_1S_TRIGGER_DiffContrast" || stimulus_type=="CL_looming_CHECKER" || stimulus_type=="CL_looming_dot_distractor"
              
           for j = 1:end_of_stim-1
            DATA_ALL_VBL(j+cum_T, 1)= {mat(j,1)};
            DATA_ALL_VBL(j+cum_T, 2)= {mat(j,2)};
            DATA_ALL_VBL(j+cum_T, 3)= log_data.stim(i); %stimulus 
            DATA_ALL_VBL(j+cum_T, 4)= {mat(j,4)}; % radius. 
            DATA_ALL_VBL(j+cum_T, 5)= {mat(j,3)}; % LOOMVALUE
            DATA_ALL_VBL(j+cum_T, 6)= {mat(j,5)}; % mouse x
            DATA_ALL_VBL(j+cum_T, 7)= {mat(j,6)}; % mouse y
            DATA_ALL_VBL(j+cum_T, 8)= {mat(j,7)}; % TRIGGERVALUE - 0 for NO LOOM, 1 for LOOM.
            DATA_ALL_VBL(j+cum_T, 9)= {[log_data.BkR(i), log_data.BkG(i), log_data.BkB(i)]}; % BKG RGB
            DATA_ALL_VBL(j+cum_T, 10)= {[log_data.dotR(i), log_data.DotG(i), log_data.DotB(i)]}; % DOT RGB. 
            DATA_ALL_VBL(j+cum_T, 11)= {animal_num};
            DATA_ALL_VBL(j+cum_T, 12)= {exp_date};
            DATA_ALL_VBL(j+cum_T, 13)= {log_data.ExpName(i)};
           end
           
      elseif stimulus_type=="CL_looming_1S_TRIGGER_ILI"
  
          for j = 1:end_of_stim-1
            DATA_ALL_VBL(j+cum_T, 1)= {mat(j,1)};
            DATA_ALL_VBL(j+cum_T, 2)= {mat(j,2)};
            DATA_ALL_VBL(j+cum_T, 3)= log_data.stim(i); %stimulus 
            DATA_ALL_VBL(j+cum_T, 4)= {mat(j,4)}; % radius. 
            DATA_ALL_VBL(j+cum_T, 5)= {mat(j,3)}; % LOOMVALUE
            DATA_ALL_VBL(j+cum_T, 6)= {mat(j,5)}; % mouse x
            DATA_ALL_VBL(j+cum_T, 7)= {mat(j,11)}; % ILI VALUE. 
            DATA_ALL_VBL(j+cum_T, 8)= {mat(j,7)}; % TRIGGERVALUE - 0 for NO LOOM, 1 for LOOM.
            DATA_ALL_VBL(j+cum_T, 9)= {[log_data.BkR(i), log_data.BkG(i), log_data.BkB(i)]}; % BKG RGB
            DATA_ALL_VBL(j+cum_T, 10)= {[log_data.dotR(i), log_data.DotG(i), log_data.DotB(i)]}; % DOT RGB. 
            DATA_ALL_VBL(j+cum_T, 11)= {animal_num};
            DATA_ALL_VBL(j+cum_T, 12)= {exp_date};
            DATA_ALL_VBL(j+cum_T, 13)= {log_data.ExpName(i)};
           end
          
      end
      
    cum_T = cum_T + end_of_stim-1;
    log_data.DATA_ALL_VBL_row(i) = cum_T;
    log_data.StimNum(i) = str2double(string(i));
end 

num_all_vbl = numel(DATA_ALL_VBL(:,1));

%% FIND # RED FRAMES 'SENT' BY MATLAB - - - - - - - Find out how many RF commands were sent from the vbl_array files.                *** MATLAB_RF_rows *** 

MATLAB_RF_rows = find([DATA_ALL_VBL{:,2}]==1);
MATLAB_RF_rows = MATLAB_RF_rows';
nRFs_MAT = numel(MATLAB_RF_rows);

% For TROUBLESHOOTING: Find the #rows difference between matlab commands to show RF. This should be 5. 

% for m = 1: n_matrows-1
%     MATLAB_RF_rows(m,2) = MATLAB_RF_rows(m+1,1)-MATLAB_RF_rows(m,1);
% end 

% Find out when the difference is not 5. 
% find(MATLAB_RF_rows(:,2)~=5)
% MATLAB_RF_rows([360, 738, 1103, 1469],2)

%% STAGE 3: INTERPOLATE AND COMBINE IN 'INFO' 












%% Option 1: Anchoring at each RF. 
% Combining arduino RFs and MATLAB RFs. 
% Index_RFs_2  = array of all the rows in 'ard' where RF starts. 
% MATLAB_RF_rows = array with all the row values in DATA_ALL_VBL where RF command sent. 

% 1- Add DATA_ALL_VBL row value to 'ard' for RFs. 

% for i = 1:nRFs_Ard
%     index = Index_RFs_2(i);
%     matlab_index = MATLAB_RF_rows(i);
% %     vbl_value = vbl_array(matlab_index);
%     ard(index, 5) = matlab_index; 
% end 
% 
% % 2- Between index_startRf(1) and indexstartRF(2) must *interpolate* between values MATLAB_RF_rows(1) and MATLAB_RF_rows(2).  
% 
% for j = 1: nRFs_Ard-1
%      x = [Index_RFs_2(j,1), Index_RFs_2(j+1,1)];
%      v = [MATLAB_RF_rows(j,1), MATLAB_RF_rows(j+1,1)];
%      xq = (Index_RFs_2(j,1) : 1 : Index_RFs_2(j+1,1));
%      xq_size = numel(xq);
%      vq = interp1(x,v,xq);
%      for j = 1:xq_size
%          ard_index2 = xq(j);
%         ard(ard_index2,5)=round(vq(j));
%      end 
% end 
   

%% OPTION 2 - anchor to where 'gap' in RFs are from the 0.5s without RFs. 

% log_data.DATA_ALL_VBL_row
% stim_ard_rows

% 1- Add DATA_ALL_VBL row value to 'ard' for RFs. 
% num_stimrows = numel(stim_ard_rows(:,1));
% 
% % Adding 4th column to stim_ard_rows: the last row of DATA_ALL_VBL for that
% % stimulus. 
% for k =1:num_stimrows
%     if k == num_stimrows && num_stimrows ~=1
%          stim_ard_rows(k,4)=log_data.DATA_ALL_VBL_row(num-1);
%     elseif height(log_data)==1
%         stim_ard_rows(1,4)=log_data.DATA_ALL_VBL_row(1);
%     else 
%     stim_ard_rows(k,4)=log_data.DATA_ALL_VBL_row(2*k);
%     end 
% end 
% 
% % 2- Between index_startRf(1) and indexstartRF(2) must *interpolate* between values MATLAB_RF_rows(1) and MATLAB_RF_rows(2).  
% 
% for j = 1:num_stimrows
%      if j == 1
%          x = [stim_ard_rows(j,1), stim_ard_rows(j,2)];
%          v = [1, stim_ard_rows(j,4)];
%          xq = (stim_ard_rows(j,1): 1 : stim_ard_rows(j,2));
%          xq_size = numel(xq);
%          vq = interp1(x,v,xq);
%      else
%          x = [stim_ard_rows(j,1), stim_ard_rows(j,2)];
%          v = [(stim_ard_rows(j-1,4)+1), stim_ard_rows(j,4)];
%          xq = (stim_ard_rows(j,1): 1 : stim_ard_rows(j,2));
%          xq_size = numel(xq);
%          vq = interp1(x,v,xq);
%      end
%      for k = 1:xq_size
%          ard_index = xq(k);
%          ard(ard_index,5)=round(vq(k));
%      end
% end 


% 3 - Extract all rows from 'ard' where camera took image. i.e. column 3 ==1. 
ard_camera = ard(ard(:,3)==1, :);
%Each ROW of this array corresponds to one basler image. 
% Column 5 now gives the DATA_ALL_VBL array row number for that image which
% contains stimulus information. 

%% MAKE TABLE 'INFO' WHICH CONTAINS BASLER IMAGE WITH SYNCHRONISED VBL_ARRAY INFORMATION. 

% save(strcat('LOGDATA_', exp_name, '.mat'), 'log_data');
% num_all_vbl = numel(DATA_ALL_VBL(:,1));

% Load the basler folder
baslerpath = fullfile(inpath, '/basler_*');
myBaslerFiles = dir(fullfile(baslerpath,'*.jpg'));
num_bas = numel(myBaslerFiles);
num_ard_cam = numel(ard_camera(:,1));

if num_bas > num_ard_cam
    num_to_use = num_ard_cam;
else 
    num_to_use = num_bas;
end 

% Define cell 'Info'
% Info = cell(num_bas-1,13);
Info = cell(num_to_use,13);

for i = 1:num_to_use %num_bas-1
     Info{i,1} = {myBaslerFiles(i).name}; 
                                             % Col 1 = basler image
%      if i == num_bas-1 || i == num_bas 
%          Info{i,2} = ard_camera(num_bas-2,5);
%         else 
%      Info{i,2} = ard_camera(i,5);% Col 2 = DATA_ALL_VBL array row value from ard_camera. 
%      end 
      if i <= 50 || i > num_all_vbl
          Info{i,2} = 0;
          Info{i,3} = 0; 
          Info{i,4} = 0; 
          Info{i,5} = 0; 
          Info{i,6} = 0; 
          Info{i,7} = 0;  
          Info{i,8} = 0;  
          Info{i,9} = 0;  
          Info{i,10} = 0;  
          Info{i,11} = 0; 
          Info{i,12} = 0;  
          Info{i,13} = 0;
      else
          row = i-50; 
          Info{i,2} = 1;
          Info{i,3} = DATA_ALL_VBL{row,3}; % Stimulus
          Info{i,4} = DATA_ALL_VBL{row,4}; % dot x / radius
          Info{i,5} = DATA_ALL_VBL{row,5}; % dot y / loom values
          Info{i,6} = DATA_ALL_VBL{row,6}; % mouse x - from python - OR ColG
          Info{i,7} = DATA_ALL_VBL{row,7}; % mouse y - from python - OR ColB OR ILI
          Info{i,8} = DATA_ALL_VBL{row,8}; % triggervalue
          Info{i,9} = DATA_ALL_VBL{row,9}; % bkg RGB
          Info{i,10} = DATA_ALL_VBL{row,10}; % dot RGB
          Info{i,11} = DATA_ALL_VBL{row,11}; % animal num
          Info{i,12} = DATA_ALL_VBL{row,12}; % exp_date 
          Info{i,13} = DATA_ALL_VBL{row,13}; % exp_name
      end      
end 
% 

%% STAGE 5: STIMROWS

% Info(:, [3:13]) = DATA_ALL_VBL(:, [3:13]);










%% 2 - Generate double array 'STIMROWS' (size = num of stim: 'num')
% STIMROWS is an array size (num x 2). Where 'num' is the number of stimuli for that 'experiment' - i.e. run of start_stimuli. 
% Therefore, each row is a stimulus. 

% Col 1: The last row value for that stimulus in DATA_ALL_VBL. 
% Col 2: Stimulus Type
% Col 3: Exp_name
% Col 4: Start row of 'Info' 
% Col 5: End row of 'Info' 
% Col 6 = # of rows for that stimulus in 'Info'

VBL_ROW = log_data.DATA_ALL_VBL_row;
num = numel(VBL_ROW); %Number of stimuli.

%Generate cell array STIMROWS. 
STIMROWS = cell(num+2,1);

% Add in information about DATA_ALL_VBL array row, stimulus type and exp_name. 
STIMROWS(1,1)={0};
STIMROWS(num+2,1)={0};
STIMROWS(1,2)={"no stimulus"};
STIMROWS(num+2,2)={"no stimulus"};
STIMROWS(1,3)={"no stimulus"};
STIMROWS(num+2,3)={"no stimulus"};

for i = 2:num+1
STIMROWS(i,1)= {VBL_ROW(i-1)}; 
STIMROWS(i,2) = log_data{i-1,1};
STIMROWS{i,3} = log_data.ExpName{i-1};
end 

% Add information about start and end rows in 'Info' 
stimuli_names = [STIMROWS{:,3}]; %'no stimulus' is the last. 
num_stim = length(stimuli_names);

for i = 2:num_stim-1
    stimulus = stimuli_names(i);
    rows_stim = find([Info{:,13}] == stimulus);
    STIMROWS{i,4}= rows_stim(1);
    STIMROWS{i,5}= rows_stim(end);
end 

% Add row values for "no stimulus"
STIMROWS{1,4} = 1; 
STIMROWS{1,5} = STIMROWS{2,4}-1; 
STIMROWS{num_stim,4} = STIMROWS{num_stim-1,5}+1; 
STIMROWS{num_stim,5} = num_bas; 
% num_stim = 10;
% Add # of rows per stimulus. 
for i = 1:num_stim
    STIMROWS{i,6}= STIMROWS{i,5}- STIMROWS{i,4}+1;
end
 
%% SAVE THE ARRAYS 

% Save STIMROWS
 save(strcat('STIMROWS_', exp_name, '.mat'), 'STIMROWS');

% If need to track the position of the mouse again,save Info as a cell
% array. Otherwise, make it into a table already. 

% if track_req == 1         
%Save the 'Info' cell array as a .mat file. 
save(strcat('INFO_', exp_name, '.mat'),'Info');

% elseif track_req == 0 
%     
%     %Make Info into table. 
%     Info_Table = cell2table(Info);
%     Info_Table.Properties.VariableNames = {'JPG', 'DATA_ALL_VBL_ROW', 'StimulusName', 'Radius', 'Loom', 'Mouse_x', 'Mouse_y', 'Trigger','MouseNum', 'Date', 'ExpName'}; 
% 
%     % SaveInfo as both a cell and a table.  
%     save(strcat('INFOTABLE_', exp_name, '.mat'), 'Info_Table');
% end 

clearvars

end

