function process_dome_data()
%% Process experiment data
% Created by Burnett 
% Last modified 26/03/20

% To be used with data acquired from the PCF dome setup - with IR light
% background.

% Start in the DATE FOLDER. e.g. '191012'
%% 
% 'track_req' refers to whether the experiment occurred PRE 2020 or POST
% 2020. This determines whether post-processing tracking of the mouse is
% required. 

%%
clear all

global inpath exp_name currD animal_num exp_date exp_name_short dir3_path dir2_path folderpath track_req new_setup output_folder ILI diff_contrasts image_size no_shelter checker distractor

% To Specify Beforehand. 

% If exp done in new box without legs then new_setup == 1
new_setup =1; 
% Specify whether the experiments need post-processign tracking or not. 
track_req = 1; 
% Diff contrasts 
diff_contrasts = 0;
% No shelter 
 no_shelter = 0; 
 image_size = 416; % 416; %518; %524; 
 ILI = 0; 
 distractor = 0; 
 checker = 0; 

dir4_path = cd; % PROJECT FOLDER i.e. C:\Data_analysis\DATA\2001_Setd5_C2\
directory4 = dir; % . .. DATE FOLDERS

for ci = 24 % WHICH PROJECT
  D3=directory4(ci).name; %animal number 
  dir3_path = strcat(dir4_path, '\', D3); % DAY FOLDER i.e. C:\Data_analysis\DATA\2001_Setd5_C2\190926
  cd(dir3_path)
%   dir3_path = cd; % PROJECT FOLDER i.e. C:\Data_analysis\DATA\2001_Setd5_C2\
   directory3 = dir('2*'); % . .. DATE FOLDERS

for q = 1 %1:length(directory3)% WHICH DATES 
  D2=directory3(q).name; %animal number 
  dir2_path = strcat(dir3_path, '\', D2); % DAY FOLDER i.e. C:\Data_analysis\DATA\2001_Setd5_C2\190926
  output_folder = dir2_path(1:end-6);
  cd(dir2_path)
  directory2 = dir(); % . .. 19F01 19F02 19RB01

for p = 3:length(directory2) %9:12 %3:12               %For which animals.  % For 210914 from animal 4!   
    
    D=directory2(p).name; %animal number airb
    folderpath = strcat(dir2_path, '\', D);% ANIMAL PATH -  Windows
    cd(folderpath)
   
%     dir_path = cd;
    directory = dir(); %All the files in animal directory. 
    
for k = 3:length(directory) % Specify which experiments.    %25 / 27        
% FOR GAIA'S ANIMALS:

    currD = directory(k).name; %Experiment folder NAME 
    inpath = strcat(folderpath, '\', currD);% Entire path of experiment folder - Windows
    animal_num = D; %dir_path(end-4:end);
    exp_date = folderpath(end-12:end-7);
    exp_name_short = currD;
    exp_name = strcat(exp_date, '_', animal_num, '_', currD);

  cd(inpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % SORT LOG FILESs
% % % % % % % % 
%                        fprintf("Sorting logs...\n");
%                          make_log_data_LOOM_BOX();
% % % % % %                         make_log_data_DISTRACTORS()
%            
%                         fprintf("Making LOGDATA, ROWS, STIMROWS and DATA...\n");
%                          make_DATA_RF_Sync();         
% % % % % % % %                           make_DATA_NO_RF()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% % % % TRACKING 
% % % % 
%                    fprintf("Making TRACKINGDATA...\n"); 
% % %                     Track_Mouse_Position_Box();
%                   Track_Mouse_Position_Box2(); % Using different method -  hopefully faster and more robust
% % %                    Track_Mouse_AVI();
% 
%                   fprintf("Merging Info and TRACKING...\n");
%                   merge_Info_TRACKING()

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% XY ARRAY 
% % % % % % % 
%        fprintf("Making xy_array.. \n");
%        make_xy()
% % % % % % % % % % % % % 
% % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % %% ANALYSIS
% % % % % % % % % % %  
%       fprintf("Making loom array..\n");
%        make_loom_array()
% % % % % % % % % % 
% % % % % % % % % % % % %       REQUIRES NEW XY_ARRAY WITH DIST FROM CENRTRE OF ARENA. 
%       fprintf("Making activity array.. \n");
%       make_ACTIVITY_ARRAY()
% % % % % % % % % % %  
%      fprintf("Making xy loom array..\n");
%      make_xy_loom_array()

% % % fprintf("Making dlc loom array..\n");
% % %  make_DLC_LOOM_array()
% % 
% % % %      fprintf("MAKE RESPONSE AND LOCO...\n");
% % % %      make_RESPONSE_TYPE_and_LOCO_INDEX()
% % 
% % % % % % % % % % % % %  
% 
%         fprintf("Run exit_analysis..\n");
%         make_exit_analysis() % 1 for acclim or banana, 0 for loom 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % LOOM VIDEOS
% %
% % % % % if k ~= 4
%     fprintf("Making stim video...\n");
%     make_stim_plot_video_LOOM()
% elseif k == 4
%     fprintf("Make full video with loom...\n")
%     make_stim_plot_video_FULLVID()
%     
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % REMOVE JPGS - MAKE VIDEOS 
% % 
%                     fprintf("Starting to make video...\n");
%                     cleanup_jpgs() %This makes mp4 video files from JPG images. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %  ONCE ALL DAYS COMPLETE - RUN THIS.  - MUST ADD GENOTYPES. 
% 
% fprintf("Making xy analysis...\n");
% make_xy_analysis()
% % % % % % % % clear
% fprintf("Making full xy array ...\n");
% make_FULL_xy_loom_array()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTO 
% % % % % % SORT LOG FILES - OPTO SETUP 
% 
%                        fprintf("Sorting logs...\n");
%                          make_log_data_LOOM_BOX();
%                          
%                         fprintf("Making LOGDATA, ROWS, STIMROWS and DATA...\n");
%                          make_DATA_OPTO();
%                          
%                           fprintf("Making TRACKINGDATA...\n"); 
%                          Track_Mouse_OPTO();

% % %                               % USE THIS SCRIPT IF YOU HAVE
% % %                               % ACCIDENTALLY ALREADY CONVERTED
% % %                               % THE IMAGES TO A VIDEO. 
% % % %                             Track_Mouse_OPTO_AVI();
%                     
%                        fprintf("Analysing...\n"); 
%                          make_xy_analysis_OPTO2()
%                          make_xy_loom_array_OPTO() 

%                          fprintf("Making video...\n");
%                          make_stim_plot_video_OPTO()


end
 
%FOR ALL EXPERIMENTS FOR THIS MOUSE ON THIS DAY DO: 
%   cd(folderpath)

fprintf("Finished with mouse... \n");
end 

% FOR ALL EXPERIMENTS FOR ALL MICE ON THIS DAY DO: 
cd(dir4_path)
 fprintf("Done! \n");
 
% close 
% cd(folderpath);
end 

end 



% FOR EACH OF THE FOLLOWING EXPERIMENTS DO: 

% 1 - Generate log_data table. This is not saved but is used by make_DATA.m
%                 fprintf("Sorting logs...\n");
%                 make_log_data_LOOM_BOX();
   

% 2 - Generate: 'LOGDATA...', 'ROWS...', 'STIMROWS..' and 'DATA...' 
%                fprintf("Making LOGDATA, ROWS, STIMROWS and DATA...\n");
%                 make_DATA_RF_Sync();
%                 make_DATA_samefps();
   %  make_DATA_noARD();
%       make_loom_distance_heatmap()
%       make_stim_plot_video_BOX()

 % 3 - Track mouse - with light background. This generates 'ORIENTDATA'
 
%    fprintf("Making TRACKINGDATA...\n"); 
%    Track_Mouse_Position_Box();
%   trackmouselight()
  % OR 
%    fprintf("Making TRACKINGDATA...\n");
    % Track_Mouse_Position_fromJPG()
 
  % 4- Create table with All the data for the entire video 
%         fprintf("Merging Info and TRACKING...\n");
%        make_ALLDATA()

%   % 5- Make EXP_SUMMARY_TABLE, ALL_DIST and two plots. Line of speed and
%   % bar of distance per stimulus. 
%      fprintf("Making EXP_SUMMARY_TABLE, ALL_DIST and two plots...\n");
%      make_EXP_SUM_TABLE()

%      fprintf("Making stim video...\n");
   %   make_stim_plot_video()
%     make_lineplotwithrect()

        %fprintf("Making Video...\n");
       % make_maskedAVI()
   
%    fprintf("Adding DLC... \n");
%     Add_DLC()
         
%    fprintf("Plotting DLC... \n");
%    make_DLC_Plots()

%    fprintf("Making Clock/AntiClock... \n");
%    make_clock_anticlock()

%   fprintf("EXP DONE \n");

%   fprintf("Making HeatMap \n");
%   make_heatmap()

%      fprintf("Making response array and loom rows.. \n");
%      make_response_array_and_loom_rows()

% fprintf("Adjust coords...\n");
% adjust_box_coords()



% To move EXP_SUM_TAB  and ALLDATA to day folder. 
% load(strcat('EXP_SUM_TAB_DLC_', exp_name, '.mat'), 'EXP_SUMMARY_TABLE');
% save(fullfile(dir2_path,strcat('EXP_SUM_TAB_DLC_', exp_name, '.mat')), 'EXP_SUMMARY_TABLE'); 
%load(strcat('ALLDATA_', exp_name, '.mat'), 'ALLDATA');
%save(fullfile(dir2_path,strcat('ALLDATA_', exp_name, '.mat')), 'ALLDATA');


%For trajectory plotting - take a look at patch for interpolated colours 
% https://uk.mathworks.com/help/matlab/ref/patch.html

%   fprintf("Making all day summaries \n");
%   make_allday()

%     fprintf("Analysing trajectory...\n");
%     analyse_traj_comp()

%  fprintf("Moving XY...\n");
% move_xy_array()


%% 

%     if D == "19RB01" || D == "19RB02" 
%     currD = directory(k).name;
%     inpath = strcat(dir_path, '\', currD);% Windows        
%     animal_num = D; %dir_path(end-5:end);
%     exp_date = dir_path(end-12:end-7);
%     exp_name = strcat(exp_date, '_', animal_num, '_', currD);
%     exp_name_short = currD;    
%        
%     else
%     currD = directory(k).name;
%     inpath = strcat(dir_path, '\', currD);% Windows
% %     inpath = strcat(dir_path, '/', currD);% Linux
%     animal_num = D; %dir_path(end-4:end);
%     exp_date = dir_path(end-11:end-6);
%     exp_name_short = currD;
%     exp_name = strcat(exp_date, '_', animal_num, '_', currD);
%     end



%% 

% GAIA BOX EXPERIMENTS
% % % 
% fprintf("Making xy analysis...\n");
% make_xy_analysis()
% make_xy_analysis_distractor_checker()

%  fprintf("Making dlc loom array...\n");
% make_DLC_LOOM_array()
% MAKE SURE THE DLC/LOOMS/POS AND SPEED FOLDERS ARE CREATED ANDDDDD ADDED
% TO PATH!!! 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % SORT LOG FILESs
% % % % % 
%                        fprintf("Sorting logs...\n");
%                          make_log_data_LOOM_BOX();
% % % % % %                         make_log_data_DISTRACTORS()
%             
%  
%                         fprintf("Making LOGDATA, ROWS, STIMROWS and DATA...\n");
%                          make_DATA_RF_Sync();
% % %                          
%                         fprintf("Adding ILI...\n");
%                         add_ILI(); 
% % % % % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % SORT LOG FILES - OPTO SETUP 
% % % % 
%                        fprintf("Sorting logs...\n");
%                          make_log_data_LOOM_BOX();
%                          
%                         fprintf("Making LOGDATA, ROWS, STIMROWS and DATA...\n");
%                          make_DATA_OPTO();
%                          
%                        fprintf("Making TRACKINGDATA...\n"); 
%                          Track_Mouse_OPTO();
                    
%                        fprintf("Analysing...\n"); 
% %                          make_xy_analysis_OPTO2()
%                          make_xy_loom_array_OPTO() 

% 
%                          fprintf("Making video...\n");
%                          make_stim_plot_video_OPTO()
%                 

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % TRACKING  - option 1 0
% % % % % 
%                    fprintf("Making TRACKINGDATA...\n"); 
% %                   Track_Mouse_Position_Box();
%                   Track_Mouse_Position_Box2(); % Using different method -  hopefully faster and more robust
% % 
%                   fprintf("Merging Info and TRACKING...\n");
%                   merge_Info_TRACKING()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % TRACKING  -
%                     fprintf("Starting to make video...\n");
%                     cleanup_jpgs() %This makes mp4 video files from JPG images. 
% % % % % %                      option 2

%                     fprintf("Cleaning DLC\n");
%                     clean_dlc_csv(346) %Variable is image size. 518/346 
% 
%                     fprintf("Merging Info and DLC...\n");
%                     merge_Info_DLC()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % ANALYSIS 
% % % 
%        fprintf("Making xy_array.. \n");
%        make_xy()
% % % % % % % % % % % %     
%       fprintf("Making loom array..\n");
%        make_loom_array()
% % 
% % %       REQUIRES NEW XY_ARRAY WITH DIST FROM CENRTRE OF ARENA. 
%       fprintf("Making activity array.. \n");
%       make_ACTIVITY_ARRAY()
% % % 
% %      fprintf("MAKE RESPONSE AND LOCO...\n");
% %      make_RESPONSE_TYPE_and_LOCO_INDEX()
% %  
%      fprintf("Making xy loom array..\n");
%      make_xy_loom_array()
% % % % % % % % % %  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % VIDEOS AND PLOTS      

%     fprintf("Making stim video...\n");
%       make_stim_plot_video_LOOM()

%         fprintf("Making speed/dist plot\n");
%         make_PLOT_Speed_DFS_Looms()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % IF CRASHEDDDD. 

%   fprintf("Sorting logs...\n");
%   make_log_data_LOOM_BOX();
  
%      fprintf("Making TRACKINGDATA...\n"); 
%      Track_Mouse_Position_Box();

% fprintf("Making xy_array.. \n");
%        make_xy_CRASH()

%      fprintf("Making xy loom array..\n");
%      make_xy_loom_array_CRASH()
