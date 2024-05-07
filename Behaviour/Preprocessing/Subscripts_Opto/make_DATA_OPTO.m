function make_DATA_OPTO()
% Created by Burnett 14/10/20

% Script for OPTO experiments where RF only presented when 'opto' should be
% switched on. 

% Much simplified form of 'make_DATA_RF_Sync.m'

% This code:
% 1 - Finds the number of 'triggers' sent to the camera to take images. 
%       These are the # of '1's in the 3rd Column of 'ard'
% 2 - Find the number of images taken by Basler. 
%        THESE TWO NUMBERS SHOULD MATCH
% 3 - Find the row values in 'ard' where RFs were recorded by the photodiode. (Index_RFs) - this is when OPTO was stimulated. 
% 4 - Find the number of images BEFORE each of these values. 
% 5 - Generate 'Info' table with Basler image names in the first column.
% Column 2 is 0 normally and 1 when the opto light is on. 

global inpath exp_name 

%  inpath = '/Users/lauraburnett/Data_Analysis_Mac/MATLAB/TroubleShooting/RedFrames/02_Loom';
% exp_date = '191206';
% animal_num = '19S09';
% exp_name = '02_Loom';

%% RED FRAMES RECORDED BY PHOTODIODE - - - - - - - Find out how many RFs the photodiode recorded from the arduino file. 
% inpath = 'C:\Data_analysis\DATA\Setd5\191216\GN3959\01_Loom';
ardfile = dir(fullfile(inpath,'*.txt')); %gets txt files in folder
ard_name = ardfile(1).name;
ard = dlmread(ard_name, ',', 60,0); %Ignore the first 50 lines. 

%% FILTER ARDUINO FILE - DORIC SHOULD ALWAYS START BEFORE MATLAB. 

% Make any value <1000 in the 4th column == 0
ard(ard(:,4)<600,4) = 0;
ard(ard(:,4)>=600,4) = 100;

% %Generate sequence of zeroes and a 100 and isolate where this is found in the fourth column of the 'ard' file. 
seq1 = zeros(30,1);
seq1(31, 1) = 100;

col4_trans = ard(:,4)'; %Transpose column 4
seq1 = seq1'; %Transpose seq
Index = strfind(col4_trans, seq1);
Index_RFs = Index';
Index_RFs = Index_RFs+30; %Array of ‘ard’ row values of all RFs. - these are when the RED LIGHT IS ON And OPTO switched on.  

if isempty(Index_RFs)
    n_RFs = 0;
    OPTO_STIMULATED = 0; 
else
    n_RFs = numel(Index_RFs(:,1));
    OPTO_STIMULATED = 1; 
    disp(['OPTO', num2str(n_RFs)])
end

%% Find the # of images to be used. 

% Find number of images 'triggered' by arduino    
n_ims_triggered = numel(ard(ard(:,3)==1));

% Load the basler folder
baslerpath = fullfile(inpath, '/basler_*');
myBaslerFiles = dir(fullfile(baslerpath,'*.jpg'));

[~,index] = sortrows({myBaslerFiles.date}.'); 
myBaslerFiles = myBaslerFiles(index); 
clear index %% THIS IS IMPORTANT. SOMETIMES TIME of IMAGE e.e. 'img_903305.jpg' reaches 10 then starts low again. This means images are not in chron order.

num_bas = numel(myBaslerFiles);

if num_bas > n_ims_triggered
    num_to_use = n_ims_triggered;
    excess = 1;
else 
     num_to_use = num_bas;
     excess = 0;
end


%% MAKE TABLE 'INFO'   
% Define cell 'Info'
Info = cell(num_to_use,2);

for i = 1:num_to_use %num_bas-1
     Info{i,1} = {myBaslerFiles(i+excess).name}; % (i+1) because the first image is taken before. 
     Info{i,2} = 0;     
end 


%% If the light was stimulated at some point. 
if OPTO_STIMULATED == 1
    
    for q = 1:n_RFs 
        Index_RFs(q,2) = numel(find(ard(1:Index_RFs(q,1),3)==1)); %Number of Images triggered BEFORE each RF.
        
        for w = 1:n_RFs
            row = Index_RFs(w,2)+1; % Plus one because the frame with the trigger is the one after the # of RFs before.
            Info{row,2} = 1;
        end
    end
    
end


if Info{1,2}==1
    Info{1,2}=0; 
end 
%% SAVE THE ARRAYS 
     
%Save the 'Info' cell array as a .mat file. 
save(strcat('INFO_', exp_name, '.mat'),'Info');

clearvars

end

