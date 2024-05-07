function merge_Info_TRACKING()
%% Merge information from 'Info' and 'TRACKING_DATA'.
% Updates the cell array 'Info' by adding the info from tracking. 
% Also creates table with the same data called 'Info_Table'. 
% Created by Burnett 07/01/20. 

global exp_name 

%Load Info - contains jpg, dot pos, mouse 'live' pos... 
load(strcat('INFO_', exp_name, '.mat'), 'Info');
load(strcat('TRACKINGDATA_', exp_name, '.mat'), 'TRACKING_DATA');

%Load TRACKING_DATA - contains the variables acquired from 'Track_Mouse..."
% Col 2: Major Axis
% Col 3: Centroid x 
% Col 4: Centroid y 
% Col 5: Phi 
% Col 6: Orientation 

% If camera has stopped before the end of the stimulus presentation. 
% for j = 1:6
%     for i = 7367:12883
%         TRACKING_DATA(i,j)={0};
%     end 
% end 

%Might not need this line. Only if want to specify exactly which columns to
%choose. 
 Info = Info(:, 1:13); 
 l_info = length(Info); 
%  Info = Info(:, 1:11); 

for i = 2:6 
A = TRACKING_DATA(1:l_info,i);
Info = [Info A];
end

%Make Info into table. 
    Info_Table = cell2table(Info);
%     Info_Table.Properties.VariableNames = {'JPG', 'DATA_ALL_VBL_ROW', 'StimulusName', 'Radius', 'Loom', 'Mouse_x', 'Mouse_y', 'Trigger', 'MouseNum', 'Date', 'Stimulus', 'Major', 'CentroidX', 'CentroidY', 'Phi', 'Orient'}; 

      Info_Table.Properties.VariableNames = {'JPG', 'DATA_ALL_VBL_ROW', 'StimulusName', 'Radius', 'Loom', 'Mouse_x', 'Mouse_y', 'Trigger','BKG_RGB', 'DOT_RGB', 'MouseNum', 'Date', 'Stimulus', 'Major', 'CentroidX', 'CentroidY', 'Phi', 'Orient'}; 

% SaveInfo as both a cell and a table.  
save(strcat('INFOTABLE_', exp_name, '.mat'), 'Info_Table');
save(strcat('INFO_', exp_name, '.mat'), 'Info');
clearvars
end 

% Info is cell 
% TRACKING_DATA is cell