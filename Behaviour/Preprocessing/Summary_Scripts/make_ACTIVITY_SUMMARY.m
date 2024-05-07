function make_ACTIVITY_SUMMARY()

% Make sure to start this script in 'C:\Data_analysis\DATA\Setd5\ALL_DAY_SUMMARIES\ACTIVITY\activity_arrays\'

%% Activity files
activity_files = dir('activity_array*'); 
ALL_ACTIVITY = {}; 
n_act_files = numel(activity_files); 

for i = 1:n_act_files

    name = activity_files(i).name; 
    load(name, 'activity_array')
%     
%     animal = name(end-17:end-12); 
%     date = name(end-24:end-19);
%     exp = name(end-10:end-4);

    %Banana or Acclim
    animal = name(end-19:end-14); 
    date = name(end-26:end-21);
    exp = name(end-12:end-4);

    % FOR MULTILOOM 
%     animal = name(end-22:end-17); 
%     date = name(end-29:end-24);
%     exp = name(end-15:end-4);

   
    act_info ={};

    act_info(1,1)= {date};
    act_info(1,2)= {animal};
    act_info(1,3)={exp};
    act_info(1,4)={activity_array(1,1)};
    act_info(1,5)={activity_array(1,2)};
    act_info(1,6)={activity_array(1,3)};
    act_info(1,7)={activity_array(1,4)};
    act_info(1,8)={activity_array(1,5)};
    act_info(1,9)={activity_array(1,6)}; %Slow 
    act_info(1,10)={activity_array(1,7)}; % OUT BOX. 
    act_info(1,11)={activity_array(1,8)}; % CENTRE 
    act_info(1,12)={activity_array(1,9)}; % EDGE 
    act_info(1,13)={activity_array(1,10)};
    act_info(1,14)={activity_array(1,11)};
    act_info(1,15)={activity_array(1,12)};
    act_info(1,16)={activity_array(1,13)}; % SLOW
    act_info(1,17)={activity_array(1,14)}; %OUT BOX 
    act_info(1,18)={activity_array(1,15)}; %CENTRE 
    act_info(1,19)={activity_array(1,16)}; %EDGE
   
    ALL_ACTIVITY = vertcat(ALL_ACTIVITY, act_info);
end 

% Add genotype to column 5. 

h =numel(ALL_ACTIVITY(:,1));
% het_animals = ["MJ2437", "MJ2440", "MJ1934", "MJ1931", "MJ2167", "MJ2168", "MJ2169"]; 
% het_animals = ["GN7375", "GN7376", "GN7390", "GN7394", "GN7398", "GN7382"]; 
% het_animals = ["MJ0086", "MJ0088", "MJ0090", "MJ0091", "MJ0093"]; 
%  het_animals = ["GN1385", "GN1386", "GN1388", "GN1394"]; 
%  het_animals= ["GN2593", "GN2594", "MJ1861", "MJ1864", "MJ1867"]; 
% het_animals = ["GN2868", "GN3903","GN3959", "GN4244", "GN4473", "GN2708", "GN2709"];
% het_animals = ["MJ0853", "MJ0593", "MJ0595", "GN4369", "GN4373"];
% het_animals = ["GN5148", "GN5596"];
het_animals = ["GN4369", "GN4373"];

for j = 1:h
          if contains(ALL_ACTIVITY{j,2}, het_animals)
               ALL_ACTIVITY{j,20} = "het";
          else 
               ALL_ACTIVITY{j,20} = "wt";
          end 
end 


% for j = 1:h
%     if ALL_ACTIVITY{j,2}=="GN3959" || ALL_ACTIVITY{j,2}=="GN4244" || ALL_ACTIVITY{j,2}=="GN4473" ||  ALL_ACTIVITY{j,2}=="GN6560" || ALL_ACTIVITY{j,2}=="GN6562" || ALL_ACTIVITY{j,2}=="GN6610" || ALL_ACTIVITY{j,2}=="GN7269" || ALL_ACTIVITY{j,2}=="GN7476" || ALL_ACTIVITY{j,2}=="GN7614" || ALL_ACTIVITY{j,2}=="GN7790"
%         ALL_ACTIVITY{j,20} = {'het'};
%     else
%         ALL_ACTIVITY{j,20} = {'wt'};
%     end 
% end


%% Make table from array

ALL_ACTIVITY_TABLE = cell2table(ALL_ACTIVITY, 'VariableNames', {'Date', 'Animal', 'Exp', 'RowsExp', 'TimeExp', 'DistExp', 'SpeedExp', 'MaxSpeed', 'PerSlow', 'PerOutBox', 'PerCentre', 'PerEdge', 'D10', 'Sp10', 'Max10', 'PerSlow10','PerOut10', 'PerCentre10', 'PerEdge10','Geno'});

save(strcat('ALL_ACTIVITY_TABLE_Ptchd1_C5_Acclim_Preloom.mat'), 'ALL_ACTIVITY_TABLE')

clear

% save(strcat('ALL_ACTIVITY_TABLE_COHORT2_ACCLIM.mat'), 'ALL_ACTIVITY_TABLE')
end 


%% BELOW IS NO LONGER NEEDED SINCE DONT NEED TO MAKE DAY SUMMARIES FIRST -  CAN MAKE ALL IN ONE. 

