function make_FULL_xy_loom_array()

%% Created by Burnett - 06/10/21

% Makes XY_FULL_TABLE for 3s before the loom starts til 10s after the loom
% starts. Makes 99 ms moving mean average of speed across this time. 

% XY_TABLE for JUST time period peri stimulus. 
% Each row is a trial. 
% Col 1 = Animal 

% Col 1 = X
% Col 2 = Y 
% Col 3 = Speed 
% Col 4 = Acc
% Col 5 = Dist C
% Col 6 = Dist Shelter. 
% Col 

global exp_name dir3_path 

if exist(strcat(dir3_path,'\', 'FULL_XY_LOOM_TABLE.mat'), 'file')
    load(strcat(dir3_path,'\', 'FULL_XY_LOOM_TABLE.mat'),'FULL_XY_TABLE')
end 

% Make empty cell array in which to enter data. 
if ~exist('FULL_XY_TABLE', 'var')
    FULL_XY_TABLE = table(); 
    current_size = 0;
     % Make longer and smoothed 'xy_loom' array as well
     xy_loom = []; 
else 
     current_size = length(FULL_XY_TABLE.Date); 
end 


if exist(strcat('ALL_LOOM_ROWS_', exp_name,'.mat'), 'file') 

    load((strcat('ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS'); 
    load(strcat('XY_table_', exp_name, '.mat'), 'xy_table');
    
    date = exp_name(1:6);
    animal = exp_name(8:13);
    exp = exp_name(15:end);
    
     loom_rows = ALL_LOOM_ROWS(1,:)'; 
     num_looms = numel(loom_rows);

     
     for i = 1:num_looms
         if isempty(FULL_XY_TABLE)
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
         
         if row3 > height(xy_table)
             row3 = height(xy_table);
         end
         
         % Start filling table
         
         FULL_XY_TABLE.Date{p} = date;
         FULL_XY_TABLE.Animal{p} = animal;
         FULL_XY_TABLE.Exp{p} = exp;
         FULL_XY_TABLE.Loom{p} = loom;
         FULL_XY_TABLE.LoomRow{p} = loom_rows(i,1);
         
         
          X = xy_table.X(row:row3);
          Y = xy_table.Y(row:row3);
          SPEED = xy_table.Speed(row:row3);
          ACC = xy_table.Acc(row:row3);
          DIST_C = xy_table.DistCentre(row:row3);
          DIST_SH = xy_table.DistShelter(row:row3);
          
         FULL_XY_TABLE.X{p} = X;
         FULL_XY_TABLE.Y{p} = Y;
         FULL_XY_TABLE.SPEED{p} = SPEED;
         FULL_XY_TABLE.ACC{p} = ACC;
         FULL_XY_TABLE.DIST_C{p} = DIST_C;
         FULL_XY_TABLE.DIST_SH{p} = DIST_SH;
         
%          het_animals = ["GN4369", "GN4373", "MJ1964", "MJ1970", "GN5148", "GN5596", "GN5542", "GN5597"];
%          het_animals = ["MJ2479", "MJ2482", "MJ2688", "MJ2690", "MJ0310", "MJ0311"]; 
%             het_animals = ["MJ0344", "MJ0345", "MJ0347", "MJ0349"];
%            het_animals = ["MJ0593", "MJ0595"];
%            het_animals = ["MJ0853", "MJ0593", "MJ0595"]; 
%             het_animals = ["M21125", "M21120", "M21121", "M21282"]; 
% Cul3
%             het_animals = ["GN2375", "GN2377", "GN2382", "GN2628","GN2637", "GN2829", "GN2832", "GN2833", "GN2901", "GN2902"]; %Cul3
% Ptchd1
%             het_animals = ["GN1385", "GN1386", "GN1388", "GN1394", "GN2593", "GN2594", "GN2708", "GN2709", "GN4369", "GN4373", "GN3903","MJ2479", "MJ2482", "MJ2688", "MJ2690", "MJ0310", "MJ0311"];
             het_animals = ["M21125", "M21120", "M21121", "M21282", "M21505", "M21507"];

         if contains(FULL_XY_TABLE.Animal{p}, het_animals)
             FULL_XY_TABLE.Geno{p} = "het";
         else
             FULL_XY_TABLE.Geno{p} = "wt";
         end
        
     end
     
     save(strcat(dir3_path,'\', 'FULL_XY_LOOM_TABLE.mat'), 'FULL_XY_TABLE')

end 

end        
