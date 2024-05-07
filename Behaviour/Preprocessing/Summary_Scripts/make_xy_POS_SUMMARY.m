function make_xy_POS_SUMMARY()

% Make sure to start this script in 'C:\Data_analysis\DATA\Setd5\ALL_DAY_SUMMARIES\LOOMS\XY_LOOM_ALL'
% Combines information from XY_LOOM - (x,y) positions over loom period. 

xyloom_files = dir('XY_LOOM_*'); 
n_xyl_files = numel(xyloom_files); 
ALL_XYLOOM = {}; 
% 
for i = 1:n_xyl_files

    name = xyloom_files(i).name; 
    load(name, 'XY_LOOM')
    
    animal = name(end-17:end-12); 
    date = name(end-24:end-19);
    exp = name(end-10:end-4);
    
%         % FOR MULTILOOM 
%     animal = name(end-22:end-17); 
%     date = name(end-29:end-24);
%     exp = name(end-15:end-4);
   
    numrows = numel(XY_LOOM(:,1))/2; 
    
    xyl_info =cell(numrows, 1564);

    for q = 1:numrows
    xyl_info(q,1)= {date};
    xyl_info(q,2)= {animal};
    xyl_info(q,3)={exp};
        for p = 1:780
        xyl_info(q,4+p)={XY_LOOM(q,p)};
        xyl_info(q,784+p)={XY_LOOM(q+1,p)};
        end 
    end 
  
    ALL_XYLOOM = vertcat(ALL_XYLOOM, xyl_info);
end 

% Add genotype to column 5. 

h =numel(ALL_XYLOOM(:,1));

% het_animals = ["MJ2437", "MJ2440", "MJ1934", "MJ1931", "MJ2167", "MJ2168", "MJ2169"]; 
% het_animals = ["GN7375", "GN7376", "GN7390", "GN7394", "GN7398", "GN7382"]; 
%  het_animals= ["GN2593", "GN2594", "MJ1861", "MJ1864", "MJ1867"];
%  het_animals = ["GN2375", "GN2377", "GN2382", "GN2628", "GN2637", "GN2754", "GN2829", "GN2832", "GN2833", "GN2900", "GN2901", "GN2902"]; 
% het_animals = ["GN2593", "GN2594", "GN1385", "GN1386", "GN1388", "GN1394"];
% het_animals = ["GN2708", "GN2709", "GN2868", "GN4369", "GN4373", "MJ0310", "MJ0311"];
het_animals = ["GN5148", "GN5596"];


for j = 1:h
          if contains(ALL_XYLOOM{j,2}, het_animals)
               ALL_XYLOOM{j,4} = "het";
          else 
               ALL_XYLOOM{j,4} = "wt";
          end 
end 
 
% for j = 1:h
%       if ALL_XYLOOM{j,2}=="GN3959" || ALL_XYLOOM{j,2}=="GN4244" || ALL_XYLOOM{j,2}=="GN4473" ||  ALL_XYLOOM{j,2}=="GN6560" || ALL_XYLOOM{j,2}=="GN6562" || ALL_XYLOOM{j,2}=="GN6610" || ALL_XYLOOM{j,2}=="GN7269" || ALL_XYLOOM{j,2}=="GN7476" || ALL_XYLOOM{j,2}=="GN7614" || ALL_XYLOOM{j,2}=="GN7790"
%         ALL_XYLOOM{j,4} = {'het'};
%     else
%         ALL_XYLOOM{j,4} = {'wt'};
%     end 
% end
% % 
%  
% for j = 1:h 
%         ALL_XYLOOM{j,4} = {'PV'}; 
% end
% % 
% %% Make table from array
% 
Date = ALL_XYLOOM(:,1); 
Animal = ALL_XYLOOM(:,2);
Exp = ALL_XYLOOM(:,3);
Geno = ALL_XYLOOM(:,4);
X = ALL_XYLOOM(:,5:784);
Y = ALL_XYLOOM(:,785:1564);
% 
ALL_XYLOOM_POS_TABLE = table(Date, Animal, Exp, Geno, X, Y);
% 
save(strcat('ALL_XYLOOM_POS_TABLE_Setd5_Julie2.mat'), 'ALL_XYLOOM_POS_TABLE')
clear

% save(strcat('ALL_XYLOOM_POS_SERT.mat'), 'ALL_XYLOOM')



















%% FOR ENTIRE VIDEOS




%% Activity files
xytable_files = dir('XY_TABLE_*'); 
% xytable_files =dir;
n_xyl_files = numel(xytable_files); 
ALL_XY = {}; 

for i = 1:n_xyl_files

    name = xytable_files(i).name; 
    load(name, 'xy_table')
%     
%    animal = name(end-17:end-12); 
%     date = name(end-24:end-19);
%     exp = name(end-10:end-4);
%     row_end = 12000;


    % For banana / acclim 
%     animal = name(end-19:end-14); 
%     date = name(end-26:end-21);
%     exp = name(end-12:end-4);
%     row_end = 12000;
%     
%         % FOR MULTILOOM 
    animal = name(end-22:end-17); 
    date = name(end-29:end-24);
    exp = name(end-15:end-4);
    row_end = 14700;

%     numrows = numel(xy_table(:,1)); 
    
    xyl_info =cell(1,5);

    xyl_info(1,1)= {date};
    xyl_info(1,2)= {animal};
    xyl_info(1,3)={exp};
    xyl_info(1,5)={xy_table.X(1:row_end)'};
    xyl_info(1,6)={xy_table.Y(1:row_end)'};
  
    ALL_XY = vertcat(ALL_XY, xyl_info);
end 

% Add genotype to column 5. 


 het_animals= ["GN2593", "GN2594", "MJ1861", "MJ1864", "MJ1867"]; 

for j = 1:n_xyl_files
          if contains(ALL_XY{j,2}, het_animals)
               ALL_XY{j,4} = "het";
          else 
               ALL_XY{j,4} = "wt";
          end 
end 


% for j = 1:n_xyl_files
%     if ALL_XY{j,2}=="GN3959" || ALL_XY{j,2}=="GN4244" || ALL_XY{j,2}=="GN4473" ||  ALL_XY{j,2}=="GN6560" || ALL_XY{j,2}=="GN6562" || ALL_XY{j,2}=="GN6610" || ALL_XY{j,2}=="GN7269" || ALL_XY{j,2}=="GN7476" || ALL_XY{j,2}=="GN7614" || ALL_XY{j,2}=="GN7790"
%         ALL_XY{j,4} = {'het'};
%     else
%         ALL_XY{j,4} = {'wt'};
%     end 
% end

% for j = 1:n_xyl_files
%     if ALL_XY{j,2}=="GN7184" || ALL_XY{j,2}=="GN7751" || ALL_XY{j,2}=="GN7154" 
%         ALL_XY{j,4} = {'STOP'};
%     else
%         ALL_XY{j,4} = {'wt'};
%     end 
% end

% for j = 1:n_xyl_files
%         ALL_XYLOOM{j,4} = {'SERT'}; 
% end

%% Make table from array

Date = ALL_XY(:,1); 
Animal = ALL_XY(:,2);
Exp = ALL_XY(:,3);
Geno = ALL_XY(:,4);
X = ALL_XY(:,5);
Y = ALL_XY(:,6);

ALL_XY_POS = table(Date, Animal, Exp, Geno, X, Y);

save(strcat('ALL_XY_POS_LOOMS_COHORT3_MULTI.mat'), 'ALL_XY_POS')

clear all 
end