function make_LOOM_SUMMARY()

% In folder "C:\Data_analysis\DATA\Setd5\ResponseArrays\LoomArrays\"
% Can already make one giant table. Do not need to make individual day
% summaries. 

%% Loom files

loom_files = dir('LOOM_ARRAY*'); 
ALL_LOOMS = {}; 
n_loom_files = numel(loom_files); 

for i = 1:n_loom_files

    name = loom_files(i).name; 
    load(name, 'loom_array')
    
    animal = name(end-17:end-12); 
    date = name(end-24:end-19);
    exp = name(end-10:end-4);
 
    loom_info ={};

    loom_info(1,1)= {date};
    loom_info(1,2)= {animal};
    loom_info(1,3)={exp};
    loom_info(1,4)={loom_array(1,1)};
    loom_info(1,5)={loom_array(1,2)};
    loom_info(1,6)={loom_array(1,3)};
    loom_info(1,7)={loom_array(1,4)};
    loom_info(1,8)={loom_array(1,5)};
    loom_info(1,9)={loom_array(1,6)};
   
    ALL_LOOMS = vertcat(ALL_LOOMS, loom_info);
    
    %% If MULTILOOM 
    
%      name = loom_files(i).name; 
%     load(name, 'loom_array')
%     
%     animal = name(end-22:end-17); 
%     date = name(end-29:end-24);
%     exp = name(end-15:end-4);
%     
%     loom_info ={};
% 
%     loom_info(1,1)= {date};
%     loom_info(1,2)= {animal};
%     loom_info(1,3)={exp};
%     loom_info(1,4)={loom_array(1,1)}; %bout total 
%     loom_info(1,5)={loom_array(1,2)}; % Bouts Q1
%     loom_info(1,6)={loom_array(1,3)}; % Bouts Q2
%     loom_info(1,7)={loom_array(1,4)}; % Bouts Q3
%     loom_info(1,8)={loom_array(1,5)}; % Bouts Q4
%     loom_info(1,9)={loom_array(1,6)}; % Looms Total
%     loom_info(1,10)={loom_array(1,7)}; % Looms Q1
%     loom_info(1,11)={loom_array(1,8)}; % Looms Q2
%     loom_info(1,12)={loom_array(1,9)}; % Looms Q3
%     loom_info(1,13)={loom_array(1,10)}; % Looms Q4
%    
%     ALL_LOOMS = vertcat(ALL_LOOMS, loom_info);
%     
end 

% Add genotype to column 5. 

h =numel(ALL_LOOMS(:,1));

% for j = 1:h
%     if ALL_LOOMS{j,2}=="GN3959" || ALL_LOOMS{j,2}=="GN4244" || ALL_LOOMS{j,2}=="GN4473" ||  ALL_LOOMS{j,2}=="GN6560" || ALL_LOOMS{j,2}=="GN6562" || ALL_LOOMS{j,2}=="GN6610" || ALL_LOOMS{j,2}=="GN7269" || ALL_LOOMS{j,2}=="GN7476" || ALL_LOOMS{j,2}=="GN7614" || ALL_LOOMS{j,2}=="GN7790"
%         ALL_LOOMS{j,10} = {'het'};
%     else
%         ALL_LOOMS{j,10} = {'wt'};
%     end 
% end

% If Multiloom 
% for j = 1:h
%     if ALL_LOOMS{j,2}=="GN3959" || ALL_LOOMS{j,2}=="GN4244" || ALL_LOOMS{j,2}=="GN4473" ||  ALL_LOOMS{j,2}=="GN6560" || ALL_LOOMS{j,2}=="GN6562" || ALL_LOOMS{j,2}=="GN6610" || ALL_LOOMS{j,2}=="GN7269" || ALL_LOOMS{j,2}=="GN7476" || ALL_LOOMS{j,2}=="GN7614" || ALL_LOOMS{j,2}=="GN7790"
%         ALL_LOOMS{j,14} = {'het'};
%     else
%         ALL_LOOMS{j,14} = {'wt'};
%     end 
% end

% FOR SERTS /WTS
% for j = 1:h
%     ALL_LOOMS{j,10} = "C57";  %{"STOP"}; %{"SERT"};  %
% end 

% Cul3

% for j = 1:h
%     if ALL_LOOMS{j,2}=="GN2375" || ALL_LOOMS{j,2}=="GN2377" || ALL_LOOMS{j,2}=="GN2382" ||  ALL_LOOMS{j,2}=="GN2628" || ALL_LOOMS{j,2}=="GN2637" 
%         ALL_LOOMS{j,14} = {'het'};
%     else
%         ALL_LOOMS{j,14} = {'wt'};
%     end 
% end 


% %Cul3_C2
% for j = 1:h
%     if ALL_LOOMS{j,2}=="GN2754" || ALL_LOOMS{j,2}=="GN2829" || ALL_LOOMS{j,2}=="GN2832" ||  ALL_LOOMS{j,2}=="GN2900" || ALL_LOOMS{j,2}=="GN2901" || ALL_LOOMS{j,2}=="GN2902" 
%         ALL_LOOMS{j,10} = {'het'};
%     else
%         ALL_LOOMS{j,10} = {'wt'};
%     end 
% end 

for j = 1:h
    if ALL_LOOMS{j,2}=="20PV01" || ALL_LOOMS{j,2}=="20PV03" || ALL_LOOMS{j,2}=="20PV04" || ALL_LOOMS{j,2}=="20PV05" || ALL_LOOMS{j,2}=="20PV08" ||  ALL_LOOMS{j,2}=="20PV09" ||  ALL_LOOMS{j,2}=="20PV13" ||  ALL_LOOMS{j,2}=="20PV14" ||  ALL_LOOMS{j,2}=="20PV16" ||  ALL_LOOMS{j,2}=="20PV17"  
        ALL_LOOMS{j,10} = {'DREADD'};
    else
        ALL_LOOMS{j,10} = {'CHERRY'};
    end 
end 


% for j = 1:h
%     if ALL_LOOMS{j,2}=="20GP02" || ALL_LOOMS{j,2}=="20GP03" || ALL_LOOMS{j,2}=="20GP05"
%         ALL_LOOMS{j,10} = {'DREADD'};
%     else
%         ALL_LOOMS{j,10} = {'CHERRY'};
%     end 
% end 


%% Make table from array

 ALL_LOOMS_TABLE = cell2table(ALL_LOOMS, 'VariableNames', {'Date', 'Animal', 'Exp', 'NumBouts', 'NumBoutsL1', 'NumBoutsL2', 'NumLooms', 'NumLoomsL1', 'NumLoomsL2','Geno'});

% IF MULTILOOM 
 ALL_LOOMS_TABLE = cell2table(ALL_LOOMS, 'VariableNames', {'Date', 'Animal', 'Exp', 'NumBouts', 'NumBoutsQ1', 'NumBoutsQ2','NumBoutsQ3', 'NumBoutsQ4', 'NumLooms', 'NumLoomsQ1', 'NumLoomsQ2','NumLoomsQ3', 'NumLoomsQ4','Geno'});

save(strcat('ALL_LOOMS_TABLE_PV_DREADD_C4.mat'), 'ALL_LOOMS_TABLE')
clear all

end 