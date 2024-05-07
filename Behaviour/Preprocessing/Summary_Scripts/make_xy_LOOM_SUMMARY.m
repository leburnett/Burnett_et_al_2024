function make_xy_LOOM_SUMMARY()

% Make sure to start this script in 'C:\Data_analysis\DATA\Setd5\ALL_DAY_SUMMARIES\LOOMS\xy_loom_arrays'
%Combines information from xy_loom - response speed over time 3s before the
%loom til 10s after the loom. 
% details = ALL_XYLOOM_TABLE(:, [1,2,3,4,6,7])
%% Activity files

% Set diff_contrasts 
diff_contrasts = 0; 

xyloom_files = dir('xy_loom_*'); 
n_xyl_files = numel(xyloom_files); 
ALL_XYLOOM = {}; 

for i = 1:n_xyl_files

    name = xyloom_files(i).name; 
    load(name, 'xy_loom')
    
    animal = name(end-17:end-12); 
    date = name(end-24:end-19);
    exp = name(end-10:end-4);
   
%     % FOR MULTILOOM 
%     animal = name(end-22:end-17); 
%     date = name(end-29:end-24);
%     exp = name(end-15:end-4);


    numrows = numel(xy_loom(:,1)); 
    
    if diff_contrasts ==0
        xyl_info =cell(numrows, 784);
        
        for q = 1:numrows
            xyl_info(q,1)= {date};
            xyl_info(q,2)= {animal};
            xyl_info(q,3)={exp};
            for p = 1:780
                xyl_info(q,4+p)={xy_loom(q,p)};
            end
        end
        
    elseif diff_contrasts == 1
        xyl_info =cell(numrows, 785);
        
        for q = 1:numrows
            xyl_info(q,1)= {date};
            xyl_info(q,2)= {animal};
            xyl_info(q,3)={exp};
            for p = 1:781
                xyl_info(q,4+p)={xy_loom(q,p)};
            end
        end
    end
        
  
    ALL_XYLOOM = vertcat(ALL_XYLOOM, xyl_info);
end 

% Add genotype to column 5. 

h =numel(ALL_XYLOOM(:,1));
% het_animals = ["MJ2159", "MJ2162", "MJ2172", "MJ2179"]; 
% het_animals = ["GN8604", "GN8650", "GN8738"]; 
% het_animals = ["MJ2437", "MJ2440"]; 
% het_animals = ["MJ0086", "MJ0088", "MJ0090", "MJ0091", "MJ0093"]; 
%  het_animals = ["GN1385", "GN1386", "GN1388", "GN1394"]; 
%  het_animals = ["MJ0423", "MJ0424", "MJ0427", "MJ0428", "MJ0429", "MJ0430"]; 
%  het_animals = ["MJ0374", "MJ0380", "MJ0381", "MJ0382"]; 
%   het_animals = ["MJ0380", "MJ0382"]; % wrongly labelled animals 
%   het_animals= ["GN2593", "GN2594", "MJ1861", "MJ1864", "MJ1867", "GN1388"]; %ptchd1 C2
% het_animals = ["MJ0853", "MJ0593", "MJ0595", "M21125",  "M21120", "M21121", "M21282"];
% het_animals = ["GN2708", "GN2709", "GN2868", "GN4369", "GN4373"];
% het_animals = ["GN5148", "GN5596", "GN5542", "GN5597"];
% het_animals = ["MJ2479", "MJ2482", "MJ2688", "MJ2690", "MJ0310", "MJ0311"]; 
%     het_animals = ["MJ0991,", "MJ0993", "MJ0995", "MJ0997", "GN3903", "MJ1970"];
%  het_animals = ["MJ0344", "MJ0345", "MJ0347", "MJ0349"];
% het_animals = ["220300", "220301", "220580"];
%  het_animals = ["220777", "220820", "220821", "220822"];
het_animals = ["M22080", "M22081", "M22085", "M22341"];

for j = 1:h
          if contains(ALL_XYLOOM{j,2}, het_animals)
               ALL_XYLOOM{j,4} = "het";
          else 
               ALL_XYLOOM{j,4} = "wt";
          end 
end 
% 
% %Setd5 Animals 
% het_animals = ["GN3959", "GN4244", "GN4473", "GN6560", "GN6562", "GN6610","GN7269", "GN7476", "GN7614", "GN7790"]; 
% 
% % STOP animals
% het_animals = ["GN7375", "GN7376", "GN7390", "GN7394", "GN7398", "GN7382"]; 
% 
% het_animals = ["MJ1931", "MJ1934"]; 
% 
% 
% % %Cul3 C1 +C2
% het_animals = ["GN2375", "GN7277", "GN2382", "GN2628", "GN2637", "GN2754", "GN2829", "GN2832", "GN2900", "GN2901", "GN2902"]; 

% 
% for j = 1:h
%     if ALL_XYLOOM{j,2}=="20PV18" || ALL_XYLOOM{j,2}=="20PV21" 
%         ALL_XYLOOM{j,785} = {'CHERRY'};
%     else
%         ALL_XYLOOM{j,785} = {'GQ'};
%     end 
% end 
% 
% for j = 1:h
%     if ALL_XYLOOM{j,1}=="201113" || ALL_XYLOOM{j,1}=="201115" 
%         ALL_XYLOOM{j,786} = {'Saline'};
%     else
%         ALL_XYLOOM{j,786} = {'CNO'};
%     end 
% end 
% for j = 1:h
%     if ALL_XYLOOM{j,1}=="210120" || ALL_XYLOOM{j,1}=="210122" 
%         ALL_XYLOOM{j,785} = {'DC'};
%     else
%         ALL_XYLOOM{j,785} = {'ILI'};
%     end 
% end 

%% Make table from array

Date = ALL_XYLOOM(:,1); 
Animal = ALL_XYLOOM(:,2);
Exp = ALL_XYLOOM(:,3);
Geno = ALL_XYLOOM(:,4);
Speed = ALL_XYLOOM(:,5:784);

if diff_contrasts == 1
    Contrast = ALL_XYLOOM(:, 785);
    ALL_XYLOOM_TABLE = table(Date, Animal, Exp, Geno, Speed, Contrast); %, Virus, CNO, ExpType
else
   ALL_XYLOOM_TABLE = table(Date, Animal, Exp, Geno, Speed);  
end 
% ExpType = ALL_XYLOOM(:,785);
% Virus = ALL_XYLOOM(:,785);
% CNO = ALL_XYLOOM(:,786);


% Add MaxSpeed and TimeToMax 

for k = 1:h
	G = cell2mat(ALL_XYLOOM_TABLE{k,5}); 
    G = G(180:end); %exclude time before loom.
	maxsp = max(G);
	t2m = find(G == max(G)); %time in frames (/60 to get in seconds)
	ALL_XYLOOM_TABLE.MaxSp{k} = maxsp;
	ALL_XYLOOM_TABLE.T2M{k} = t2m;
end 

save(strcat('ALL_XYLOOM_TABLE_Dafna_D1_.mat'), 'ALL_XYLOOM_TABLE');
% save(strcat('ALL_XYLOOM_STOP.mat'), 'ALL_XYLOOM')

clear 

end