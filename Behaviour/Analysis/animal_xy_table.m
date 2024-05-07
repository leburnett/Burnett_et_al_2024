% make animal_xy_table

all_animals = unique(all_xy_analysis.Animal);
n_animals = numel(all_animals);

xy_analysis_animals = zeros(n_animals, 6);

for j = 1:n_animals
    ani = all_animals{j};
    all_ani = find(string(all_xy_analysis.Animal) == ani & cell2mat(all_xy_analysis.MaxSpEscape)>20 & all_xy_analysis.Trial==1) ; %cell2mat(all_xy_analysis.ReturnToShelter)==1 
    
    if numel(all_ani)>0
        
        g = all_xy_analysis.Geno{all_ani(1)};
        
        if g == "wt"
            geno = 1;
        else 
            geno = 2; 
        end 
        
        % Assign animal number to  first column of table.
        xy_analysis_animals(j,1) = str2double(ani(3:end));
        xy_analysis_animals(j,2) = geno;
        xy_analysis_animals(j,3) = numel(all_ani);
        cv = 4;
        
        for k = [9, 12, 13, 17]
            all_vals = cell2mat(table2array(all_xy_analysis(all_ani, k)));
            xy_analysis_animals(j,cv) = nanmean(all_vals);
            cv = cv+1;
        end
    else
        xy_analysis_animals(j,1) = ani;
        xy_analysis_animals(j,2) = geno;
        xy_analysis_animals(j,3) = numel(all_ani);
        cv = 4;
    end
    
end

Ani = xy_analysis_animals(:,1);
Geno = xy_analysis_animals(:,2);
Return2Shelter = xy_analysis_animals(:,4);
NumTrials = xy_analysis_animals(:,3);
MaxSp = xy_analysis_animals(:,5);
T2M = xy_analysis_animals(:,6);
T2A = xy_analysis_animals(:,7);

animal_xy_table = table(Ani, Geno, Return2Shelter, NumTrials, MaxSp, T2M, T2A);
