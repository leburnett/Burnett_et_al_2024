%% Find info for alllll cells


% Open all files '*GRF.mat'
files = dir('*info.mat'); 
num_files = numel(files);

all_info = [];

for i = 1:num_files
    filename = files(i).name;
    load(filename)
    n_cells = numel(ids_info(:,1)); 
    depth_probe = str2num(filename(end-16:end-13)); 
    
    %Flipping depth back to -800. 
    for jj = 1:n_cells
        d1 = ids_info(jj, 4);
        diff1 = 800 - d1;
        d2 = depth_probe - diff1;
        ids_info(jj,11) = d2;
    end
    
    all_info = vertcat(all_info, ids_info);
end

save('210125_ALL_RF_INFO_FLIP.mat', 'all_info') 



% Add WT/HET
for j = 1:2066
    if all_info(j,9)==7270 || all_info(j,9)==7788 || all_info(j,9)==7475 ||all_info(j,9)==7616
        all_info(j,12) = 1; %WT
    else 
        all_info(j,12) = 0; %HET
    end 
end 
     
allWT = find(all_info(:,12)==1); 
allHET = find(all_info(:,12)==0); 


figure
histogram(all_info(allWT,11))
hold on 
histogram(all_info(allHET,11))
xlabel('Depth')
ylabel('#VR Cells')

histogram(all_info(allWT,11), 'Normalization', 'pdf')
hold on 
histogram(all_info(allHET,11), 'Normalization', 'pdf')
xlabel('Depth')
ylabel('PDF')



a1 = mean(all_info(allWT,11))
a2 = mean(all_info(allHET,11))

[p,h] = ranksum(a1,a2)



