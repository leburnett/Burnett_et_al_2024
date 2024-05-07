%% Script for making full arrays of good temporal and spatial receptive fields. 
% Created by Burnett 14/12/20

% Make GTRF_full

% Use good_RF
% GRF_r
% GRF_info
% tempRFs

% Extract row from TempRFs which are in 'good_RF' add the columns from
% GRF_info. 

% Open all files '*GRF.mat'
files = dir('*GRF.mat'); 
num_files = numel(files);


%% TEMPORAL ARRAY 

% Make array with ALL The TEMPORAL ARRAYS: 

all_GTRF = [];

for i = 1:num_files
    filename = files(i).name;
    load(filename)
    
    %Make GTRF. 
    GTRF = tempRFs(good_RF,:);
    GTRF(:,38) = GRF_info(:,2); % Amplitude
    GTRF(:,39) = GRF_info(:,3); % Channel
    GTRF(:,40) = GRF_info(:,5); % Num Spikes
    GTRF(:,41) = GRF_info(:,10); % Trial
    
    new_filename = filename(1:end-7);
    filename_to_save = strcat(new_filename,'GTRF.mat');
    save(filename_to_save, 'GTRF')
    
    all_GTRF = vertcat(all_GTRF, GTRF);
end

save('201214_ALL_GTRF.mat', 'all_GTRF', 'gtrf_data', 'gtrf_table')

% Add WT/HET
for j = 1:226
    if all_GTRF(j,35)==7270 || all_GTRF(j,35)==7788 || all_GTRF(j,35)==7475 ||all_GTRF(j,35)==7616
        all_GTRF(j,42) = 1; %WT
    else 
        all_GTRF(j,42) = 0; %HET
    end 
end 
     

%% Looking at the data: 
   
allWT = find(all_GTRF(:,42)==1);
allHET = find(all_GTRF(:,42)==0);

RFwt = all_GTRF(allWT,:);
RFhet = all_GTRF(allHET,:);

figure
imagesc(RFwt(:,1:30))
colorbar

figure
imagesc(RFhet(:,1:30))
colorbar

a1 = sortrows(RFwt, 40);
a2 = sortrows(RFhet, 40);

figure
imagesc(a1(:,1:30))
colorbar

figure
imagesc(a2(:,1:30))
colorbar


for j = 1:226
   av_end = mean(all_GTRF(j,25:28)); 
   colmax = find(all_GTRF(j, 1:30) == max(all_GTRF(j, 1:30)));
   colmin = find(all_GTRF(j, 1:30) == min(all_GTRF(j, 1:30)));
   
   all_GTRF(j, 43) = av_end;
   
   if av_end <= 0
       all_GTRF(j, 44) = 0;
   elseif av_end > 0
       all_GTRF(j, 44) = 1;
   end
   
   if numel(colmax)>1
       colmax = colmax(1);
   end
   all_GTRF(j, 45) = colmax;
  
    if numel(colmin)>1
       colmin = colmin(1);
   end
   all_GTRF(j, 46) = colmin;
 
end 

%%
  gtrf_data = all_GTRF(:, 1:30);
  gtrf_table = all_GTRF(:, 31:46); 
  gtrf_table = array2table(gtrf_table, 'VariableNames', {'Diff', 'cx', 'cy', 'Date', 'Animal','Depth', 'ID', 'Amp','Channel', 'Spikes', 'Trial', 'Geno', 'AvEnd', 'OnOff', 'MaxCol', 'MinCol'});

%%

for i = 1:226
    if gtrf_table{i,13}<=0 && gtrf_table{i,12}==1
        gtrf_table.Group(i) = 1; 
    elseif gtrf_table{i,13}>0 && gtrf_table{i,12}==1
        gtrf_table.Group(i) = 2; 
    elseif gtrf_table{i,13}<=0 && gtrf_table{i,12}==0
        gtrf_table.Group(i) = 3; 
    elseif gtrf_table{i,13}>0 && gtrf_table{i,12}==0
        gtrf_table.Group(i) = 4;
    end 
end 
  
all1 = find(gtrf_table.Group ==1);
all2 = find(gtrf_table.Group ==2);
all3 = find(gtrf_table.Group ==3);
all4 = find(gtrf_table.Group ==4);

allwt = find(gtrf_table.Geno ==1);
allhet = find(gtrf_table.Geno ==0);

var = gtrf_table.Depth; 
group = gtrf_table.Group; 

boxplot(var, group)
xticks([1,2,3,4])
xticklabels({'WT-ON', 'WT-OFF', 'HET-ON', 'HET-OFF'})
xtickangle(45)
ylabel('Depth from surface (um)')

[p,t,stats] = anova1(var, group)
[c,m,h] = multcompare(stats)
        
vals1= gtrf_table{all1,16}; 
vals2 = gtrf_table{all3,16}; 
[h,p] = kstest2(vals1, vals2)










%% SPATIAL ARRAY 


% Open all files '*GRF.mat'
files = dir('*GRF.mat'); 
num_files = numel(files);

% Make array with ALL The SPATIAL ARRAYS: 

all_GSRF = [];
gsrf_info = [];

for i = 1:num_files
    filename = files(i).name;
    load(filename)
    
    %Remove extra depth column and column for 'grf'
    GRF_info(:,4) = [];
    GRF_info(:,7) = [];
    
    % Normalise each cell.
    num = numel(GRF_r(:,1));
    GRF_r2 = []; 
        
    for j = 1:num
        GRF_r2(j, :) = zscore(GRF_r(j,:));
        
        average_val = mean(GRF_r2(j,:));
        [max_val, maxi] = max(GRF_r2(j,:));
        [min_val, mini] = min(GRF_r2(j,:));
        
        GRF_info(j, 9) = average_val;
        GRF_info(j, 10) = max_val;
        GRF_info(j, 11) = maxi;
        GRF_info(j, 12) = min_val;
        GRF_info(j, 13) = mini;
        
        average60to100 = mean(GRF_r2(j, 70:90));
%         averageall = mean(GRF_r2(j, :));
        if average60to100 >= 0
            onoff = 1;
            GRF_info(j, 14) = 1;
        elseif average60to100 < 0
            onoff = 0;
            GRF_info(j, 14) = 0;
        end
        
        GRF_info(j, 15) = average60to100; 
        
        if onoff == 1
            maxon = max(GRF_r2(j,70:90));
            data = GRF_r2(j,40:120);
            data = data-(maxon/2);
            vals = sign(data);
            vals2 = diff(vals);
            vall = find(vals2 ~=0);
            if numel(vall) == 2
                difval = diff(vall);
            elseif numel(vall)~=2
                difval = NaN;
            end
            
        elseif onoff == 0
            minoff = min(GRF_r2(j,70:90));
            data = GRF_r2(j,40:120);
            data = data+(abs(minoff)/2);
            vals = sign(data);
            vals2 = diff(vals);
            vall = find(vals2 ~=0);
            if numel(vall) == 2
                difval = diff(vall);
            elseif numel(vall)~=2
                difval = NaN;
            end
        end
        
        GRF_info(j, 16) = difval;
        
    end
    
    all_GSRF = vertcat(all_GSRF, GRF_r2);
    gsrf_info = vertcat(gsrf_info, GRF_info);
    
end



% Add WT/HET
for j = 1:226
    if gsrf_info(j,7)==7270 || gsrf_info(j,7)==7788 || gsrf_info(j,7)==7475 ||gsrf_info(j,7)==7616
        gsrf_info(j,17) = 1; %WT
    else 
        gsrf_info(j,17) = 0; %HET
    end 
end 

save('201214_ALL_GSRF.mat', 'all_GSRF', 'gsrf_data', 'gsrf_table')

allWT = find(gsrf_info(:,17)==1);
allHET = find(gsrf_info(:,17)==0);

allWT_on = find(gsrf_info(:,17)==1 & gsrf_info(:,14)==1);
allHET_on = find(gsrf_info(:,17)==0 & gsrf_info(:,14)==1);
allWT_off = find(gsrf_info(:,17)==1 & gsrf_info(:,14)==0);
allHET_off = find(gsrf_info(:,17)==0 & gsrf_info(:,14)==0);

all_GSRF_WT = all_GSRF(allWT,:);
all_GSRF_HET = all_GSRF(allHET,:);

gsrf_info_WT = gsrf_info(allWT,:);
gsrf_info_HET = gsrf_info(allHET,:);

% ADD column of depth to values. 
all_GSRF_WT(:, 162) = gsrf_info_WT(:,5); 
all_GSRF_HET(:, 162) = gsrf_info_HET(:,5); 

all_GSRF_WT = sortrows(all_GSRF_WT, 162);
all_GSRF_HET = sortrows(all_GSRF_HET, 162);

figure
subplot(1,2,1)
imagesc(all_GSRF_WT(:, 1:161))
subplot(1,2,2)
imagesc(all_GSRF_HET(:, 1:161))


allON = find(gsrf_info(:,14)==1);
allOFF = find(gsrf_info(:,14)==0);

figure
subplot(1,2,1)
imagesc(all_GSRF(allON,:))
subplot(1,2,2)
imagesc(all_GSRF(allOFF,:))
colorbar

% PLOT 
figure
subplot(2,2,1)
imagesc(all_GSRF(allWT_on,:))
colorbar
title('WT - ON')
colormap(redblue)

subplot(2,2,2)
imagesc(all_GSRF(allWT_off,:))
colorbar
title('WT - OFF')
colormap(redblue)

subplot(2,2,3)
imagesc(all_GSRF(allHET_on,:))
colorbar
title('HET - ON')
colormap(redblue)

subplot(2,2,4)
imagesc(all_GSRF(allHET_off,:))
colorbar
title('HET - OFF')
colormap(redblue)


% GRF info 
% Col 1 = ID
% COl 2 = Amp
% Col 3 = Channel
% Col 4 = spikes 
% Col 5 = depth 
% Col 6 = Date
% Col 7 = Animal
% Col 8 = Trial 


for i = 1:12
    GRF_r2(i,:) = zscore(GRF_r(i,:));
end 


imagesc(GRF_r2)
figure 
imagesc(GRF_r)

colorbar


[p,h] =ranksum(gsrf_info(allWT_off,15), gsrf_info(allHET_off, 15))


save('201214_ALL_GSRF.mat', 'all_GSRF', 'gsrf_info');


for k = 1:226
  if gsrf_info(k,17)==1 && gsrf_info(k,14)==1
      gsrf_info(k,18) = 1; 
  elseif gsrf_info(k,17)==1 && gsrf_info(k,14)==0
      gsrf_info(k,18) = 2; 
  elseif gsrf_info(k,17)==0 && gsrf_info(k,14)==1
      gsrf_info(k,18) = 3; 
  elseif gsrf_info(k,17)==0 && gsrf_info(k,14)==0
      gsrf_info(k,18) = 4; 
  end 
end 


var = gsrf_info(:,16);
group = gsrf_info(:,18);

boxplot(var, group)

[p,h] =ranksum(gsrf_info(allWT_on,16), gsrf_info(allHET_on, 16))
[p,h] =ranksum(gsrf_info(allWT_off,16), gsrf_info(allHET_off, 16))











