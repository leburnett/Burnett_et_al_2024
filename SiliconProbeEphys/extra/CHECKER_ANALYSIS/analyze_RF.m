%% Analyse RF properties. 

% Script generated by Burnett 29/10/20, adapted from Olga's script
% 'analyze_checkers.m'. 

% The purpose of the script is to take the 'RF' 4D array - of all the
% values of the pixels of the screen over time for each cluster and analyse
% the temporal and spatial Rf properties. 

% All files in directory. 
files = dir('*.mat'); 
num_files = numel(files);

% First, find properties of checker. 
ncol = 30; %This was always 30. 
texRect_size=[608, 342];
sq_w=single(floor(texRect_size(1)/ncol)); %size of the square
nrow=max(1,floor(texRect_size(2)/sq_w)); % number of checker rows
% npatt=ceil(nfr/nframes)+10;
% frame_array=zeros(nrow*sq_w, ncol*sq_w,npatt,'uint8');
imheight = nrow*sq_w;
imwidth = ncol*sq_w; 
checker_size=sq_w;

%make crop size be twice the size of the checker
crop_size_half=checker_size*4; %
full_crop_size=crop_size_half*2+1;

for i = 1:num_files
    fname = files(i).name;
    load(fname);
    RFsum = RFsummary; 
    n_gcl = numel(RFsum(:,1,1));
%     RFsum = RFsum2;
        
%% RF - spatial analysis

%% PLOT ALL RF - with ID
% % figure
% % for j = 1:n_gcl
% %     A = squeeze(RFsum(j,:,:)); % 321 * 321 square - 103,041 pixels. 
% %     subplot(10,7,j)
% %     imagesc(A)
% %     colormap(redblue)
% %     title(j)
% % end 

%% Sort RFs manually into 'good' and 'bad'
bad_RF = []; 
good_RF = [];
prompt = 'Press y for GCL, n for BCL: ';
bkg_resp = squeeze(nanmean(squeeze(RFsum(:,:,:)))); %Mean of allll activity over all cells. 

for j = 1:6 %n_gcl
    A = squeeze(RFsum(j,:,:)); % 321 * 321 square - 103,041 pixels. 
    A2= A-bkg_resp; 
    figure
    imagesc(A2)
    colorbar
    colormap(redblue)
%     m2 = mean(nanmean(A));
%     caxis([m2-(m2/2), m2+(m2/2)]) 
    x = input(prompt, 's'); 
    if x == 'n'
        bad_RF = [bad_RF, j];
    elseif x == 'y'
        good_RF = [good_RF, j];
    end
    close
end 
   
% Subtracting the mean? 
% meanA = mean(A);
% A2 = A-meanA;
% imagesc(23)
% colorbar
% colormap(redblue)
%% Plot Good and Bad RFs to check. 

% % Plot Good RF to check. 
%  n_grf = numel(good_RF); 
% figure
% for j2 = 1:n_grf
%     cl_id = good_RF(j2);
%     subplot(7,3,j2)
%     A = squeeze(RFsum(cl_id,:,:));
%     imagesc(A)
%     colormap(redblue)
%     title(string(cl_id))
%     xticks([])
%     xticklabels({})
%     yticks([])
%     yticklabels({})
% end 
% sgtitle('Good RF')
% pause 
% close

% 
% % Plot Bad RF to check. 
% n_brf = numel(bad_RF); 
% figure
% for j3 = 1:n_brf
%     cl_id = bad_RF(j3);
%     subplot(7,6,j3)
%     A = squeeze(RFsum(cl_id,:,:));
%     imagesc(A)
%     colormap(redblue)
%     title(string(cl_id))
%     xticks([])
%     xticklabels({})
%     yticks([])
%     yticklabels({})
% end 
% sgtitle('Bad RF')


% If there are any that should be in different groups - manually add/remove them
% from 'good_RF'. 

%% Then make new 'RFsum' with ONLY the good clusters (GRF).  - SPATIAL RFs
n_grf = numel(good_RF); 

GRF_full = zeros(n_grf, 161, 161); 
GRF_r = zeros(n_grf, 161); 

for j2 = 1:n_grf
    cl_id = good_RF(j2);
    A = squeeze(RFsum(cl_id,:,:));
    GRF_full(j2, :, :)= A; %full RF. 
    
%     figure
%     subplot(n_grf,2,j2*2-1)
%     imagesc(A)
%     colormap(redblue)
%     title(string(cl_id))
%     xticks([])
%     xticklabels({})
%     yticks([])
%     yticklabels({})
%     subplot(n_grf,2,j2*2)
    Am = mean(A,1);
    Am2 = mean(A,2)';
    Am3 = [Am;Am2];
    Am3 = mean(Am3); %radial average. 
%     plot(smooth(Am3,6), 'k')
    GRF_r(j2, :)= Am3;
end 


%PLOT RF LINES 
% figure
% for j2 = 1:n_grf
%     cl_id = good_RF(j2);
%     subplot(7,3,j2)
%     A = squeeze(RFsum(cl_id,:,:));
%     Am = mean(A,1);
%     Am2 = mean(A,2)';
%     Am3 = [Am;Am2];
%     Am3 = mean(Am3);
%     hold on 
%     plot(smooth(Am3,6), 'k')
%     xticks([])
%     xticklabels({})
% end 
% sgtitle('Good RF')



exp_name = fname(1:end-13);
save_folder = '/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/RF/GRF/GRF_LONG';
save(fullfile(strcat(save_folder,'\', exp_name, '_GRF.mat')), 'GRF_full', 'GRF_r', 'good_RF', 'bad_RF');
% save(strcat(exp_name, '_GRF.mat'), 'GRF_full', 'GRF_r');

end 
%% 





%% GRF_r

n = numel(GRF_r(:,1));
figure
for i = 1:n
%     figure
plot(smooth(GRF_r(i,:)), 'k')
hold on 
plot([0 180], [127.5 127.5], 'k:')
end
title('GN7614-2000')









%%  IF STARTING FROM 'RF' 

% % % % % % If starting from ALL RF - ALL TIME: 
% % % % % %array to visualize the summary of RF (cropped and only 1 frame)
% 
% ncol = 30; %This was always 30. 
% texRect_size=[608, 342];
% sq_w=single(floor(texRect_size(1)/ncol)); %size of the square
% nrow=max(1,floor(texRect_size(2)/sq_w)); % number of checker rows
% % npatt=ceil(nfr/nframes)+10;
% % frame_array=zeros(nrow*sq_w, ncol*sq_w,npatt,'uint8');
% imheight = nrow*sq_w;
% imwidth = ncol*sq_w; 
% checker_size=sq_w;
% 
% %make crop size be twice the size of the checker
% crop_size_half=checker_size*4; %
% full_crop_size=crop_size_half*2+1;
% 
% n_gcl = numel(goodclusters);
% RFsum=zeros(n_gcl,full_crop_size, full_crop_size);
% 
% for ci=1:n_gcl
%     RFi=squeeze(RF(ci,:,:,:)); % selecting all pixels over time for an individual cluster. 
%     vari=var(RFi,0,3); % Find the variance over time. 
%     [mval,cx]=max(max(vari,[],1));
%     [mval,cy]=max(max(vari,[],2));    
%     %disp([ci,mval]);
%     meani=mean(RFi(:));
%     %find the largest deviation at the RF in time
%     [mval,ct]=max(abs(squeeze(RFi(cy,cx,:))-meani)); % - here do NOT squeeze to keep time component. 
%     %Crop the RF
%     RFi_tocrop=zeros(crop_size_half*2+imheight,crop_size_half*2+imwidth);
%     RFi_tocrop(crop_size_half+1:crop_size_half+imheight,crop_size_half+1:crop_size_half+imwidth)=RFi(:,:,ct);
%     RFsum(ci,:,:)=RFi_tocrop(cy:cy+2*crop_size_half,cx:cx+2*crop_size_half);
% end
% 
% save(fullfile(strcat('C:\Data_analysis\Silicon_Probe_Analysis\RF_files\201006_GN7614_Test4_2000.mat')), 'RFsum')
% 
% clear all

% % % % % 
% % % % % %sort clusters by depth
% % % % % ids_depth_sorted=sortrows(ids_depth,2);


%%  Stage 2: Once you have all the GRF. 
% Reading through all the GRFs - to make COMBINED WT/HET arrays. 
% Move into directory with the 'GRP.mat' files for HET/WT. These have been
% separated manually. 

GRF_WT = [];

files = dir('*.mat'); 
num_files = numel(files);

for i = 1:num_files
    fname = files(i).name;
    load(fname);
    n_gcl = numel(GRF_r(:,1));
    GRF_WT = vertcat(GRF_WT, GRF_r); 
end 

%% GRF_full - WT

GRF_WT_ALL = zeros(10,161,161);

files = dir('*.mat'); 
num_files = numel(files);
v = 1; 

for i = 1:num_files
    fname = files(i).name;
    load(fname);
    n_gcl = numel(GRF_r(:,1));
    GRF_WT_ALL(v:v+n_gcl-1,:,:) = GRF_full; 
    v = v+n_gcl;  
end 

%% SAVE
save('201030_GRF_WT.mat', 'GRF_WT', 'GRF_WT_ALL')

%% Do the same for HET. 

GRF_HET = [];

files = dir('*.mat'); 
num_files = numel(files);

for i = 1:num_files
    fname = files(i).name;
    load(fname);
    n_gcl = numel(GRF_r(:,1));
    GRF_HET = vertcat(GRF_HET, GRF_r); 
end 

%% GRF_full

GRF_HET_ALL = zeros(10,161,161);

files = dir('*.mat'); 
num_files = numel(files);
v = 1; 

for i = 1:num_files
    fname = files(i).name;
    load(fname);
    n_gcl = numel(GRF_r(:,1));
    GRF_HET_ALL(v:v+n_gcl-1,:,:) = GRF_full; 
    v = v+n_gcl;  
end 


%% SAVE
 
save('201030_GRF_HET.mat', 'GRF_HET', 'GRF_HET_ALL')

%% If you open one 'GRF.mat' file - can plot all the 'selected RFs' for this recording. 

% To plot ALL of the RFs.  
 n_grf = numel(GRF_r(:,1)); 
figure
for j2 = 1:n_grf
    subplot(2,5,j2)
    A = squeeze(GRF_full(j2,:,:));
    imagesc(A)
    colormap(redblue)
%     title(string(j2))
    xticks([])
    xticklabels({})
    yticks([])
    yticklabels({})
end 
sgtitle('7476 - HET - 2000')


%% Once you have the COMBINED WT/HET arrays, the following code will plot them. 

%% Plots of traces - WT
% All traces - shifted down by 'offset' so that they are not overlapping. 

offset = 0; 
figure
for i = 1:225
plot(GRF_WT(i,:)-offset, 'k')
hold on 
offset = offset+10; 
end 

[idx, c] = kmeans(GRF_WT, 10); 
GRF_WT2 = [GRF_WT,idx]; 
GRF_WT2 = sortrows(GRF_WT2, 162);

offset = 0; 
figure
for i = 1:225
plot(smooth(GRF_WT2(i,1:160),6)-offset, 'k')
hold on 
offset = offset+15; 
end 

%% RF + Traces subplot - WT
% THis is a two-part subplot, with the RF heatmap on the LHS and the trace
% on the RHS. 
% Choose the values in 'vals' to choose which RFs to show. 

% vals = [1:1:10]; 
vals = [2,3,4,5,11,12,21,22,23,24];

figure
for i = 1:10
    j = vals(i);
    A = GRF_WT(j,1:160);
    AB = smooth(A, 6);
    A2 = squeeze(GRF_WT_ALL(j,:,:)); 
  
    subplot(10,3,[(i*3)-2])
    imagesc(A2)
    colormap(redblue)
    axis off
    
    subplot(10,3,[(i*3)-1, (i*3)])
    plot(AB, 'k')
    xticks([])
    xticklabels({})
end 


%% RF + Traces subplot - HET 

% vals = [70:1:79]; 
vals = [25,22,3,4,5,6,7,8,9,10];

figure
for i = 1:10
    j = vals(i);
    A = GRF_HET(j,1:160);
    AB = smooth(A, 6);
    A2 = squeeze(GRF_HET_ALL(j,:,:)); 
  
    subplot(10,3,[(i*3)-2])
    imagesc(A2)
    colormap(redblue)
    axis off
    
    subplot(10,3,[(i*3)-1, (i*3)])
    plot(AB, 'k')
    xticks([])
    xticklabels({})
end 





%%  TRYING TO QUANTIFY RFs.  %%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% FIND MAX / MIN - current example = HET
num_HET = numel(GRF_HET(:,1));



for i = 1:num_HET
    data = GRF_HET(i,1:160); 
    maxD = max(data);
    minD = min(data); 
    meanD = mean(data); 
    ONd = maxD - meanD; 
    OFFd = meanD-minD;
    GRF_HET(i, 162) = maxD; 
    GRF_HET(i, 163) = minD; 
    GRF_HET(i, 164) = ONd; 
    GRF_HET(i, 165) = OFFd;
    GRF_HET(i, 166) = ONd -OFFd; 
    
    %on or off cell
    if GRF_HET(i, 166) >0 
        GRF_HET(i, 167) = 1; 
    elseif GRF_HET(i, 166) <=0
        GRF_HET(i, 167) = 0;
    end 
    
%     if GRF_HET(i, 167) ==1
%         yval_halfwidth = ((maxD-meanD)/2)+meanD;
%         Q = find(data == round(yval_halfwidth)); 
%         
%     end  
        
       % FIND WHERE CURVE INTERSECTS WITH STRAIGHT LINE.      
end 

offc = find(GRF_HET(:,167)==0);
onc = find(GRF_HET(:,167)==1);
   



%% Finding width at half peak. 

num_HET = numel(GRF_HET(:,1)); 
voi_HET = zeros(num_HET, 8);

for j = 1:num_HET
    
    data = GRF_HET(j,1:160); 
    data2 = squeeze(GRF_HET_ALL(j,1:160, 1:160));
    
    figure
    imagesc(data2)
    colorbar
    colormap(redblue)
    m2 = mean(data);
    caxis([m2-(m2/2), m2+(m2/2)]) % - normalise so that white is the mean of all the pixels.
    
    prompt = 'Press y if good, else b for skip: ';
    x = input(prompt, 's'); 
    
    if x == 'n' % If the RF does not look adequate. 
        close
        voi_HET(j,1) = 0; 
        voi_HET(j,2:8) = NaN; 
        continue
        
    elseif x == 'y' % Else continue 
        voi_HET(j,1) = 1; 
        close
        
        %Find the min and max values of data. 
        minD = min(data);
        voi_HET(j,2) = minD;
        
        maxD = max(data);
        voi_HET(j,3) = maxD;
        
        voi_HET(j,4) = m2;
        
        gain_min = abs(minD - m2);
        voi_HET(j,5) = gain_min;
        
        gain_max = abs(maxD - m2);
        voi_HET(j,6) = gain_max;
        
        gain_s = gain_max-gain_min; % if negative then off, if positive then on.
        voi_HET(j,7) = gain_s;
        
        if gain_s <=0 
             halfmax = m2-(gain_min/2);
        elseif gain_s >0 
            halfmax = m2+(gain_max/2);
        end 
            
        data3 = data-halfmax;
        sdata = sign(data3);
%         plot(data3)
%         plot([0 160], [0 0], 'k')
        sd_data = diff(sdata);
        xvals_halfpeak  = find(sd_data ~=0);
        numx = numel(xvals_halfpeak);
        if numx ==2
            whp = xvals_halfpeak(2)-xvals_halfpeak(1);
            voi_HET(j,8) = whp;
        elseif numx ~= 2
            voi_HET(j,8) = 0;
            disp('Warning: more than 2 points of intersection. ')
        end
    end
 
end 

save('VOI_RF_HET.mat', 'voi_HET')

% figure
% plot(data)
% hold on 
% plot([0 160], [m2 m2], 'k')



%% WT - VOI

num_WT = numel(GRF_WT(:,1)); 
voi_WT = zeros(num_WT, 8);

for j = 1:num_WT
    
    data = GRF_WT(j,1:160); 
    data2 = squeeze(GRF_WT_ALL(j,1:160, 1:160));
    
    figure
    imagesc(data2)
    colorbar
    colormap(redblue)
    m2 = mean(data);
    caxis([m2-(m2/2), m2+(m2/2)]) % - normalise so that white is the mean of all the pixels.
    
    prompt = 'Press y if good, else n for skip: ';
    x = input(prompt, 's'); 
    
    if x == 'n' % If the RF does not look adequate. 
        close
        voi_WT(j,1) = 0; 
        voi_WT(j,2:8) = NaN; 
        continue
        
    elseif x == 'y' % Else continue 
        close
        voi_WT(j,1) = 1; 
        
        %Find the min and max values of data. 
        minD = min(data);
        voi_WT(j,2) = minD;
        
        maxD = max(data);
        voi_WT(j,3) = maxD;
        
        voi_WT(j,4) = m2;
        
        gain_min = abs(minD - m2);
        voi_WT(j,5) = gain_min;
        
        gain_max = abs(maxD - m2);
        voi_WT(j,6) = gain_max;
        
        gain_s = gain_max-gain_min; % if negative then off, if positive then on.
        voi_WT(j,7) = gain_s;
        
        if gain_s <=0 
             halfmax = m2-(gain_min/2);
        elseif gain_s >0 
            halfmax = m2+(gain_max/2);
        end 
            
        data3 = data-halfmax;
        sdata = sign(data3);
%         plot(data3)
%         plot([0 160], [0 0], 'k')
        sd_data = diff(sdata);
        xvals_halfpeak  = find(sd_data ~=0);
        numx = numel(xvals_halfpeak);
        if numx ==2
            whp = xvals_halfpeak(2)-xvals_halfpeak(1);
            voi_WT(j,8) = whp;
        elseif numx ~= 2
            voi_WT(j,8) = 0;
            disp('Warning: more than 2 points of intersection. ')
        end
    end
 
end 

save('VOI_RF_WT.mat', 'voi_WT')


%% Quantifying values of interest. 

numRF_WT = numel(find(voi_WT(:,1)==1)); 
numRF_HET = numel(find(voi_HET(:,1)==1)); 

rowsnzWT = find(voi_WT(:,8)~=0 & ~isnan(voi_WT(:,8)));
meanwhpWT = nanmean(voi_WT(rowsnzWT, 8));

rowsnzHET = find(voi_HET(:,8)~=0 & ~isnan(voi_HET(:,8)));
meanwhpHET = nanmean(voi_HET(rowsnzHET, 8));


whp_WT = []; 
for i = 1:72
    val  = rowsnzWT(i); 
    whp1 = voi_WT(val,8); 
    whp_WT = [whp_WT, whp1]; 
end 

whp_HET = []; 
for i = 1:70
    val  = rowsnzHET(i); 
    whp1 = voi_HET(val,8); 
    whp_HET = [whp_HET, whp1]; 
end 

figure
plot(whp_WT, 'k.', 'MarkerSize', 20)
axis([0 75 0 50])
hold on 
plot(whp_HET, 'r.', 'MarkerSize', 20)

[p,h] = ranksum(whp_WT, whp_HET)
% [h,p] = kstest2(whp_WT, whp_HET)

figure
boxplot([whp_WT])

whp_WT(2,:) = 0; 
whp_HET(2,:) = 1; 
whp = horzcat(whp_WT, whp_HET); 

figure
bp = boxplot(whp(1,:), whp(2,:), 'ColorGroup', whp(2,:), 'Color', 'k');
title('WHP - all RFs')
xticks([1,2])
xticklabels({'WT', 'HET'})
hold on 
x = ones(1,72);
% plot(x,whp_WT(1,:), 'k.', 'MarkerSize', 18)
scatter(x,whp_WT(1,:),150, 'k.' ,'jitter', 'on', 'jitteramount', 0.05)
x2 = ones(1,70)*2;
% plot(x2,whp_HET(1,:), 'r.', 'MarkerSize', 18)
scatter(x2,whp_HET(1,:),150, 'r.' ,'jitter', 'on', 'jitteramount', 0.05)

savefig(gcf, 'WHP_WT_HET_ALL.fig')
close






%%  Making FULL plots of all - 'good' RFs. 

rowszero_WT = find(voi_WT(:,8)==0);
rowszero_HET =  find(voi_HET(:,8)==0);

figure
for j = 1:72
    val  = rowsnzWT(j); 
    data2 = squeeze(GRF_WT_ALL(val,1:160, 1:160));
    subplot(18,4,j)
    imagesc(data2)
    xticks([])
    yticks([])
    xticklabels({})
    yticklabels({})
%     colorbar
    colormap(redblue)
    m2 = mean(mean(data2));
    caxis([m2-(m2/2), m2+(m2/2)]) % - normalise so that white is the mean of all the pixels.
%     pause 
%     close
end 
sgtitle('RF WT')


figure
for j = 1:70
    val  = rowsnzHET(j); 
    data2 = squeeze(GRF_HET_ALL(val,1:160, 1:160));
    subplot(18,4,j)
%     figure
    imagesc(data2)
%     colorbar
    xticks([])
    yticks([])
    xticklabels({})
    yticklabels({})
    colormap(redblue)
    m2 = mean(mean(data2));
    caxis([m2-(m2/2), m2+(m2/2)]) % - normalise so that white is the mean of all the pixels.
%     pause 
%     close
end 
sgtitle('RF HET')

savefig(gcf, '201030_RF_WT.fig')
close











%% Temporal analysis. 
% Using the 'good RFs' already found. 

% see: 'Make_TemporalRF_Array.m'







% 
% files = dir('*.mat'); 
% num_files = numel(files);
% 
% % First, find properties of checker. 
% ncol = 30; %This was always 30. 
% texRect_size=[608, 342];
% sq_w=single(floor(texRect_size(1)/ncol)); %size of the square
% nrow=max(1,floor(texRect_size(2)/sq_w)); % number of checker rows
% % npatt=ceil(nfr/nframes)+10;
% % frame_array=zeros(nrow*sq_w, ncol*sq_w,npatt,'uint8');
% imheight = nrow*sq_w;
% imwidth = ncol*sq_w; 
% checker_size=sq_w;
% 
% %make crop size be twice the size of the checker
% crop_size_half=checker_size*4; %
% full_crop_size=crop_size_half*2+1;
% 
% for i = 1:num_files
%     fname = files(i).name;
%     load(fname);
%     n_gcl = numel(RFsum(:,1,1));
% 
% 
% 


% % % % % 
% % % % % 
% % % % % nframes = 30; 
% % % % % %array of croppped RF over 30 frames. 
% % % % % RFsum=zeros(n_gcl,full_crop_size, full_crop_size, nframes);
% % % % % 
% % % % % for ci=1:n_gcl
% % % % %     RFi=squeeze(RF(ci,:,:,:)); % selecting all pixels over time for an individual cluster. 
% % % % %     vari=var(RFi,0,3); % Find the variance over time. 
% % % % %     [mval,cx]=max(max(vari,[],1));
% % % % %     [mval,cy]=max(max(vari,[],2));    
% % % % %     %disp([ci,mval]);
% % % % %     meani=mean(RFi(:));
% % % % %     %find the largest deviation at the RF in time
% % % % %     [mval,ct]=max(abs(squeeze(RFi(cy,cx,:))-meani)); % - here do NOT squeeze to keep time component. 
% % % % %     %Crop the RF
% % % % %     RFi_tocrop=zeros(crop_size_half*2+imheight,crop_size_half*2+imwidth);
% % % % %     RFi_tocrop(crop_size_half+1:crop_size_half+imheight,crop_size_half+1:crop_size_half+imwidth)=RFi(:,:,ct);
% % % % %     RFsum(ci,:,:)=RFi_tocrop(cy:cy+2*crop_size_half,cx:cx+2*crop_size_half);
% % % % % end
% % % % % 
% % % % % GRF_full = zeros(n_grf, 161, 161); 
% % % % % GRF_r = zeros(n_grf, 161); 
% % % % % 
% % % % % n_grf = numel(good_RF); 
% % % % % for j2 = 1:n_grf
% % % % %     cl_id = good_RF(j2);
% % % % %     A = squeeze(RFsum(cl_id,:,:));
% % % % %     GRF_full(j2, :, :)= A; %full RF. 
% % % % %     
% % % % %     figure
% % % % %     subplot(n_grf,2,j2*2-1)
% % % % %     imagesc(A)
% % % % %     colormap(redblue)
% % % % %     title(string(cl_id))
% % % % %     xticks([])
% % % % %     xticklabels({})
% % % % %     yticks([])
% % % % %     yticklabels({})
% % % % %     subplot(n_grf,2,j2*2)
% % % % %     Am = mean(A,1);
% % % % %     Am2 = mean(A,2)';
% % % % %     Am3 = [Am;Am2];
% % % % %     Am3 = mean(Am3); %radial average. 
% % % % %     plot(smooth(Am3,6), 'k')
% % % % %     GRF_r(j2, :)= Am3;
% % % % % end 
% % % % % 





