% Create Spatial RF array. 
% Created by Burnett 19/11/20
% % % % 

% % % % % % If starting from ALL RF - ALL TIME: 
% % % % % %array to visualize the summary of RF (cropped and only 1 frame)

% Properties of the checker. 

ncol = 30; 
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


% 
% for ci=1:n_gcl
%     RFi=squeeze(RF(ci,:,:,:)); % selecting all pixels over time for an individual cluster. 
%     vari=var(RFi,0,3); % Find the variance over time. 
%     [mval,cx]=max(max(vari,[],1));
%     [mval,cy]=max(max(vari,[],2));    
%     %disp([ci,mval]);
%     meani=mean(RFi(:));
%     %find the largest deviation in the RF across time
%     [mval,ct]=max(abs(squeeze(RFi(cy,cx,:))-meani)); % - here do NOT squeeze to keep time component. 
%     %Crop the RF
%     RFi_tocrop=zeros(crop_size_half*2+imheight,crop_size_half*2+imwidth);
%     RFi_tocrop(crop_size_half+1:crop_size_half+imheight,crop_size_half+1:crop_size_half+imwidth)=RFi(:,:,ct);
%     RFsum(ci,:,:)=RFi_tocrop(cy:cy+2*crop_size_half,cx:cx+2*crop_size_half);
% end
% 
% 


%% MEAN ACROSS PIXEL WITH HIGEST VAR + 1 PIXEL EITHER SIDE (both x and y)

% Need to find CENTRE of RF and plot the spatial RFs only in these places.

% Only generate this once!! It will be added to each time. 
SpatRF = struct('RF', {}, 'Animal', {}, 'Depth', {}, 'cx', {}, 'cy', {}, 'ct', {});  

%%% THINGS TO DO EACH TIME %%%
files = dir; 
numfiles = length(files);

for j = 3:numfiles
% 1 - Load 'RF' 
   fname = files(j).name;
   load(fname);

% 2 - add depth of recording and animal # %%%%%%%%%%%%%
depthprobe = str2num(fname(1:4)); 
animal = str2num(fname(13:16));

%sort clusters by depth
n_gcl = numel(goodclusters);
ids_depth(1:n_gcl, 3) = 1:1:n_gcl;
ids_depth_sorted=sortrows(ids_depth,2);

spatialrf = struct('RF', {}, 'Animal', {}, 'Depth', {}, 'cx', {}, 'cy', {}, 'ct', {});  

% % % % % % Spatial plot for the entire FOV - not cropped. 

for i = 1: n_gcl
    
    RFsum=zeros(1,full_crop_size, full_crop_size);
    
    clid = ids_depth_sorted(i,3); %1 does not work.
    if clid ~=0
        depth1 =  ids_depth_sorted(i,2);
        depth1 =  800 - depth1;
        realdepth = depthprobe - depth1;
        
        cluster1 = squeeze(RF(clid,:,:,:)); % all info for this cluster
        
        vari=var(cluster1,0,3);  % Find variance within this cluster.
        [mval,cx]=max(max(vari,[],1));  % Find the x/y values of highest variance.
        [mval,cy]=max(max(vari,[],2));
        
        if cx > 1 && cx <599 && cy >1 && cy < 339 %ensure that max dev is inside FOV.
            %disp([ci,mval]);
            meani=mean(cluster1(:)); %Find mean over entire time for this cell.
            % Find 'ct' the index of the time when the variance is largest.
            [mval,ct]=max(abs(squeeze(cluster1(cy,cx,:))-meani)); % - here do NOT squeeze to keep time component.
            %     %Crop the RF
            RFi_tocrop=zeros(crop_size_half*2+imheight,crop_size_half*2+imwidth);
            RFi_tocrop(crop_size_half+1:crop_size_half+imheight,crop_size_half+1:crop_size_half+imwidth)=cluster1(:,:,ct);
            RFsum=RFi_tocrop(cy:cy+2*crop_size_half,cx:cx+2*crop_size_half);
            %  else
            %      meani=mean(cluster1(:)); %Find mean over entire time for this cell.
            % %     tempRFCROP(ci,:)
            %     tpp = squeeze(cluster1(cy:cy,cx:cx,:))-meani;
            % %     tpp = squeeze(tin(1,1,:));
        end
        
        spatialrf(i).RF = RFsum; 
        spatialrf(i).Animal = animal; 
        spatialrf(i).Depth = realdepth; 
        spatialrf(i).cx = cx; 
        spatialrf(i).cy = cy;
        spatialrf(i).ct = ct; 
    end
end 

 % Vertically concatenate
SpatRF = [SpatRF; spatialrf]; 
end 


save('201112_TemporalRF_CROP_HET.mat', 'TempRF')

% % % 
% % % for ci = 1:20 %n_gcl
% % %     RFi=squeeze(RF(ci,:,:,:)); % RF per cell. 
% % %     vari=var(RFi,0,3);  % Find variance within this cluster. 
% % %     [mval,cx]=max(max(vari,[],1));  % Find the x/y values of highest variance. 
% % %     [mval,cy]=max(max(vari,[],2));   
% % %     
% % %     if cx > 1 && cx <599 && cy >1 && cy < 339
% % %     %disp([ci,mval]);
% % %     meani=mean(RFi(:)); %Find mean over entire time for this cell. 
% % % %     tempRFCROP(ci,:) 
% % %     tin = squeeze(RFi(cy-1:cy+1,cx-1:cx+1,:))-meani;
% % %     tpl = squeeze(tin(1,1,:));
% % %     tpl1 = squeeze(tin(2,2,:));
% % %     tpl2 = squeeze(tin(3,3,:));
% % %     tpp = mean(horzcat(tpl, tpl1, tpl2),2);
% % %     end 
% % % %     figure
% % % %     plot(smooth(tpl), 'Color', 'r', 'LineWidth', 1.2) 
% % % %     hold on 
% % % %     plot(smooth(tpl1), 'Color', 'c', 'LineWidth', 1.2) 
% % % %     plot(smooth(tpl2), 'Color', 'm', 'LineWidth', 1.2) 
% % % %     plot(smooth(tpp), 'Color', 'k', 'LineWidth', 1.2) 
% % % %     title(string(ci))
% % %     
% % % 
% % %     end 
% % % end 






%%

% figure 
% offset = 0; 
% for i = 1:890 
%     plot(smooth(TempRF(i,1:30)-offset), 'LineWidth', 1.2)
%     hold on 
%     offset = offset+10;
% end 

T2  = sortrows(SpatRF, 31);

%  T2 = T2(1:859, :);

% T2 = T2(1:1088, :);


% Remove rows with Rf outside of normal VF place. 
T2B  = T2; 
%rowstorm = find(T2(:, 34)<200 | T2(:, 34)>500 | T2(:, 35)<50 | T2(:, 35)>300); 
 rowstorm = find(T2(:, 34)<250 | T2(:, 34)>450 | T2(:, 35)<100 | T2(:, 35)>250); 
T2B(rowstorm, :) = []; 

n = numel(T2B(:,1));

%zscore
for i = 1:n
    T2B(i, 1:30) = zscore(T2B(i, 1:30)); 
end 


for i = 1:n%310
    [mv, mi] = max(T2B(i,1:30));
    T2B(i,36) = mv;
    T2B(i,37) = mi;
    [mv2, mi2] = min(T2B(i,1:30));
    T2B(i,38) = mv2;
    T2B(i,39) = mi2;
    
    av_val = mean(T2B(i, 25:29));
    T2B(i,40) = av_val;
    if av_val <=0
        T2B(i, 41) = 0; 
    else 
        T2B(i, 41) = 1; 
    end 
    
    if T2B(i,30)<T2B(i,28)
        T2B(i,42) = 1; 
    else 
        T2B(i,42) = 0; 
    end 
    
end 
    

allON = find(T2B(:, 42) == 1); 
allOFF = find(T2B(:, 42) == 0); 

% Make new arrays of Off/On cells
TON = T2B(allON, :); 
TOFF = T2B(allOFF, :);


figure
plot((TOFF(:, 1:30))', 'Color', [0.75, 0.75, 0.75])
hold on 
plot(mean(TOFF(:, 1:30))', 'Color', 'k', 'LineWidth', 1.2)
title('OFF')

figure
plot((TON(:, 1:30))', 'Color', [0.75, 0.75, 0.75])
hold on 
plot(mean(TON(:, 1:30))', 'Color', 'k', 'LineWidth', 1.2)
title('ON')


% save('201116_T2B_HET.mat', 'T2B', 'TON', 'TOFF');

TON2 = sortrows(TON, 32);
TOFF2 = sortrows(TOFF, 32);

% offset = 0; 
% figure
% for j = 1:147
% plot((TOFF2(j, 1:30))'-offset, 'Color', [0.75, 0.75, 0.75])
% hold on 
% offset = offset+3; 
% end 

figure
imagesc(TOFF2(:, 1:30))
colormap(redblue)
    depthsoff = string(TOFF2(:,32)); 
    noff = numel(TOFF2(:,1));
    yticks([1:1:noff])
    yticklabels(depthsoff)



figure
imagesc(TON2(:, 1:30))
colormap(redblue)
    depthson = string(TON2(:,32)); 
    non = numel(TON2(:,1));
    yticks([1:1:non])
    yticklabels(depthson)





%% 
TON_sorted = sortrows(TON, 32);
figure 
imagesc(TON_sorted(:, 1:30))
colormap(redblue)



%% Run through the animals. 
% Make a plot - line plot of all traces + Mean
% Heatmap by depth - add yticklabels as DEPTH. 

all_animals = unique(T2B(:, 33));


for i = 1:4 
    figure
    animal = all_animals(i); 
    allANI = find(T2B(:, 33) == animal); 
    TANI = T2B(allANI, :); 
    
    allON = find(TANI(:, 42) == 1);  
    allOFF = find(TANI(:, 42) == 0); 

    TON = TANI(allON, :); 
    TOFF = TANI(allOFF, :);
    
    subplot(3,2,1)
plot((TOFF(:, 1:30))', 'Color', [0.75, 0.75, 0.75])
hold on 
plot(mean(TOFF(:, 1:30))', 'Color', 'k', 'LineWidth', 1.2)
title('OFF')

    subplot(3,2,2)
plot((TON(:, 1:30))', 'Color', [0.75, 0.75, 0.75])
hold on 
plot(mean(TON(:, 1:30))', 'Color', 'k', 'LineWidth', 1.2)
title('ON')

    TON2 = sortrows(TON, 32);
    TOFF2 = sortrows(TOFF, 32);

    subplot(3,2,[3,5])
    imagesc(TOFF2(:, 1:30))
    depthsoff = string(TOFF2(:,32)); 
    noff = numel(TOFF2(:,1));
    yticks([1:1:noff])
    yticklabels(depthsoff)
    colormap(redblue)
    
    subplot(3,2,[4,6])
    imagesc(TON2(:, 1:30))
    depthson = string(TON2(:,32)); 
    non = numel(TON2(:,1));
    yticks([1:1:non])
    yticklabels(depthson)
    colormap(redblue)
    
    sgtitle(string(animal))
end 



%% Position of RF. on versus off 


for i = 1:4 
    figure
    animal = all_animals(i); 
    allANI = find(T2B(:, 33) == animal); 
    TANI = T2B(allANI, :); 
    
    allON = find(TANI(:, 42) == 1);  
    allOFF = find(TANI(:, 42) == 0);
    
%     nON = numel(allON);
%     nOFF = numel(allOFF); 

    TON = TANI(allON, :); 
    TOFF = TANI(allOFF, :);
    
    subplot(2,2,1)
plot((TOFF(:, 1:30))', 'Color', [0.75, 0.75, 0.75])
hold on 
plot(mean(TOFF(:, 1:30))', 'Color', 'k', 'LineWidth', 1.2)
title('OFF')

    subplot(2,2,2)
plot((TON(:, 1:30))', 'Color', [0.75, 0.75, 0.75])
hold on 
plot(mean(TON(:, 1:30))', 'Color', 'k', 'LineWidth', 1.2)
title('ON')

    subplot(2,2,3)
%     plot(T2(:,34), T2(:,35),'k.', 'MarkerSize', 4)
%     hold on 
    plot(TOFF(:,34), TOFF(:,35), 'k.', 'MarkerSize', 15)
    axis([0 600 0 340])
    
    subplot(2,2,4)
%     plot(T2(:,34), T2(:,35),'k.', 'MarkerSize', 4)
%     hold on 
    plot(TON(:,34), TON(:,35), 'k.', 'MarkerSize', 15)
    axis([0 600 0 340])
    
    sgtitle(string(animal))
end 
















% T2B(:,36) = maxvals(:,2);
% T2B(:,37) = minvals(:,2);





% 
% figure 
% offset = 0; 
% for i = 830:859 
%     plot(smooth(T2(i,1:30)-offset), 'LineWidth', 1.2)
%     hold on 
%     offset = offset+10;
% end 
% 
% 
% figure 
% for i = 1:1088 
%     plot(smooth(T2(1000,1:30)), 'LineWidth', 1.2)
%     hold on 
% end 
% 
% figure
% imagesc(smooth(T2(1:50, 1:30))')
% imagesc((T2(800:1088, 1:30)))


% imagesc((T3(1000:1080, 1:30)))


% Find time to MAX and time to MIN - Add columns. Then sort by these
% values. 

% for i = 1:1088
%  [mv, mi] = max(T2(i,1:30));
% maxvals(i,1) = mv; 
% maxvals(i,2) = mi; 
% [mv2, mi2] = min(T2(i,1:30));
% minvals(i,1) = mv2; 
% minvals(i,2) = mi2; 
% end 
% 
% T2(:,36) = maxvals(:,2);
% T2(:,37) = minvals(:,2);
% 
% 
% T3  = sortrows(T2, 36);
% T4 =  sortrows(T2, 37);
% 
% imagesc((T3(1:1088, 1:30)))
% colorbar
% 
% imagesc((T4(1:1088, 1:30)))

%%
figure

for i = 1:587 
    xval = TON_sorted(i, 34);
    yval = TON_sorted(i, 35);
    plot(xval, yval, 'k.', 'MarkerSize', 20)
    hold on 
end 


for i = 1:501 
    xval = TOFF(i, 34);
    yval = TOFF(i, 35);
    plot(xval, yval, 'm.', 'MarkerSize', 20)
    hold on 
end 

% Remove rows when RF is found in weird area of VF. 
rowstorm = find(TON(:, 34)<250 | TON(:, 34)>450 | TON(:, 35)<100 | TON(:, 35)>250); 



