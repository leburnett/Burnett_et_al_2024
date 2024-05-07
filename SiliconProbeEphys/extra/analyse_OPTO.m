% Analyse OPTO

% Burnett - 23/07/21
% To be used after 'Core_Loom_Processing.m'

%% 
files = dir('*Spikes_Opto*'); 
num_files = numel(files);
het_animals = ["7614", "7476", "7790", "7269","1970", "2833", "1385", "1386", "1388", "1394"]; 

%

W= 595; 

L1 = [];
R1 = []; 
L2 = [];
R2 = [];



for i = 1:num_files
    fname = files(i).name;
    load(fname);
    
    ani = str2double(fname(10:13));
    probe_depth = str2double(fname(21:24));
    date =  str2double(fname(1:6)); 
    
    if contains(string(ani), het_animals)
        geno = 0; 
    else 
        geno = 1;
    end 
    
      n_cl = numel(all_spikes_L1(:,1));
    
    %% Speed1 
    
%     W = W1; 
   
    w = numel(all_spikes_L1(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            all_spikes_L1(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        all_spikes_L1(di, W+1) = date; 
        all_spikes_L1(di, W+2) = ani; 
        all_spikes_L1(di, W+3) = real_depth;
        all_spikes_L1(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        L1 = vertcat(L1, all_spikes_L1); 
        
        
        %% Speed2
        
%         W = W2; 
   
    w = numel(all_spikes_L2(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            all_spikes_L2(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        all_spikes_L2(di, W+1) = date; 
        all_spikes_L2(di, W+2) = ani; 
        all_spikes_L2(di, W+3) = real_depth;
        all_spikes_L2(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        R1 = vertcat(R1, all_spikes_L2); 
        
        
        %% Speed 3
        
%         W = W3; 
        
         w = numel(all_spikes_L3(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            all_spikes_L3(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        all_spikes_L3(di, W+1) = date; 
        all_spikes_L3(di, W+2) = ani; 
        all_spikes_L3(di, W+3) = real_depth;
        all_spikes_L3(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        L2 = vertcat(L2, all_spikes_L3); 
        
        
        
        %% Speed4
        
%         W = W4; 
        
         w = numel(all_spikes_L4(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            all_spikes_L4(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        all_spikes_L4(di, W+1) = date; 
        all_spikes_L4(di, W+2) = ani; 
        all_spikes_L4(di, W+3) = real_depth;
        all_spikes_L4(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        R2 = vertcat(R2, all_spikes_L4); 
          
        
end 
    

L1(:, 1:W) = L1(:, 1:W)*60; 
R1(:, 1:W) = R1(:, 1:W)*60; 
L2(:, 1:W) = L2(:, 1:W)*60; 
R2(:, 1:W) = R2(:, 1:W)*60; 

 save('210723_Spiking_Ptchd1_OPTO.mat', 'L1', 'L2', 'R1', 'R2');
 
%%
allWT = find(L1(:, W1+4) == 1); 
allHET = find(L1(:, W1+4) == 0); 
 

% Colour 
v1 = 0.1; 
col = [1-v1 (114/255)-v1 (32/255)-v1]; % orange ptchd1
% col = 'r'; 

%% L1

figure
subplot(1,2,1)
imagesc(L1(allWT, 1:W1))
title('WT')
caxis([0 200])
subplot(1,2,2)
imagesc(L1(allHET, 1:W1))
caxis([0 200])
title('HET')


figure
subplot(2,1,1)
plot(mean(L1(allWT, 1:W1)), 'k')
hold on 
plot(mean(L1(allHET, 1:W1)), 'Color', col)
box off
hold off
subplot(2,1,2)
plot(smooth(mean(L1(allWT, 1:W1))), 'k')
hold on 
plot(smooth(mean(L1(allHET, 1:W1))), 'Color', col)
hold off
box off

%% Speed 2


figure
subplot(1,2,1)
imagesc(R1(allWT, 1:W2))
title('WT')
caxis([0 200])
subplot(1,2,2)
imagesc(R1(allHET, 1:W2))
caxis([0 200])
title('HET')



figure
subplot(2,1,1)
plot(mean(R1(allWT, 1:W2)), 'k')
hold on 
plot(mean(R1(allHET, 1:W2)), 'Color', col)
hold off
box off
subplot(2,1,2)
plot(smooth(mean(R1(allWT, 1:W2))), 'k')
hold on 
plot(smooth(mean(R1(allHET, 1:W2))), 'Color', col)
hold off
box off

%%  Speed 3 

figure
subplot(1,2,1)
imagesc(L2(allWT, 1:W3))
title('WT')
caxis([0 200])
subplot(1,2,2)
imagesc(L2(allHET, 1:W3))
caxis([0 200])
title('HET')



figure
subplot(2,1,1)
plot(mean(L2(allWT, 1:W3)), 'k')
hold on 
plot(mean(L2(allHET, 1:W3)), 'Color', col)
hold off
box off
subplot(2,1,2)
plot(smooth(mean(L2(allWT, 1:W3))), 'k')
hold on 
plot(smooth(mean(L2(allHET, 1:W3))), 'Color', col)
hold off
box off


%% Speed 4


figure
subplot(1,2,1)
imagesc(R2(allWT, 1:W4))
caxis([0 200])
title('WT')
subplot(1,2,2)
imagesc(R2(allHET, 1:W4))
caxis([0 200])
title('HET')



figure
subplot(2,1,1)
plot(mean(R2(allWT, 1:W4)), 'k')
hold on 
plot(mean(R2(allHET, 1:W4)), 'Color', col)
hold off
box off
subplot(2,1,2)
plot(smooth(mean(R2(allWT, 1:W4))), 'k')
hold on 
plot(smooth(mean(R2(allHET, 1:W4))), 'Color', col)
hold off
box off


%% Speed 5


figure
subplot(1,2,1)
imagesc(L5(allWT, 1:W5))
caxis([0 200])
title('WT')
subplot(1,2,2)
imagesc(L5(allHET, 1:W5))
caxis([0 200])
title('HET')



figure
subplot(2,1,1)
plot(mean(L5(allWT, 1:W5)), 'k')
hold on 
plot(mean(L5(allHET, 1:W5)), 'Color', col)
hold off
box off
subplot(2,1,2)
plot(smooth(mean(L5(allWT, 1:W5))), 'k')
hold on 
plot(smooth(mean(L5(allHET, 1:W5))), 'Color', col)
hold off
box off

%% KMEANS - finding clusters

 n = 10; 
 
 data = L5(:, 1:W5); 
 vals = kmeans(data, n);
%  data2 = sortrows(data, vals);
 
 for i = 1:n
     v1 = find(vals == i); 
     
    subplot(n,1,i)
    plot(mean(data(v1, :)));

 end 
 
 %% Remove cells with very low firing rates
 
 