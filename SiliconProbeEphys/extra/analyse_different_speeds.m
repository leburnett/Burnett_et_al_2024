% Analyse DIFFERENT SPEEDS

% Burnett - 22/07/21
% To be used after 'Core_Loom_Processing.m'

%% 

files = dir('*Spikes_5LoomsDS*'); 
num_files = numel(files);
het_animals = ["7614", "7476", "7269", "2833","3557", "4124", "1394", "2709", "4369"]; 

% 5 SPEEDS 

% L1 = 200 
% L2 = 275
% L3 = 350
% L4 = 425
% L5 = 500

W1 = 200; 
W2 = 275; 
W3 = 350; 
W4 = 425; 
W5 = 500; 

L1 = [];
L2 = []; 
L3 = [];
L4 = [];
L5 = []; 



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
    
      n_cl = numel(allsp_L1(:,1));
    
    %% Speed1 
    
    W = W1; 
   
    w = numel(allsp_L1(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            allsp_L1(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        allsp_L1(di, W+1) = date; 
        allsp_L1(di, W+2) = ani; 
        allsp_L1(di, W+3) = real_depth;
        allsp_L1(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        L1 = vertcat(L1, allsp_L1); 
        
        
        %% Speed2
        
        W = W2; 
   
    w = numel(allsp_L2(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            allsp_L2(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        allsp_L2(di, W+1) = date; 
        allsp_L2(di, W+2) = ani; 
        allsp_L2(di, W+3) = real_depth;
        allsp_L2(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        L2 = vertcat(L2, allsp_L2); 
        
        
        %% Speed 3
        
        W = W3; 
        
         w = numel(allsp_L3(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            allsp_L3(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        allsp_L3(di, W+1) = date; 
        allsp_L3(di, W+2) = ani; 
        allsp_L3(di, W+3) = real_depth;
        allsp_L3(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        L3 = vertcat(L3, allsp_L3); 
        
        
        
        %% Speed4
        
        W = W4; 
        
         w = numel(allsp_L4(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            allsp_L4(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        allsp_L4(di, W+1) = date; 
        allsp_L4(di, W+2) = ani; 
        allsp_L4(di, W+3) = real_depth;
        allsp_L4(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        L4 = vertcat(L4, allsp_L4); 
        
        
              %% Speed5 
        
        W = W5; 
        
         w = numel(allsp_L5(1,:));
    if w~= W
        num_missing = W-w;
        added = zeros(1, num_missing);
        for j = 1:n_cl
            allsp_L5(j,(w+1:W)) =  added;
        end 
    end 
        
    for di = 1:n_cl
        real_depth = probe_depth - (800-ids_depth(di,2)); % 0 is TIP OF PROBE - NOT - TOMAS FLIPPED - 0 is TOP of probe. 
        allsp_L5(di, W+1) = date; 
        allsp_L5(di, W+2) = ani; 
        allsp_L5(di, W+3) = real_depth;
        allsp_L5(di, W+4) = geno; 
    end 

     %Vertically concatenate all the cells for each 
        L5 = vertcat(L5, allsp_L5); 
        
        
end 
    

L1(:, 1:W1) = L1(:, 1:W1)*60; 
L2(:, 1:W2) = L2(:, 1:W2)*60; 
L3(:, 1:W3) = L3(:, 1:W3)*60; 
L4(:, 1:W4) = L4(:, 1:W4)*60; 
L5(:, 1:W5) = L5(:, 1:W5)*60; 

 save('210827_Spiking_Ptchd1_DifferentSpeeds_N2.mat', 'L1', 'L2', 'L3', 'L4', 'L5');
 
%%
allWT = find(L1(:, W1+4) == 1); 
allHET = find(L1(:, W1+4) == 0); 
 

% Colour 
v1 = 0.1; 
col = [1-v1 (114/255)-v1 (32/255)-v1]; % orange ptchd1
% col = 'r'; 
% col = 'm';
% 
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
imagesc(L2(allWT, 1:W2))
title('WT')
caxis([0 200])
subplot(1,2,2)
imagesc(L2(allHET, 1:W2))
caxis([0 200])
title('HET')



figure
subplot(2,1,1)
plot(mean(L2(allWT, 1:W2)), 'k')
hold on 
plot(mean(L2(allHET, 1:W2)), 'Color', col)
hold off
box off
subplot(2,1,2)
plot(smooth(mean(L2(allWT, 1:W2))), 'k')
hold on 
plot(smooth(mean(L2(allHET, 1:W2))), 'Color', col)
hold off
box off

%%  Speed 3 

figure
subplot(1,2,1)
imagesc(L3(allWT, 1:W3))
title('WT')
caxis([0 200])
subplot(1,2,2)
imagesc(L3(allHET, 1:W3))
caxis([0 200])
title('HET')



figure
subplot(2,1,1)
plot(mean(L3(allWT, 1:W3)), 'k')
hold on 
plot(mean(L3(allHET, 1:W3)), 'Color', col)
hold off
box off
subplot(2,1,2)
plot(smooth(mean(L3(allWT, 1:W3))), 'k')
hold on 
plot(smooth(mean(L3(allHET, 1:W3))), 'Color', col)
hold off
box off


%% Speed 4


figure
subplot(1,2,1)
imagesc(L4(allWT, 1:W4))
caxis([0 200])
title('WT')
subplot(1,2,2)
imagesc(L4(allHET, 1:W4))
caxis([0 200])
title('HET')



figure
subplot(2,1,1)
plot(mean(L4(allWT, 1:W4)), 'k')
hold on 
plot(mean(L4(allHET, 1:W4)), 'Color', col)
hold off
box off
subplot(2,1,2)
plot(smooth(mean(L4(allWT, 1:W4))), 'k')
hold on 
plot(smooth(mean(L4(allHET, 1:W4))), 'Color', col)
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
 
 data = L5(allWT, 1:W5); 
 vals = kmeans(data, n);
%  data2 = sortrows(data, vals);
 figure
 for i = 1:n
     v1 = find(vals == i); 
     
    subplot(n,1,i)
    plot(mean(data(v1, :)), 'k');

 end 
 
 
  n = 10; 
 
 data = L5(allHET, 1:W5); 
 vals = kmeans(data, n);
%  data2 = sortrows(data, vals);
 figure
 for i = 1:n
     v1 = find(vals == i); 
     
    subplot(n,1,i)
    plot(mean(data(v1, :)), 'r');

 end 
 %% Remove cells with very low firing rates
 
 
 
 
 %%
 
 
 
figure
subplot(5,1,1)
plot(mean(L1(allWT, 1:W1)), 'k')
hold on 
plot(mean(L1(allHET, 1:W1)), 'Color', col)
box off
hold off
% xlim([45 105])

subplot(5,1,2)
plot(mean(L2(allWT, 1:W2)), 'k')
hold on 
plot(mean(L2(allHET, 1:W2)), 'Color', col)
box off
hold off
% xlim([45 105])

subplot(5,1,3)
plot(mean(L3(allWT, 1:W3)), 'k')
hold on 
plot(mean(L3(allHET, 1:W3)), 'Color', col)
box off
hold off
% xlim([45 105])

subplot(5,1,4)
plot(mean(L4(allWT, 1:W4)), 'k')
hold on 
plot(mean(L4(allHET, 1:W4)), 'Color', col)
box off
hold off
% xlim([45 105])

subplot(5,1,5)
plot(mean(L5(allWT, 1:W5)), 'k')
hold on 
plot(mean(L5(allHET, 1:W5)), 'Color', col)
box off
hold off
% xlim([45 105])


 