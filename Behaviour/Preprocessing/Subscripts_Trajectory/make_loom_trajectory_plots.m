% Make Trajectory Plots for post loom til shelter. 

% Make spatial heatmaps  - hotter colour for speed assocaited with that
% space. Do mice react faster when seeing loom from certain position in the
% box? 

% Find all the het and the wt rows.
% Make two arrays for each - one x and one y. 

allWT = find(string(ALL_XYLOOM_POS_TABLE.Geno) == "wt"); 
allHET = find(string(ALL_XYLOOM_POS_TABLE.Geno) == "het"); 

n = height(ALL_XYLOOM_POS_TABLE);
       
xWT = [];
yWT = []; 
xHET = [];
yHET = []; 

for i = 1:n
    if string(ALL_XYLOOM_POS_TABLE.Geno{i}) == "wt"
        G = cell2mat(ALL_XYLOOM_POS_TABLE{i,5});
        xWT = vertcat(xWT, G); 
        G2 = cell2mat(ALL_XYLOOM_POS_TABLE{i,6});
        yWT = vertcat(yWT, G2); 
        
    elseif string(ALL_XYLOOM_POS_TABLE.Geno{i}) == "het"
        F = cell2mat(ALL_XYLOOM_POS_TABLE{i,5});
        xHET = vertcat(xHET, F);
        F2 = cell2mat(ALL_XYLOOM_POS_TABLE{i,6});
        yHET = vertcat(yHET, F2);        
    end
end

%% Plot Trajectories

%Info of the shelter. 
boxC = [430 78]; 
radius = 90; 

% WT 

figure
for q = 1:50 
    for j = 180:779  
        x = xWT(q,j);
        y = 518 - yWT(q,j); 
        x2 = xWT(q,j+1);
        y2 = 518 - yWT(q,j+1);
        plot([x, x2],[y,y2],'k')
        hold on 
    end 
plot(xWT(q,180), (518 - yWT(q,180)), 'r.', 'MarkerSize', 15)
end 
axis([0 520 0 520])
hold on 
viscircles(boxC, radius, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 0.5)

title('Escape Trajectories - WT')
savefig(gcf, 'EscapeTrajectories_WT.fig')
close

%%
figure
for q = 1:20 
    for j = 180:779  
        x = xHET(q,j);
        y = 518 - yHET(q,j); 
        x2 = xHET(q,j+1);
        y2 = 518 - yHET(q,j+1);
        plot([x, x2],[y,y2],'r')
        hold on 
    end 
plot(xHET(q,180), (518 - yHET(q,180)), 'k.', 'MarkerSize', 15)
end 
axis([0 520 0 520])
hold on 
viscircles(boxC, radius, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
% title('First 20 Escape Trajectories - HET')
title('Escape Trajectories - HET')
savefig(gcf, 'EscapeTrajectories_HET.fig')
close




%% PER ANIMAL
    
x_all = [];
y_all = []; 
for i = 1:77
    if string(ALL_XYLOOM_POS_TABLE.Animal{i}) == "GN6620"
        G = cell2mat(ALL_XYLOOM_POS_TABLE{i,5});
        x_all = vertcat(x_all, G); 
        G2 = cell2mat(ALL_XYLOOM_POS_TABLE{i,6});
        y_all = vertcat(y_all, G2); 
    end     
end

n = numel(x_all(:,1)); 
figure
for q = 1:n
    for j = 180:779  
        x = x_all(q,j);
        y = 518 - y_all(q,j); 
        x2 = x_all(q,j+1);
        y2 = 518 - y_all(q,j+1);
        plot([x, x2],[y,y2],'r')
        hold on 
    end 
plot(x_all(q,180), (518 - y_all(q,180)), 'k.', 'MarkerSize', 15);
end 
axis([0 520 0 520]);
hold on 
viscircles(boxC, radius, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5);
title('Escape Trajectories - GN6620');


savefig(gcf, 'EscapeTrajectories_6620.fig')
close 



%% Response Speed versus Distance from Shelter. 

% ADd 7th column to ALL_XYLOOM_POS_TABLE - distance from box centre when
% loom takes place. 

boxCentre = [430 440]; 

for r = 1:77
    x_array = cell2mat(ALL_XYLOOM_POS_TABLE{r,5});
    y_array = cell2mat(ALL_XYLOOM_POS_TABLE{r,6});
    x_loom = x_array(1,180); 
    y_loom = y_array(1,180); 
    dist_centre = pdist([boxCentre(1), boxCentre(2) ; x_loom y_loom])*(32/518); 
    
    ALL_XYLOOM_POS_TABLE.DistCentre{r} = dist_centre; 
end 


% Adding 8th column to ALL_XYLOOM_POS_TABLE - the AVERAGE RESPONSE SPEED per
% row. 

for d = 1:93
    speed_array = cell2mat(ALL_XYLOOM_TABLE{d,5});
    av = mean(speed_array); 
    ALL_XYLOOM_POS_TABLE.Speed{d} = av; 
end 


allWT = find(string(ALL_XYLOOM_POS_TABLE.Geno) == "wt"); 
allHET = find(string(ALL_XYLOOM_POS_TABLE.Geno) == "het"); 

figure
for i = 1:93
    if string(ALL_XYLOOM_POS_TABLE.Geno{i}) == "wt"
        marker = 'k.'; 
    elseif string(ALL_XYLOOM_POS_TABLE.Geno{i}) == "het"
        marker = 'r.'; 
    end 
    Dis = ALL_XYLOOM_POS_TABLE.DistCentre{i}; 
    Sp = ALL_XYLOOM_POS_TABLE.Speed{i}; 
    plot(Dis, Sp, marker, 'MarkerSize', 20)
    hold on 
end 

axis([10 35 2 14])
title('Speed of Response versus Distance from Shelter when loom happened')
xlabel('Distance from shelter - cm')
ylabel('Speed of response - cm/s')

savefig(gcf, 'ResponseSpeedvsDistFromShelter.fig')
close


%% Individual analysis of Acclim Day. - Trajectories and Speed analysis

      
n = 5023; 
figure
 for j = 1:n-1  
        x = xy_array(j,3);
        y = 518 - xy_array(j,4); 
        x2 = xy_array(j+1,3);
        y2 = 518 - xy_array(j+1,4);
        plot([x, x2],[y,y2],'r')
        hold on 
end 
axis([0 520 0 520]);
title('Acclim Trajectory - GN4473');


savefig(gcf, 'Acclim_Traj_4473.fig')
close


% acclim_array2 = zeros(8,5); 

%Each ROW = 1 animal. 
animal = 8; 

% Col 1  = total D travelled - cm. 
total_distance_travlled = sum(xy_array(:,5)); 
acclim_array2(animal,1) = total_distance_travlled;

% Col 2  = average speed 
acclim_array2(animal,2) = mean(xy_array(:,5));

% Col 3  = Max Speed 
fast = find(xy_array(:,5)<180); 
fastrows = xy_array(fast,5);
acclim_array2(animal,3) = max(fastrows);
% max(xy_array(:,5));

% Col 4 = Average speed (threshold for speeds >2cm/s)
faster = find(xy_array(:,5)>2); 
av_fast = mean(xy_array(faster,5)); 
acclim_array2(animal,4) = av_fast;

% Col 5 = Time spent outside of box. 
outside = find(xy_array(:,6)>0); 
acclim_array2(animal,5) = numel(outside)/60;

acclim_array2(:,1) = acclim_array2(:,1)/100; 

TotalD = acclim_array2(:,1); % metres!!! 
Speed = acclim_array2(:,2);
Max = acclim_array2(:,3);
AverageFast = acclim_array2(:,4);
Outside = acclim_array2(:,5);

Indiv_Animals_Acclim_Table = table(TotalD, Speed, Max, AverageFast, Outside, 'RowNames', {'3956', '3959', '4242', '4244', '4468', '4469', '4472', '4473'}); 
save('Indiv_Animals_Acclim_Table.mat', 'Indiv_Animals_Acclim_Table'); 

%  WT vs HET 

Acclim_Array3 = zeros(2,5); 

%D 
Acclim_Array3(1,1) = mean(acclim_array2([1,3,5,6], 1)); 
Acclim_Array3(2,1) = mean(acclim_array2([2,4,7,8], 1)); 

for p = 2:5
Acclim_Array3(1,p) = mean(acclim_array2([1,3,5,6], p)); 
Acclim_Array3(2,p) = mean(acclim_array2([2,4,7,8], p)); 
end











%% WT vs HET spatial heatmap.

% Vertcat x  - all WTs. 
% vertcat y - all WTs. 

xy_WT = []; 
xy_HET = []; 

xy_animal = xy_array(:,3:4); 

% xy_WT = vertcat(xy_WT, xy_animal); 
% save('xy_WT.mat', 'xy_WT')

% xy_HET = vertcat(xy_HET, xy_animal); 
% save('xy_HET.mat', 'xy_HET')

%% Spatial heatmap for WTs

Xdata=xy_WT(:,1);
Ydata=xy_WT(:,2);

edges=0:74:518; %largest value of x is 348 and largest value of y is 349 so going to 360 not 362 is sufficient.

n=size(Xdata);

Xbindata=discretize(Xdata, edges, 'IncludedEdge', 'Right');
xy_WT(:,3)=Xbindata;
Ybindata=discretize(Ydata, edges, 'IncludedEdge', 'Right');
xy_WT(:,4)=Ybindata;

%% Adding DeltaF/F data into centroids_array
data = xy_WT;

%% Create array of Spatial Activity
NumF = zeros(7,7);

for i = 1:n  %number of rows/frames in the video
NumF(data(i,3), data(i,4)) = NumF(data(i,3),data(i,4))+1; %Sum the number of frames 
end

NumF_log = log(NumF);
figure
imagesc(NumF_log)
colorbar
title('Log of the time spent in each location - WT')
colormap('parula')
savefig(gcf, 'Log_spatialHEATMAP_WT_HOT.fig')










%% Spatial heatmap for HETs

Xdata=xy_HET(:,1);
Ydata=xy_HET(:,2);

edges=0:74:518; %largest value of x is 348 and largest value of y is 349 so going to 360 not 362 is sufficient.

n=size(Xdata);

Xbindata=discretize(Xdata, edges, 'IncludedEdge', 'Right');
xy_HET(:,3)=Xbindata;
Ybindata=discretize(Ydata, edges, 'IncludedEdge', 'Right');
xy_HET(:,4)=Ybindata;

%% Adding DeltaF/F data into centroids_array
data = xy_HET;

%% Create array of Spatial Activity
NumF = zeros(7,7);

for i = 1:n  %number of rows/frames in the video
NumF(data(i,3), data(i,4)) = NumF(data(i,3),data(i,4))+1; %Sum the number of frames 
end

%% Save all the figures
NumF_log = log(NumF);
figure
imagesc(NumF_log)
colorbar
title('Log of the time spent in each location - HET')
colormap('hot')
savefig(gcf, 'Log_spatialHEATMAP_HET_HOT.fig')


