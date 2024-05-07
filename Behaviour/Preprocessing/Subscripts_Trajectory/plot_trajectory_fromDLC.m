% Plot the trajectory of the mouse from DLC results. 

%% Read in the values using csvread. 
filename  = 'C:\Data_analysis\DLC\WhiteBKG\06_Loom_19S09_191023DeepCut_resnet50_WhiteBKGOct25shuffle1_500000.csv';
csvDATA = csvread(filename,3);
% Col 2 is Nose X
% Col 3 is Nose y
% Col 4 is Nose likelihood.

%% Making a table - less useful but can see the labels. 
%csvTABLE = readtable(filename, 'ReadVariableNames', false, 'HeaderLines', 1);

%% If likelihood is below 0.8 - make the coordinate the same as the coordinate before. 
n = length(csvDATA);

for i = 1:n
    if csvDATA(i,4)< 0.8
        csvDATA(i,3)= csvDATA(i-1,3);
        csvDATA(i,2)= csvDATA(i-1,2);
    end 
end 

%% Plot the trajectory of the nose. From Make_EXP_SUM.m

figure 
hold on 
cmap = hsv(n);

for j = 1: n-1
    x1 = csvDATA(j,2);
    x2 = csvDATA(j+1,2); 
    y1 = csvDATA(j,3); 
    y2 = csvDATA(j+1,3);
 plot([x1 x2],[y1 y2], 'color', cmap(j,:))
 %drawnow %using drawnow makes the plot animated. 
end 


%% From dynamic_traj.m 

h = animatedline('Color', 'k', 'LineStyle', ':');
h2 = animatedline('LineStyle', 'none', 'Marker', '.', 'Color', 'r');

for i = 1:n
    x = mouse_coords(i,1);
    y = mouse_coords(i,2);
    addpoints(h, mouse_coords(i,1), mouse_coords(i,2));
    drawnow
    x2 = dot_coords(i,1);
    y2 = dot_coords(i,2);
    addpoints(h2, dot_coords(i,1), dot_coords(i,2));
    drawnow
end 

%% Making gif. 

% Make gif video. 
% % gif('track_darkbkg.gif');
% % gif. 