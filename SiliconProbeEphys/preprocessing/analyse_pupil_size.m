%% Analyse pupil size from DLC
% Created by Burnett 25/06/21

csv_files = dir('*.csv');
n_files = length(csv_files);

for ii = 1:n_files
filename = csv_files(ii).name;
data_all = readmatrix(filename);

nchr = numel(filename);
if nchr == 89
exp_name = filename(1:38);
elseif nchr == 92
   exp_name = filename(1:41); 
elseif nchr == 84
    exp_name = filename(1:31); 
end 
    
% TOp , right, bottom, left. 3 columns each 

% Col 1 = frame#
% Col 2 = top x
% Col 3 = top y
% Col4 = prob

% x values = 2, 5, 8, 11
% y values = 3, 6, 9, 12

% figure
% plot(data_all(1,2), data_all(1,3), 'k.');
% hold on 
% plot(data_all(1,5), data_all(1,6), 'k.');
% plot(data_all(1,8), data_all(1,9), 'k.');
% plot(data_all(1,11), data_all(1,12), 'k.');
% plot(centre(1), centre(2), 'r.')

for p = [2,3,5,6,8,9,11,12]
data_all(:, p) = movmean(data_all(:, p), 5);
data_all(:,p) = filloutliers(data_all(:,p), 'center');
end 

%%
nframes = numel(data_all(:,1)); 

d = [];
d1 = []; 
d2 = []; 
d3 = []; 
d4= []; 
top = []; 
right = []; 
bottom = []; 
left = []; 

data = []; 

% Want centre and radius for each frame! 
for i = 1: nframes

top = [data_all(i,2), data_all(i,3)];
right = [data_all(i,5), data_all(i,6)];
bottom = [data_all(i,8), data_all(i,9)];
left = [data_all(i,11), data_all(i,12)];

mid_tb = [(top(1)+bottom(1))/2, (top(2)+bottom(2))/2];
mid_lr = [(left(1)+right(1))/2, (left(2)+right(2))/2];

centre = [(mid_tb(1)+mid_lr(1))/2, (mid_tb(2)+mid_lr(2))/2];

d1 = pdist(vertcat(centre,top));
d2 = pdist(vertcat(centre,left));
d3 = pdist(vertcat(centre,bottom));
d4 = pdist(vertcat(centre,right));

d = mean([d1,d2,d3,d4]);

%
data.topX(i) = data_all(i,2);
data.topY(i) = data_all(i,3);
data.topP(i) = data_all(i,4);

data.leftX(i) = data_all(i,5);
data.leftY(i) = data_all(i,6);
data.leftP(i) = data_all(i,7);

data.bottomX(i) = data_all(i,8);
data.bottomY(i) = data_all(i,9);
data.bottomP(i) = data_all(i,10);

data.rightX(i) = data_all(i,11);
data.rightY(i) = data_all(i,12);
data.rightP(i) = data_all(i,13);

data.CentreX(i) = centre(1);
data.CentreY(i) = centre(2);
data.dist(i) = d; 

end 

av_centreX = mean(data.CentreX);
av_centreY = mean(data.CentreY); 

for j = 1:nframes
    
    movX = av_centreX - data.CentreX(j);
    movY = av_centreY - data.CentreY(j);
    
    data.movX(j) = movX;
    data.movY(j) = movY;
end 

figure
ax1 = subplot(2,1,1);
plot(data.dist)
title('Radius')
ax2 = subplot(2,1,2);
plot(data.movX, 'r');
hold on 
plot(data.movY, 'k'); 
title('Centroid x,y')
linkaxes([ax1, ax2], 'x')

savefig(gcf, strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Eye_DLC/eye_files/Eye_movement_', exp_name, '.fig'));

save(strcat('/Users/lauraburnett/Data_Analysis_Mac/SiliconProbe/Eye_DLC/eye_files/', exp_name, '_eye_data.mat'), 'data');

close
data_all = []; 
end 

















