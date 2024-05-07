%% Have 'ACTIVITY\xy_arrays' folder.  

% Load XY_POS_TABLE

allWT = find(string(ALL_XY_POS_TABLE.Geno) =="wt");
allHET = find(string(ALL_XY_POS_TABLE.Geno) =="het");

%% Plot trajectories during the Acclim period.  

xWT = ALL_XY_POS_TABLE.X(allWT);
yWT = ALL_XY_POS_TABLE.Y(allWT);

xvalues = cell2mat(TRACKING_DATA(:,3));
yvalues = cell2mat(TRACKING_DATA(:,4));


xHET = ALL_XY_POS_TABLE.X(allHET);
yHET = ALL_XY_POS_TABLE.Y(allHET);

% boxC = [430 78]; 
boxC = [410 104]; 
radius = 90; 

% WT trajectories 
for j = 1:12
    xvalues = xWT{j,1};
    yvalues = yWT{j,1};
    for p = 1:12798
    x = xvalues(p);
    y = 518 - yvalues(p); 
    x2 = xvalues(p+1);
    y2 = 518 - yvalues(p+1);
    plot([x, x2],[y,y2],'k')
    hold on 
    end  
end 
axis([0 520 0 520])
hold on 
viscircles(boxC, radius, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);

% HET trajectories 
for j = 1:12
    xvalues = xHET{j,1};
    yvalues = yHET{j,1};
    for p = 1:12799
    x = xvalues(p);
    y = 518 - yvalues(p); 
    x2 = xvalues(p+1);
    y2 = 518 - yvalues(p+1);
    plot([x, x2],[y,y2],'r')
    hold on 
    end  
end 
axis([0 520 0 520])
hold on 
viscircles(boxC, radius, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);

%% Per ANIMAL 


allANIMAL = find(string(ALL_XY_POS_TABLE.Animal) =="GN6620");
xANIMAL = ALL_XY_POS_TABLE.X(allANIMAL);
yANIMAL = ALL_XY_POS_TABLE.Y(allANIMAL);

for j = 1:3
    xvalues = xANIMAL{j,1};
    yvalues = yANIMAL{j,1};
    for p = 1:12799
    x = xvalues(p);
    y = 518 - yvalues(p); 
    x2 = xvalues(p+1);
    y2 = 518 - yvalues(p+1);
    plot([x, x2],[y,y2],'r')
    hold on 
    end  
end 
axis([0 520 0 520])
hold on 
viscircles(boxC, radius, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
title('Trajectory during Acclim - GN6620')
savefig(gcf, 'Traj_Acclim_GN6620.fig')
close

%%


xvalues = cell2mat(TRACKING_DATA(:,3)); %xy_table.X; 
yvalues = cell2mat(TRACKING_DATA(:,4)); %xy_table.Y;

n = 12800;
% n = height(xy_table); 
hold on 
for j = 1:n-1  
        x = xvalues(j);
        y = 518 - yvalues(j); 
        x2 = xvalues(j+1);
        y2 = 518 - yvalues(j+1);
        plot([x, x2],[y,y2],'r')
        hold on 
end 
axis([0 520 0 520])
hold on 
viscircles(boxC, radius, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);

title('Acclim - all HET - 191218')
savefig(gcf, 'Acclim_HET_191218.fig')
close



% 
% figure
% for j = 1:20223  
%         x = xy_HET(j,1);
%         y = 518 - xy_HET(j,2); 
%         x2 = xy_HET(j+1,1);
%         y2 = 518 - xy_HET(j+1,2);
%         plot([x, x2],[y,y2],'r')
%         hold on 
% end 
% axis([0 520 0 520])
% % hold on 
% % viscircles(boxC, radius, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
% 
% title('Acclim Trajectories - all HET')
% savefig(gcf, 'AcclimTrajectories_WT.fig')
% close
% xvalues = cell2mat(TRACKING_DATA(:,3));
% yvalues = cell2mat(TRACKING_DATA(:,4));
% 
% xvalues = xy_array(:,1);
% yvalues = xy_array(:,2);
% 
% figure
%   for p = 1:12800
%     x = xvalues(p);
%     y = 518 - yvalues(p); 
%     x2 = xvalues(p+1);
%     y2 = 518 - yvalues(p+1);
%     plot([x, x2],[y,y2],'b')
%     hold on 
%   end 