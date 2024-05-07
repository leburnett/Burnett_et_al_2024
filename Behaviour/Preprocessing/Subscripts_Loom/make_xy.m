function make_xy()

% Makes array 'xy_array' for ALL videos regardless of loom/ no loom. 

    % "xy_array"
    % Col 1 = x position of mouse
    % Col 2 = y position of mouse
    % Col 3- x values - outliers filled
    % Col 4 - y values - outliers filled
    % Col 5 - SPEED -  distance between frames (cm/s)
    % Col 6 - ACC - change in speed between frames (cm/s/s))
    % Col 7 - Distance of mouse from the centre of the arena. (cm)
    % Col 8 - Distance from edge of shelter (cm)
    
global exp_name track_req new_setup output_folder image_size

if track_req == 1
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    num_bas = numel(Info(:,1));
elseif track_req == 0
    load(strcat('INFOTABLE_', exp_name, '.mat'), 'Info_Table');
    num_bas = height(Info_Table);
end 

if new_setup ==0 %Image size is 518 x 518 
%     xy_array = zeros(num_bas, 8);
%     boxC = [image_size/2, image_size/2];  %- centre of the arena
%     shelterC = [430, 440]; % centre of the shelter
%     box_size_cm = 32; 
% %     image_size = 518; 
%     sz_ratio = box_size_cm/image_size; %cm/pixels
%     pix_ratio = image_size/box_size_cm; % pixels/cm
%     fps = 60; 
%     shelter_size = 90; % pixels - roughly 6cm.
    
    %COHORT 2
    
    xy_array = zeros(num_bas, 8);
    boxC = [262,262]; %[258, 258];  %- centre of the arena
    shelterC = [430, 430]; % centre of the shelter
    box_size_cm = 32; 
%     image_size = 524;%518; 
    sz_ratio = box_size_cm/image_size; 
    pix_ratio = image_size/box_size_cm; % pixels/cm
    fps = 60; 
    shelter_size = 90; % pixels - roughly 6cm.
    
elseif new_setup ==1  % Image size is 400 x 400 
    xy_array = zeros(num_bas, 8);
    boxC = [image_size/2, image_size/2];  %- centre of the arena
    shelterC = [330, 325]; % centre of the shelter
    box_size_cm = 32; 
%     image_size = 416; 
    sz_ratio = box_size_cm/image_size; %cm/pixels. 
    pix_ratio = image_size/box_size_cm; % pixels/cm
    fps = 60; 
    shelter_size = 60;
end 

   % NEW BOX - TEST BOUNDARIES 
% imname = 'img_188188385780835.jpg';
% im = imread(imname);
% imshow(im); 
%  hold on 
%  viscircles(boxC, 60, 'Color', 'r')
% viscircles(boxC, 150, 'Color', 'b')
%  viscircles(shelterC, 60, 'Color', 'g')

% OLD BOX
% imname = 'img_3277485698605659.jpg';
% im = imread(imname);
% im = im(4:521, 23:540);
% imshow(im); 
% hold on 
% viscircles(boxC, 90, 'Color', 'r')
% viscircles(boxC, 190, 'Color', 'b')
% viscircles(shelterC, 90, 'Color', 'm')



%% Add xmouse, y mouse from Info/Info_T
% 
if track_req == 1
    
    if new_setup ==1 
         for i = 1:num_bas
             if i ==1
                 xy_array(i,1) = Info{2,15}; %Now with the new setup the first frame is not in time with the rest. 
                 xy_array(i,2) = Info{2,16};
             else
                 xy_array(i,1) = Info{i,15};
                 xy_array(i,2) = Info{i,16};
             end
         end 
    elseif new_setup == 0 
        
%         for i = 1:num_bas
%             xy_array(i,1) = Info{i,13}; 
%             xy_array(i,2) = Info{i,14};
%         end
        for i = 1:num_bas
            xy_array(i,1) = Info{i,15}; 
            xy_array(i,2) = Info{i,16};
        end

    end 
   
elseif track_req == 0
     for i = 1:num_bas
        xy_array(i,1) = Info_Table.Mouse_x(i); 
        xy_array(i,2) = Info_Table.Mouse_y(i);
     end 
end 
    
%%  Clean up x,y values. REmove outliers. 

% 1 - Find points that are OUTSIDE of the box. 
    % IF X  > 460 
    % or X < 48 
    % or Y > 480
    % or Y < 38

%     points_out_box = find(xy_array(:,1)< 35 | xy_array(:,1)> 480 | xy_array(:,2)< 35 | xy_array(:,2)> 480); 
%     if ~isempty(points_out_box)
%         fprintf('Number of points out of range: \n')
%         disp(numel(points_out_box))
%     end 

    % Find when new blocks appear. 
%     for i = 2:numel(points_out_box)
%         if points_out_box(i,1)-points_out_box(i-1,1) == 1
%             points_out_box(i,2) = 0;
%         else
%             points_out_box(i,2) = 1;
%         end 
%     end 
%     
%     blocks = find(points_out_box(:,2)~=0); 
%     
%     for j = 1:numel(blocks)+1
%         if j ==1 
%             row = points_out_box(1,1); 
%         else 
%             points_row = blocks(j-1); 
%             row = points_out_box(points_row,1); 
%         end 
%         row_before = row-1; 
%         row_after = points_out_box((blocks(j)-1),1)+1;
%         row_diff = row_after-row_before;
%              
%             %Linearly interpolate between the values of the rows found. 
%            for s = [1,2] %For col 1 and 2 (x,y)
%                 value1 = xy_array(row_before,s);
%                 value2 = xy_array(row_after,s);
%                 vals = linspace(value1, value2, row_diff+1);
%                 row_vals = linspace(row_before, row_after, row_diff+1);
% 
%                 for p = 1:row_diff +1
%                     row_to_change = row_vals(p);
%                     val = vals(p);
%                     xy_array(row_to_change,s)=val; 
%                 end 
%            end
%     end 
            
    
    %%
    d_array = [];
    
    for k = 2:num_bas
        A = pdist([xy_array(k-1,1), xy_array(k-1,2); xy_array((k),1), xy_array((k),2)]);
        d_array(k,1) = A*sz_ratio*fps; % (cm/s)
        % If wanting to categorise each row add the lines below
%         d_array(k,2)= A; 
%         if A>50
%          d_array(k,3)= 0; 
%         else
%          d_array(k,3)= 1;
%         end 
    end

    % FIND values of speed >160cm/s
    speed_threshold = 160; %cm/s 
    over_thr = find(d_array>speed_threshold);
    n_over_thr = numel(over_thr); 
    pixels_per_frame = (speed_threshold/fps)*pix_ratio; 
 
    % Change the first value in over_250 then recalculate the distance.
    % Keep doing this until all values in over_250 have been changed. 
    while ~isempty(over_thr)
        
        row = over_thr(1); % This is the row with the large leap! 
        row_before = over_thr(1)-1; % This row is assumed to be ok! 
        
        if  pdist([xy_array(row+1,1), xy_array(row+1,2); xy_array(row_before,1), xy_array(row_before,2)])< pixels_per_frame
            row_after = row+1; 
            xy_array(row,1) = mean([xy_array(row_before,1), xy_array(row_after,1)]);
            xy_array(row,2) = mean([xy_array(row_before,2), xy_array(row_after,2)]);
        else
            row_before = row-5;
            row_after = row+5;
            row_diff = row_after - row_before;
            
            for s = [1,2] %For col 1 and 2 (x,y)
                value1 = xy_array(row_before,s);
                value2 = xy_array(row_after,s);
                vals = linspace(value1, value2, row_diff+1);
                row_vals = linspace(row_before, row_after, row_diff+1);
                
                for p = 1:row_diff +1
                    row_to_change = row_vals(p);
                    val = vals(p);
                    xy_array(row_to_change,s)=val;
                end
            end 
        end
        
        d_array = [];
        for k = 2:num_bas
            A = pdist([xy_array(k-1,1), xy_array(k-1,2); xy_array((k),1), xy_array((k),2)]);
            d_array(k,1) = A*sz_ratio*fps; %Speed(cm/s)
        end
        
        over_thr = find(d_array>160);
        n_over_thr = numel(over_thr);
        
%         if new_setup == 0 
% %         if xy_array(row, 1)> 420 && xy_array(row,2)<90 || xy_array(row-1, 1)> 420 && xy_array(row-1,2)<90
% %             disp(exp_name)
% %             disp(row)
% % %             break
% %             error('Mouse hidden by leg check row')
% %         end 
%         end 
        
    end 
   
% If only a few rows to change and nthey are quite complex:
%
% 
%  row_before = 3975;
%  row_after = 3980; 
%  row_diff = row_after - row_before; 
% 
%    for s = [1,2] %For col 1 and 2 (x,y)
%                 value1 = xy_array(row_before,s);
%                 value2 = xy_array(row_after,s);
%                 vals = linspace(value1, value2, row_diff+1);
%                 row_vals = linspace(row_before, row_after, row_diff+1);
% 
%                 for p = 1:row_diff +1
%                     row_to_change = row_vals(p);
%                     val = vals(p);
%                     xy_array(row_to_change,s)=val; 
%                 end 
%    end
% % % %    
   %%
% % % % % % % 
% % % % 
% for i = 1:10253
%     if xy_array(i,1)<130
%         xy_array(i,1) = xy_array(i-1,1);
%         xy_array(i,2) = xy_array(i-1,2);
%     end 
% end 
    
    
% % % % % % % 
%     else 
%             current_row  = row+1; 
%             %This while loop will keep adding to the current row until it
%             %finds one which is less than the required distance from the
%             %last 'good' frame. 
%        
%             while pdist([xy_array(current_row,1), xy_array(current_row,2); xy_array(row_before,1), xy_array(row_before,2)])> 35
%                 current_row = current_row +1; 
%                 if current_row > over_160(2)
%                     break
%                 end 
%             end 
%           
%             
%             row_after = current_row;
%             row_diff = row_after-row_before;
%              
%             %Linearly interpolate between the values of the rows found. 
%              for s = [1,2] %For col 1 and 2 (x,y)
%                 value1 = xy_array(row_before,s);
%                 value2 = xy_array(row_after,s);
%                 vals = linspace(value1, value2, row_diff+1);
%                 row_vals = linspace(row_before, row_after, row_diff+1);
% 
%                 for p = 1:row_diff +1
%                     row_to_change = row_vals(p);
%                     val = vals(p);
%                     xy_array(row_to_change,s)=val; 
%                 end 
%              end
%         end 


% OLD CODE

%      xy_array(row,1) = mean([xy_array(row_before,1), xy_array(row_after,1)]);
%             xy_array(row,2) = mean([xy_array(row_before,2), xy_array(row_after,2)]);
                   
%         row_after = over_250(1)+1;
% 
%         if size(over_250)>=2
%             if row_after == over_250(2) && pdist([xy_array(over_250(2),1), xy_array(over_250(2),2); xy_array(over_250(1)-1,1), xy_array(over_250(1)-1,2)])> 5 && over_250(3) == over_250(2)+1 
%                 row_after = over_250(3);
%             end
%         end 

%     
%% Add updated x,y values to xy_array. 
    
    xy_array(:,3) = filloutliers(xy_array(:,1), 'linear', 'movmean', 15);
    xy_array(:,4) = filloutliers(xy_array(:,2), 'linear', 'movmean', 15);

%% Update next columns with speed, acc and distance from edge of box. 

    % Distance between frames - eg speed. 
    for k = 2:num_bas
        A = pdist([xy_array(k-1,3), xy_array(k-1,4); xy_array((k),3), xy_array((k),4)]);
        xy_array(k,5) = A*sz_ratio*fps; %Distance between frames - in cm!! 
    end
    
    %% UPDATE SPEED - make it a moving mean over 83ms window. 
    
     xy_array(:,5) = movmean(xy_array(:,5), 5);
    
     %%
    
     % Difference in speed between frames - eg ACC. 
    for k = 2:num_bas
        B = xy_array(k,5)- xy_array(k-1,5);
        xy_array(k,6) = B*sz_ratio*fps; % cm/s^2
    end
    
    % Distance from the centre of the box (arena). 
    for i = 1:num_bas 
        C = pdist([xy_array(i,3), xy_array(i,4); boxC(1),boxC(2)]);
        xy_array(i,7)= C*sz_ratio; %distance from centre of the box. - cm.  
    end
  
    % Distance from the EDGE of the shelter. 
    dbox_thresh = shelter_size*sz_ratio; % 90pixels from centre ~6cm.
     for i = 1:num_bas 
        D = pdist([xy_array(i,3), xy_array(i,4); shelterC(1),shelterC(2)]);
        xy_array(i,8)= (D*sz_ratio)-dbox_thresh;  
     end

    % pdist works by finding the distance between PAIRS of coords. The
    % first variable is a matrix where each ROW is an X,Y coordinate and it
    % compares between them. 
    % [(x1,y1);(x2,y2)]; 

%%  Make Table
    XRaw = xy_array(:,1);
    YRaw = xy_array(:,2);
    X = xy_array(:,3);
    Y = xy_array(:,4);
    Speed = xy_array(:,5);
    Acc = xy_array(:,6);
    DistCentre = xy_array(:,7);
    DistShelter = xy_array(:,8);
 
    xy_table = table(XRaw, YRaw,X, Y, Speed, Acc, DistCentre, DistShelter); 
    
%% Save table and array. 
    save(strcat('XY_table_', exp_name, '.mat'), 'xy_table');
    save(strcat('XY_array_', exp_name, '.mat'), 'xy_array');
    
    % Save in grouped folder. 
      % Save in 'Exp Project' folder. 
    xy_folder = strcat(output_folder, 'SUMMARIES\ACTIVITY\XY_ARRAYS\');
    
    if ~exist(xy_folder,'dir')
        mkdir(xy_folder)
    end

    save(fullfile(strcat(xy_folder,'XY_ARRAY_', exp_name,'.mat')), 'xy_array'); 
    
    % 
    xy_table_folder = strcat(output_folder, 'SUMMARIES\ACTIVITY\XY_ARRAYS\');
    
    if ~exist(xy_table_folder,'dir')
        mkdir(xy_table_folder)
    end
    
    save(fullfile(strcat(xy_table_folder,'XY_TABLE_', exp_name,'.mat')), 'xy_table');    

    clearvars
    
% exp_name = '191216_GN4472_01_Loom'
end 

   % Cohort2 
  %   save(fullfile(strcat('C:\Data_analysis\DATA\Setd5-Cohort2\SUMMARIES\ACTIVITY\xy_arrays\XY_table_', exp_name, '.mat')), 'xy_table');

  % Cohort 1 - No shelter
  %     save(fullfile(strcat('C:\Data_analysis\DATA\1912_Setd5\ALL_DAY_SUMMARIES\ACTIVITY\xy_array\NoShelter\XY_table_', exp_name, '.mat')), 'xy_table'); 
   
  % Cohort 1 
  %   save(fullfile(strcat('C:\Data_analysis\DATA\Setd5\ALL_DAY_SUMMARIES\ACTIVITY\xy_array\XY_table_', exp_name, '.mat')), 'xy_table');
   
%     % Cohort 2 - NoShelter
%     save(fullfile(strcat('C:\Data_analysis\DATA\BOTH_SETD5_COHORTS_\NoShelter\XY_ARRAYS\XY_table_', exp_name, '.mat')), 'xy_table'); 

     % SERT
%     save(fullfile(strcat('C:\Data_analysis\DATA\2003_SERT\SUMMARIES\XY_ARRAYS\XY_table_', exp_name, '.mat')), 'xy_table'); 
    