%% Calculating HEAD ANGLE. 
function make_dlc_analysis()

% Using DLC data. 

% # 1 - CLEAN THE DLC DATA with 'clean_dlc_csv.m' 
% This creates the file "DLC_POS...." which contains 'dlc_table'. 
% Also creates the file "DLC_SPEED...." which contains 'dlc_speed'. 

% # 2 - CREATE 'CROPPPED' VERSIONS OF THIS TABLE FOR THE INDIVIDUAL LOOMS. 
%  'make_DLC_LOOM_array.m' 
%   Requires full length vidoes to be DLC processed. 

% # 3 - ANALYSE THESE ARRAYS. 

% For the videos which only contain the LOOM time - they are 900 rows (15s). 
% Loom happens at row 300. (5s) and it records the position of the animal
% for 10s after. 

% Centroid position.
% BACK LEFT PAW - - - - (Col 10,11)
% BACK RIGHT PAW - - - - (Col 13,14)
% NOSE - - - - - (Col 1,2)
% TailBase- - - - (Col 16,17)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% THESE THINGS WILL CHANGE:

% Coordinates of centre of shelter - COHORT 1. 
% shelter = [430, 440]; 
% loom_pos = [256 256]; % Centre of image. 
% image_size = 512; 
% shelter_size = 90; 


% NEW SETUP 
shelter = [330, 325]; 
loom_pos = [200 200]; % Centre of image. 
image_size = 400; 
shelter_size = 60; 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Should be in folder with 'DLC_LOOM_POS....'
% In other folder should be 'DLC_LOOM_SPEED....' 

% #1 - Read in DLC data.  
DLC_tables = dir('DLC_LOOM*');
n_files = length(DLC_tables);

%% Make Table

sz = [n_files, 19]; 
varNames = {'Date', 'Animal', 'Exp', 'Loom', 'AngAtLoom',  'LoomX', 'LoomY', 'OrientToShelter', 'LatToOrient2Shelter', 'AngVel', 'OrientToLoom','LatToOrient2Loom', 'AngVelLoom', 'Return2Shelter', 'row_in_shelter', 'dshelt_loom', 'time_loom_2_shelter', 'speed_2_shelter', 'Geno'}; 
varTypes = {'string','string', 'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double','string'};
dlc_analysis = table('Size', sz, 'VariableNames', varNames, 'VariableTypes', varTypes);

for j = 1:n_files 
    dlc_name = DLC_tables(j).name;
    
    date = dlc_name(14:19);
    animal = dlc_name(21:26);
    exp = dlc_name(28:34);
    loom = dlc_name(36:37);
    
    % Load the file.
    load(dlc_name, 'dlc_pos_array');
    dlc_tbl = dlc_pos_array; 
    
    %Create 'dlc' which is only the points we will want for the time being.
    %  'dlc' has 8 COLUMNS.
    dlc = dlc_tbl(:, [1,2,10,11,13,14,16,17]);
    
    
    %% 1 - Distance between NOSE and Shelter Centre
    % Column 9 will be the distance between the nose and the centre of the shelter.
    dshelt = zeros(900,1);
    for i = 1:900
        dshelt(i,1) = pdist([dlc{i,1}, dlc{i,2}; shelter(1), shelter(2)])*32/image_size;
    end
    
    dlc.Nose_Shelter = dshelt;
    
%     plot(table2array(dlc(:,9)), 'r')
%     hold on
%     plot([300 300], [0 30], 'k:')
    
    
    %% 2 - Angle between (Mid of Left and Right Paws)-NOSE and Shelter.
    % in 360 degrees.
     array = zeros(900,1);
     centroid = zeros(900,1);
     
    for i = 1:900
        
        %angle between nose and centre of the shelter. 
        centre =  ([dlc{i,3}, dlc{i,5}]+[dlc{i,4}, dlc{i,6}])/2; %B - centre of the back limbs.
        AB = shelter-centre;
        CB = [dlc{i,1}, dlc{i,2}]-centre;
        ang = atan2d((det([AB;CB])),dot(AB,CB)); %angle in degrees. Can use 'atan2d' or 'atan2 + *180/pi'.
        array(i,1) = ang;
        centroid(i,1) = centre(1);
        centroid(i,2) = centre(2);
        
        % Angle between nose and centre of loom. 
        DB = loom_pos - centre; 
        ang_loom = atan2d((det([DB;CB])),dot(DB,CB)); 
        array(i,2) = ang_loom;
        
        % ang = atan2((det([AB;CB])),dot(AB,CB))*180/pi; %angle in degrees.
        % ang_360 = wrapTo360(ang); %converts +/- degrees to 360.
    end
    
    % smooth data - 1s linear regression. 
    array2 = smoothdata(array(:,1),'loess', 6); %ang from shelter
    array3 = smoothdata(array(:,2),'loess', 6); % ang from loom
    
    %Add columns of orientation towards shelter and loom into 'dlc'
    dlc.AngShelter = array2;
    dlc.AngLoom = array3;
    
    %% 3 - Speed and Distances. 

    % SPEED of mouse
    for i = 2:900
        centroid(i,4) = pdist([centroid(i,1), centroid(i,2);centroid(i-1,1), centroid(i-1,2)]);
    end 
    % Add speed of point between the back legs of the animal to 'dlc'
    centroid(:,4) = smoothdata(centroid(:,3), 'rlowess', 5); 
    dlc.Speed = centroid(:,4); 
    
    acc = zeros(900,1);
    % ACC of mouse 
    for i = 2:900
        acc(i, 1) = dlc.Speed(i)- dlc.Speed(i-1);
    end
%     acc(:,2) = smoothdata(acc(:,1), 'rlowess', 4);
    acc(:,2) = smooth(acc(:,1), 4);
    dlc.Acc = acc(:,2); 
    
    
    
    % Distance of centroid (between back of legs) from edge of shelter 
    size_ratio = 32/image_size; 
    dbox_thresh = shelter_size*size_ratio; 
    
    for i = 1:900
        D = pdist([centroid(i,1), centroid(i,2); shelter(1), shelter(2)])*size_ratio;
        centroid(i,5) = D-dbox_thresh; 
    end
    
    dshelt2 = centroid(:,5); 
    dlc.DShelt = dshelt2; 
    
%     plot(table2array(dlc(:,10)), 'k')
%     hold on
%     plot([300 300], [0 30], 'k:')
%     

%% 

    % Add 360 degrees column.
    ang_360 = wrapTo360(array2);
    dlc.Ang360 = ang_360;
    
    angLoom_360 = wrapTo360(array3);
    dlc.AngLoom360 = angLoom_360;
    
    
    % Add column of change in angle.
    array2(1,2) = 0;
    
    for i = 2:900
        array2(i,2) = array2(i,1)-array2(i-1,1);
    end
    
    dlc.DeltaDeg = abs(array2(:,2));
    
    % Add column of CUMULATIVE change in 360 angle.
    % ang_360(1,2) = 0;
    %
    % for i = 2:900
    %     ang_360(i,3) = ang_360(i,2)+ang_360(i-1,3);
    % end
    % dlc.CumAngle = array(:,2);
    
%% Make TABLE 

    % Add values to table. 
    dlc_analysis.Date(j) = date;
    dlc_analysis.Animal(j) = animal;
    dlc_analysis.Exp(j) = exp;
    dlc_analysis.Loom(j) = loom;
    
    % Heading Angle wrt SHELTER when loom starts
    ang_loom =  dlc.AngShelter(300);
    dlc_analysis.AngAtLoom(j) = ang_loom; 
    
    ang_loom360 =  dlc.Ang360(300);
    dlc_analysis.AngAtLoom360(j) = ang_loom360; 
   
    
    % Heading angle wrt LOOM when loom starts 
    ang_loom_loom = dlc.AngLoom(300);
    dlc_analysis.AngLoomAtLoom(j) = ang_loom_loom;
    
    ang_loom_loom360 =  dlc.AngLoom360(300);
    dlc_analysis.AngLoomAtLoom360(j) = ang_loom_loom360; 
    
    
     % Add position when loom happened 
           position_loom_x = centroid(300,1); 
           position_loom_y = centroid(300,2);
           
    dlc_analysis.LoomX(j) = position_loom_x;
    dlc_analysis.LoomY(j) = position_loom_y;
    
    
    % Does the mouse orient towards within 10 degrees of the shelter during
    % the 10s after the loom happens? 
        data = dlc.AngShelter(300:end);
        data = abs(data);
        
        vals = find((data < 0.5)); % Find when the mouse is turning to face shelter - angle < 0.5 degrees. 
        
        if isempty(vals)
            vals = find((data < 1));
        end 
        
            if isempty(vals)
                vals = find((data < 2.5));
            end
        
                if isempty(vals)
                    vals = find((data < 5));
                end 
        
                    if isempty(vals)
                        vals = find((data < 10));
                    end 
        
                         if isempty(vals) %If from loom to 10s post-loom never orients towards shelter. 
                            dlc_analysis.OrientToShelter(j) = 0; 
                            dlc_analysis.LatToOrient2Shelter(j) = 0; 
                            dlc_analysis.AngVel(j) = 0; 
%                             ang_values.DegreesTot(j) = 0;  
                         else
                             
                            dlc_analysis.OrientToShelter(j) = 1;
                            
                            rows_to_orient_from_loom_start = vals(1)-1; 
                            row_when_orient_to_shelter = rows_to_orient_from_loom_start+300; % Minus one row. 
                            lat_sec = rows_to_orient_from_loom_start/60; % latency in seconds. 
                            dlc_analysis.LatToOrient2Shelter(j) = lat_sec; %Number of frames until head is oriented to shelter. 
    
                            angleto_shelter = dlc.AngShelter(row_when_orient_to_shelter);  
                            ang_velocity = abs(ang_loom - angleto_shelter)*lat_sec;
                            dlc_analysis.AngVel(j) = ang_velocity; 
    
%                             % Total degrees moved from loom start til heading towards shelter. 
%                             degrees_turned = sum(abs(dlc.DeltaDeg(300:300+val)));
%                             ang_values.DegreesTot(j) = degrees_turned;   
             
                         end 
        
                        
     % Does the mouse orient towards within 10 degrees of the LOOM WHILE ITS HAPPENING?
     
        data2 = dlc.AngLoom(300:525);
        data2 = abs(data2);
        
        vals2 = find((data2 < 0.5)); % Find when the mouse is turning to face shelter - angle < 0.5 degrees. 
        
        if isempty(vals2)
            vals2 = find((data2 < 1));
        end 
        
            if isempty(vals2)
                vals2 = find((data2 < 2.5));
            end
        
                if isempty(vals2)
                    vals2 = find((data2 < 5));
                end 
        
                    if isempty(vals2)
                        vals2 = find((data2 < 10));
                    end 
        
                         if isempty(vals2) %If from loom to 10s post-loom never orients towards shelter. 
                            dlc_analysis.OrientToLoom(j) = 0; 
                            dlc_analysis.LatToOrient2Loom(j) = 0; 
                            dlc_analysis.AngVelLoom(j) = 0; 
%                             ang_values.DegreesTot(j) = 0;  
                         else
                             
                            dlc_analysis.OrientToLoom(j) = 1;
                            
                            rows_to_orient_to_loom = vals2(1)-1; 
                            row_when_orient_to_loom = rows_to_orient_to_loom+300; % Minus one row. 
                            lat_sec2 = rows_to_orient_to_loom/60; % latency in seconds. 
                            dlc_analysis.LatToOrient2Loom(j) = lat_sec2; %Number of frames until head is oriented to shelter. 
    
                            angleto_loom = dlc.AngLoom(row_when_orient_to_loom);  
                            ang_velocity2 = abs(ang_loom - angleto_loom)*lat_sec2;
                            dlc_analysis.AngVelLoom(j) = ang_velocity2; 
    
%                             % Total degrees moved from loom start til heading towards shelter. 
%                             degrees_turned = sum(abs(dlc.DeltaDeg(300:300+val)));
%                             ang_values.DegreesTot(j) = degrees_turned;   
             
                         end 

        % Does the mouse return to the shelter within 10s of loom starting?                
    
        shelter_rows = find(dshelt2(300:end)<0) + 300; % row value when mouse back in shelter.
        dshelt_loom = dlc.DShelt(300);
        
        if isempty(shelter_rows)
            dlc_analysis.Return2Shelter(j) = 0;
            dlc_analysis.row_in_shelter(j) = NaN; 
            dlc_analysis.dshelt_loom(j) = dshelt_loom;
            dlc_analysis.time_loom_2_shelter(j) = NaN;
            dlc_analysis.speed_2_shelter(j) = NaN;
        
        else 
            dlc_analysis.Return2Shelter(j) = 1;
            row_in_shelter = shelter_rows(1);
            
            time_loom_2_shelter = (row_in_shelter-300)/60;
            speed_2_shelter = dshelt_loom/time_loom_2_shelter;
            dlc_analysis.row_in_shelter(j) = row_in_shelter;
            dlc_analysis.dshelt_loom(j) = dshelt_loom;
            dlc_analysis.time_loom_2_shelter(j) = time_loom_2_shelter;
            dlc_analysis.speed_2_shelter(j) = speed_2_shelter;
        end 
        
    
    %beginning of head rotation or beginning of acceleration. 
%  if animal == "MJ0593" || animal == "MJ0595" || animal == "GN3959" || animal == "GN4244" || animal == "GN4473" ||animal == "GN6560" || animal == "GN6562" || animal == "GN6610" || animal == "GN7269" ||animal == "GN7476" || animal == "GN7790" || animal == "GN7614"
%      dlc_analysis.Geno(j) = "HET";
%  else 
%      dlc_analysis.Geno(j) = "WT";
%  end 
 
%   if  animal == "GN2375" || animal == "GN2377" || animal == "GN2382" || animal == "GN2628"  ||  animal == "GN2637"  ||  animal == "GN2754"||  animal == "GN2902"||  animal == "GN2901"||  animal == "GN2900"|| animal == "GN2832" ||  animal == "GN2829" ||  animal == "GN2833"
%               dlc_analysis.Geno(j) = "het";
%    else 
%               dlc_analysis.Geno(j) = "wt";
%    end 

% Ptchd1
 if animal == "GN1385" || animal == "GN1386" || animal == "GN1387" || animal == "GN1388" || animal == "GN1394" ||animal == "GN2593" || animal == "GN2594" || animal == "GN2708" || animal == "GN2709" ||animal == "GN3903" || animal == "GN4369" || animal == "GN4373"
     dlc_analysis.Geno(j) = "HET";
 else 
     dlc_analysis.Geno(j) = "WT";
 end 

 % Save the dlc file with angles etc. 
 save(strcat(date, '_', animal, '_', exp, '_', loom, '_HEADING.mat'), 'dlc')
 
end 

save('221031_DLC_ANALYSIS_Ptchd1.mat', 'dlc_analysis'); 

% save('200828_ANG_TABLE_Setd5_C1_C2.mat', 'ang_values')


% PLOT - heading angle at the time of loom start vs speed to reach shelter.

% figure
% for i = 1:63
%     if ang_values.Geno(i) == "WT" 
%         marker = 'k.'; 
%     else
%         marker = 'r.';
%     end 
%         x = abs(ang_values.AngAtLoom(i)); 
%         y = ang_values.LatToShelterAng(i); 
%         plot(x,y,marker, 'MarkerSize', 25);
%         hold on 
% end 
% title('Distance from shelter vs time to shelter. ')
% ylabel('Speed - cm/s')
% xlabel('Distance - cm')


end 




% if isempty(shelter_rows) %%%%%% IF THE ANIMAL DOES NOT GO BACK TO THE SHELTER. - Either freeze or no response. 
%         
% %             shelter_rows = 0; % row value when mouse back in shelter.
% %             row_shelter = 0; 
%         
%          % Add position when loom happened 
%                 position_loom_x = centroid(300,1); 
%                 position_loom_y = centroid(300,2);
%     ang_values.LoomX(j) = position_loom_x;
%     ang_values.LoomY(j) = position_loom_y;
%     
%                     d_shelter_loom = pdist([position_loom_x, position_loom_y; position_shelter_x, position_shelter_y])*32/image_size; 
%                     d_shelter_loom_nose = dlc.Nose_Shelter(300); 
%     ang_values.d_shelter_loom(j) = d_shelter_loom; 
%     ang_values.d_shelter_loom_nose(j) = d_shelter_loom_nose; 
%     ang_values.t_shelter_loom(j) = 0; 
%     ang_values.speed_loom_shelter(j) = 0; 
% 
%         else 
%             row_shelter = shelter_rows(1); 
%     
%             % Add position when loom happened 
%                 position_loom_x = centroid(300,1); 
%                 position_loom_y = centroid(300,2);
%     ang_values.LoomX(j) = position_loom_x;
%     ang_values.LoomY(j) = position_loom_y;
%    
%                 position_shelter_x = centroid(row_shelter,1);
%                 position_shelter_y = centroid(row_shelter,2);
%     
%                 d_shelter_loom = pdist([position_loom_x, position_loom_y; position_shelter_x, position_shelter_y])*32/image_size; 
% %     d_shelter_loom = pdist([position_loom_x, position_loom_y; shelter(1), shelter(2)])*32/image_size; 
%                 t_shelter_loom = (row_shelter-300)/60;
%                 speed_loom_shelter = d_shelter_loom/t_shelter_loom; 
%                 d_shelter_loom_nose = dlc.Nose_Shelter(300); 
%     
%     ang_values.d_shelter_loom(j) = d_shelter_loom; 
%     ang_values.d_shelter_loom_nose(j) = d_shelter_loom_nose; 
%     ang_values.t_shelter_loom(j) = t_shelter_loom; 
%     ang_values.speed_loom_shelter(j) = speed_loom_shelter; 
%         end 


%% Make future script of 'Gait Analysis' - using information about the paw movements. 
