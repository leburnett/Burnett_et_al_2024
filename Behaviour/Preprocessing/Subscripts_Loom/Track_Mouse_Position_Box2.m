function Track_Mouse_Position_Box2() %% TRACKING MOUSE - CENTROID LOCATER - FROM JPG
 
% Created by Burnett 11/12/19

 % To be run after 'make_DATA.m', incombination with process_dome_data.m 
 % This function uses a folder of jpg images and tracks the mouse position.
 
 % It outputs a table 'TRACKING_DATA' which contains:
% Col 1: Name of JPG
% Col 2: Major Axis
% Col 3: Centroid x 
% Col 4: Centroid y 
% Col 5: Phi 
% Col 6: Orientation 
 
 global inpath exp_name 

%% Load Basler
baslerpath = (fullfile(inpath, '/basler_*'));
%  baslerpath = 'C:\Data_analysis\DATA\2005_PV_DREADDS\200608\20PV01\02_Loom\basler_15_56_37_200608';
myBaslerFiles = dir(fullfile(baslerpath,'*.jpg'));
[~,index] = sortrows({myBaslerFiles.date}.'); myBaslerFiles = myBaslerFiles(index); clear index %% THIS IS IMPORTANT. SOMETIMES TIME of IMAGE e.e. 'img_903305.jpg' reaches 10 then starts low again. This means images are not in chron order.
num_bas = numel(myBaslerFiles);

    
    %% Make Table
    TRACKING_DATA = cell(num_bas, 6);
    
    %% Add jpg names to table.
    for i = 1:num_bas
        TRACKING_DATA{i,1} = string(myBaslerFiles(i).name);
    end
     
    % imshow(bkg)
    %% For making the TRACKING_DATA table

     se = strel('octagon',6);
     se2 = strel('square',5);
    
    for k = 1:num_bas

        %Read in the video frames.
        imname = string(TRACKING_DATA{k,1});
        imfolder = string(myBaslerFiles(k).folder);
        imnamefull = strcat(imfolder,'\', imname);
        im = imread(imnamefull);
        
        im2 = im; 
        im2(im2>48)=255; 
        im3 = imdilate(im2, se2);
        im3b = imerode(im3, se);
        im4 = imbinarize(im3b);
        im5 = imcomplement(im4); % Need ths to swap b/w 
%         im6 = imclearborder(im5,8);
        CC = bwconncomp(im5);
        stats = regionprops(CC, 'Area','MajorAxisLength','Centroid','MajorAxisLength','Orientation');
        [val, idx] = max([stats.Area]);
        
        if isempty(stats)
            a = 0;
            b = 0;
            Xc = 0;
            Yc = 0;
            phi = 0;
            orient = 0;
            
            if k ==1
                TRACKING_DATA(k,2) = {a};
                TRACKING_DATA(k,3)= {Xc};
                TRACKING_DATA(k,4)= {Yc};
                TRACKING_DATA(k,5)= {phi};
                TRACKING_DATA(k,6) = {orient};
            else
                TRACKING_DATA(k,2) = TRACKING_DATA(k-1,2);
                TRACKING_DATA(k,3)= TRACKING_DATA(k-1,3);
                TRACKING_DATA(k,4)= TRACKING_DATA(k-1,4);
                TRACKING_DATA(k,5)= TRACKING_DATA(k-1,5);
                TRACKING_DATA(k,6) = TRACKING_DATA(k-1,6);
            end
            
        else
            
%                   L = labelmatrix(CC);
%                   x= stats(idx).Centroid(1);
%                   y= stats(idx).Centroid(2);
%                   imshow(im)
%                    hold on
%                   plot(x,y,'r.', 'MarkerSize', 20)
%                   drawnow
%                   pause(0.5)
            
            a = stats(idx).MajorAxisLength; 
            Xc = stats(idx).Centroid(1);
            Yc = stats(idx).Centroid(2);
            phi = deg2rad(-stats(idx).Orientation);
            orient = stats(idx).Orientation;
            
            TRACKING_DATA(k,2) = {a};
            TRACKING_DATA(k,3)= {Xc};
            TRACKING_DATA(k,4)= {Yc};
            TRACKING_DATA(k,5)= {phi};
            TRACKING_DATA(k,6) = {orient};
        end
        
    end
    
    save(strcat('TRACKINGDATA_', exp_name, '.mat'),'TRACKING_DATA');
    
clearvars
end 
