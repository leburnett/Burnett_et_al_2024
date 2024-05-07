function Track_Mouse_OPTO() %% TRACKING MOUSE - CENTROID LOCATER - FROM JPG
 % Created by Burnett 14/10/20

 % To be run after 'make_DATA_OPTO.m', incombination with process_dome_data.m 
 % This function uses a folder of jpg images and tracks the mouse position.

 % It loads 'Info' and adds the following columns.
 
% Col 1: Image Name 
% Col 2: RF  
% Col 3: Centroid x 
% Col 4: Centroid y 
% Col 5: Major Axis
% Col 6: Phi 
% Col 7: Orientation 
 
 global inpath exp_name  

 
if exist(strcat('INFO_', exp_name, '.mat'), 'file')
    load(strcat('INFO_', exp_name, '.mat'),'Info')

num_rows = length(Info); 
    
    %% Read in first image.
    
    %%
    % if new_setup == 0
        % bkg = imread('C:\Data_analysis\Core_Scripts\Loom_Analysis_Pipeline\JPGS\BOX2b.jpg');
        % bkg = rgb2gray(bkg);
        % bkg = bkg(4:521, 23:540); %CROP
    % else
        % bkg = imread('C:\Data_analysis\Core_Scripts\Loom_Analysis_Pipeline\JPGS\BOX6.jpg');
        % bkg = rgb2gray(bkg);
        %% bkg = bkg(5:404, 10:409);
        % bkg = bkg(10:409, 1:400);
    % end
    
    % imshow(bkg)
    %Set BaslerPath.
    baslerpath = fullfile(inpath, '/basler_*');
    myBaslerFiles = dir(fullfile(baslerpath,'*.jpg'));
%     num_bas = length(myBaslerFiles);
    
    [~,index] = sortrows({myBaslerFiles.date}.'); 
    myBaslerFiles = myBaslerFiles(index); 
    clear index 
    
    se = strel('octagon',6);
    se2 = strel('square',5);
     
    for k = 1:num_rows

        %Read in the video frames.
        imname = string(Info{k,1});
        imfolder = string(myBaslerFiles(k).folder);
        imnamefull = strcat(imfolder,'\', imname);
        im = imread(imnamefull);
        
%         if new_setup ==0
%             im = im(4:521, 23:540); %CROP %%%%%%%%%%%%%%%% For NEW SETUP COMMENT OUT CROPPING.
%         else
%             im = im(10:409, 1:400); 
%             im = imresize(im, [400 400]);
%         end

        % Heavily threshold the mouse. 
        im2 = im; 
        im2(im2>75)=255; 
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
                Info(k,3) = {Xc};
                Info(k,4)= {Yc};
                Info(k,5)= {a};
                Info(k,6)= {phi};
                Info(k,7) = {orient};
            else
                Info(k,3) = Info(k-1,3);
                Info(k,4)= Info(k-1,4);
                Info(k,5)= Info(k-1,5);
                Info(k,6)= Info(k-1,6);
                Info(k,7) = Info(k-1,7);
            end
            
        else
            
%                   x= stats(idx).Centroid(1);
%                   y= stats(idx).Centroid(2);
%                   imshow(im4)
%                   hold on
%                   plot(x,y,'r.', 'MarkerSize', 20)
%                   drawnow
%                   pause(0.1)
            
            a = stats(idx).MajorAxisLength; %Used to be /2 - not sure why....
            Xc = stats(idx).Centroid(1);
            Yc = stats(idx).Centroid(2);
            phi = deg2rad(-stats(idx).Orientation);
            orient = stats(idx).Orientation;
            
            Info(k,3) = {Xc};
            Info(k,4)= {Yc};
            Info(k,5)= {a};
            Info(k,6)= {phi};
            Info(k,7) = {orient};
        end
        
    end
    
    save(strcat('Info_', exp_name, '.mat'),'Info');
    
clearvars

end

end 

%% For visualising the tracking
% 
% 
% se = strel('octagon',9);
% se2 = offsetstrel('ball',40, 20);
% centre = [260 254];
% radius = 232;
% % 
% % gif('track_darkbkg.gif');
% % 
% for k = 700:1200 %num_bas
%     
%     %Read in the video frames. 
%     imname = string(TRACKING_DATA{k,1});
%     im = imread(imname);
%     
%     % Cropping image to circle of dome floor.  
%     im2 = uint8(zeros(imageSize));
%     im2 = im.*mask; %im2 is uint8 (512 x 512)
% 
%     %Subtract the background from the image. 
%     im3 = im2 - bkg_sub;
% 
%     %Remove pixels with value <50. Removes noise. 
%     im3(find(im3<60))=0; 
%     im4 = imadjust(im3, [0.1, 0.35], []);
%         
%     %Remove the tail 
%     CC = bwconncomp(im4);
%     stats = regionprops(CC, 'MajorAxisLength', 'Area');
%     L = labelmatrix(CC);
%     mask2 = ismember(L, find([stats.MajorAxisLength] <=40)); 
%     mask3 = ismember(L,find([stats.MajorAxisLength]>10));
%     mask4 = ismember(L,find([stats.Area]>50));
%     mask5 = mask2.*mask3.*mask4;
%     im4b = double(im4);
%     v4 = im4b.*mask5; %V4 now has the tail removed. 
%     v4 = uint8(v4);
% 
%     %Then dilate and close the image. 
%     v5 = imdilate(v4, se);
%     v5 = imclose(v5, se2);
%     %v5 = imadjust(v5, [0.1, 0.25], []);
% 
%     %Remove the small blobs that appear due to the tail. 
%     CC2 = bwconncomp(v5);
%     stats2 = regionprops(CC2,'Area', 'Centroid'); %creates array with 'stats' values about each component.
% %     L2 = labelmatrix(CC2);
% %     mask3 = ismember(L2, find([stats2.Area] >=3000)); 
% %     v6 = double(v5);
% %     v6 = v6.*mask3;
% %     v6 = uint8(v6);
%   
%      imshow(v5);
%      hold on 
%      
% %      %Draw outline of the mouse
%      B = bwboundaries(v5);
%      visboundaries(B)
%      
%     % Plot perimeter of dome.
%     viscircles(centre, radius, 'Color', 'c', 'LineWidth', 0.5);
%    
%  [val, idx]=max([stats2.Area]);
%  
%     %Draw elipse and plot centroid. 
%     t = linspace(0,2*pi,50);
%     a = stats2(idx).MajorAxisLength/2;
%     b = stats2(idx).MinorAxisLength/2;
%     Xc = stats2(idx).Centroid(1);
%     Yc = stats2(idx).Centroid(2);
%     phi = deg2rad(-stats2(idx).Orientation);
%     x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%     y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%     
% % Plot ellipse
% %     plot(x,y,'y','Linewidth',2)
% 
% %Plot Centroid 
%     plot(Xc,Yc, 'r.', 'MarkerSize', 15);
%     drawnow
% %     gif
%     pause(0.1);
% %   
%  end 

% 


