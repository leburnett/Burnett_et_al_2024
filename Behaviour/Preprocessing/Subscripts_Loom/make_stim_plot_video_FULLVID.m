function make_stim_plot_video_FULLVID()
%% Making the reconstruction video
%Created by Burnett. Modified 26/04/22 for multiloom experiments.
% Read in the entire video - not just the periods where the loom happens. 

global inpath exp_name 

index_files = dir('ALL_LOOM_ROWS*');
if length(index_files) ==1 
%     load(index_files.name, 'ALL_LOOM_ROWS');
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
%     load(strcat('XY_array_', exp_name, '.mat'), 'xy_array');
    
    baslerpath = (fullfile(inpath, '\basler_*'));
    myBaslerFiles = dir(fullfile(baslerpath,'*.jpg'));
    num_images = length(myBaslerFiles);
    
    basler_avi_name = fullfile(inpath, strcat(exp_name, '_LOOM.avi'));
    stimvidobj_Bas = VideoWriter(basler_avi_name, 'Motion JPEG AVI');
    stimvidobj_Bas.FrameRate = 60;
    stimvidobj_Bas.Quality = 100; 
    open(stimvidobj_Bas)
    start_row = 2; % first frame is snapshot before recording
    end_row = num_images-1;
    
    for k = start_row:end_row
        
        dot_colour = Info{k, 10};
        if dot_colour == 0
            dot_colour = [0,0,0];
        elseif dot_colour == [0, 40, 80]
            dot_colour = [0.8 0.8 0.8];
        elseif dot_colour == [0, 35, 70]
            dot_colour = [0.7 0.7 0.7];
        elseif dot_colour == [0, 25, 50]
            dot_colour = [0.5 0.5 0.5];
        end
        
        imname = char(Info{k,1});
        a1 = imread(imname); %a1 is uint8;
        imshow(a1)
        hold on
        
        %plot dot position - loom should be in centre of image.
        dot_x = 200; % was 259
        dot_y = 208;
        
        if Info{k,5} ==1 %Looming ON
            radius = Info{k,4};
            if radius == 0
                radius = 0.01;
            end
            plot(dot_x, dot_y, 'Marker', 'o', 'MarkerSize', radius, 'MarkerFaceColor', dot_colour, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
        
        %plot mouse position
        %             mouse_x = xy_array(k,3); %Uses tracking position found post-processing.
        %             mouse_y = xy_array(k,4);
        %             plot(mouse_x, mouse_y,'w.','MarkerSize',12);
        %
        currF = getframe();
        currF = imresize(uint8(currF.cdata),[416 416]);
        writeVideo(stimvidobj_Bas, currF);
        hold off
        pause(0.01)
    end
    
    %Must close video in order for it to save!
    close(stimvidobj_Bas);

elseif isempty(index_files)
    fprintf('No Looms. \n');
end 

close
end 

