function make_stim_plot_video_LOOM()
%% Making the reconstruction video
%Created by Burnett. Modified 02/09/20

global inpath exp_name ILI

index_files = dir('ALL_LOOM_ROWS*');
if length(index_files) ==1 
    load(index_files.name, 'ALL_LOOM_ROWS');
    
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    load(strcat('XY_array_', exp_name, '.mat'), 'xy_array');
    
    baslerpath = (fullfile(inpath, '\basler_*'));
    myBaslerFiles = dir(fullfile(baslerpath,'*.jpg'));
    num_looms = numel(ALL_LOOM_ROWS(1,:));
    
%     if ILI ==1 
%         frames_before = 299;
%         frames_after = 900; 
%     elseif ILI == 0  
        frames_before = 299;
        frames_after = 600; 
%     end 
    
    
    for j = 1:num_looms

         basler_avi_name = fullfile(inpath, strcat(exp_name, '_0', string(j), '.avi'));
         stimvidobj_Bas = VideoWriter(basler_avi_name, 'Motion JPEG AVI');
         stimvidobj_Bas.FrameRate = 60;
         stimvidobj_Bas.Quality = 100; 
         open(stimvidobj_Bas)
         
        start_loom_row = ALL_LOOM_ROWS(1,j)-5;
        start_row = start_loom_row-frames_before;
        end_row = start_loom_row + frames_after;
       
        for k = start_row:end_row
            
            dot_colour = Info{k, 10}; 
%             if dot_colour == 0 
%                 dot_colour = [0,0,0];
%             elseif dot_colour == [0, 40, 80]
%                 dot_colour = [0.8 0.8 0.8]; 
%             elseif dot_colour == [0, 35, 70]
%                 dot_colour = [0.7 0.7 0.7];
%             elseif dot_colour == [0, 25, 50]
%                 dot_colour = [0.5 0.5 0.5];
%             end 
            
            if dot_colour == 0 | dot_colour == [0,0,0]
                dot_colour = [0,0,0];
            elseif dot_colour == [0, 40, 80] | dot_colour == [0,37,75]
                dot_colour = [0.8 0.8 0.8]; 
            elseif dot_colour == [0, 35, 70] | dot_colour == [0 25,50]
                dot_colour = [0.7 0.7 0.7];
            elseif dot_colour == [0, 25, 50] | dot_colour == [0,12,25]
                dot_colour = [0.5 0.5 0.5];
            end 
            
            imname = char(Info{k,1});
            a1 = imread(imname); %a1 is uint8;
%             a1 = a1(4:521, 23:540); % now 518 x 518 uint8 - THIS WILL DEPEND ON NEW/OLD SETUP. 
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
    end 

elseif isempty(index_files)
    fprintf('No Looms. \n');
end 

close 

end 

