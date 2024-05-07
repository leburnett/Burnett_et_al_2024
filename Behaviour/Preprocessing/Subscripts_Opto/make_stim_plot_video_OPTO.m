function make_stim_plot_video_OPTO()
%% Making the reconstruction video
%Created by Burnett. Modified 02/09/20

global inpath exp_name 

if exist(strcat('INFO_', exp_name,'.mat'), 'file') 
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    
    opto_rows = find(cell2mat(Info(:,2))==1);
    num_looms = numel(opto_rows);

    baslerpath = (fullfile(inpath, '\basler_*'));
    myBaslerFiles = dir(fullfile(baslerpath,'*.jpg'));
    
    load(strcat('LOGDATA_', exp_name, '.mat'), 'log_data');
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
    
    for j = 1:num_looms

         basler_avi_name = fullfile(inpath, strcat(exp_name, '_0', string(j), '.avi'));
         stimvidobj_Bas = VideoWriter(basler_avi_name, 'Motion JPEG AVI');
         stimvidobj_Bas.FrameRate = 60;
         stimvidobj_Bas.Quality = 100; 
         open(stimvidobj_Bas)
         
        start_loom_row = opto_rows(j,1);
        start_row = start_loom_row - 300;
        end_row = start_loom_row + 480;
        framen = 1; 
        
        for k = start_row:end_row
            
            laser_power = log_data.EC(2);
%             laser_power = log_data.T_pulse(2);
%             laser_power = log_data.FreqPulse(2);
            
            imname = char(Info{k,1});
            a1 = imread(imname); %a1 is uint8;
%             a1 = a1(4:521, 23:540); % now 518 x 518 uint8 - THIS WILL DEPEND ON NEW/OLD SETUP. 
            imshow(a1)
            hold on
            
            x = 5; 
            y = 5 ;
            w = 410; 
            h = 410; 
            
            if framen >300 && framen <361
                rectangle('Position', [x,y,w,h], 'EdgeColor', 'c', 'LineWidth', 6);
                text(340,390, num2str(laser_power), 'FontSize', 20, 'Color','w')
            end
            
            %plot mouse position
%              mouse_x = Info{k,3}; %Uses tracking position found post-processing. 
%              mouse_y = Info{k,4};
%             plot(mouse_x, mouse_y,'w.','MarkerSize',12);
%             
            currF = getframe();
            currF = imresize(uint8(currF.cdata),[416 416]);
            writeVideo(stimvidobj_Bas, currF);
            hold off
            pause(0.01)
            framen = framen +1; 
        end

    %Must close video in order for it to save!
    close(stimvidobj_Bas);
    end 

elseif isempty(index_files)
    fprintf('No OPTO. \n');
end 

close
end 

