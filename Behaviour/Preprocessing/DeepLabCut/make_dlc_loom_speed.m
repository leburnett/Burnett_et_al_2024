function make_dlc_loom_speed()

% Makes table: dlc_loom_speed 
% Each column of the table corresponds to the speed of one limb in the time
% 3s before the loom til 10s after the loom. 

global exp_name output_folder image_size

if exist(strcat('ALL_LOOM_ROWS_', exp_name,'.mat'), 'file') 

    load((strcat('ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS'); 
    load(strcat('DLC_',exp_name, '.mat'), 'dlc_tbl');

    loom_rows = ALL_LOOM_ROWS(1,:)'; 
    num_looms = numel(loom_rows);

    num = 780; 
    frames_before = 179; 
    frames_after = 600; 
       
 %  Make folder in summaries for dlc-speed.  
    loom_sp_folder = strcat(output_folder, '\SUMMARIES\DLC\Loom_Speed\');

    if ~exist(loom_sp_folder,'dir')
        mkdir(loom_sp_folder)
    end
    
  for i = 1:num_looms
        
         dlc_loom_speed = table('Size', [num, 8], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'Nose', 'FL', 'FR', 'BL', 'BR', 'TB', 'TM', 'TT'});  
      
         row = loom_rows(i,1) - frames_before; % Start 3s (180 frames) before the loom stimulus starts 
         row3 = loom_rows(i,1) + frames_after;  % end 10s after loom starts. 
         
         if row<0 
             row = 1;
         end 
         
         if row3 > height(dlc_tbl)
            row3 = height(dlc_tbl);
         end 
         
         
         for q = 1:8
             if q ==1
                 % Nose
                 for j = row:row3
                     dlc_loom_speed.Nose(j-(row-1)) = pdist([dlc_tbl.N_x(j), dlc_tbl.N_y(j) ; dlc_tbl.N_x(j-1), dlc_tbl.N_y(j-1)])*(32/image_size)*60;
                 end
             elseif q ==2
                 % FL
                 for j = row:row3
                     dlc_loom_speed.FL(j-(row-1)) = pdist([dlc_tbl.FL_x(j), dlc_tbl.FL_y(j) ; dlc_tbl.FL_x(j-1), dlc_tbl.FL_y(j-1)])*(32/image_size)*60;
                 end
             elseif q ==3
                 % FR
                 for j = row:row3
                     dlc_loom_speed.FR(j-(row-1)) = pdist([dlc_tbl.FR_x(j), dlc_tbl.FR_y(j) ; dlc_tbl.FR_x(j-1), dlc_tbl.FR_y(j-1)])*(32/image_size)*60;
                 end
             elseif q ==4
                 % BL
                 for j = row:row3
                     dlc_loom_speed.BL(j-(row-1)) = pdist([dlc_tbl.BL_x(j), dlc_tbl.BL_y(j) ; dlc_tbl.BL_x(j-1), dlc_tbl.BL_y(j-1)])*(32/image_size)*60;
                 end
             elseif q ==5
                 % BR
                 for j = row:row3
                     dlc_loom_speed.BR(j-(row-1)) = pdist([dlc_tbl.BR_x(j), dlc_tbl.BR_y(j) ; dlc_tbl.BR_x(j-1), dlc_tbl.BR_y(j-1)])*(32/image_size)*60;
                 end
             elseif q ==6
                 % TB
                 for j = row:row3
                     dlc_loom_speed.TB(j-(row-1)) = pdist([dlc_tbl.TB_x(j), dlc_tbl.TB_y(j) ; dlc_tbl.TB_x(j-1), dlc_tbl.TB_y(j-1)])*(32/image_size)*60;
                 end
             elseif q ==7
                 % TM
                 for j = row:row3
                     dlc_loom_speed.TM(j-(row-1)) = pdist([dlc_tbl.TM_x(j), dlc_tbl.TM_y(j) ; dlc_tbl.TM_x(j-1), dlc_tbl.TM_y(j-1)])*(32/image_size)*60;
                 end
             elseif q ==8
                 % TT
                 for j = row:row3
                     dlc_loom_speed.TT(j-(row-1)) = pdist([dlc_tbl.TT_x(j), dlc_tbl.TT_y(j) ; dlc_tbl.TT_x(j-1), dlc_tbl.TT_y(j-1)])*(32/image_size)*60;
                 end
             end
         end
         
        dlc_loom_speed = smoothdata(dlc_loom_speed, 'rloess');
        dlc_loom_speed_array = movmean(table2array(dlc_loom_speed), 5); %83ms 
        
        loom_number = sprintf('%.f',i);
        loom_file_name = strcat(exp_name, '_L',loom_number);
        save((strcat('DLC_LOOM_SPEED_', loom_file_name,'.mat')), 'dlc_loom_speed');  
        save(fullfile(strcat(loom_sp_folder,'DLC_LOOM_SPEED_', loom_file_name,'.mat')), 'dlc_loom_speed'); 
        save(fullfile(strcat(loom_sp_folder,'DLC_LOOM_SPEED_ARRAY_', loom_file_name,'.mat')), 'dlc_loom_speed_array');  
  end
  
else
end 

end   