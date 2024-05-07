function make_xy_loom_array()

%% Created by Burnett - 30/01/20 

% Makes xy_array for 3s before the loom starts til 10s after the loom
% starts. Makes 99 ms moving mean average of speed across this time. 

% Save as an array with a single row - where each column is the
% instantaneous speed at this point. 

% Also makes XY_LOOM - the (x,y) position during the same time period (3s
% pre to 10s post loom). 

global exp_name output_folder diff_contrasts

if exist(strcat('ALL_LOOM_ROWS_', exp_name,'.mat')) 

    load((strcat('ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS'); 
    load(strcat('XY_array_', exp_name, '.mat'), 'xy_array');
    load(strcat('INFOTABLE_', exp_name, '.mat'), 'Info_Table');

     bin_window = 5; %average over 83 ms. 
     loom_rows = ALL_LOOM_ROWS(1,:)'; 
     num_looms = numel(loom_rows);

     xy_loom = [];
     XY_LOOM = [];
     
  for i = 1:num_looms
        
         row = loom_rows(i,1) - 179; % Start 3s (180 frames) before the loom stimulus starts 
         row3 = loom_rows(i,1) + 600;  % end 10s after loom starts. 
         
         if row<0 
             row = 1;
         end 
         
         if row3 > length(xy_array)
            row3 = length(xy_array);
         end 
         
         % Speed array 
          XY_loom = [];
          XY_loom = xy_array(row:row3,:); 
          
          M = movmean(XY_loom(:,5), bin_window)'; %col 5 = speed! Movmean over 83ms. 
          
          if numel(M) < 780 && row ~=1
              diff = 780- numel(M); 
              M(numel(M)+1:1081)= zeros(1,diff);
          elseif numel(M) < 780 && row ==1
              diff = 780 - numel(M); 
              added_zeros = zeros(1,diff); 
              M = [added_zeros, M];
          end 
          
          xy_loom = vertcat(xy_loom, M); %speed array 
          
          % Position array 
          
          XY_loom2 = XY_loom(:,3:4)'; % col 3 and 4 = adjusted x and y. 
          if numel(XY_loom2(1,:)) < 1081 && row ~=1
              diff = 780- numel(XY_loom2(1,:)); 
              num_orig = numel(XY_loom2(1,:));
              XY_loom2(1,num_orig+1:780)= zeros(1,diff); 
              XY_loom2(2,num_orig+1:780)= zeros(1,diff);
          elseif numel(XY_loom2(1,:)) < 1081 && row ==1
              diff = 780- numel(XY_loom2(1,:)); 
              added_zeros2 = zeros(2,diff);
              XY_loom2= [added_zeros2, XY_loom2]; 
          end 
          
          XY_LOOM = vertcat(XY_LOOM, XY_loom2); % x,y position array. 
  end
  
  if diff_contrasts == 1
      % Add Contrast of LOOM
      for jj = 1:num_looms
           loom_row = loom_rows(jj,1);
           contrast = Info_Table.DOT_RGB{loom_row};
           if contrast == [0,0,0]
               cval = 4; 
           elseif contrast == [0 12 25]
               cval = 3;
           elseif contrast == [0 25 50]
               cval = 2;
           elseif contrast == [0 37 75]
               cval = 1; 
           end 
           xy_loom(jj, 781) = cval;
      end
  end
  
     % Save in 'Exp Project' folder. 
    xyloom_folder = strcat(output_folder, 'SUMMARIES\LOOMS\XYLOOM_ARRAYS\');

    if ~exist(xyloom_folder,'dir')
        mkdir(xyloom_folder)
    end

    save(fullfile(strcat(xyloom_folder,'XY_LOOM_ARRAY_', exp_name,'.mat')), 'xy_loom');  
    
    % SAVE 
    xyloom_pos_folder = strcat(output_folder, 'SUMMARIES\LOOMS\XYLOOM_POS_ARRAYS\');
  
    if ~exist(xyloom_pos_folder,'dir')
        mkdir(xyloom_pos_folder)
    end
    
    save(fullfile(strcat(xyloom_pos_folder,'XY_LOOM_POS_ARRAY_', exp_name,'.mat')), 'XY_LOOM');  
  
else
end 

end        


