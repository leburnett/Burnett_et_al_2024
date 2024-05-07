function make_loom_array()

%% Created by Burnett - 18/03/20
% Last modified 26/03/20

% FIND THE NUMBER OF LOOMS/BOUTS AND THE ROWS WHICH CORRESPOND TO THEIR START.

% To be used after 'make_xy_array.m'
% Finds - loom rows and loom bout rows.
% Generates arrays: RESPONSE_TYPE and LOCOINDEX

global exp_name output_folder ILI

if ILI == 0
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
%     n = numel(Info(:,1));
    
%     rows_looming = find(cell2mat(Info(:,5))==1);
     rows_looming = find(cell2mat(Info(:,4))==5);
     
    % Find rows in 'Info' - Col5 that == 1. corresponds to when the stimulus
    % was looming. If this == 0. Then code stops here. NO LOOMS OCCURRED.
    % % % % % % % % % %  % % % % % % % % % % % % % % % % % % % % % %  % % %
    
    if ~isempty(rows_looming)
        
        % New method of finding looms and loom bouts - 25/03/20
       n_looms = numel(rows_looming);
        
        if n_looms==1 % ONLY ONE LOOM OCCURRED.
            NLooms = 1;
            rows_looms = rows_looming(1);
            NBouts = 1;
            rows_bouts = rows_looming(1);
            ALL_LOOM_ROWS(1,1)= rows_looms(1);
        else
            
            % Find gap betwee rows when looming.
            for i = 2:numel(rows_looming)
                rows_looming(i,2) = rows_looming(i,1)-rows_looming(i-1,1);
            end
            
            find_looms = find(rows_looming(:,2)>1); % Find gaps between individual looms.
            
            
            NLooms = numel(find_looms(:,1))+1;
            
            for j = 1:numel(find_looms)
                r = find_looms(j,1);
                find_looms(j,2)= rows_looming(r,1);
            end
            
            first_loom = rows_looming(1,1);
            rows_looms = [first_loom, find_looms(:,2)']'; %Vertical array of rows in 'Info' that correspond to the start of a loom.
            
            % Finding Bouts
            find_bouts = find(rows_looming(:,2)>70); % Gaps between loom bouts.
            
            if isempty(find_bouts)
                NBouts = 1;
                rows_bouts = rows_looms(1,1);
            else
                NBouts = numel(find_bouts)+1;
                
                for k = 1:numel(find_bouts)
                    r = find_bouts(k,1);
                    find_bouts(k,2)= rows_looming(r,1);
                end
                
                rows_bouts = [first_loom, find_bouts(:,2)']';
            end
            
            % 'ALL_LOOM_ROWS' is an array where each COLUMN is a new LOOM BOUT. Each
            % ROW contains the row_value in 'Info' for the beginning of that individual loom
            % presentation.
            
            ALL_LOOM_ROWS(1,1)= rows_looms(1);
            
            for j = 2:NLooms
                row_diff = rows_looms(j,1)- rows_looms(j-1,1);
                cols = size(ALL_LOOM_ROWS,2);
                nrows = nnz(ALL_LOOM_ROWS(:,cols));
                if row_diff >70 %If the number of rows difference between the two looms is <49 then they are part of the same 'loom_bout'
                    ALL_LOOM_ROWS(1,cols+1)=rows_looms(j,1);
                else %If they were >49 rows apart, they belong to different loom_bouts.
                    ALL_LOOM_ROWS(nrows+1,cols)=rows_looms(j,1);
                end
            end
        end
        
        %% SAVE ALL_LOOM_ROWS
        save((strcat('ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS');
        
        
        %% Save in 'SUMMARIES' folder.
        all_looms_folder = strcat(output_folder, 'SUMMARIES\LOOMS\ALL_LOOM_ROWS\');
        
        if ~exist(all_looms_folder,'dir')
            mkdir(all_looms_folder)
        end
        
        save(fullfile(strcat(all_looms_folder,'ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS');
        
        %%
%         ALL_BOUT_ROWS = rows_bouts; % Each column is a new bout of looms.
        
%         %% Then make LOOM_ARRAY
%         if exp_name(end-4:end) == "_Loom" % NORMAL LOOM EXP
%             
%             % load(strcat('STIMROWS_', exp_name,'.mat'), 'STIMROWS');
%             start_row_L2 = 1200 ; %STIMROWS{8,4};
%             
%             loom_array = zeros(1,6);
%             % Col 1 = # Bouts in total.
%             % col 2 = # bouts L1
%             % Col 3 = # Bouts L2
%             % Col 4 = # Looms
%             % Col 5 = # Looms L1
%             % Col 6 = # Looms L2
%             
%             loom_array(1,1) = NBouts;
%             loom_array(1,4) = NLooms;
%             
%             for k = 1:NBouts
%                 start_row = ALL_LOOM_ROWS(1,k); % row when loom bout starts
%                 n_looms_in_bout = nnz(ALL_LOOM_ROWS(:,k));
%                 if start_row < start_row_L2
%                     loom_array(1,5)= loom_array(1,5) + n_looms_in_bout;
%                     loom_array(1,2)= loom_array(1,2) + 1;
%                 elseif start_row >= start_row_L2
%                     loom_array(1,6)= loom_array(1,6) + n_looms_in_bout;
%                     loom_array(1,3)= loom_array(1,3) + 1;
%                 end
%             end
%             
%         elseif exp_name(end-4:end) == "iLoom" % MULTILOOM LOOM EXP
%             
%             loom_array = zeros(1,10);
%             % Col 1 = # Bouts in total.
%             % Col 2 = # Bouts - First Quarter
%             % Col 3 = # Bouts - Second Quarter
%             % Col 4 = # Bouts - Third Quarter
%             % Col 5 = # Bouts - Fourth Quarter
%             % Col 6 = # Looms
%             % Col 7 = # Looms - First Quarter
%             % Col 8 = # Looms - Second Quarter
%             % Col 9 = # Looms - Third Quarter
%             % Col 10 = # Looms - Fourth Quarter
%             Q1 = round(n*0.25);
%             Q2 = round(n*0.25);
%             Q3 = round(n*0.25);
%             
%             loom_array(1,1) = NBouts;
%             loom_array(1,6) = NLooms;
%             
%             for k = 1:NBouts
%                 start_row = ALL_LOOM_ROWS(1,k); % row when loom bout starts
%                 n_looms_in_bout = nnz(ALL_LOOM_ROWS(:,k));
%                 if start_row < Q1
%                     loom_array(1,7)= loom_array(1,7) + n_looms_in_bout;
%                     loom_array(1,2)= loom_array(1,2) + 1;
%                 elseif start_row >= Q1 && start_row < Q2
%                     loom_array(1,8)= loom_array(1,8) + n_looms_in_bout;
%                     loom_array(1,3)= loom_array(1,3) + 1;
%                 elseif start_row >= Q2 && start_row < Q3
%                     loom_array(1,9)= loom_array(1,9) + n_looms_in_bout;
%                     loom_array(1,4)= loom_array(1,4) + 1;
%                 elseif start_row >= Q3 && start_row < n
%                     loom_array(1,10)= loom_array(1,10) + n_looms_in_bout;
%                     loom_array(1,5)= loom_array(1,5) + 1;
%                 end
%             end
%         end
%         
%         %Save locally
%         save((strcat('LOOM_ARRAY_', exp_name,'.mat')), 'loom_array');
%         
%         % Save in 'Exp Project' LOOM folder.
%         loom_folder = strcat(output_folder, 'SUMMARIES\LOOMS\LOOM_ARRAYS\');
%         
%         if ~exist(loom_folder,'dir')
%             mkdir(loom_folder)
%         end
%         
%         save(fullfile(strcat(loom_folder,'LOOM_ARRAY_', exp_name,'.mat')), 'loom_array');
%         
%     else
        % % %    If later interested in the time between looms:
        % % %      % open loom bout rows
        % % %      % Change orientation to vertical
        % % %      % Make column 2 the difference in rowsbetween the loom bouts
        % % %      % Make column 3 the differnce in time between the loom bouts.
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
elseif ILI == 1 % If it is a different ILI trial - then will think looms are part of a different bout when they're not!! 
    
    load(strcat('INFO_', exp_name, '.mat'), 'Info');
%     n = numel(Info(:,1));
    
%     rows_looming = find(cell2mat(Info(:,5))==1);
    rows_looming = find(cell2mat(Info(:,4))==5);
    % Find rows in 'Info' - Col5 that == 1. corresponds to when the stimulus
    % was looming. If this == 0. Then code stops here. NO LOOMS OCCURRED.
    % % % % % % % % % %  % % % % % % % % % % % % % % % % % % % % % %  % % %
    
    if ~isempty(rows_looming)
        
        % New method of finding looms and loom bouts - 25/03/20
        if numel(rows_looming)==1 % ONLY ONE LOOM OCCURRED.
            NLooms = 1;
            rows_looms = rows_looming(1);
            NBouts = 1;
            rows_bouts = rows_looming(1);
            ALL_LOOM_ROWS(1,1)= rows_looms(1);
        else
            
            % Find gap betwee rows when looming.
            for i = 2:numel(rows_looming)
                rows_looming(i,2) = rows_looming(i,1)-rows_looming(i-1,1);
            end
            
            find_looms = find(rows_looming(:,2)>1); % Find gaps between individual looms.
            
            NLooms = numel(find_looms(:,1))+1;
            
            for j = 1:numel(find_looms)
                r = find_looms(j,1);
                find_looms(j,2)= rows_looming(r,1);
            end
            
            first_loom = rows_looming(1,1);
            rows_looms = [first_loom, find_looms(:,2)']'; %Vertical array of rows in 'Info' that correspond to the start of a loom.
           
            
            % Add a column  with the ILI of each of these looms.
            for jj = 1:NLooms
                row = rows_looms(jj,1); 
                ili_val = Info{row,7}; 
                rows_looms(jj,2) = ili_val;
            end 
            
            rows_looms(2:end, 3)=  diff(rows_looms(:,1));  

            % Finding Bouts - THIS WILL DEPEND ON ILI!! 
            for i2 =1:NLooms
                if rows_looms(i2,3)> (60*rows_looms(i2,2))+50
                    rows_looms(i2,4) = 1;
                else 
                    rows_looms(i2,4) = 0;
                end 
            end 
                    
            find_bouts = find(rows_looms(:,4)==1); % Gaps between loom bouts.
            
            if isempty(find_bouts)
                NBouts = 1;
                rows_bouts = rows_looms(1,1);
            else
                NBouts = numel(find_bouts)+1;
                
                for k = 1:numel(find_bouts)
                    r = find_bouts(k,1);
                    find_bouts(k,2)= rows_looms(r,1);
                end
                
                rows_bouts = [first_loom, find_bouts(:,2)']';
            end
            
            % 'ALL_LOOM_ROWS' is an array where each COLUMN is a new LOOM BOUT. Each
            % ROW contains the row_value in 'Info' for the beginning of that individual loom
            % presentation.
            
            ALL_LOOM_ROWS(1,1)= rows_looms(1);
            
            for j = 2:NLooms
                diffbout =  rows_looms(j,4);
                cols = size(ALL_LOOM_ROWS,2);
                nrows = nnz(ALL_LOOM_ROWS(:,cols));
                if diffbout == 1 %If the number of rows difference between the two looms is <49 then they are part of the same 'loom_bout'
                    ALL_LOOM_ROWS(1,cols+1)=rows_looms(j,1);
                elseif diffbout == 0 %If they were >49 rows apart, they belong to different loom_bouts.
                    ALL_LOOM_ROWS(nrows+1,cols)=rows_looms(j,1);
                end
            end
        end
        
        %% SAVE ALL_LOOM_ROWS
        save((strcat('ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS', 'rows_looms');
        
        
        %% Save in 'SUMMARIES' folder.
        all_looms_folder = strcat(output_folder, 'SUMMARIES\LOOMS\ALL_LOOM_ROWS\');
        
        if ~exist(all_looms_folder,'dir')
            mkdir(all_looms_folder)
        end
        
        save(fullfile(strcat(all_looms_folder,'ALL_LOOM_ROWS_', exp_name,'.mat')), 'ALL_LOOM_ROWS', 'rows_looms');
        
        %%
%         ALL_BOUT_ROWS = rows_bouts; % Each column is a new bout of looms.
%         
%         %% Then make LOOM_ARRAY
%         if exp_name(end-4:end) == "_Loom" % NORMAL LOOM EXP
%             
%             % load(strcat('STIMROWS_', exp_name,'.mat'), 'STIMROWS');
%             start_row_L2 = 1200 ; %STIMROWS{8,4};
%             
%             loom_array = zeros(1,6);
%             % Col 1 = # Bouts in total.
%             % col 2 = # bouts L1
%             % Col 3 = # Bouts L2
%             % Col 4 = # Looms
%             % Col 5 = # Looms L1
%             % Col 6 = # Looms L2
%             
%             loom_array(1,1) = NBouts;
%             loom_array(1,4) = NLooms;
%             
%             for k = 1:NBouts
%                 start_row = ALL_LOOM_ROWS(1,k); % row when loom bout starts
%                 n_looms_in_bout = nnz(ALL_LOOM_ROWS(:,k));
%                 if start_row < start_row_L2
%                     loom_array(1,5)= loom_array(1,5) + n_looms_in_bout;
%                     loom_array(1,2)= loom_array(1,2) + 1;
%                 elseif start_row >= start_row_L2
%                     loom_array(1,6)= loom_array(1,6) + n_looms_in_bout;
%                     loom_array(1,3)= loom_array(1,3) + 1;
%                 end
%             end
%             
%         elseif exp_name(end-4:end) == "iLoom" % MULTILOOM LOOM EXP
%             
%             loom_array = zeros(1,10);
%             % Col 1 = # Bouts in total.
%             % Col 2 = # Bouts - First Quarter
%             % Col 3 = # Bouts - Second Quarter
%             % Col 4 = # Bouts - Third Quarter
%             % Col 5 = # Bouts - Fourth Quarter
%             % Col 6 = # Looms
%             % Col 7 = # Looms - First Quarter
%             % Col 8 = # Looms - Second Quarter
%             % Col 9 = # Looms - Third Quarter
%             % Col 10 = # Looms - Fourth Quarter
%             Q1 = round(n*0.25);
%             Q2 = round(n*0.25);
%             Q3 = round(n*0.25);
%             
%             loom_array(1,1) = NBouts;
%             loom_array(1,6) = NLooms;
%             
%             for k = 1:NBouts
%                 start_row = ALL_LOOM_ROWS(1,k); % row when loom bout starts
%                 n_looms_in_bout = nnz(ALL_LOOM_ROWS(:,k));
%                 if start_row < Q1
%                     loom_array(1,7)= loom_array(1,7) + n_looms_in_bout;
%                     loom_array(1,2)= loom_array(1,2) + 1;
%                 elseif start_row >= Q1 && start_row < Q2
%                     loom_array(1,8)= loom_array(1,8) + n_looms_in_bout;
%                     loom_array(1,3)= loom_array(1,3) + 1;
%                 elseif start_row >= Q2 && start_row < Q3
%                     loom_array(1,9)= loom_array(1,9) + n_looms_in_bout;
%                     loom_array(1,4)= loom_array(1,4) + 1;
%                 elseif start_row >= Q3 && start_row < n
%                     loom_array(1,10)= loom_array(1,10) + n_looms_in_bout;
%                     loom_array(1,5)= loom_array(1,5) + 1;
%                 end
%             end
%         end
%         
%         %Save locally
%         save((strcat('LOOM_ARRAY_', exp_name,'.mat')), 'loom_array');
%         
%         % Save in 'Exp Project' LOOM folder.
%         loom_folder = strcat(output_folder, 'SUMMARIES\LOOMS\LOOM_ARRAYS\');
%         
%         if ~exist(loom_folder,'dir')
%             mkdir(loom_folder)
%         end
%         
%         save(fullfile(strcat(loom_folder,'LOOM_ARRAY_', exp_name,'.mat')), 'loom_array');
%         
%     else
        % % %    If later interested in the time between looms:
        % % %      % open loom bout rows
        % % %      % Change orientation to vertical
        % % %      % Make column 2 the difference in rowsbetween the loom bouts
        % % %      % Make column 3 the differnce in time between the loom bouts.
    end
    
    
    
end
%%

% SAVE ALL 'loom_array' files in grouped folder.

%       save(fullfile(strcat('C:\Data_analysis\DATA\Setd5-Cohort2\SUMMARIES\LOOMS\num_loom_arrays\','loom_array_', exp_name,'.mat')), 'loom_array');
%       save(fullfile(strcat('C:\Data_analysis\DATA\Setd5\ALL_DAY_SUMMARIES\LOOMS\loom_arrays\','loom_array_', exp_name,'.mat')), 'loom_array');
%            save(fullfile(strcat('C:\Data_analysis\DATA\1912_Setd5\ALL_DAY_SUMMARIES\LOOMS\loom_arrays\NoShelter\','loom_array_', exp_name,'.mat')), 'loom_array');

% Cohort 2 - NoShelter
%         save(fullfile(strcat('C:\Data_analysis\DATA\2001_Setd5-Cohort2\SUMMARIES\NOSHELTER\LOOMS\loom_arrays\','loom_array_', exp_name,'.mat')), 'loom_array');

% DiffContrast
%         save(fullfile(strcat('C:\Data_analysis\DATA\BOTH_SETD5_COHORTS_\DiffContrasts\LOOM_ARRAYS\','loom_array_', exp_name,'.mat')), 'loom_array');

% NO Shelter