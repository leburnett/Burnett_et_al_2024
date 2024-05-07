function make_exit_analysis()
%% BANANA ANALYSIS

% Created by Burnett
% 19/01/22
global dir3_path currD
% Call this script in each folder with BANANA/MULTILOOM 
% The only variable that needs adding is the 'BANANA' variable. 
% BANANA = 1 when it is a BANANA acclim experiment, if it is a multiloom
% experiment then BANANA = 0; 

core_path = dir3_path;

if exist(strcat(core_path,'\', 'EXIT_ANALYSIS.mat'), 'file')
    load(strcat(core_path,'\', 'EXIT_ANALYSIS.mat'),'exit_analysis')
end 

if exist(strcat(core_path,'\', 'EXITS.mat'), 'file')
    load(strcat(core_path,'\', 'EXITS.mat'),'exits')
end 

% Make empty cell array in which to enter data. 
if ~exist('exit_analysis', 'var')
    exit_analysis = table(); 
    current_size = 0;
%     exits = table();
else 
     current_size = length(exit_analysis.Date); 
end

if ~exist('exits', 'var')
    exits = table();
end


if isempty(exit_analysis)
    p = 1;
else
    p =  current_size + 1; %which row to enter data into in the table.
end

% DEFAULTS - bear in mind difference in image size between C1 and C2. 
box_size_cm = 32; 
image_size = 416;
% image_size = 521; %518/524 
sz_ratio = box_size_cm/image_size; 
fps = 60; 
         
%% LOAD RELEVANT FILES. 

xy_file = dir("XY_array*");
xy_file = xy_file.name;
load(xy_file, 'xy_array')

% tr_file = dir("TRACK*");
% tr_file = tr_file.name;
% load(tr_file, 'TRACKING_DATA')

if contains(currD, "Acclim")
    % BANANA OR ACCLIM.
    date = str2double(xy_file(10:15));
    ani_ch = xy_file(end-19:end-14); 
    ani_num =  str2double(xy_file(end-17:end-14));
%     ani_num =  str2double(xy_file(end-14)); % FOR MATTEO - SHORTER
    exp = xy_file(end-12:end-4);
elseif contains(currD, "Loom")|| contains(currD, "Mult") 
    % Multiloom - 9 
%     date = str2double(xy_file(10:15));
%     ani_ch = xy_file(end-22:end-17);
%     ani_num =  str2double(xy_file(end-20:end-17));
%     exp = xy_file(end-15:end-4);

    % Mult / Loom  - 4 
    date = str2double(xy_file(10:15));
    ani_ch = xy_file(end-15:end-12);
    ani_num =  str2double(xy_file(end-15:end-12));
%     ani_num =  str2double(xy_file(end-12)); % FOR MATTEO - SHORTER
    exp = xy_file(end-10:end-4); 
end

exit_analysis.Date(p) = date;
exit_analysis.Animal(p) = ani_num;
exit_analysis.Exp{p} = exp;


%% Make PLOT 
% Assess whether TRACKING or XY is best. 
% 
% xWT = cell2mat(TRACKING_DATA(:,3));
% yWT = cell2mat(TRACKING_DATA(:,4));
% 
% nWT = numel(xWT);
% 
% figure;
% subplot(1,2,1)
% for q = 1:nWT-1
%         x = xWT(q);
%         y = image_size - yWT(q); 
%         x2 = xWT(q+1);
%         y2 = image_size - yWT(q+1);
%         plot([x, x2],[y,y2],'k')
%         hold on 
% end 
% box off
% axis([0 512 0 512])
% axis square 
% 
% % XY ARRAY 
% 
xWT = (xy_array(:,3));
yWT = (xy_array(:,4));

nWT = numel(xWT);
% 
% subplot(1,2,2)
% for q = 1:nWT2-1
%         x = xWT2(q);
%         y = image_size - yWT2(q); 
%         x2 = xWT2(q+1);
%         y2 = image_size - yWT2(q+1);
%         plot([x, x2],[y,y2],'r')
%         hold on 
% end 
% box off
% axis([0 512 0 512])
% axis square 

%% Choose which to use. 
% val = input('Type 1 for tracking_data (left) or 2 for xy_array (right)');
% close

% if val ==1  %TRACKING DATA
    Xvals = xWT; 
    Yvals = yWT; 
    
    % Calculate speed from X / Y 
    sp = zeros(nWT, 1);
    for k = 2:nWT
        A = pdist([xWT(k-1), yWT(k-1); xWT(k), yWT(k)]);
        sp(k,1) = A*sz_ratio*fps; % (cm/s)
    end
%     
% elseif val ==2  % XY ARRAY
%     Xvals = xWT2; 
%     Yvals = yWT2; 
%     
%     % Calculate speed from X / Y 
%     sp = zeros(nWT2, 1);
%     for k = 2:nWT2
%         A = pdist([xWT2(k-1), yWT2(k-1); xWT2(k), yWT2(k)]);
%         sp(k,1) = A*sz_ratio*fps; % (cm/s)
%     end
%     
% end 
    
exit_analysis.X{p} = Xvals;
exit_analysis.Y{p} = Yvals;
exit_analysis.Speed{p} = sp;


%% Use xy_array/ xy_table to look at DistShelter. 

dshelt = xy_array(:, 8);
n = numel(dshelt);

exit_analysis.DShelt{p} = dshelt;

% Threshold distance to trigger the loom:
if image_size == 518
    dist_loom_trig = 7.8;
elseif image_size == 524 
    dist_loom_trig = 9.2;
elseif image_size == 416
    dist_loom_trig = 14;
end
% % % % % %

% % Uncomment if you want to plot the figure:
% figure; plot(smooth(smooth(dshelt)), 'k');
% hold on; plot([0 n], [0 0], 'r');
% box off
% axis([0 n -10 30])
% if BANANA == 0
%     plot([0 n], [dist_loom_trig dist_loom_trig], 'b:')
% end
% savefig(gcf, fullfile('/Users/lauraburnett/Data_Analysis_Mac/Loom_Behaviour/Setd5/BANANA/PLOTS/DSHELT', strcat(string(date),'_',ani_ch, '_', exp,'.fig')))
% close
% % % % % %

% Find times out of shelter.
xx = find(diff(sign(smooth(dshelt))))+1; % This finds when 'sign' changes ie when goes from negative to positive.

if dshelt(1)>0 && numel(xx)>0 % If mouse starts outside the shelter then don't look at the first entrance back into the shelter.
    xx(1)= [];
end 

if dshelt(end)>0 && numel(xx)>0% If mouse ends outside the shelter remove that value. 
    xx(end)= [];
end

if ~isempty(xx)
    nallout = numel(xx);
    numout = nallout/2;
    vv = 1:2:nallout;
    vv2 = 2:2:nallout;
    
    xx_OUT = xx(vv); % Frame OUT
    xx_OUT(:,2) = xx(vv2); % Frame IN
    % xx_OUT(2:numout,3) = diff(xx_OUT(:,1)); % Difference between EXITs in frames
    % xx_OUT(:,4) = xx_OUT(:,3)/60; % Difference between EXITs in seconds
    
    % Inter Bout Interval
    for l = 1:numout-1
        IBI = xx_OUT(l+1,1)-xx_OUT(l,2);
        xx_OUT(l,3) = IBI; % in frames
        xx_OUT(l,4) = IBI/60; % in seconds not frames.
        
    end
    
    
    for i = 1:numout
        out = xx_OUT(i,1);
        in = xx_OUT(i,2);
        
        % Exit duration (in - out)
        exit_dur = in - out;
        xx_OUT(i,5) = exit_dur; % in frames
        xx_OUT(i,6) = exit_dur/60; % in seconds not frames.
        
        % Maximum distance travelled during exit.
        maxdist = max(dshelt(out:in));
        maxdist_f = find(dshelt(out:in) == max(dshelt(out:in)))+out;
        if numel(maxdist_f)>1
            maxdist_f = maxdist_f(1);
        end 
        xx_OUT(i,7) = maxdist;
        
        % Col 6 = Did the exit pass the trigger threshold?
        if maxdist > dist_loom_trig
            xx_OUT(i,8) = 1;
        else
            xx_OUT(i,8) = 0;
        end
        
        % Speed from 'out' to max distance.
        frames_out2max = maxdist_f-out;
        speed_out2max = maxdist/(frames_out2max/60);
        xx_OUT(i,9) = speed_out2max;
        
        % Speed from max distance to 'in'
        frames_max_to_in =  in - (maxdist_f);
        speed_max_to_in = maxdist/(frames_max_to_in/60);
        xx_OUT(i,10) = speed_max_to_in;
        
        if xx_OUT(i,8) ==1
            % Speed from trigger distance til 'in'
            trigdist_f = find(dshelt(out:in)>dist_loom_trig);
            trigdist_f = trigdist_f(1)+out;
            frames_trig2in =  in - (trigdist_f);
            time_trig2in = (frames_trig2in/60);
            xx_OUT(i,11) = time_trig2in;
        else
            xx_OUT(i, 11)= NaN;
        end
        
    end
    
    %% Remove any 'exits' where the mouse does not go more than 2cm - might just
    % be a head poke.
    xx_OUT(xx_OUT(:,7)<2, :) = [];
    
    % Recalculate the number of times out having removed the 'head pokes'.
    numout = numel(xx_OUT(:,1));
    numout_trig = numel(find(xx_OUT(:,8)==1));
    
    %% Add the summary of these values to the 'banana_analysis' table.
    if numout > 0
        % Find mean D
        mean_d = nanmean(xx_OUT(:,7));
        
        % Find Max d
        max_d = max(xx_OUT(:,7));
        if numel(max_d)>1
            max_d = max_d(1);
        end
        
        % Find mean interbout interval
        if numout>1
        IBI = nanmean(xx_OUT(:,4));
        else 
            IBI = NaN;
        end 
        
        % Find mean exit duration
        m_exit_duration = nanmean(xx_OUT(:,6));
        
        % Find mean out speed
        m_outspeed = nanmean(xx_OUT(:,9));
        
        % Find mean speed max2in
        m_spmax2in = nanmean(xx_OUT(:,10));
        
        % Find mean speed trg2in
        m_trig2in = nanmean(xx_OUT(:,11));
        
        % Add to banana_analysis:
        exit_analysis.NumOut(p) = numout;
        exit_analysis.NumOutTrig(p) = numout_trig;
        
        exit_analysis.MeanD(p) = mean_d;
        exit_analysis.MaxD(p) = max_d;
        exit_analysis.InterBoutInterval(p) = IBI;
        exit_analysis.ExitDur(p) = m_exit_duration;
        
        exit_analysis.outspeed(p) = m_outspeed;
        exit_analysis.sp_max2in(p) = m_spmax2in;
        exit_analysis.trig2in(p) = m_trig2in;
        
    elseif numout == 0
        exit_analysis.NumOut(p) = 0;
        exit_analysis.NumOutTrig(p) = 0;
        exit_analysis.MeanD(p) = NaN;
        exit_analysis.MaxD(p) = NaN;
        exit_analysis.InterBoutInterval(p) = NaN;
        exit_analysis.ExitDur(p) = NaN;
        exit_analysis.outspeed(p) = NaN;
        exit_analysis.sp_max2in(p) = NaN;
        exit_analysis.trig2in(p) = NaN;
    end
    
else
    numout = 0;
    exit_analysis.NumOut(p) = 0;
    exit_analysis.NumOutTrig(p) = 0;
    exit_analysis.MeanD(p) = NaN;
    exit_analysis.MaxD(p) = NaN;
    exit_analysis.InterBoutInterval(p) = NaN;
    exit_analysis.ExitDur(p) = NaN;
    exit_analysis.outspeed(p) = NaN;
    exit_analysis.sp_max2in(p) = NaN;
    exit_analysis.trig2in(p) = NaN;
    
end
%% Need to see if ALL_LOOM_ROWS exists
% Find % Looms triggered. 

loom_file = dir("ALL_LOOM*");
if ~isempty(loom_file)
    loom_file = loom_file.name;
    load(loom_file, 'ALL_LOOM_ROWS')
    
    n_bouts = numel(ALL_LOOM_ROWS(1,:));
    n_looms = nnz(ALL_LOOM_ROWS);
    
    exit_analysis.NumBOUTS(p) = n_bouts;
    exit_analysis.NumLOOMS(p) = n_looms;
    exit_analysis.AvLOOMS(p) = n_looms/n_bouts;
    exit_analysis.LoomRows{p} = ALL_LOOM_ROWS(1,:);
    
else
    exit_analysis.NumBOUTS(p) = 0;
    exit_analysis.NumLOOMS(p) = 0;
    exit_analysis.AvLOOMS(p) = 0;
    exit_analysis.LoomRows{p} = NaN;
end

% het_animals = [3959, 4244, 4473, 6558, 6560, 6562, 6610, 7269, 7476, 7614, 7790, 5148, 5596, 5542, 5597, 0593, 0595];
% het_animals = [2375, 2377, 2382, 2628, 2637, 2754, 2829, 2832, 2833, 2900, 2901, 2902];
% het_animals = [1385, 1386, 1388, 1394, 2593, 2594, 2708, 2709, 3903, 4369, 4373];
% het_animals = [1125, 1120, 1121, 1282];
 het_animals = [2080, 2081, 2085, 2341];
% het_animals = [2479, 2482, 2688, 2690, 0310, 0311]; 
%   het_animals = [0344, 0345, 0347, 0349];
% het_animals = [0300, 0301, 0580, 0622];

if ismember(ani_num, het_animals)
    geno = 2;
else
    geno = 1;
end
exit_analysis.geno(p) = geno;

if numout>0
    % Convert xx_OUT into table: banana_exits
    Animal = repmat(ani_num, numout, 1);
    Date = repmat(date, numout, 1);
    Exp = repmat({exp}, numout, 1);
    ExitNum = (1:1:numout)';
    FrameOut = xx_OUT(:,1);
    FrameIn = xx_OUT(:,2);
    FramesIBI = xx_OUT(:,3);
    IBI = xx_OUT(:,4);
    FramesExit = xx_OUT(:,5);
    ExitDur = xx_OUT(:,6);
    MaxDistance = xx_OUT(:,7);
    LoomTrig = xx_OUT(:,8);
    OutSpeed = xx_OUT(:,9);
    Max2InSpeed = xx_OUT(:,10);
    Trig2In = xx_OUT(:,11);
    Geno = repmat(geno, numout, 1);
    
    b = table(Animal, Date, Exp, ExitNum, FrameOut, FrameIn, FramesIBI, IBI, FramesExit, ExitDur, MaxDistance, LoomTrig, OutSpeed, Max2InSpeed, Trig2In, Geno);
    exits = vertcat(exits, b);
    save(strcat(core_path,'\','EXITS.mat'),'exits')
end

save(strcat(core_path,'\','EXIT_ANALYSIS.mat'),'exit_analysis')


end 