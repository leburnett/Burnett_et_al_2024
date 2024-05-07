% OPTOGENETICS - PRE-PROCESSING
% Combine all the arrays together
% Burnett - 28/4/22

%% Step 1 - LOAD THE FILES WITH ALL THE DATA

load('Setd5VGlut2_C1_N4_all.mat');

% Check all the files have the same size:
sz1 = height(xy_analysis);
sz2 = numel(xy_opto(:,1));
sz3 = height(xy_opto_info);
sz4 = numel(XY_OPTO(:,1));

if sz1 == sz2 && sz1 == sz3
    sizes_correct = 1;
else
    sizes_correct = 0;
end

if sz1*2 == sz4
    sz_double = 1
else 
    sz_double = 0
end 

% Check how many days/ animals there are!
all_animals1 = unique(string(xy_analysis.Animal));
all_animals2 = unique(string(xy_opto_info.Animal));

all_dates1 = unique(xy_analysis.Date);
all_dates2 = unique(xy_opto_info.Date);


%%  C1

a1 = xy_analysis;
a2 = xy_opto;
a3 = XY_OPTO;
a4 = xy_opto_info; 
szA = sz1; 
%%  C2

b1 = xy_analysis;
b2 = xy_opto;
b3 = XY_OPTO;
b4 = xy_opto_info; 
szB = sz1; 

%%  C3

c1 = xy_analysis;
c2 = xy_opto;
c3 = XY_OPTO;
c4 = xy_opto_info; 
szC = sz1; 
%%  C4

d1 = xy_analysis;
d2 = xy_opto;
d3 = XY_OPTO;
d4 = xy_opto_info; 
szD = sz1; 
%%  C5

e1 = xy_analysis;
e2 = xy_opto;
e3 = XY_OPTO;
e4 = xy_opto_info; 
szE = sz1; 

%% COMBINE THE ARRAYS

xy_analysis = vertcat(a1,b1,c1,d1,e1);
xy_opto = vertcat(a2,b2,c2,d2,e2);
XY_OPTO = vertcat(a3,b3,c3,d3,e3);
xy_opto_info = vertcat(a4,b4,c4,d4,e4);

save('220428_Setd5_VG_OPTO_DATA_C1-C5.mat', 'xy_analysis', 'xy_opto', 'xy_opto_info', 'XY_OPTO');
% saved here: '/Users/lauraburnett/Data_Analysis_Mac/OPTO/Setd5/DATA'


%% DON'T USE XY_ANALYSIS - somehow not the same... xy_opto and xy_opto_info should be the same!!!

% xy_opto_info

% Date
% Animal
% Exp
% Loom
% Loom row
% Num Pulses - sort of irrelevant
% T_pulse
% FreqPulse
% EC
% Geno

% Columns to ADD:

% 1 - Distance from shelter edge when light triggered!

% 2 - Return to shelter? - Does the mouse return to the shelter within the 10s

% 3 - Speed at the time of the light turning on
% 4 - Speed in the 300-800ms after the light is on.
% 5 - Speed in the first ~6s before light on.
% 6 - Log(Speed immed / speed at loom)

% 7 - Maximum speed from light on til end of 10s after light. 
% 8 - Time to maximum speed
% 9 - Log(Max Speed/ speed at loom)

% 10 - Maximum acceleration from light on - 10s after
% 11 - Time to the maximum acc 


n_all = height(xy_opto_info);
%% ADD VARIABLES TO XY_OPTO_INFO

for i = 1:n_all 
    
    speed = xy_opto(i, :);
    acc = diff(speed);
    
    X = XY_OPTO((i*2)-1, :);
    Y = 416 - XY_OPTO((i*2), :); % Take value away from 416 to get correct Y value - flipped. 
    
    light_on_row = 600; 
    
    %% DISTANCE 
  
    % Distance from shelter. 
    shelter_size = 90;
    shelterC = [360, 60]; %[330, 325];
    sz_ratio = 32/416;
    
    dbox_thresh = shelter_size*sz_ratio; % 90pixels from centre ~6cm.
    
    Q = NaN(1200,1);
    
    for j = 1:1200
        D = pdist([X(j), Y(j); shelterC(1),shelterC(2)]);
        Q(j,1) = (D*sz_ratio)-dbox_thresh; 
    end
    
    % (1) - Distance from shelter when light ON
    d_shelt_light_on = Q(light_on_row);
    xy_opto_info.DShelt_LIGHT{i} = d_shelt_light_on;
    
    % Find when first returns back to shelter after the light is on:
    in_shelter = find(Q(light_on_row+5:end)<0);
    
    % (2) - Return to shelter? 
    if ~isempty(in_shelter)
        first_row_in_shelter = in_shelter(1); %number of rows after when the mouse is in the shelter. Time in frames. 
        row_shelter = light_on_row+5+first_row_in_shelter; % The actual row number when mouse is back in shelter. 
        
        time_to_shelter = first_row_in_shelter/60; % time in seconds til mouse back in shelter. 
        speed_to_shelter = d_shelt_light_on/time_to_shelter; 
        
        return_to_shelter = 1;
    else 
        first_row_in_shelter = NaN;
        return_to_shelter = 0;
    end 
    
    xy_opto_info.Return2Shelter{i} = return_to_shelter;
    xy_opto_info.ROW_Shelter{i} = first_row_in_shelter;
    
    %% SPEED 
    
    % (3) - Speed at LIGHT ON
    sp_at = nanmean(speed(light_on_row-5:light_on_row+5));
    xy_opto_info.SpAt{i} = sp_at;
    
     % (4) - Speed 300-800ms after light on. 
    sp_immed = nanmean(speed(light_on_row+18:light_on_row+48));
    xy_opto_info.SpImmed{i} = sp_immed;
    
    % (5)- Speed in the first ~6s before light on.
    sp_before = nanmean(speed(1:500));
    xy_opto_info.SpB4{i} = sp_before;
    
    % (6) - Log(Speed immed / speed at loom)
    delta_immed = sp_immed/sp_at;
    log_delta_imm = log(delta_immed);
    xy_opto_info.DeltaImmed{i} = delta_immed;
    xy_opto_info.LogDeltaImmed{i} = log_delta_imm;
    
    %% MAX SPEED 
    
    % (7) - Maximum speed from light on til end of 10s after light. 
    maxsp = max(speed(light_on_row:end));
    maxsp_light = max(speed(light_on_row:light_on_row+180));
    
    xy_opto_info.MaxSp{i} = maxsp;
    xy_opto_info.MaxSpLIGHT{i} = maxsp_light;
    
    rows_maxsp = find(speed(light_on_row:end) == maxsp); % How many rows til the max speed is reached. 
    rows_maxsp_LIGHT = find(speed(light_on_row:light_on_row+180) == maxsp_light); % How many rows til the max speed is reached. 

    if numel(rows_maxsp)>1
        rows_maxsp = rows_maxsp(1);
    end 
    
    if numel(rows_maxsp_LIGHT)>1
        rows_maxsp_LIGHT = rows_maxsp_LIGHT(1);
    end 
    
    row_MAXSP = rows_maxsp+light_on_row;  % Row where max speed reached. 
    
    % (8) - Time to maximum speed
    time_to_maxsp = rows_maxsp/60; % Time in seconds
    xy_opto_info.T2M{i} = time_to_maxsp;
    
    time_to_maxsp_Light = rows_maxsp_LIGHT/60; % Time in seconds
    xy_opto_info.T2M_LIGHT{i} = time_to_maxsp_Light;
    
    % (9) - Log(Max Speed/ speed at loom)
    delta_max = maxsp/sp_at;
    log_delta_max = log(delta_max);
    
    xy_opto_info.DeltaMax{i} = delta_max;
    xy_opto_info.LogDeltaMax{i} = log_delta_max;
    
    delta_max_light = maxsp_light/sp_at;
    log_delta_max_light = log(delta_max_light);
    
    xy_opto_info.DeltaMaxLight{i} = delta_max_light;
    xy_opto_info.LogDeltaMaxLight{i} = log_delta_max_light;
    
    %% ACCELERATION

     % (10) - Maximum acceleration from light on - 10s after
     maxacc = max(acc(light_on_row:end));
     maxacc_light = max(acc(light_on_row:light_on_row+180)); %Within light on 
     
     xy_opto_info.MaxAcc{i} = maxacc;
     xy_opto_info.MaxAccLight{i} = maxacc_light;

     rows_maxacc = find(acc(light_on_row:end) == maxacc); % How many rows til the max acc is reached. 
     rows_maxacc_light = find(acc(light_on_row:light_on_row+180) == maxacc_light); % How many rows til the max acc is reached. 

    if numel(rows_maxacc)>1
        rows_maxacc = rows_maxacc(1);
    end 
    
    if numel(rows_maxacc_light)>1
        rows_maxacc_light = rows_maxacc_light(1);
    end 
    
    row_ACC = rows_maxacc+light_on_row;  % Row where max speed reached. 
    
    % (11) - Time to the maximum acc 
     time_to_maxacc = rows_maxacc/60; % Time in seconds
     xy_opto_info.T2MaxAcc{i} = time_to_maxacc;
     
     time_to_maxacc_light = rows_maxacc_light/60; % Time in seconds
     xy_opto_info.T2MaxAccLIGHT{i} = time_to_maxacc_light;
     
     %% Group 
     geno = xy_opto_info.Geno{i};
     if geno == "wt"
         gp = 1;
     else
         gp = 2;
     end 
     
     if xy_opto_info.Animal{i} == "MJ0580"
         gp = 3; % CONTROl
     end 
    
     xy_opto_info.Group{i} = gp; 
end 


%% Visualise trajectory
figure; 
subplot(2,1,1)
for j = 1:1199
    plot([X(j), X(j+1)], [Y(j), Y(j+1)], 'k')
    hold on 
end
box off
axis([0 416 0 416])
plot(360, 60, 'b.', 'MarkerSize', 50)
viscircles(shelterC, shelter_size, 'Color', 'k')  
subplot(2,1,2)
plot(Q, 'k')       
         

%% SPLIT BASED ON LP/PW/FREQ
% Rows of xy_opto_info from '/Users/lauraburnett/Data_Analysis_Mac/OPTO/Setd5/DATA/220428_Setd5_VG_OPTO_DATA_C1-C5.mat'. 

LP = [1:88, 109:232, 237:556, 916:1119, 1968:2424, 2676:2782, 2820:2855];
PW = [557:673, 1120:1361, 1485:1730, 2425:2543, 2783:2801, 2856:2863];
FR = [674:915, 1362:1484, 1731:1967, 2544:2675, 2802:2819];

Hz = [89:108, 233:236]; %20Hz

% Laser Power
xy_optoA = xy_opto(LP, :);
xy_opto_infoA= xy_opto_info(LP, :);
% XY_OPTOA = XY_OPTO([LP*2, (LP*2)-1])

% Pulse Width
xy_optoB = xy_opto(PW, :);
xy_opto_infoB= xy_opto_info(PW, :);

% Frequency
xy_optoC = xy_opto(FR, :);
xy_opto_infoC= xy_opto_info(FR, :);

save('220429_Setd5_OPTO_C1-C5_Split.mat', 'xy_optoA', 'xy_optoB', 'xy_optoC', 'xy_opto_infoA', 'xy_opto_infoB', 'xy_opto_infoC');


%% TRY AND FIGURE OUT HOW TO GET XY_OPTO AS WELL. 

% % FOR XY_OPTO
% 
% nLP = numel(LP);
% LP2 = NaN(nLP*2,1);
% 
% for jj = 1:nLP
%     LP2(jj*2-1) = LP(jj);
%     LP2(jj*2) = LP(jj+1)
%     
%     
% end 


%% THEN WE'RE READY TO GOOOOO... 
% Move to 'analyse_optogenetics'



































%% VISUALISE RESPONSES OF ONE ANIMAL ON ONE DAY:

DAY_TO_ANA = "210906";

n = numel(xy_opto(:,1));
col = 'r';
frames_b4 = 595;
total_lighton_frames = 60;  % light was on for 1s. - 60 frames.

all_animals = unique(string(xy_opto_info.Animal));
n_animals = numel(all_animals);

for i = 1:n_animals
    ani = string(all_animals{i}); 
    all_ani = find(xy_analysis.Animal == ani); 
    n_trials_ani = numel(all_ani);
    
    all_EC = unique(cell2mat(xy_analysis.EC)); 
    n_EC = numel(all_EC); 
    figure
    
    for j = 1:n_EC
        EC_str = all_EC(j,1); 
        all_ani_EC = find(cell2mat(xy_analysis.EC) == EC_str & xy_analysis.Animal == ani & cell2mat(xy_analysis.DShelt_Start)>0 & string(xy_analysis.Date) == DAY_TO_ANA);

        col = 1 -[j/n_EC j/n_EC j/n_EC];

        if ~isempty(all_ani_EC)
            if numel(all_ani_EC)>1
            av_resp = mean(xy_opto(all_ani_EC, :));
            elseif numel(all_ani_EC)==1
                av_resp = (xy_opto(all_ani_EC, :));
            end 
%             av_resp = smooth(av_resp);
            plot(av_resp, 'Color', col, 'LineWidth', 1.2);
            hold on 
        end 
        plot([frames_b4 frames_b4], [0 90], 'k:', 'LineWidth', 1.2)
        plot([frames_b4+total_lighton_frames frames_b4+total_lighton_frames], [0 90], 'k:', 'LineWidth', 1)
        xticklabels({''})
    end 
    title(ani)
    axis([500 800 0 80])
    box off
    ylabel('Speed - cm/s')
    xticks([510, 540, 570, 600, 630, 660, 690, 720, 750, 780])
    xticklabels({'-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'})
    ax = gca;
    ax.FontSize = 14; 
end 


%% OLD:

% 
% %% COMBINE LP DATA 
% % 
% % xy_analysis = vertcat(xy_analysis, xy_analysisA);
% % xy_opto = vertcat(xy_opto, xy_optoA);
% % xy_opto_info = vertcat(xy_opto_info, xy_opto_infoA);
% % 
% % save('220426_LaserPower_C1-C5_Mod.mat', 'xy_analysis', 'xy_opto', 'xy_opto_info');
% % 
% 
% %% ROWS IN XY_OPTO BUT NOT XY_ANALYSIS
% 
% rows_not_found = [];
% 
% for i = 1:634
%     
%     dat = xy_opto_infoA.Date{i};
%     ani = xy_opto_infoA.Animal{i};
%     expp = xy_opto_infoA.Exp{i};
%     Looom = xy_opto_infoA.Loom{i};
%     
%     roww = find(string(xy_analysisA.Date) == dat & string(xy_analysisA.Animal) == ani & string(xy_analysisA.Exp) == expp & string(xy_analysisA.Loom) == Looom);
%     
%     if isempty(roww)
%         rows_not_found = [rows_not_found, i];
%     end 
%       
% end 
% 
% 
% 
% %% ROWS IN XY_ANALYSIS BUT NOT XY_OPTO
% 
% rows_not_found2 = [];
% 
% for i = 1:587
%     
%     dat = xy_analysisA.Date{i};
%     ani = xy_analysisA.Animal{i};
%     expp = xy_analysisA.Exp{i};
%     Looom = xy_analysisA.Loom{i};
%     
%     roww = find(string(xy_opto_infoA.Date) == dat & string(xy_opto_infoA.Animal) == ani & string(xy_opto_infoA.Exp) == expp & string(xy_opto_infoA.Loom) == Looom);
%     
%     if isempty(roww)
%         rows_not_found2 = [rows_not_found2, i];
%     end 
%       
% end 
% 
% 
% %%
% 
% % rows_not_found = [];
% % 
% % for i = 1:358
% %     
% %     dat = xy_opto_info.Date{i};
% %     ani = xy_opto_info.Animal{i};
% %     expp = xy_opto_info.Exp{i};
% %     Looom = xy_opto_info.Loom{i};
% %     
% %     roww = find(string(xy_analysis.Date) == dat & string(xy_analysis.Animal) == ani & string(xy_analysis.Exp) == expp & string(xy_analysis.Loom) == Looom);
% %     
% %     if isempty(roww)
% %         rows_not_found = [rows_not_found, i];
% %     end 
% %       
% % end 
% 
% save('Setd5VGlut2_C5_N2_all.mat', 'xy_analysis', 'xy_opto', 'xy_opto_info');
% 
%