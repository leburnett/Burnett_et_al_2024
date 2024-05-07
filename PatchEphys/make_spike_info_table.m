%% make spike_info_table
% Burnett - 08/08/22

% Script to generate a table with:

% TABLE CONTENTS
% Ani
% Cell
% Geno
% Cohort
% Sweep
% Current Inj
% Spike #
% Voltage (over time)
% dv/ dt
% Peak amp
% AP thresh
% 1/2 width
% Rise time 
% Decay Time
clear all
close all

%% 
files = dir('*.abf');
nfiles = length(files);

ani = 9; 
genoo = 2; 

cohort = 1; 
DTX = 1; 

%% 1 - SET THIS INFO FIRST TO ENSURE DON'T FORGET. 

for i = 6
    
    celll = i;
    % Set the i value for the abf file to open!
    % i = celll;
    
    % Load the file
    fn = files(i).name;
    [d,si,h]=abfload(fn);
    d = squeeze(d);
    n = numel(d(1,:));
    
    % d = data where each column = 1 sweep
    d2 = d';
    % Cut data around point where current is injected.
    d2 = d2(:, 7501:60000);
    n_sweeps = numel(d2(:, 1));
    wid = 1000;
    
    % Initialise table
    spike_info_table = table();
    
    % Go through all the sweeps - find the timing of the peak of the spikes.
    for p = 1:n_sweeps
        
        d3 = d2(p, :);
        % Find when the voltage is > 0mV
        spiks = find(d3>0);
        dsp = diff(spiks);
        % Find when there is a diffference > 1 between the values where the voltage is > 0mV
        dsp_rows = find(dsp>1);
        
        if ~isempty(dsp_rows)
            
            % THe number of spikes - based on these blocks of > 0mV
            n_sp = numel(dsp_rows)+1;
            
            % Find the timing of the beginning and end of the V > 0mV
            % Col 1 = START
            % Col 2 = END
            spike_info = cell(n_sp, 19);
            for jj = 1:n_sp
                
                % Col 1 - Animal
                spike_info(jj,1) = {ani};
                % Col 2 - Cell
                spike_info(jj,2) = {celll};
                % Col 3 - Geno
                spike_info(jj,3) = {genoo};
                % Col 4 - Cohort
                spike_info(jj,4) = {cohort};
                % Col 5 - Sweep
                spike_info(jj,5) = {p};
                % Col 6 - SpikeNum
                spike_info(jj,6) = {jj};
                
                % Col 7 = time when voltage > 0mV
                if jj == 1
                    spike_info(1,7) = {spiks(1)};
                elseif jj == n_sp
                    spike_info(n_sp,7) = {spiks(dsp_rows(n_sp-1)+1)};
                else
                    spike_info(jj, 7) = {spiks(dsp_rows(jj-1)+1)};
                end
                
                % Col 8 - time when voltage < 0mV
                if jj == n_sp
                    spike_info(jj,8) = {spiks(end)};
                else
                    spike_info(jj, 8) = {spiks(dsp_rows(jj))};
                end
                
                % Col 9 - How many 0.02ms units are above 0mV?
                spike_info(jj, 9) = {(diff([spike_info{jj,7}, spike_info{jj,8}]))*0.02};
                
                % Col 10 = MAX voltage
                spike_info(jj,10) = {max(d3(spike_info{jj,7}:spike_info{jj,8}))};
                
                % Col 11 - Time point when voltage was max during this spike
                time_during_spike = find(d3(spike_info{jj,7}:spike_info{jj,8}) == max(d3(spike_info{jj,7}:spike_info{jj,8})));
                if numel(time_during_spike)>1
                    time_during_spike = time_during_spike(1);
                end
                t_full = spike_info{jj,7}+time_during_spike;
                spike_info(jj,11) = {t_full};
                
                if (t_full-wid) <0 || (t_full+wid)>52500
                    spike_info(jj,12) = {zeros(1, 2001)};
                    spike_info(jj,13) = {zeros(1, 2000)};
                    spike_info(jj,14) = {NaN};
                    spike_info(jj,15) = {NaN};
                    spike_info(jj,16) = {NaN};
                    spike_info(jj,17) = {NaN};
                    spike_info(jj,18) = {DTX};
                    spike_info(jj, 19) = {NaN};
                else
                    
                    volt_data = d3(t_full-wid:t_full+wid); % 20ms either side (1000*0.02)
                    dv_data = diff(volt_data)/0.02;
                    
                    % Col 12 = Voltage data - around 1 spike
                    spike_info(jj,12) = {volt_data};
                    
                    % Col 13 = dvdt data - around 1 spike
                    spike_info(jj,13) = {dv_data};
                    
                    % % % % VARIABLES
                    
                    % Col 14 - AP THRESHOLD - mV when dv/dt > 20mV/S.
                    t_more20 = find(dv_data(500:1000)>20);
                    
                    if isempty(t_more20) % Sometimes small spikes don't reach this threshold
                        spike_info(jj, 14) = {NaN};
                        spike_info(jj, 15) = {NaN};
                        spike_info(jj, 16) = {NaN};
                        spike_info(jj, 17) = {NaN};
                        spike_info(jj, 18) = {DTX};
                        spike_info(jj, 19) = {NaN};
                    else
                        t20 = t_more20(1)+500;
                        APthresh = volt_data(t20);
                        spike_info(jj, 14) = {APthresh};
                        
                        % Col 15 - Rise time - time from threshold mV to peak mV
                        peakt = wid;
                        riset = (peakt-t20)*0.02;
                        spike_info(jj, 15) = {riset};
                        
                        % Col 16 - Width half peak
                        hp_val = ((spike_info{jj,10} - APthresh)/2);
                        hp = APthresh + hp_val;
                        vd2 = volt_data- hp;
                        d_sign = sign(vd2(500:1500));
                        
                        tUP = find(diff(d_sign)==2);
                        if isempty(tUP)
                            tt0 = find(diff(d_sign)==1);
                            if isempty(tt0)
                                spike_info(jj, 16) = {NaN};
                                spike_info(jj, 17) = {NaN};
                                spike_info(jj, 18) = {DTX};
                                spike_info(jj, 19) = {NaN};
                                return
                            else
                                tUP = tt0(1);
                            end
                            
                        elseif numel(tUP)>1
                            tUP = tUP(1);
                        end
                        
                        tDOWN = find(diff(d_sign)==-2);
                        if isempty(tDOWN)
                            tt0 = find(diff(d_sign)==-1);
                            tDOWN = tt0(1);
                        elseif numel(tDOWN)>1
                            tDOWN = tDOWN(end);
                        end
                        whp = (tDOWN-tUP)*0.02;
                        spike_info(jj, 16) = {whp};
                        
                        % Col 17 - Decay Time - from peak to momst negative value.
                        % From peak to most negative value.
                        
                        mostneg = find(volt_data(1000:1500)==min(volt_data(1000:1500)));
                        if numel(mostneg)>1
                            mostneg = mostneg(1);
                        end
                        decayt = mostneg*0.02;
                        spike_info(jj, 17) = {decayt};
                        
                        % Col 18 - DTX?
                        spike_info(jj,18) = {DTX};
                        
                        % Col 19 - Min aplitude (1000-1500) - how
                        % hyperpolarised does the cell get?
                        % 1000 = peak
                        minAmp = min(volt_data(1000:1500));
                        spike_info(jj, 19) = {minAmp};
                        
                    end
                end
                
                
            end
            
            sp_table = cell2table(spike_info, "VariableNames", ["Ani", "Cell", "Geno", "Cohort", "Sweep", "SpikeN", "T0St", "T0End", "Tabove0", "MaxAmp", "TMaxFull","V", "dvdt", "APThresh", "RiseT", "whp", "DecayT", "DTX", "MinAmp"]);
            
            % Add the info from this sweep to the table
            spike_info_table = vertcat(spike_info_table, sp_table);
            
        end
    end
    
    % savefold = '/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/SPIKE_INFO_TABLES';
    savefold = '/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/Ptchd1/Analysis/Spike_Info_Tables';
%     savefold = '/Users/lauraburnett/Data_Analysis_Mac/PATCH/DATA/Cul3/Spike_Info_Tables';
    save(fullfile(savefold, strcat('Spike_Info_Table_M', string(ani),'_Cell', string(celll), '_Cohort', string(cohort),'_Geno', string(genoo),'_DTX', string(DTX), '.mat')), 'spike_info_table');
    
    
end
%%

% 
% 
% %% PHASE PLOTS
% Vdata = spike_info_table.V;
% DVDT = spike_info_table.dvdt;
% figure; plot(nanmean(Vdata(:, 1:2000)), nanmean(DVDT))
% 
% %% AVERAGE AP 
% figure; plot(1:1:2001, nanmean(Vdata))
% 
% 

%% Combine all individual files into one giant table.

spike_table = table();

files = dir('Spike_Info_Table_*');
n_files = numel(files);

for k = 1:n_files
   fname = files(k).name;
   load(fname, 'spike_info_table')
    
   spike_table = vertcat(spike_table, spike_info_table);
    
end 

save('221107_Spike_Table_Ptchd1_N8_withDTX.mat', 'spike_table');



%% ADDING TO THE TABLE AFTERWARDS

% 220902 - Adjusting APT - from > 20mV/ms to > 10mV/ms
% Also adding peak d2v/d2t rise and falls. 
wid = 1000;

for jj = 1:height(spike_table)
 
                % Col 12 = Voltage data - around 1 spike
               volt_data =  spike_table{jj,12};
               
                 % Col 13 = dvdt data - around 1 spike
                dv_data = spike_table{jj,13};
                
                % % % % VARIABLES
                
                % Col 14 - AP THRESHOLD - mV when dv/dt > 20mV/S.
                t_more10 = find(dv_data(500:1000)>10);
                
                if isempty(t_more10) % Sometimes small spikes don't reach this threshold
                  
                    % IGNORE CASE 
                    
                else
                    t10 = t_more10(1)+500;
                    APthresh = volt_data(t10);
                    spike_table.APThresh(jj) = APthresh;
                   
                    % Col 15 - Rise time - time from threshold mV to peak mV
                    peakt = wid; % The data is already aligned to the peak of the AP at 1000. 
                    riset = (peakt-t10)*0.02;
                    spike_table.RiseT(jj) = riset;
                    
                    % Col 16 - Width half peak
                    hp_val = ((spike_table.MaxAmp(jj) - APthresh)/2);
                    hp = APthresh + hp_val;
                    vd2 = volt_data- hp;
                    d_sign = sign(vd2(500:1500));
                    
                    tUP = find(diff(d_sign)==2);
                    if isempty(tUP)
                        tt0 = find(diff(d_sign)==1);
                        if isempty(tt0)
                                % IGNORE
                            return
                        else
                            tUP = tt0(1);
                        end
                        
                    elseif numel(tUP)>1
                        tUP = tUP(1);
                    end
                    
                    tDOWN = find(diff(d_sign)==-2);
                    if isempty(tDOWN)
                        tt0 = find(diff(d_sign)==-1);
                        tDOWN = tt0(1);
                    elseif numel(tDOWN)>1
                        tDOWN = tDOWN(end);
                    end
                    whp = (tDOWN-tUP)*0.02;
                    spike_table.whp(jj) = whp;
                    % PEAK d2v/d2t
                    
                    % Look at only the RISE - from APT to MAX
                    d2vd2t_RISE = diff(dv_data(t10:1000));
                    maxd2_rise = max(d2vd2t_RISE);
                    spike_table.maxd2_rise(jj) = maxd2_rise;
                    
                    d2vd2t_FALL = diff(dv_data(900:1500));
                    maxd2_fall = min(d2vd2t_FALL);
                    spike_table.maxd2_fall(jj) = maxd2_fall;
                    
                     % Col 17 - Decay Time - from peak to momst negative value.
                    % From peak to most negative value.
                    
                    mostneg = find(volt_data(1000:1250)==min(volt_data(1000:1250)));
                    if numel(mostneg)>1
                        mostneg = mostneg(1);
                    end
                    decayt = mostneg*0.02;
                    spike_table.DecayT(jj) = decayt;

                end 
end 

save('221104_Cul3_SPIKE_TABLE_NEW_RHEOBASE.mat', 'spike_table');


















