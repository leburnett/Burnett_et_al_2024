
%% Analyse LOOM stimulus
% Generated by Burnett - 25/02/22
% Modified by Burnett 10/09/22 - for 2021 data 

% This time - only ONE rep of HC looms (5 x 5 Looms), then 5 x 5 LC looms,
% the 1 x 5 looms at 5 different speeds. 
% Does start with 15 x looms ( Multiloom) 


% Analysis script for generating spiking information from ephys output data.

% data found here:
% '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/5LOOMS/HC/SPIKES';

% Set variables
recsp = 20000; % 20 kHz!
% save_path = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/SPIKESS/5Looms/Loom_Spike_Res';
save_path = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/5LOOMS/HC/SPIKES/RESULTS';
% save_path = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Cul3/5LOOMS/SPIKES/LC/RESULTS'; 

%% % FIND FILES, FIELNAME, MOUSE, DATE, DEPTH

files = dir('*_5Looms*');
n_tot = numel(files);

for jjjj = 64 %65:67 %1:n_tot
        
        filename = files(jjjj).name;
        load(filename)
        
        date = str2num(filename(1:6));
        ani = str2num(filename(10:13));
        test = str2num(filename(19));
        probe_depth = str2num(filename(21:24));
        
        n_gcl = length(loom_spike_t);
        
        %% TEST to see what the spiking looks like - real timing.
        % Combine all 5 repetitions of the looming stimulus at their actual timing.
        
%                 all_loom_spiket = cell(n_gcl,1);
%         
%                 for jjj = 1:n_gcl
%                    s1 = cell2mat(loom_spike_t(jjj, 1));
%                    s2 = cell2mat(loom_spike_t(jjj, 2));
%                    s3 = cell2mat(loom_spike_t(jjj, 3));
%                    s4 = cell2mat(loom_spike_t(jjj, 4));
%                    s5 = cell2mat(loom_spike_t(jjj, 5));
%         
%                    spikz = [s1',s2',s3',s4',s5'];
%                    all_loom_spiket(jjj, 1) = {spikz};
%                 end
        
%                  figure
%                 for i = 1:n_gcl
%                     spiket = cell2mat(all_loom_spiket(i, 1));
%                     nspikes = numel(spiket);
%                     y = ones(nspikes,1)*i;
%                     plot(spiket, y, 'k.');
%                     hold on
%                 end
        %
        
        %% 1 - Process data - get info on spiket for all 5 looms.
        % TIME IS 5.3 SECONDS FOR THE LOOM STIMULUS
        % Get timing of beginning / end of each stimulus and create
        % variable of spike times for each loom bout rep. (5 * Looms)
        
%         n_gcl = numel(loom_spike_t(:,1));
        
        for kk = 1:5
            
            figure
            for i = 1:n_gcl
                spiket = cell2mat(loom_spike_t(i, kk));
                nspikes = numel(spiket);
                y = ones(nspikes,1)*i;
                plot(spiket, y, 'k.');
                hold on
            end
            
            axis tight
            rng = xlim;
            if kk ==1
                x1 = rng(1);
                x1b = x1+recsp*5.3;
                
            elseif kk==2
                x2 = rng(1);
                x2b = x2+recsp*5.3;
                
            elseif kk ==3
                x3 = rng(1);
                x3b = x3+recsp*5.3;
                
            elseif kk ==4
                x4 = rng(1);
                x4b = x4+recsp*5.3;
                
            elseif kk== 5
                x5 = rng(1);
                x5b = x5+recsp*5.3;
                
            end
            
            close
            
        end
        
        
        % Spike timing - for RASTER PLOT
        L1SP =  loom_spike_t(:,1);
        L2SP =  loom_spike_t(:,2);
        L3SP =  loom_spike_t(:,3);
        L4SP =  loom_spike_t(:,4);
        L5SP =  loom_spike_t(:,5);
        
        %         rng_all(1,1) = (x1c-x1)/recsp;
        %         rng_all(2,1) = (x2c-x2)/recsp;
        %         rng_all(3,1) = (x3c-x3)/recsp;
        %         rng_all(4,1) = (x4c-x4)/recsp;
        %         rng_all(5,1) = (x5c-x5)/recsp;
        
        
%         %% Check total time of recording - should be 40s - NEED THE PREVIOUS FIGURE TO BE OPEN!
%         axis tight
%         rng = xlim;
%         x1 = rng(1);
%         x2 = x1+recsp*5.3;
%         
%         % Add vertical lines for where bar changes.
%         xvls = x1+(recsp):recsp*0.73:x2-(recsp*0.27);
%         for i = 1:6
%             plot([xvls(i) xvls(i)], [0 n_gcl], 'r')
%             hold on
%         end
        
        % 19 BINS!!
%         nbins = numel(xvls)-1;

        
        %% Find the spikes per frame
        
        % LOOM1
        txvals = x1:recsp/60:x1b;
        nbinst = numel(txvals)-1;
        
        spikes_per_time1 = zeros(n_gcl, nbinst);
        
        for i = 1:n_gcl
            spiket = cell2mat(loom_spike_t(i, 1));
            spikes_per_time1(i, :) = histcounts(spiket, txvals);
        end
        
        % LOOM2
        txvals = x2:recsp/60:x2b;
        
        spikes_per_time2 = zeros(n_gcl, nbinst);
        
        for i = 1:n_gcl
            spiket = cell2mat(loom_spike_t(i, 2));
            spikes_per_time2(i, :) = histcounts(spiket, txvals);
        end
        
        % LOOM3
        txvals = x3:recsp/60:x3b;
        
        spikes_per_time3 = zeros(n_gcl, nbinst);
        
        for i = 1:n_gcl
            spiket = cell2mat(loom_spike_t(i, 3));
            spikes_per_time3(i, :) = histcounts(spiket, txvals);
        end
        
        % LOOM4
        txvals = x4:recsp/60:x4b;
        
        spikes_per_time4 = zeros(n_gcl, nbinst);
        
        for i = 1:n_gcl
            spiket = cell2mat(loom_spike_t(i, 4));
            spikes_per_time4(i, :) = histcounts(spiket, txvals);
        end
        
        % LOOM5
        txvals = x5:recsp/60:x5b;
        nbinst = numel(txvals)-1;
        
        spikes_per_time5 = zeros(n_gcl, nbinst);
        
        for i = 1:n_gcl
            spiket = cell2mat(loom_spike_t(i, 5));
            spikes_per_time5(i, :) = histcounts(spiket, txvals);
        end
        
        %% Spikes per time - all 25 presentations of the loom
        % L25 - Spikes per frame across  all 25 reps.
        % L25_ST = spike times across all 25 loom reps.
        % L5REP = spikes per frame, average over all 25 reps for 5 REPS.
        
        xvls = 1*60:0.73*60:(5.3-0.27)*60;
        
        for ii = 1:n_gcl
            
            baseline = nanmean(horzcat(spikes_per_time1(ii,1:60), spikes_per_time2(ii,1:60),spikes_per_time3(ii,1:60),spikes_per_time4(ii,1:60),spikes_per_time5(ii,1:60)));
            
            L25(ii, :) = horzcat(spikes_per_time1(ii,:), spikes_per_time2(ii,:),spikes_per_time3(ii,:),spikes_per_time4(ii,:),spikes_per_time5(ii,:));
            L25_norm(ii,:) = L25(ii, :)-baseline;
            
            L25_ST(ii, :) = {horzcat(cell2mat(L1SP(ii))', cell2mat(L2SP(ii))', cell2mat(L3SP(ii))', cell2mat(L4SP(ii))', cell2mat(L5SP(ii))')};
            
            all_reps = vertcat(spikes_per_time1(ii,:), spikes_per_time2(ii,:),spikes_per_time3(ii,:),spikes_per_time4(ii,:),spikes_per_time5(ii,:));
            %             all_reps_norm = vertcat(zscore(spikes_per_time1(ii,:)), zscore(spikes_per_time2(ii,:)), zscore(spikes_per_time3(ii,:)), zscore(spikes_per_time4(ii,:)) ,zscore(spikes_per_time5(ii,:)));
            all_reps_norm = all_reps-baseline;
            
            L5REP(ii, :) = nanmean(all_reps);
            L5REP_norm(ii, :) = nanmean(all_reps_norm);
            
            
            % Find the average across ONE Loom Rep from the average of 5 * 5REPS
            all_looms = vertcat(all_reps(:, xvls(1):xvls(2)-1), all_reps(:, xvls(2):xvls(3)-1), all_reps(:, xvls(3):xvls(4)-1), all_reps(:, xvls(4):xvls(5)-1), all_reps(:, xvls(5):xvls(6)-1));
            all_looms_norm = vertcat(all_reps_norm(:, xvls(1):xvls(2)-1), all_reps_norm(:, xvls(2):xvls(3)-1), all_reps_norm(:, xvls(3):xvls(4)-1), all_reps_norm(:, xvls(4):xvls(5)-1), all_reps_norm(:, xvls(5):xvls(6)-1));
            
            L1REP(ii, :) = nanmean(all_looms);
            L1REP_norm(ii, :) = nanmean(all_looms_norm);
        end
        
        %         figure; imagesc(L1REP); figure; imagesc(L1REP_norm)
        %         figure; imagesc(L5REP); figure; imagesc(L5REP_norm)
        %         figure; imagesc(L25); figure; imagesc(L25_norm)
        
        %% 
      
        dates = repmat(date, n_gcl, 1);
        anis = repmat(ani, n_gcl, 1);
        tests = repmat(test, n_gcl, 1);
        probe_depths = repmat(probe_depth, n_gcl, 1);
        
        
        %% 
                
        save(fullfile(save_path, strcat(string(date), '_', string(ani), '_', string(test), '_', string(probe_depth), '_LOOM.mat')), 'L25_ST', 'dates', 'anis', 'tests', 'probe_depths', 'L25', 'L25_norm', 'L5REP', 'L5REP_norm', 'L1REP', 'L1REP_norm', 'ids_depth');

        clearvars -except recsp save_path files n_tot 
        
        close all
    
end


clear