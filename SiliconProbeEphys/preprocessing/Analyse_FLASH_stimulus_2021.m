% FLASH ANALYSIS from SPIKING data. 
% Generated by Burnett - 18/02/22
% Modified by Burnett 10/09/22 - for 2021 data 

% For the analysis of the flash stimulus in the ephys setup.

% From "spike_timing" data - save a file for each recording - i.e. one  depth for one animal with:

% 2020: 
% 10 REPS
% 1s flash
% 0.5 before and 0.5 after. 
% 20s in total!
% STARTS - 0.5s [0 0 0] then 10* (1s [0 255 255]) then 0.5s [0 0 0]. 

% 2021
% 10 x flashs - of 1.5s, then 1s, then 0.5s, then 0.25s

%% SAVES - PER EXPERIMENT AT ONE DEPTH FOR EACH ANIMAL. 
% 
% spikes_per_stim  - Spikes binned by the ON/OFF times. 
% 
% norm_spikes_per_stim  - as above - zscored. 
% 
% spikes_per_time - spikes per frame
% 
% spikes_rep_av  - spikes per frame averaged over the 10 reps. 
% 
% norm_spikes_per_time -  spikes per frame - zscored
% 
% norm_spikes_rep_av - spikes per frame averaged over the 10 reps - zscored. 

%%
clear
close all

FLASHTIME = "1-5"; 

recsp = 20000; % 20 kHz!

% save_path = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/SPIKESS/Flashes/FLASH_OFO';
% save_path = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/FLASH/SPIKES/RESULTS2'; 
save_path = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Cul3/FLASH/SPIKES/RESULTS'; 

%% % FIND FILES, FIELNAME, MOUSE, DATE, DEPTH

files = dir('*.mat');
n_tot = numel(files);

for jjjj = 1:n_tot
    
        
        filename = files(jjjj).name;
        load(filename)
        
        date = str2num(filename(1:6));
        ani = str2num(filename(10:13));
        test = str2num(filename(19));
        probe_depth = str2num(filename(21:24));
        
        
        %% 1 - Raster plot 
        n_gcl = numel(flash_spike_t);
        n_ids = numel(ids_depth(:,1));
        
        if n_gcl == n_ids
        
        figure
        for i = 1:n_gcl
            spiket = cell2mat(flash_spike_t(i, 1));
            nspikes = numel(spiket);
            y = ones(nspikes,1)*i;
            plot(spiket, y, 'k.');
            hold on
        end
        
        %% Set parameters for different flash lengths
        if FLASHTIME == "1-5"
            % for 1.5s - 5 reps - 0.75s OFF then 1.5s ON - 1.5s off etc.
            % 905 frames, 15.083s
            tb4 = 0.75;
            tfl = 1.5;
            freps = 10;
        elseif FLASHTIME == "1-0"
            tb4 = 0.5;
            tfl = 1;
            freps = 20;
        elseif FLASHTIME == "0-5"
            tb4 = 0.25;
            tfl = 0.5;
            freps = 40;
        elseif FLASHTIME == "0-25"
            tb4 = 0.75;
            tfl = 1.5;
            freps = 80;
        end
        
        %% Check total time of recording - should be 40s - NEED THE PREVIOUS FIGURE TO BE OPEN!
        axis tight
        rng = xlim;
        x1 = rng(1);
        x2 = x1+recsp*(tfl*freps);
        
%         % Add vertical lines for where bar changes.
%         xvls = x1+(recsp*tb4):recsp*tfl:x2-(recsp*tb4);
%         for i = 1:freps
%             plot([xvls(i) xvls(i)], [0 n_gcl], 'b')
%             hold on
%         end
%         
         xvls = x1:recsp*tfl*2:x2;
%         for i = 1:freps
%             plot([xvls(i) xvls(i)], [0 n_gcl], 'g')
%             hold on
%         end
        
        % 19 BINS!!
        nbins = numel(xvls)-1;
        close
        
        %% Bins values by each bar "stimulus" and count the number of spikes per direction - 8 BINS.
        
        spikes_per_stim = zeros(n_gcl, nbins);
        for i = 1:n_gcl
            spiket = cell2mat(flash_spike_t(i, 1));
            spikes_per_stim(i, 1:nbins) = histcounts(spiket, xvls);
        end
        
        % NORMALSE - PER STIMULUS
        norm_spikes_per_stim = zeros(n_gcl, nbins);
        
        for jj = 1:n_gcl
            rho = spikes_per_stim(jj,:);
            rho2 = zscore(rho);
            norm_spikes_per_stim(jj, :) = rho2;
        end
        
        %% Plot HEATMAP of the normalised activity during each bar stimulus.
        %  figure; imagesc(norm_spikes_per_stim)
        
        %% Create AVERAGE response to OFF-ON-OFF for each cell over 10 reps!
        
        % SUM THE NUMBER OF SPIKES IN EACH FRAME - 0.0167s
        txvals = x1:recsp/60:x2;
        nbinst = numel(txvals)-1;
        
        spikes_per_time = zeros(n_gcl, nbinst);
        
        for i = 1:n_gcl
            spiket = cell2mat(flash_spike_t(i, 1));
            spikes_per_time(i, :) = histcounts(spiket, txvals);
        end
        
        %% THEN CHOP INTO CHUNKS for each REP - Flash OFF = in the middle of array. 
        
        nt = numel(spikes_per_time(1,:));
        rept = nt/(freps/2);
        
        spikes_rep_av = zeros(n_gcl, rept); 
        
        for j = 1:n_gcl
            data = spikes_per_time(j, 1:nt);
            data2 = reshape(data, [rept,(freps/2)])'; % Reshape to 10 reps - 120 frames each
            md = mean(data2);
            spikes_rep_av(j, :) = md;
        end


        %% z-scored firing:
        norm_spikes_per_time = zeros(n_gcl, nbinst);
        
        for jj = 1:n_gcl
            
            rho = spikes_per_time(jj,:);
            rho2 = zscore(rho);
            
            norm_spikes_per_time(jj, :) = rho2;
        end
        
        %%
        norm_spikes_rep_av = zeros(n_gcl, rept); %nbinst/10
        
        for j = 1:n_gcl
            data = norm_spikes_per_time(j, 1:nt);
            data2 = reshape(data, [rept,(freps/2)])'; % Reshape to 10 reps - 120 frames each
            md = mean(data2);
            norm_spikes_rep_av(j, :) = md;
        end
        

        %% Combine data

        dates = repmat(date, n_gcl, 1);
        anis = repmat(ani, n_gcl, 1);
        tests = repmat(test, n_gcl, 1);
        probe_depths = repmat(probe_depth, n_gcl, 1);
        
        OFO = [spikes_rep_av, dates, anis, tests, probe_depths, ids_depth(:,2)];
        OFO_NORM = [norm_spikes_rep_av, dates, anis, tests, probe_depths, ids_depth(:,2)];
        
        %% SAVE
        
        save(fullfile(save_path, strcat(string(date), '_', string(ani), '_', string(test), '_', string(probe_depth), '_FLASH_', FLASHTIME, '.mat')), 'OFO', 'OFO_NORM', 'spikes_per_stim','norm_spikes_per_stim','spikes_per_time', 'spikes_rep_av','norm_spikes_per_time','norm_spikes_rep_av');
        
        clearvars -except recsp save_path files n_tot FLASHTIME

        close all
        
        else
            
            disp(filename)
            disp(n_gcl - n_ids)
            
        end 
   
    
end 

clear
















%% GROUP RESPONSES BY KMEANS - HEATMAP

% figure; imagesc(OFO_NORM(:, 1:120)); colormap(redblue)

% n_k = 3; 
% 
% OFO_NORM(:, 126) = kmeans(OFO_NORM(:, 1:120), n_k);
% OFO2 = sortrows(OFO_NORM, 126);
% 
% figure; subplot(1,5,1:4); imagesc(OFO2(:,1:120)); colormap(redblue)
% subplot(1,5,5); imagesc(OFO2(:, 126))
% 
% figure
% for j = 1:n_k 
% all_type = find(OFO2(:, 126)== j);
% subplot(n_k,1,j); plot(nanmean(OFO2(all_type, 1:120)), 'k'); hold on; plot([30 30], [-0.5 2], 'r:'); plot([90 90], [-0.5 2], 'r:'); box off
% end 
%%

%%


% figure; plot(smooth(nanmean(norm_spikes_per_time1)), 'k'); hold on;
%  plot(smooth(nanmean(norm_spikes_per_time2)), 'r');
%  
% %
%   figure; plot((nanmean(spikes_rep_av1)), 'k'); hold on;
%  plot((nanmean(spikes_rep_av2)), 'r');
% 
%  
%  %
% 







