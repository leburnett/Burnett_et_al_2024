% Combine 'full stimulus' analysis

% LOAD table with information about the depth of the sSC in different animals. 
% load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/220628_Setd5_Rec_Depth_Info.mat', 'animal_depth_info_table')

%Setd5
% het_animals = [7269, 7476, 7614];

% Setd5
% %spiking data 
% dir1 = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Full_Stimulus'; 
% files1 = dir(fullfile(dir1, '*FULL*'));
% nfiles =numel(files1);
% 
% dir2 = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/Full_Stimulus';
% files2 =  dir(fullfile(dir2, '*FULL*'));
% 
%% 
Full_Table = table(); 

for i = 1:nfiles
   
    fname = fullfile(files1(i).folder, files1(i).name);
    load(fname)

    % Ammend the depth of each cell! 
    animal_num =fname(end-19:end-16); 
    row = find(animal_depth_info_table.Animal == str2num(animal_num));
    sSC_depth = animal_depth_info_table.SC_DEPTH_NEW(row); % depth of beginning of sSC
    probe = full_stim_spikes.probe_depths(1); %probe depth
    n_gcl = numel(full_stim_spikes(:,1)); %number of units in this recording
    
    for jj = 1:n_gcl
        d = probe - (800 - full_stim_spikes.id_depth(jj));
        d2 = sSC_depth*-1 - d;
        full_stim_spikes.Depth(jj) = d2;
        
        if ismember(str2num(animal_num), het_animals)
            full_stim_spikes.Geno = ones(n_gcl, 1)*2;
        else
            full_stim_spikes.Geno = ones(n_gcl, 1);
        end 
        
        max_val = max(full_stim_spikes.SPT(jj, :));
        min_val = min(full_stim_spikes.SPT(jj, :));
        rng = max_val-min_val;
        
        full_stim_spikes.max(jj) = max_val;
        full_stim_spikes.min(jj) = min_val;
        full_stim_spikes.rng(jj) = rng;
        
    end 
    
    Date = full_stim_spikes.dates;
    Ani = full_stim_spikes.anis;
    ID_NUM = full_stim_spikes.id_num;
    Depth = full_stim_spikes.Depth;
    MeanSPF = full_stim_spikes.mean_spf;
    MeanSPT = full_stim_spikes.mean_spt;
    Geno = full_stim_spikes.Geno;
    Test = full_stim_spikes.tests;
    Max = full_stim_spikes.max;
    Min = full_stim_spikes.min;
    Rng = full_stim_spikes.rng;
    
    % Load vals over Acclim stimulus at the beginning
    fname2 = fullfile(files2(i).folder, files2(i).name);
    load(fname2)
    
    MeanSPF_ACC = full_stim_spikes.mean_spf;
    MeanSPT_ACC = full_stim_spikes.mean_spt;
    
    tbl = table(ID_NUM, Date, Ani, Test, Geno, Depth, MeanSPF, MeanSPT, MeanSPF_ACC, MeanSPT_ACC, Max, Min, Rng);
   
    Full_Table = vertcat(Full_Table, tbl);

end 

save('220803_Full_Table.mat', 'Full_Table')

%% Ptchd1
% Spike timing data
% Modfiied by Burentt - 13/09/2022 for PTchd1/ Cul3 data

load('/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/220912_AniDepth_Ptchd1.mat', 'animal_depth_info_table')

% Ptchd1
het_animals = [1385, 1386, 1394, 2709, 4369];

dir1 = '/Users/lauraburnett/Data_Analysis_Mac/Silicon_Probe/RESULTS_NewCMAP/Ptchd1/FULL';
files1 = dir(fullfile(dir1, '*Full_Acclim_Start*'));
nfiles =numel(files1);

recsp = 20000; 

Full_Table = table();

for i = 1:nfiles
    
    fname = fullfile(files1(i).folder, files1(i).name);
    load(fname)
    expname = files1(i).name;
    
    % Ptchd1
    animal_num = expname(10:13);
    date = str2num(expname(1:6));
    test = str2num(expname(19));
    row = find(animal_depth_info_table.Ani == str2num(animal_num) & animal_depth_info_table.Date == date);
    if ~isempty(row)
        
        sSC_depth = animal_depth_info_table.Depth(row); % depth of beginning of sSC
        probe = str2num(expname(21:24));
        
        if ismember(str2num(animal_num), het_animals)
            genoo = 2;
        else
            genoo = 1;
        end
        
        n_gcl = length(full_spike_t(:,1)); %number of units in this recording
        
        figure
        for ii = 1:n_gcl
            spiket = cell2mat(full_spike_t(ii, 1));
            nspikes = numel(spiket);
            y = ones(nspikes,1)*ii;
            plot(spiket, y, 'k.');
            hold on
        end
        
        axis tight
        rng = xlim;
        x1 = rng(1);
        x2 = rng(2);
        
        xvls = x1:recsp/60:x2;
        nbins = numel(xvls)-1;
        close
        
        spikes_full = NaN(n_gcl, 1805);
        DATES = [];
        TESTT =[];
        ANINUM = [];
        DEPTH = [];
        TOTALSP = [];
        GENOO = [];
        MAXVAL = [];
        AVVAL= [];
        RANGEVAL = [];
        
        for jj = 1:n_gcl % Run through units
            
            % Count spikes in bins  - 1 frame (1/60s) bins.
            spiket = cell2mat(full_spike_t(jj, 1));
            data = histcounts(spiket, xvls);
            
            spikes_full(jj, 1:nbins) = data;
            
            d = probe - (800 - ids_depth(jj,2));
            d2 = sSC_depth - d;
            
            DATES(jj,1) = date;
            TESTT(jj,1) = test;
            ANINUM(jj,1) = str2num(animal_num);
            DEPTH(jj,1) = d2;
            TOTALSP(jj,1) = sum(spikes_full(jj, :));
            GENOO(jj,1) = genoo;
            
            max_val = max(data)*60;
            av_val = nanmean(data)*60;
            rng = range(data)*60;
            
            MAXVAL(jj,1) = max_val;
            AVVAL(jj,1) = av_val;
            RANGEVAL(jj,1) = rng;
        end
        
        
        Date = DATES;
        Test = TESTT;
        Ani = ANINUM;
        ID_NUM = ids_depth(:,1);
        Depth =  DEPTH;
        SPT = spikes_full;
        Geno =  GENOO;
        Max = MAXVAL;
        Av = AVVAL;
        Rng = RANGEVAL;
        
        tbl = table(ID_NUM, Date, Ani, Test, Geno, Depth, SPT, Max, Av, Rng);
        
        Full_Table = vertcat(Full_Table, tbl);
    end
    
end

    
   save('220913_Full_Table_Ptchd1.mat', 'Full_Table', 'animal_depth_info_table');
 
    
    
    
    
    