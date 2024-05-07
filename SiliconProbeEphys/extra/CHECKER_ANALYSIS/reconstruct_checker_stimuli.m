function [frame_array, nframes, sq_w]=reconstruct_checker_stimuli(stimuli_log_file)
% checker stimuli presented for the duration of the time t
% each random pattern will be presented for the nframes number of frames
% the random pattern can be black-and-white 'bw', or greyscale: 'uniform' or 'normal'.
% if the random pattern is 'bw', the value to each square will be assigned
% using a uniform distribution. It will be binarized by comparing the value against the threshold.
% the threshold value will be defined by 'fraction',
% which is the number_of_wite_pixels/number_of_wite_pixels
% if the distribution is defined as 'uniform', the value of pixels will be
% drawn from uniform distribution and the fraction value does not apply.
% if the distribution is defined as 'normal', the value of pixels will be
% drawn from normal distribution with the mean shifted to the value
% fraction
% if distribution is defined 'sparse' then the whole screen will be gray
% and the number of random pixels = fraction will be white and the same number
% will be black.
% shift_fraction is the fraction [0..1] of the size of a checker square that
% will be use to shift the checker patter. If 0 no shift will be done
% the idea to move checkers was originally implmented by Jan.

% parameters fadein_t and fadeout_t specify the amount of time spent in
% addition to t, in order to slowly increase and decrease the intensity of
% the stimulus.
%
%06.08.2020
%O. Symonova



[tstart,tend, ifi, t, nframes, ncol, distr, fraction, shift_fraction, channels,   fadeintime, fadeouttime,fadeinframes, stimulusframes, fadeoutframes, s] = ...
    read_checker_stimuli_general_log(stimuli_log_file);
   

nfr=stimulusframes;

texRect_size=[608, 342];
       

sq_w=single(floor(texRect_size(1)/ncol)); %size of the square
shift_size=round(shift_fraction*sq_w);
% nrow=floor(texRect_size(2)/sq_w); % number of checker rows
nrow=max(1,floor(texRect_size(2)/sq_w)); % number of checker rows

rng(s); 

ii=0;
npatt=ceil(nfr/nframes)+10;

frame_array=zeros(nrow*sq_w, ncol*sq_w,npatt,'uint8');
% Loop until a key is pressedScreen('DrawTexture', texid(2), checkerMaskTex,[],[],0,0);    
for i=0:nfr  
    % generate a new checker image every nframes    
    if mod(i,nframes)==0  
        %generate random pattern
        Cgen = checker_image(nrow, ncol, distr, fraction);
        %generate random shift
        if shift_fraction>0
             shift = randi([-shift_size ,shift_size],1,2);
             Cgen=imresize(Cgen,sq_w,'nearest');
             Cgen = circshift(Cgen,shift);             
        end
        frame_array(:,:,ii+1)=uint8(Cgen);
        ii=ii+1;
    end       
end
debug=0;
if debug
    v=VideoWriter('rec_test.avi');
    open(v);
    for i=1:ii
        for j=1:nframes
            writeVideo(v,squeeze(frame_array(:,:,i)));
        end
    end
    close(v);
end
end




function [tstart,tend, ifi, t, nframes, ncol, distr, fraction, shift_fraction, channels,   fadeintime, fadeouttime,fadeinframes, stimulusframes, fadeoutframes, s] = read_checker_stimuli_general_log(stimuli_log_file)
    %% read the file with stimuli parameters
    f=fopen(stimuli_log_file);
    tline=fgetl(f);
    while  ~feof(f) && tline==""
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'checker'))
        disp('The file does not contain information about scanning dot stimulus');
        fclose(f);
        return;
    end

    tline=fgetl(f);
    tstr=tline(length('Start time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    tend=datevec(tstr,formatIn);

    tline=fgetl(f);
    ifi_str=tline(length('Freq: ')+1:end);
    ifi=str2double(ifi_str); %flip rate : frames per second

    tline=fgetl(f);
    t_str=tline(length('t\t')+1:end);
    t=str2double(t_str); %requested time

    tline=fgetl(f);
    res=strfind(tline, 'nframes');
    if isempty(res)
        nfr_str=tline(length('Frames presented:\t')+1:end);
        frames_presented=str2num(nfr_str); %number of frames each random pattern was presented
        tline=fgetl(f);
    end


    nfr_str=tline(length('nframes\t')+1:end);
    nframes=str2num(nfr_str); %number of frames each random pattern was presented

    tline=fgetl(f);
    ncol_str=tline(length('ncol\t')+1:end);
    ncol=str2num(ncol_str); %number of columns

    tline=fgetl(f);
    distr=tline(length('distr\t')+1:end); %name of the distribution

    tline=fgetl(f);
    fraction_str=tline(length('fraction\t')+1:end);
    fraction=str2double(fraction_str); %used in uniform distribution
    
    tline=fgetl(f);
    shift_fraction_str=tline(length('shift fraction\t')+1:end);
    shift_fraction=str2double(shift_fraction_str); %for moving checkers    

    tline=fgetl(f);
    ri_tr=tline(length('maxRedIntensity\t')+1:end);
    maxRedIntensity=str2num(ri_tr); %maximum red intencity

    tline=fgetl(f);
    ch_str=tline(length('channels\t[')+1:end-1);
    ch1=str2num(ch_str(1)); %used in uniform distribution
    ch2=str2num(ch_str(3));
    ch3=str2num(ch_str(5));
    channels=[ch1 ch2 ch3];

    %% for the newer version with fadein and fadeout options
    tline=fgetl(f); %this is either seed structure or fadein time
    fres=strfind(tline, 'Seed');
    if isempty(fres) %this is not seed structure, hence fadein
        fit_str=tline(length('fadein time: ')+1:end);
        fadeintime=str2double(fit_str); %fadein time
        tline=fgetl(f); %fadeout time
        fot_str=tline(length('fadeout time: ')+1:end);
        fadeouttime=str2double(fot_str); 
        tline=fgetl(f); % fadinframes
        fif_str=tline(length('fadein frames presented: ')+1:end);
        fadeinframes=str2num(fif_str); 
        tline=fgetl(f); % stimulus frames
        stf_str=tline(length('stimulus frames presented: ')+1:end);
        stimulusframes=str2num(stf_str);
        tline=fgetl(f); % fadeout frames
        fof_str=tline(length('fadeout frames presented: ')+1:end);
        fadeoutframes=str2num(fof_str);    
    else    
        fadeintime=0;   
        fadeouttime=0; 
        fadeinframes=0;     
        stimulusframes=0;     
        fadeoutframes=0; 
    end 

    %now read info about the random generator seed:
    if isempty(fres)
        tline=fgetl(f);
    end
    tline=fgetl(f);
    seed_type=tline(length('Type\t')+1:end);
    tline=fgetl(f);
    seed_id= str2num(tline(length('Seed\t')+1:end));
    tline=fgetl(f);
    seed_size= str2num(tline(length('State (size)\t')+1:end));
    tline=fgetl(f);
    nn = strsplit(tline, {' '});
    nn= nn(~cellfun('isempty',nn));  
    nnum=uint32(cellfun(@str2num, nn))';

    fclose(f);

    %restore the random seed
    s.Type=seed_type;
    s.Seed = seed_id;
    s.State = nnum;
end
