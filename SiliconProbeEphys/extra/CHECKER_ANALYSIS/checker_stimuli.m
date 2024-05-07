function checker_stimuli(t, nframes, ncol, distr, fraction, shift_fraction, red_intensity, channels, fadein_t, fadeout_t, outputfile)
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

% parameters fadein_t and fadeout_t specify the amount of time spent in
% addition to t, in order to slowly increase and decrease the intensity of
% the stimulus.


global win texid redmasktex;

%%initialization of default values
if ~exist('t', 'var') || isempty(t)
    t = 120;
end

if ~exist('nframes', 'var') || isempty(nframes)
    nframes = 5;
end

if ~exist('red_intensity', 'var') || isempty(red_intensity)
    red_intensity = 127;
end
if ~exist('ncol', 'var') || isempty(ncol)
    ncol = 20;
end

if ~exist('distr', 'var') || isempty(distr)
    distr = 'uniform';
end

switch distr    
    case 'normal'
        if ~exist('fraction', 'var') || isempty(fraction) || (fraction<1 || fraction>255)
            fraction = 127;
        end
    case 'sparse'
        if ~exist('fraction', 'var') || isempty(fraction) || (fraction<0 || fraction>50)
            fraction = 3;
        end    
    otherwise % black and white with the given fraction
        if ~exist('fraction', 'var') || isempty(fraction) || (fraction<0 || fraction>1)
            fraction = 0.5;
        end        
end

if ~exist('channels', 'var') || isempty(channels)
%     channels=[0 0 1]; %blue
     channels=[0 1 1];  %green and blue     
 elseif length(channels)~=3
    disp('channels should be a three component vector (RGB). Exitting...');
    return;
end

foldername=[];
tstamp =datestr(clock,'YYmmdd_HH_MM_FFF');
if ~exist('outputfile', 'var') || ~ischar(outputfile)
    tstamp_folder =datestr(clock,'YYmmdd');
    foldername=['F:\Data\Tomas\',tstamp_folder];
%     foldername=['C:\DATA\EPhys\stimuli_info\',tstamp_folder];
    if exist(foldername,'dir')==0
        mkdir(foldername)
    end      
    outputfile = [foldername,'\stimuli_',mfilename(), '_',tstamp,'.log'];   
else
    k=strfind(outputfile,'\');
    if isempty(k)
       tstamp_folder =datestr(clock,'YYmmdd');
       foldername=['C:\Data\Tomas\',tstamp_folder];
       if exist(foldername,'dir')==0
           mkdir(foldername)
       end
    else
       foldername=outputfile(1:k(end));
       if exist(foldername,'dir')==0
           mkdir(foldername)
       end
    end
end

if ~exist('shift_fraction', 'var') || isempty(shift_fraction)
    shift_fraction = 0;
end


if ~exist('fadein_t', 'var') || isempty(fadein_t)
    fadein_t = 0;
end
if ~exist('fadeout_t', 'var') || isempty(fadeout_t)
    fadeout_t = 0;
end

ifi = Screen('GetFlipInterval',win);

fadeinscale=1*ifi/fadein_t;
fadeoutscale=1*ifi/fadeout_t;

start_fadeout = t+fadein_t;

fadeinframes=0;
stimulusframes=0;
fadeoutframes=0;
total_t=t+fadein_t+fadeout_t;
%%
winRect = Screen('Rect', win);
texRect = Screen('Rect', texid(2));
texRect_size=[texRect(3)-texRect(1) texRect(4)-texRect(2)]; 



% Sync us to the vertical retrace
vbl = Screen('Flip', win);

% Draw flash texid(2) with black
Screen('FillRect', texid(2), [0 0 0]);

%create the complement of the area where the stimulus will be drawn. We
%will flash red syncing signal in this area
% lut_image=Screen('GetImage', texid(1));
% mask=~(lut_image(:,:,1)&lut_image(:,:,2));
% redMask=zeros(texRect_size(2)*2,texRect_size(1), 4,'uint8');
% redMask(:,:,1)=repelem(mask,2,1);
% redMask(:,:,4)=redMask(:,:,1)*255;
% redMask(:,:,1)=redMask(:,:,1)*red_intensity;
% redMaskTex = Screen('MakeTexture', win, redMask);                 

% We set PTB to wait one frame before re-drawing
waitframes = 1;
frame_count=0;
tstart = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');

sq_w=single(floor(texRect_size(1)/ncol)); %size of the square
shift_size=round(shift_fraction*sq_w);
% nrow=floor(texRect_size(2)/sq_w); % number of checker rows
nrow=max(1,floor(texRect_size(2)/sq_w)); % number of checker rows
ch_w=single(sq_w*ncol); % width of the checher image
ch_h=single(sq_w*nrow); % height of the checher image
% checker image is centered on the full screen, i.e. [608 342]
istart=max(floor((texRect_size(1)-ch_w)/2),1); 
jstart=max(floor((texRect_size(2)-ch_h)/2),1);
% check_image=zeros(texRect_size(2),texRect_size(1),1, 'uint8');
%pre-allocate the image
C=zeros(texRect_size(2),texRect_size(1),3, 'uint8');
%index of the channel where to draw the picture: red, green, or blue
chind = find(channels);

rseed=rng;
vbl_array=zeros([round((t+10)/ifi) 2]);
startfadeout=0;
startstimulus=0;
rf=1;
%farray=[];
t0=GetSecs;
% Loop until a key is pressedScreen('DrawTexture', texid(2), checkerMaskTex,[],[],0,0);    
while GetSecs -t0 < total_t
    %determine how much to dimm checker pattern at fadein and fadeout
    if GetSecs -t0 < fadein_t
        imscale=frame_count*fadeinscale;
    elseif GetSecs -t0 < start_fadeout
            imscale=1;
            if startstimulus==0
                startstimulus=1;
                fadeinframes=frame_count;
            end
    else        
        if startfadeout==0
            startfadeout=1;
            stimulusframes=frame_count-fadeinframes;
        end        
        imscale=max(0,1-(frame_count-stimulusframes-fadeinframes)*fadeoutscale);
    end
    
    % generate a new checker image every nframes    
    if mod(frame_count,nframes)==0   
        if frame_count~=0
            Screen('Close', checkerMaskTex);
        end
        C(:)=0;
        %generate random pattern
        Cgen = checker_image(nrow, ncol, distr, fraction);
        %generate random shift
        if shift_fraction>0
             shift = randi([-shift_size ,shift_size],1,2);
             Cgen=imresize(Cgen,sq_w,'nearest');
             Cgen = circshift(Cgen,shift);
             for ci=1:numel(chind)
                 C(jstart:jstart+ch_h-1,istart:istart+ch_w-1,chind(ci))=Cgen;  
             end
        else
            for ci=1:numel(chind)
                C(jstart:jstart+ch_h-1,istart:istart+ch_w-1,chind(ci))=imresize(Cgen,sq_w,'nearest');  
            end       
        end
    
%         C=cat(3,zeros(size(check_image),'uint8'),zeros(size(check_image),'uint8'),check_image);        
        checkerMaskTex = Screen('MakeTexture', win, uint8(double(C)*imscale));
    end
%    farray = cat(3,farray,C(:,:,2));
    % morph and draw the image
    Screen('DrawTexture', texid(2), checkerMaskTex,[],[],0,0);
    morphtexid = correct_distortion(win, texid);
    Screen('DrawTexture', win, morphtexid,texRect,winRect,[],0);
    
    % draw red surrounding area every 5th frame
    if mod(frame_count,5)==0  
        %activate the blending mode in order to combine two textures: red
        %frames and stimulus itself
        Screen(win,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('DrawTexture', win, redmasktex);      
        Screen(win,'BlendFunction',GL_ONE, GL_ZERO);
    end
         
%         Screen('DrawTexture', win, morphtexid);
%         Screen('DrawTexture', win,texid(2));
    Screen('DrawingFinished', win);
        
    % Flip to the screen on the next vertical retrace
    vbl = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi, 2);
    vbl_array(frame_count+1,1)=vbl;    
    if mod(frame_count,5)==0
        vbl_array(frame_count+1,2)=1;    
    end
    % Now increment the frame counter for the next loop            
    frame_count=frame_count+1;
    Screen('Close', morphtexid);      
end
tend = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');
Screen('Close', checkerMaskTex);

Screen('FillRect', win, [0 0 0]);
vbl = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi);
%Stops LAbView Aquisition To add Everywhere
data_to_send=['|',num2str(2)];

if fadeout_t>0
    fadeoutframes=frame_count-fadeinframes-stimulusframes;
else
    stimulusframes=frame_count-fadeinframes;
end

%%save('farray.mat','farray');

% write stats to the stimuli.log file
fileID = fopen(outputfile,'a');
fprintf(fileID, "\n");
fprintf(fileID, "Stimuli: %s\n", mfilename());
fprintf(fileID, "Start time: %s\n", tstart);
fprintf(fileID, "End time: %s\n", tend);
fprintf(fileID, "Freq: %.4f\n", ifi);
fprintf(fileID, "t\t %d\n", t);
fprintf(fileID, "nframes\t %d\n", nframes);
fprintf(fileID, "ncol\t %d\n", ncol);
fprintf(fileID, "distr\t %s\n", distr);
fprintf(fileID, "fraction\t %.2f\n", fraction);
fprintf(fileID, "shift fraction\t %.2f\n", shift_fraction);
fprintf(fileID, "maxRedIntensity\t %d\n", red_intensity);
fprintf(fileID, "channels\t [%d %d %d]\n", channels(1:3));
fprintf(fileID, "fadein time: %.2f\n", fadein_t);
fprintf(fileID, "fadeout time: %.2f\n", fadeout_t);
fprintf(fileID, "fadein frames presented: %d\n", fadeinframes);
fprintf(fileID, "stimulus frames presented: %d\n", stimulusframes);
fprintf(fileID, "fadeout frames presented: %d\n", fadeoutframes);

fprintf(fileID, "Seed structure:\n");
fprintf(fileID, "Type\t %s\n", rseed.Type);
fprintf(fileID, "Seed\t %d\n", rseed.Seed);
fprintf(fileID, "State (size)\t %d\n", length(rseed.State));
for i=1:length(rseed.State)
    fprintf(fileID, "%d ", rseed.State(i));
end
fprintf(fileID, "\n");
fclose(fileID);

fredname=[foldername,'\','redframes_',mfilename(),'_',tstamp,'.mat'];
save(fredname,'vbl_array');