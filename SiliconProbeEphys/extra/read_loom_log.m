function [nrep, ifi, tloom, tbtw, Speed, MaxR, Pos, loomCol, bgCol] = read_loom_log(logfile)
 %% read the file with stimuli parameters
    f=fopen(logfile);
    tline=fgetl(f);
    while  ~feof(f) && tline==""
        tline=fgetl(f);
    end

%     if ~contains(strfind(tline,'looming'))
%         disp('The file does not contain information about loom stimuli');
%         fclose(f);
%         return;
%     end

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
    nrep_str=tline(length('Number of looms: ')+1:end);
    nrep=str2num(nrep_str); %number of frames presented

    tline=fgetl(f);
    tloom_str=tline(length('Time per loom: ')+1:end);
    tloom=str2double(tloom_str);
    
    tline=fgetl(f);
    tbtw_str=tline(length('tafter: ')+1:end);
    tbtw =str2double(tbtw_str); 
    
    tline=fgetl(f);
    speed_str=tline(length('Speed of expansion: ')+1:end);
    Speed =str2double(speed_str);
     
    tline=fgetl(f);
    MaxR_str=tline(length('Max R: ')+1:end);
    MaxR =str2double(MaxR_str);
    
    tline=fgetl(f);
    pos_str=tline(length('Position of loom(out of [10,20]): ')+1:end);
    Pos =str2double(pos_str);
    
    tline=fgetl(f);
    loomcol_str=tline(length('Loom color: ')+1:end);
    loomCol =str2double(loomcol_str);
    
    tline=fgetl(f);
    bgcol_str=tline(length('Background color: ')+1:end);
    bgCol=str2double(bgcol_str); %number of frames presented

%     tline=fgetl(f);
%     ch_str=tline(length('Channels: ')+1:end);
%     
%     %find openning brackets
%     ob=strfind(ch_str,'[');
%     %find closing brackets
%     cb=strfind(ch_str,']');
%     nchannels=length(ob);
%     channels=zeros(nchannels,3,'uint8');
%     for ci=1:nchannels
%         chistr=ch_str(ob(ci)+1:cb(ci)-1);
%         chicell=strsplit(chistr,' ');
%         channels(ci,:)= uint8(cellfun(@str2num, chicell));
%     end  
    
%     tline=fgetl(f);
%     nfr_str=tline(length('Stimulus frames: ')+1:end);
%     nframes=str2num(nfr_str); %number of frames presented  
end
