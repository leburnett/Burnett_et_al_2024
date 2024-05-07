function make_log_data_LOOM_BOX()
%Created 11/12/19 by Burnett

%Function used to convert values from log file created during LOOM experiments IN BOX into a table which can then be read out later to recreate the stimulus. 
% TO be called from 'process_dome_data.m'

global inpath exp_name log_data

% inpath = 'C:\Data_analysis\DATA\Setd5\191216\GN3956\01_Loom';

%% Save the vbl_arrays from the .mat files into a table to be accessed later
my_matfiles = dir(fullfile(inpath,'*.mat'));
outname_MAT= strcat('MATDATA_',exp_name,'.mat');
save(string(outname_MAT), 'my_matfiles');

%%
myFiles = dir(fullfile(inpath,'*.log'));
numFiles = length(myFiles);

%% Make Table
log_data = table();

%%
for k = 1:numFiles
f=fopen(myFiles(k).name);

tline=fgetl(f);
ststr=tline(length('Stimuli: ')+1:end);
stim=string(ststr);
log_data.stim(k)=num2cell(stim);

%% Conditional log entry for each stimulus type. 

%% 
if log_data.stim{k}=="CL_green_screen_1S"

tline=fgetl(f);
ts_str=tline(length('Start_green: ')+1:end);
tstart=string(ts_str);
log_data.tstart(k)=tstart;

tline=fgetl(f);
te_str=tline(length('End_green: ')+1:end);
tend=string(te_str);
log_data.tend(k)=tend;

tline=fgetl(f);
dur_str=tline(length('Duration: ')+1:end);
dur = str2double(dur_str);
log_data.duration(k)=dur;

tline=fgetl(f);
RF_str=tline(length('RedFrame: ')+1:end);
RF = str2double(RF_str);
log_data.redframe(k)=RF;
    
tline=fgetl(f);
BkR_str=tline(length('Background Colour R: ')+1:end);
Bkred=str2double(BkR_str);
log_data.BkR(k)=Bkred;

tline=fgetl(f);
BkG_str=tline(length('Background Colour G: ')+1:end);
BkGreen=str2double(BkG_str);
log_data.BkG(k)=BkGreen;

tline=fgetl(f);
BkB_str=tline(length('Background Colour B: ')+1:end);
BkBlue=str2double(BkB_str);
log_data.BkB(k)=BkBlue;

tline=fgetl(f);
exp_str=tline(length('ExpName: ')+1:end);
exp=string(exp_str);
log_data.ExpName(k)=exp;
    

%%
elseif log_data.stim{k}=="CL_looming_1S_TRIGGER_simplified" || log_data.stim{k}=="CL_looming_1S_TRIGGER_ILI"
    
tline=fgetl(f);
ts_str=tline(length('Start time: ')+1:end);
tstart=string(ts_str);
log_data.tstart(k)=tstart;

tline=fgetl(f);
te_str=tline(length('End time: ')+1:end);
tend=string(te_str);
log_data.tend(k)=tend;

tline=fgetl(f);
dur_str=tline(length('Duration: ')+1:end);
dur = str2double(dur_str);
log_data.duration(k)=dur;

tline=fgetl(f);
fr_str=tline(length('Flip rate: ')+1:end);
ifi=str2double(fr_str); 
log_data.ifi(k)=ifi;

tline=fgetl(f);
nloom_str=tline(length('NLooms: ')+1:end);
nloom=str2double(nloom_str);
log_data.NLoom(k)=nloom;

tline=fgetl(f);
tloom_str=tline(length('t_loom: ')+1:end);
tloom=str2double(tloom_str);
log_data.t_loom(k)=tloom;

tline=fgetl(f);
dpf_str=tline(length('Position: ')+1:end);
dpf=string(dpf_str); 
log_data.position(k)=dpf;

tline=fgetl(f);
rad_str=tline(length('Starting radius: ')+1:end);
radius=str2double(rad_str);
log_data.radius(k)=radius;

tline=fgetl(f);
bound_str=tline(length('Boundary: ')+1:end);
bound=str2double(bound_str);
log_data.Boundary(k)=bound;
 
tline=fgetl(f);
sp_str=tline(length('Speed: ')+1:end);
sp=str2double(sp_str);
log_data.Speed(k)=sp;

tline=fgetl(f);
red_str=tline(length('Dot Colour R: ')+1:end);
red=str2double(red_str);
log_data.dotR(k)=red;
   
tline=fgetl(f);
gr_str=tline(length('Dot Colour G: ')+1:end);
green=str2double(gr_str);
log_data.DotG(k)=green;

tline=fgetl(f);
blue_str=tline(length('Dot Colour B: ')+1:end);
blue=str2double(blue_str);
log_data.DotB(k)=blue;

tline=fgetl(f);
BkR_str=tline(length('Background Colour R: ')+1:end);
Bkred=str2double(BkR_str);
log_data.BkR(k)=Bkred;

tline=fgetl(f);
BkG_str=tline(length('Background Colour G: ')+1:end);
BkGreen=str2double(BkG_str);
log_data.BkG(k)=BkGreen;

tline=fgetl(f);
BkB_str=tline(length('Background Colour B: ')+1:end);
BkBlue=str2double(BkB_str);
log_data.BkB(k)=BkBlue;

tline=fgetl(f);
tth_str=tline(length('T_thresh: ')+1:end);  
tth=str2double(tth_str);
log_data.t_thresh(k)=tth;

%New version - like position 191204. 
tline=fgetl(f);
trig_str=tline(length('Trigger: ')+1:end);
trig=string(trig_str); 
log_data.trigger(k)=trig;
 
tline=fgetl(f);
exp_str=tline(length('ExpName: ')+1:end);
exp=string(exp_str);
log_data.ExpName(k)=exp;

elseif log_data.stim{k}=="red_frame_triggered_opto"
    
tline=fgetl(f);
ts_str=tline(length('Start time: ')+1:end);
tstart=string(ts_str);
log_data.tstart(k)=tstart;

tline=fgetl(f);
te_str=tline(length('End time: ')+1:end);
tend=string(te_str);
log_data.tend(k)=tend;

tline=fgetl(f);
dur_str=tline(length('Duration: ')+1:end);
dur = str2double(dur_str);
log_data.duration(k)=dur;

tline=fgetl(f);
wait_str=tline(length('WaitT: ')+1:end);
wt=str2double(wait_str); 
log_data.waitT(k)=wt;

tline=fgetl(f);
fr_str=tline(length('Flip rate: ')+1:end);
ifi=str2double(fr_str); 
log_data.ifi(k)=ifi;

tline=fgetl(f);
red_str=tline(length('Dot Colour R: ')+1:end);
red=str2double(red_str);
log_data.dotR(k)=red;
   
tline=fgetl(f);
gr_str=tline(length('Dot Colour G: ')+1:end);
green=str2double(gr_str);
log_data.DotG(k)=green;

tline=fgetl(f);
blue_str=tline(length('Dot Colour B: ')+1:end);
blue=str2double(blue_str);
log_data.DotB(k)=blue;

tline=fgetl(f);
BkR_str=tline(length('Background Colour R: ')+1:end);
Bkred=str2double(BkR_str);
log_data.BkR(k)=Bkred;

tline=fgetl(f);
BkG_str=tline(length('Background Colour G: ')+1:end);
BkGreen=str2double(BkG_str);
log_data.BkG(k)=BkGreen;

tline=fgetl(f);
BkB_str=tline(length('Background Colour B: ')+1:end);
BkBlue=str2double(BkB_str);
log_data.BkB(k)=BkBlue;

tline=fgetl(f);
tth_str=tline(length('T_thresh: ')+1:end);  
tth=str2double(tth_str);
log_data.t_thresh(k)=tth;

tline=fgetl(f);
trig_str=tline(length('Trigger: ')+1:end);
trig=string(trig_str); 
log_data.trigger(k)=trig;

tline=fgetl(f);
np_str=tline(length('NumPulses: ')+1:end);
np = str2double(np_str);
log_data.NumPulses(k)=np;

tline=fgetl(f);
tp_str=tline(length('T_pulse: ')+1:end);
tp = str2double(tp_str);
log_data.T_pulse(k)=tp;

tline=fgetl(f);
fp_str=tline(length('FreqPulse: ')+1:end);
fp = str2double(fp_str);
log_data.FreqPulse(k)=fp;
 
tline=fgetl(f);
ec_str=tline(length('EC: ')+1:end);
ec = str2double(ec_str);
log_data.EC(k)=ec;
 
tline=fgetl(f);
exp_str=tline(length('ExpName: ')+1:end);
exp=string(exp_str);
log_data.ExpName(k)=exp;

elseif log_data.stim{k}=="red_frame_triggered_opto_LOOM"
    
tline=fgetl(f);
ts_str=tline(length('Start time: ')+1:end);
tstart=string(ts_str);
log_data.tstart(k)=tstart;

tline=fgetl(f);
te_str=tline(length('End time: ')+1:end);
tend=string(te_str);
log_data.tend(k)=tend;

tline=fgetl(f);
dur_str=tline(length('Duration: ')+1:end);
dur = str2double(dur_str);
log_data.duration(k)=dur;

tline=fgetl(f);
wait_str=tline(length('WaitT: ')+1:end);
wt=str2double(wait_str); 
log_data.waitT(k)=wt;

tline=fgetl(f);
d2l_str=tline(length('delay2loom: ')+1:end);
d2l=str2double(d2l_str); 
log_data.delay2loom(k)=d2l;

tline=fgetl(f);
fr_str=tline(length('Flip rate: ')+1:end);
ifi=str2double(fr_str); 
log_data.ifi(k)=ifi;

tline=fgetl(f);
nloom_str=tline(length('NLooms: ')+1:end);
nloom=str2double(nloom_str);
log_data.NLoom(k)=nloom;

tline=fgetl(f);
tloom_str=tline(length('t_loom: ')+1:end);
tloom=str2double(tloom_str);
log_data.t_loom(k)=tloom;

tline=fgetl(f);
dpf_str=tline(length('Position: ')+1:end);
dpf=string(dpf_str); 
log_data.position(k)=dpf;

tline=fgetl(f);
rad_str=tline(length('Starting radius: ')+1:end);
radius=str2double(rad_str);
log_data.radius(k)=radius;

tline=fgetl(f);
bound_str=tline(length('Boundary: ')+1:end);
bound=str2double(bound_str);
log_data.Boundary(k)=bound;
 
tline=fgetl(f);
sp_str=tline(length('Speed: ')+1:end);
sp=str2double(sp_str);
log_data.Speed(k)=sp;

tline=fgetl(f);
red_str=tline(length('Dot Colour R: ')+1:end);
red=str2double(red_str);
log_data.dotR(k)=red;
   
tline=fgetl(f);
gr_str=tline(length('Dot Colour G: ')+1:end);
green=str2double(gr_str);
log_data.DotG(k)=green;

tline=fgetl(f);
blue_str=tline(length('Dot Colour B: ')+1:end);
blue=str2double(blue_str);
log_data.DotB(k)=blue;

tline=fgetl(f);
BkR_str=tline(length('Background Colour R: ')+1:end);
Bkred=str2double(BkR_str);
log_data.BkR(k)=Bkred;

tline=fgetl(f);
BkG_str=tline(length('Background Colour G: ')+1:end);
BkGreen=str2double(BkG_str);
log_data.BkG(k)=BkGreen;

tline=fgetl(f);
BkB_str=tline(length('Background Colour B: ')+1:end);
BkBlue=str2double(BkB_str);
log_data.BkB(k)=BkBlue;

tline=fgetl(f);
tth_str=tline(length('T_thresh: ')+1:end);  
tth=str2double(tth_str);
log_data.t_thresh(k)=tth;

tline=fgetl(f);
trig_str=tline(length('Trigger: ')+1:end);
trig=string(trig_str); 
log_data.trigger(k)=trig;

tline=fgetl(f);
np_str=tline(length('NumPulses: ')+1:end);
np = str2double(np_str);
log_data.NumPulses(k)=np;

tline=fgetl(f);
tp_str=tline(length('T_pulse: ')+1:end);
tp = str2double(tp_str);
log_data.T_pulse(k)=tp;

tline=fgetl(f);
fp_str=tline(length('FreqPulse: ')+1:end);
fp = str2double(fp_str);
log_data.FreqPulse(k)=fp;
 
tline=fgetl(f);
ec_str=tline(length('EC: ')+1:end);
ec = str2double(ec_str);
log_data.EC(k)=ec;
 
tline=fgetl(f);
exp_str=tline(length('ExpName: ')+1:end);
exp=string(exp_str);
log_data.ExpName(k)=exp;


elseif log_data.stim{k}=="CL_looming_1S_TRIGGER_DiffContrast" % NEW DIFF CONTRASTS STIMULUS - 17/05/22

        
tline=fgetl(f);
ts_str=tline(length('Start time: ')+1:end);
tstart=string(ts_str);
log_data.tstart(k)=tstart;

tline=fgetl(f);
te_str=tline(length('End time: ')+1:end);
tend=string(te_str);
log_data.tend(k)=tend;

tline=fgetl(f);
dur_str=tline(length('Duration: ')+1:end);
dur = str2double(dur_str);
log_data.duration(k)=dur;

tline=fgetl(f);
fr_str=tline(length('Flip rate: ')+1:end);
ifi=str2double(fr_str); 
log_data.ifi(k)=ifi;

tline=fgetl(f);
nloom_str=tline(length('NLooms: ')+1:end);
nloom=str2double(nloom_str);
log_data.NLoom(k)=nloom;

tline=fgetl(f);
tloom_str=tline(length('t_loom: ')+1:end);
tloom=str2double(tloom_str);
log_data.t_loom(k)=tloom;

tline=fgetl(f);
dpf_str=tline(length('Position: ')+1:end);
dpf=string(dpf_str); 
log_data.position(k)=dpf;

tline=fgetl(f);
rad_str=tline(length('Starting radius: ')+1:end);
radius=str2double(rad_str);
log_data.radius(k)=radius;

tline=fgetl(f);
bound_str=tline(length('Boundary: ')+1:end);
bound=str2double(bound_str);
log_data.Boundary(k)=bound;
 
tline=fgetl(f);
sp_str=tline(length('Speed: ')+1:end);
sp=str2double(sp_str);
log_data.Speed(k)=sp;

tline=fgetl(f);
red_str=tline(length('Dot Colour R: ')+1:end);
red=str2double(red_str);
log_data.dotR(k)=red;
   
tline=fgetl(f);
gr_str=tline(length('Dot Colour G: ')+1:end);
green=str2double(gr_str);
log_data.DotG(k)=green;

tline=fgetl(f);
blue_str=tline(length('Dot Colour B: ')+1:end);
blue=str2double(blue_str);
log_data.DotB(k)=blue;

tline=fgetl(f);
BkR_str=tline(length('Background Colour R: ')+1:end);
Bkred=str2double(BkR_str);
log_data.BkR(k)=Bkred;

tline=fgetl(f);
BkG_str=tline(length('Background Colour G: ')+1:end);
BkGreen=str2double(BkG_str);
log_data.BkG(k)=BkGreen;

tline=fgetl(f);
BkB_str=tline(length('Background Colour B: ')+1:end);
BkBlue=str2double(BkB_str);
log_data.BkB(k)=BkBlue;

tline=fgetl(f);
CV_str=tline(length('ContrastValue: ')+1:end);  
CV=str2double(CV_str);
log_data.ContrastValue(k)=CV;

tline=fgetl(f);
tth_str=tline(length('T_thresh: ')+1:end);  
tth=str2double(tth_str);
log_data.t_thresh(k)=tth;

%New version - like position 191204. 
tline=fgetl(f);
trig_str=tline(length('Trigger: ')+1:end);
trig=string(trig_str); 
log_data.trigger(k)=trig;
 
tline=fgetl(f);
exp_str=tline(length('ExpName: ')+1:end);
exp=string(exp_str);
log_data.ExpName(k)=exp;

end

fclose(f);

end

%% Find out how long each stimulus ran for and add to log_data
T1 = 0;
log_data.t_cumul(1)=0;

for p = 1:numFiles 
tstart = log_data.tstart(p);
formatIn='yyyy/mm/dd HH:MM:SS:FFF';
tstart=datevec(tstart,formatIn);
tstart = duration(tstart(4:6));

tend = log_data.tend(p);
formatIn='yyyy/mm/dd HH:MM:SS:FFF';
tend=datevec(tend,formatIn);
tend = duration(tend(4:6));

time_new = tend-tstart;
time_exp = seconds(time_new);
log_data.exp_dur(p)= time_exp;

T1 = T1 +time_new;
log_data.t_cumul(p)=seconds(T1);
end

% save(string(outname), 'log_data');
save(strcat('LOGDATA_', exp_name, '.mat'), 'log_data');
end

