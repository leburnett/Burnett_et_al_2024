% Opto data - speed at vs speed_im

% Made by Burnett 05/05/23

% DATA
% /Users/lauraburnett/Data_Analysis_Mac/OPTO/Setd5/DATA/220429_Setd5_Opto_LASERPOWER_N8_4WT_4HET_DATA_FOR_PLOTS.mat


%%

allWT = (find(string(xy_opto_infoA.Geno)=="wt")); %& xy_return.Trial <= ntrials));
allHET = (find(string(xy_opto_infoA.Geno)=="het")); 

n_ret = height(xy_opto_infoA);


rowsb4 = 600; %175; 

for i = 1:n_ret
    
    % G = speed data
    G = (xy_optoA(i,:));
    
    sp_at = nanmean(G(rowsb4-15:rowsb4));
    sp_immed = nanmean(G(rowsb4+5:rowsb4+35));% 45
    maxxx = max(G(rowsb4:end));

%     sp_at = nanmean(G(rowsb4-10:rowsb4+5));
%     sp_immed = nanmean(G(rowsb4-5:rowsb4+10));% 45
%     maxxx = max(G(rowsb4:end));

%     xy_opto_infoA.SpAt{i} = sp_at;
%     xy_opto_infoA.SpImmed{i} = sp_immed;
    xy_opto_infoA.maxxx{i} = maxxx;

end 
    
% xy_opto_infoA = xy_opto_info;

%% Trials sorted by L2M group 

all_EC = unique(cell2mat(xy_opto_infoA.EC));


EC_VALS = [1.9, 2, 2.25, 2.5, 2.75, 3];

%%

figure
for i = 1:n_ret
    
    g = xy_opto_infoA.Geno{i};
    %     t = xy_opto_infoA.Trial(i);
    %     d = xy_opto_infoA.Day(i);
    %     maxsp = xy_opto_infoA.MaxSp{i};
    spat = xy_opto_infoA.SpAt{i};
    EC = xy_opto_infoA.EC{i};
    if EC < 3.7 %ismember(EC, EC_VALS) || EC ==2.7 || EC ==3.25 %|| EC ==3.5
        if  spat>3  && spat<20 %&& maxsp > 35
            %     if d==5
            
            if g == "wt"
                
                col = [0.6 0.6 0.6];
                
                if EC <= 1.75 %== EC_VALS(1)
                    subplot(2,6,1);
                elseif EC > 1.75 && EC <=2.1  %EC > 1.75 && EC <2.25%== EC_VALS(2)
                    subplot(2,6,2);
                elseif EC > 2.1 && EC <=2.25 %EC == 2.25 %EC_VALS(3)
                    subplot(2,6,3);
                elseif EC == 2.5 %EC_VALS(4)
                    subplot(2,6,4);
                elseif EC >2.5 && EC<3.0 %== EC_VALS(5) || EC == 2.7  
                    subplot(2,6,5);
                elseif EC >= 3.0 && EC<3.75  %== EC_VALS(6) || EC== 3.25 %|| EC == 3.5
                    subplot(2,6,6);
                else  
                end
                
                y1 = xy_opto_infoA.SpAt{i};
                y2 = xy_opto_infoA.SpImmed{i};
                %     y3 = xy_return.maxxx{i};
                %     plot([1 2 3], [y1 y2 y3], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 1.2)
                plot([1 2], [y1 y2], '-', 'Color', col, 'LineWidth', 0.7, 'Marker', '.', 'MarkerSize', 12) % 'Marker', '.', 'MarkerSize', 30,
                hold on;
                box off
                ylim([-5 70])
                xlim([0.5 2.5])
                xticks([1 2])
                xticklabels({''})
                ax = gca;
                ax.TickDir  = 'out';
                ax.LineWidth = 1.2;
                ax.TickLength = [0.03 0.03];
                if EC >1.9
                    ax.YAxis.Visible = 'off';
                end 
                
            elseif g == "het"
                
                col = 'r'; %[1 0.45 0.13];
                
                if EC <= 1.75 %== EC_VALS(1)
                    subplot(2,6,7);
                elseif EC > 1.75 && EC <=2.1 %EC > 1.75 && EC <2.25%== EC_VALS(2)
                    subplot(2,6,8);
                elseif EC > 2.1 && EC <=2.25 %EC == 2.25 %EC_VALS(3)
                    subplot(2,6,9);
                elseif EC == 2.5 %EC_VALS(4)
                    subplot(2,6,10);
                elseif EC >2.5 && EC<3.1 %== EC_VALS(5) || EC == 2.7  
                    subplot(2,6,11);
                elseif EC >= 3.1 && EC<3.75   %== EC_VALS(6) || EC== 3.25 %|| EC == 3.5
                    subplot(2,6,12);
                else
                end
                
                y1 = xy_opto_infoA.SpAt{i};
                y2 = xy_opto_infoA.SpImmed{i};
                
                %     plot([1 2 3], [y1 y2 y3], '-', 'Marker', '.', 'MarkerSize', 30, 'Color', col, 'LineWidth', 1.2)
                plot([1 2], [y1 y2], '-', 'Color', col, 'LineWidth', 0.7, 'Marker', '.', 'MarkerSize', 12)
                hold on;
                box off
                ylim([-5 70])
                xlim([0.5 2.5])
                xticks([1 2])
                ax = gca;
                ax.TickDir  = 'out';
                ax.LineWidth = 1.2;
                ax.TickLength = [0.03 0.03];
                xticklabels({''})
                if EC >1.9
                    ax.YAxis.Visible = 'off';
                end 
                
            end
        end
    end
    
end


f = gcf;
f.Position = [259   436   720   352];%[432   111   692   549]; 



%% 6 comparisons

%   elseif EC == EC_VALS(5) || EC == 2.7
%   elseif EC == EC_VALS(6) || EC ==3.25

EC_VALS = [1.9, 2, 2.25, 2.5, 2.75, 3];

EC = EC_VALS(1);

rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)==EC & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);

% rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)>2.65 & cell2mat(xy_opto_infoA.EC)<2.8 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);

spat = cell2mat(xy_opto_infoA{rows, 14});
spim = cell2mat(xy_opto_infoA{rows, 15});

[h,p] = ttest2(spat, spim)
n = numel(spim)

nanmean(spat)
nanmean(spim)


%% Add column to xy_opto_info about whether speedIM is </>/= 1SD from speedAT

allWT = find(string(xy_opto_infoA.Geno)=="wt" ); 
avWTspeed = nanmean(cell2mat(xy_opto_infoA.SpAt(allWT)));
stdWT = nanstd(cell2mat(xy_opto_infoA.SpAt(allWT)));

minWT = avWTspeed-(stdWT/2);
maxWT = avWTspeed+(stdWT/2);

% minWT = avWTspeed-(stdWT);
% maxWT = avWTspeed+(stdWT);

allHET = find(string(xy_opto_infoA.Geno)=="het"); 
avHETspeed = nanmean(cell2mat(xy_opto_infoA.SpAt(allHET)));
stdHET = nanstd(cell2mat(xy_opto_infoA.SpAt(allHET)));

minHET = avHETspeed-(stdHET/2);
maxHET = avHETspeed+(stdHET/2);

% minHET = avHETspeed-(stdHET);
% maxHET = avHETspeed+(stdHET);

for i = 1:height(xy_opto_infoA)
    
    geno = string(xy_opto_infoA.Geno{i});
    spim = xy_opto_infoA.SpImmed{i};
    
    if geno == "wt"
        
        if spim < minWT
            xy_opto_infoA.GRP{i} = 1;
        elseif spim >= minWT && spim <= maxWT
            xy_opto_infoA.GRP{i} = 2;
        elseif spim > maxWT
            xy_opto_infoA.GRP{i} = 3;
        end 
        
    elseif geno == "het"
        
        if spim < minHET
            xy_opto_infoA.GRP{i} = 1;
        elseif spim >= minHET && spim <= maxHET
            xy_opto_infoA.GRP{i} = 2;
        elseif spim > maxHET
            xy_opto_infoA.GRP{i} = 3;
        end 
    end 
    
end 



%% STACKED BAR CHART for each EC val

EC_VALS = [1.9, 2, 2.25, 2.5, 2.75, 3];

VALSW = []; 
VALSH = []; 

for i = 1:6 
    
    for j = 1:2
        
        if j == 1
            % WT
            
            if i == 1
                rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)<1.75 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,1);
            elseif i ==2
                rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>1.9 & cell2mat(xy_opto_infoA.EC)<=2.1 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,2);
            elseif i == 3
                  rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>2.1 & cell2mat(xy_opto_infoA.EC)<2.3 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,3);
            elseif i == 4
                rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,4);
            elseif i == 5
                rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>2.5 & cell2mat(xy_opto_infoA.EC)<3.0 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,5);
            elseif i == 6
                rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>=3.0 & cell2mat(xy_opto_infoA.EC)<3.75 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,6);
            end
            
            n = numel(rows);
            all1 = find(cell2mat(xy_opto_infoA.GRP(rows))==1);
            all2 = find(cell2mat(xy_opto_infoA.GRP(rows))==2);
            all3 = find(cell2mat(xy_opto_infoA.GRP(rows))==3);
            
            vals = [numel(all1)/n, numel(all2)/n, numel(all3)/n];
            
            VALSW = vertcat(VALSW, vals);
            
        elseif j == 2
            
            if i == 1
                rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)<1.75 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,1);
            elseif i ==2
                rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)>1.8 & cell2mat(xy_opto_infoA.EC)<2.25 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,2);
            elseif i == 3
                rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)>2.1 & cell2mat(xy_opto_infoA.EC)<2.3 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,3);
            elseif i == 4
                rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,4);
            elseif i == 5
                rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)>2.5 & cell2mat(xy_opto_infoA.EC)<3.1 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,5);
            elseif i == 6
                rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)>3.1 & cell2mat(xy_opto_infoA.EC)<3.75 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
%                 subplot(2,6,6);
            end
            
            n = numel(rows);
            all1 = find(cell2mat(xy_opto_infoA.GRP(rows))==1);
            all2 = find(cell2mat(xy_opto_infoA.GRP(rows))==2);
            all3 = find(cell2mat(xy_opto_infoA.GRP(rows))==3);
            
            vals = [numel(all1)/n, numel(all2)/n, numel(all3)/n];
            
            VALSH = vertcat(VALSH, vals);
            
        end
    end 
 
end 

        figure
        subplot(2,1,1) 
        b1 = bar(VALSW,'stacked');
        b1(1).LineWidth = 1.2;
        b1(2).LineWidth = 1.2;
        b1(3).LineWidth = 1.2;
        b1(1).FaceColor = [0.2 0.2 0.2]; %[0.39 0 0.6];
        b1(2).FaceColor = [0.8 0.8 0.8];
        b1(3).FaceColor = [1 1 1]; %[0.3 0.7 0.9];
        ax = gca;
        ax.LineWidth = 1.25;
        ax.TickDir = 'out';
        ax.TickLength = [0.02 0.02];
       
        box off
        subplot(2,1,2)
        b2 = bar(VALSH,'stacked');  
        b2(1).FaceColor = [0.2 0.2 0.2 ]; %[0.39 0 0.6];
        b2(2).FaceColor = [0.8 0.8 0.8];
        b2(3).FaceColor = [1 1 1]; %[0.3 0.7 0.9];
        b2(1).LineWidth = 1.2;
        b2(2).LineWidth = 1.2;
        b2(3).LineWidth = 1.2;
        box off
        ax = gca;
        ax.LineWidth = 1.25;
        ax.TickDir = 'out';
        ax.TickLength = [0.02 0.02];
        
        f = gcf;
        f.Position = [440   104   420   694]; 
        
        
        %% Make stacked horizontal bar bar charts - per EC val - WT versus HET. 
        
        close
        i = 6; 
        
        VALS = vertcat(VALSH(i,:), VALSW(i,:)); 
        
        figure; 
        b1 = barh(VALS, 'stacked');
        b1(1).LineWidth = 1.2;
        b1(2).LineWidth = 1.2;
        b1(3).LineWidth = 1.2;
        b1(1).FaceColor = [0.2 0.2 0.2]; %[0.39 0 0.6];
        b1(2).FaceColor = [0.8 0.8 0.8];
        b1(3).FaceColor = [1 1 1]; %[0.3 0.7 0.9];
        ax = gca;
        box off
        ax.LineWidth = 1.5;
        ax.TickDir = 'out';
        ax.TickLength = [0.04 0.04];
        ax.YAxis.Visible = 'off';
        ax.XTick = [0, 0.5, 1];
        f = gcf;
        f.Position = [680   918   193   180]; 
   

%% STATS - NEW GROUPS 
% 
% rows = find(string(xy_opto_infoA.Geno)=="het" & cell2mat(xy_opto_infoA.EC)<1.75 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);
rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>=3.0 & cell2mat(xy_opto_infoA.EC)<3.7 & cell2mat(xy_opto_infoA.SpAt)>3 & cell2mat(xy_opto_infoA.SpAt)<20);

spat = cell2mat(xy_opto_infoA{rows, 14}); % 17
spim = cell2mat(xy_opto_infoA{rows, 15}); % 18

% [h,p] = ttest2(spat, spim)
[p, h] = ranksum(spat, spim)
n = numel(spim)

nanmean(spat)
nanmean(spim)

%                 rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)<1.75 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);

%                 rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>1.9 & cell2mat(xy_opto_infoA.EC)<2.35 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);

%                 rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);

%                 rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);

%                 rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>2.65 & cell2mat(xy_opto_infoA.EC)<3.0 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);

%                 rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>=3.0 & cell2mat(xy_opto_infoA.EC)<3.3 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);


%%  Stats - distribution / proportion - all LPs

allWT = find(string(xy_opto_infoA.Geno) == "wt");
allHET = find(string(xy_opto_infoA.Geno) == "het");

wtvals = cell2mat(xy_opto_infoA.GRP(allWT));
hetvals = cell2mat(xy_opto_infoA.GRP(allHET));

w1 = repmat("WT", numel(wtvals),1);
h1 = repmat("HET", numel(hetvals),1);

x1 = [w1', h1']';
x2 = [wtvals; hetvals];
[tbl, chi2stat, pval] = crosstab(x1, x2)
       
%  tbl =
% 
%    110    59   158
%    197    73    67
% 
% 
% chi2stat =
% 
%    62.8077
% 
% 
% pval =
% 
%    2.2987e-14      
%        
% df = 2
       

%% For individual LPs

allWT = find(string(xy_opto_infoA.Geno) == "wt" & cell2mat(xy_opto_infoA.EC)>3 & cell2mat(xy_opto_infoA.EC)<=3.75);
allHET = find(string(xy_opto_infoA.Geno) == "het" & cell2mat(xy_opto_infoA.EC)>3 & cell2mat(xy_opto_infoA.EC)<=3.75);

% allWT = find(string(xy_opto_infoA.Geno) == "wt" & cell2mat(xy_opto_infoA.EC)==2.5);
% allHET = find(string(xy_opto_infoA.Geno) == "het" & cell2mat(xy_opto_infoA.EC)==2.5);

wtvals = cell2mat(xy_opto_infoA.GRP(allWT));
hetvals = cell2mat(xy_opto_infoA.GRP(allHET));

w1 = repmat("WT", numel(wtvals),1);
h1 = repmat("HET", numel(hetvals),1);

x1 = [w1', h1']';
x2 = [wtvals; hetvals];
[tbl, chi2stat, pval] = crosstab(x1, x2)

       
% rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)<1.75 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);
% rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>1.9 & cell2mat(xy_opto_infoA.EC)<2.35 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);
% rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);
% rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)==EC_VALS(i) & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);
% rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>2.65 & cell2mat(xy_opto_infoA.EC)<3.0 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);
% rows = find(string(xy_opto_infoA.Geno)=="wt" & cell2mat(xy_opto_infoA.EC)>=3.0 & cell2mat(xy_opto_infoA.EC)<3.3 & cell2mat(xy_opto_infoA.SpAt)>5 & cell2mat(xy_opto_infoA.SpAt)<25);
% 
%        
       
       
       