
% Analyse the speed at (SpAt) the loom start verus the speed immediately
% (SpIm) afterwards. Give 250ms for detection. 

% Created by Burnett - April 20th 2022

%% Add column to xy_opto_info about whether speedIM is </>/= 1SD from speedAT

allWT = find(string(xy_return.Geno)=="wt" ); 
avWTspeed = nanmean(cell2mat(xy_return.sp_at(allWT)));
stdWT = nanstd(cell2mat(xy_return.sp_at(allWT)));

minWT = avWTspeed-(stdWT/2);
maxWT = avWTspeed+(stdWT/2);


allHET = find(string(xy_return.Geno)=="het"); 
avHETspeed = nanmean(cell2mat(xy_return.sp_at(allHET)));
stdHET = nanstd(cell2mat(xy_return.sp_at(allHET)));

minHET = avHETspeed-(stdHET/2);
maxHET = avHETspeed+(stdHET/2);

% minHET = avHETspeed-(stdHET);
% maxHET = avHETspeed+(stdHET);

for i = 1:height(xy_return)
    
    geno = string(xy_return.Geno{i});
    spim = xy_return.sp_immed{i};
    
    if geno == "wt"
        
        if spim < minWT
            xy_return.GRP{i} = 1;
        elseif spim >= minWT && spim <= maxWT
            xy_return.GRP{i} = 2;
        elseif spim > maxWT
            xy_return.GRP{i} = 3;
        end 
        
    elseif geno == "het"
        
        if spim < minHET
            xy_return.GRP{i} = 1;
        elseif spim >= minHET && spim <= maxHET
            xy_return.GRP{i} = 2;
        elseif spim > maxHET
            xy_return.GRP{i} = 3;
        end 
    end 
    
end 


%% STACKED BAR CHART

VALSW = []; 
VALSH = []; 

    
    for j = 1:2
        
        if j == 1
            % WT
            
            rows = find(string(xy_return.Geno)=="wt" & cell2mat(xy_return.sp_at)>3 & cell2mat(xy_return.sp_at)<25);
            
            n = numel(rows);
            all1 = find(cell2mat(xy_return.GRP(rows))==1);
            all2 = find(cell2mat(xy_return.GRP(rows))==2);
            all3 = find(cell2mat(xy_return.GRP(rows))==3);
            
            vals = [numel(all1)/n, numel(all2)/n, numel(all3)/n];
            
            VALSW = vertcat(VALSW, vals);
            
        elseif j == 2
            

            rows = find(string(xy_return.Geno)=="het" & cell2mat(xy_return.sp_at)>3 & cell2mat(xy_return.sp_at)<25);

            n = numel(rows);
            all1 = find(cell2mat(xy_return.GRP(rows))==1);
            all2 = find(cell2mat(xy_return.GRP(rows))==2);
            all3 = find(cell2mat(xy_return.GRP(rows))==3);
            
            vals = [numel(all1)/n, numel(all2)/n, numel(all3)/n];
            
            VALSH = vertcat(VALSH, vals);
            
        end
    end 
 

        figure
        b1 = bar(vertcat(VALSW, VALSH), 'BarLayout', 'stacked');
        b1(1).LineWidth = 1.2;
        b1(2).LineWidth = 1.2;
        b1(3).LineWidth = 1.2;
        b1(1).FaceColor = [0.2 0.2 0.2]; %[0.39 0 0.6];
        b1(2).FaceColor = [0.8 0.8 0.8];
        b1(3).FaceColor = [1 1 1]; %[0.3 0.7 0.9];
        ax = gca;
        ax.LineWidth = 1.55;
        ax.TickDir = 'out';
        ax.TickLength = [0.03 0.03];
        box off
        f = gcf;
        f.Position = [926   520   255   344]; 
        

%%  Stats - distribution / proportion - all LPs

allWT = find(string(xy_return.Geno) == "wt");
allHET = find(string(xy_return.Geno) == "het");

wtvals = cell2mat(xy_return.GRP(allWT));
hetvals = cell2mat(xy_return.GRP(allHET));

w1 = repmat("WT", numel(wtvals),1);
h1 = repmat("HET", numel(hetvals),1);

x1 = [w1', h1']';
x2 = [wtvals; hetvals];
[tbl, chi2stat, pval] = crosstab(x1, x2)



