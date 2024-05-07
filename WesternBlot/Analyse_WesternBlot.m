% Western calculation

wtvals = data(1:3, 3);
hetvals = data(4:6, 3);

% 2 - test normality with Ks test. - more than 30 samples.
% [h,p] =kstest(dataWT)
% [h,p] =kstest(dataHET)

dataWT = wtvals;
dataHET = hetvals;

% Less than 30 samples
[H, pValue] = swtest(dataWT)
[H, pValue] = swtest(dataHET)

% 3 - Levene test for homogenity of variance. 
nWT = numel(dataWT);
nHET = numel(dataHET);

col1 = vertcat(dataWT, dataHET);
col2 = vertcat(ones(nWT,1), ones(nHET,1)*2);

X = horzcat(col1, col2);

%Levene's test for variance:
Levenetest(X)


nanmean(wtvals)
nanmean(hetvals)
[p,h] = ttest2(wtvals, hetvals)
[p,h] = ranksum(wtvals, hetvals)

%% PLOT 

close
semwt = nanstd(wtvals)/sqrt(3); 
semhet = nanstd(hetvals)/sqrt(3); 

figure
bar(1, nanmean(wtvals), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
bar(2, nanmean(hetvals),'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% WT
for i = 1:numel(wtvals)
    xval =  1;
    yval = wtvals(i);
   
    marker = 'k.';
    scatter(xval, yval, 1000,  marker, 'jitter', 'on');
    hold on 
end 

% HET
for i = 1:numel(hetvals)
    xval =  2;
    yval = hetvals(i);
   
    marker = 'r.';
    scatter(xval, yval, 1000,  marker, 'jitter', 'on');
    hold on 
end 

% axis([0 3 -1000 0])
axis([0 3 0 1])
box off
ax = gca;
ax.XAxis.Visible = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.02, 0.02];
ax.LineWidth= 1.5;

errorbar(1, nanmean(wtvals), semwt, 'k', 'LineWidth', 1.5);
errorbar(2, nanmean(hetvals), semhet, 'r', 'LineWidth', 1.5);

f = gcf;
f.Position = [680   750   230   348]; 
 title('Ctx')



