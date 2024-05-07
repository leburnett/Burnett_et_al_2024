% Normalise Laser Power across animals. 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Remove 20Hz trials. 
xy_opto(89:108, :) = []; 
xy_opto_info(89:108, :) = []; 
xy_analysis(89:108, :) = []; 

% Remove MJ0580 - not Cre+ve
all80 = find(xy_opto_info.Animal == "MJ0580");
xy_opto(all80,:) = [];
xy_opto_info(all80, :) = []; 
xy_analysis(all80, :) = []; 

% Freezes
all80 = find(xy_opto_info.Animal == "MJ1415");
xy_opto(all80,:) = [];
xy_opto_info(all80, :) = []; 
xy_analysis(all80, :) = []; 

% Remove vals of EC powers with few trials/ animals. 
vals = find(cell2mat(xy_opto_info.EC) == 2.15 | cell2mat(xy_opto_info.EC) == 2.3 | cell2mat(xy_opto_info.EC) == 2.9 | cell2mat(xy_opto_info.EC) == 3.75 | cell2mat(xy_opto_info.EC) == 4);
xy_opto(vals, :)= []; 
xy_opto_info(vals, :) = [];
xy_analysis(vals, :) = []; 

vals = find(cell2mat(xy_opto_info.EC) == 0 | cell2mat(xy_opto_info.EC) == 1.75 | cell2mat(xy_opto_info.EC) == 1.85 | cell2mat(xy_opto_info.EC) == 2.7);
xy_opto(vals, :)= []; 
xy_opto_info(vals, :) = [];
xy_analysis(vals, :) = []; 

% vals = find(cell2mat(xy_opto_info.EC) == 2.1 | cell2mat(xy_opto_info.EC) == 2.2 | cell2mat(xy_opto_info.EC) == 2.4 | cell2mat(xy_opto_info.EC) == 3.5);
% xy_opto(vals, :)= []; 
% xy_opto_info(vals, :) = [];
% xy_analysis(vals, :) = []; 


%% PLOT GRAPHS SHOWING THE LOGISITIC CURVES FOR EACH MOUSE. 

% Maybe need to rethink how I am calculating 'Return to Shelter'  - relook at how I calculate it... 

% 25 = Return to shelter. 
% 28 - MaxSp
  column = 25; %MaxSp 
  
for i = 1:n_animals
      ani = all_animals(i);
%       all_ani = find(string(xy_analysis.Animal) == ani ); & cell2mat(xy_analysis.MaxSpEscape) >20
       all_ani = find(string(xy_analysis.Animal) == ani);
      xy_analysis2 = xy_analysis(all_ani, :);
      geno = string(xy_analysis2.Geno(1));
      
      all_EC = unique(cell2mat(xy_analysis2.EC));
      n_EC = numel(all_EC);
      
      xvals = []; 
      yvals = []; 
      
      for j = 1:n_EC
          EC_val = all_EC(j,1);
          EC_rows = find(cell2mat(xy_analysis2.EC) == EC_val);
          n_ECrows = numel(EC_rows);
          
          dt = cell2mat(xy_analysis2{EC_rows, column});
          val = mean(dt);
          sem_val = std(dt)/sqrt(numel(dt));
          all_EC(j,2) = val;
          all_EC(j,3) = sem_val;
          all_EC(j,4) = numel(EC_rows);
          
          xvals = [xvals, ones(1,n_ECrows)*EC_val];
          yvals = [yvals, dt']; 
      end
      
      [Qpre, p, sm, varcov] = fit_logistic(xvals,yvals);
      
      figure
      subplot(1,2,1)
      if geno == "wt"
          errorbar(all_EC(:,2)', all_EC(:,3)', 'k', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 28)
      elseif geno == "het"
          errorbar(all_EC(:,2)', all_EC(:,3)', 'r', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 28)
      end
      box off
      hold on
      xticks(1:1:numel(all_EC(:,1)))
      xticklabels(string(all_EC(:,1)))
      %  axis([0 5 -0.1 1.1])
      axis([0 numel(all_EC(:,1))+1 0 1.1])
%       ylabel('Speed - cm/s')
      ax = gca;
      ax.FontSize = 18;
      xlabel('Arbitrary Laser Power')

      subplot(1,2,2)
      % Add logistic curve.
      plot(xvals,yvals,'o')
       plot(xvals,Qpre)    % best fitting logistic
      % Add text values of 'p'
      text(0.5,0.9,string(p))
       ax = gca;
      ax.FontSize = 18;
      xlabel('Arbitrary Laser Power')
      axis([0 4.5 0 1.1])
      
      sgtitle((ani))
end 
  

%% CREATE AN ARRAY - animal - all three values of p. 

% p is 3 element vector containing parameters describing the logistic:
%       thalf, Qinf, and alpha

%   Q(t) = Qinf/(1 + exp(-alpha*(t-thalf)))

%       thalf is symmetric inflection point
%       Qinf is value as t --> infinity
%       alpha is time decay constant

%  sm is 3 element vector giving 1-sigma confidence limits of parameters
%       e.g., thalf = p(1) +/- sm(1)
%       simply double the values in sm to get 95% confidence limits

%   varcov is the complete 3x3 variance-covariance matrix for investigating
%       how model paramters co-vary with each other.
%       sm is sqrt(diag(varcov))
%
%   Example:
%       Qinf = 10.2; alpha = 0.33; thalf = 108.5;
%       t = 100:120;
%       Q = Qinf./(1 + exp(-alpha*(t-thalf)));
%       noise = randn(1,length(t));
%       Qd = Q+noise;
%       Qpre = fit_logistic(t,Qd);
%       figure(1)
%       clf
%       hold on
%       plot(xvals,yvals,'o')  % data
%       plot(xvals,Qpre)    % best fitting logistic


column = 36; %MaxSp

norm_array = zeros(n_animals, 7); 

for i = 1:n_animals
ani = all_animals(i);
ani2 = all_animals{i};
ani_num = str2double(ani2(3:end)); 

all_ani = find(string(xy_analysis.Animal) == ani);
xy_analysis2 = xy_analysis(all_ani, :);
geno = string(xy_analysis2.Geno(1));

all_EC = unique(cell2mat(xy_analysis2.EC));
n_EC = numel(all_EC);

xvals = [];
yvals = [];

for j = 2:n_EC
    EC_val = all_EC(j,1);
    EC_rows = find(cell2mat(xy_analysis2.EC) == EC_val);
    n_ECrows = numel(EC_rows);
    
    dt = cell2mat(xy_analysis2{EC_rows, column});
    val = mean(dt);
    sem_val = std(dt)/sqrt(numel(dt));
    all_EC(j,2) = val;
    all_EC(j,3) = sem_val;
    all_EC(j,4) = numel(EC_rows);
    
    xvals = [xvals, ones(1,n_ECrows)*EC_val];
    yvals = [yvals, dt'];
end

[Qpre, p, sm, varcov] = fit_logistic(xvals,yvals);

if isreal(sm(1))
norm_array(i, 1) = ani_num; 
norm_array(i, 2) = p(1);
norm_array(i, 3) = sm(1);
norm_array(i, 4) = p(2);
norm_array(i, 5) = sm(2);
norm_array(i, 6) = p(3);
norm_array(i, 7) = sm(3);
else
 norm_array(i, 1) = ani_num; 
norm_array(i, 2) = p(1);
norm_array(i, 3) = NaN;
norm_array(i, 4) = p(2);
norm_array(i, 5) = NaN;
norm_array(i, 6) = p(3);
norm_array(i, 7) = NaN;  
end 

end 

%%

Y_VALS = []; 

figure
for i = 1:n_animals
x0 = norm_array(i,2); 
Qinf = norm_array(i,4); 
k = norm_array(i,6); 
x_vals = [1:0.5:3.5];
n = 1; 

for j = x_vals
    x_vals(n) = j;
    y_vals(n) = (Qinf/(1+exp(-k*(j-x0))));
    n = n+1;
end 

if i == 3 || i == 5 || i == 7 || i == 9 || i ==11
    col = 'r';
else 
    col = 'k';
end 

plot(x_vals, y_vals, col)
hold on 

Y_VALS = vertcat(Y_VALS, y_vals);
end 

speed_WT = Y_VALS([1,2,4,6,8], :);
speed_HET = Y_VALS([3,5,7,11], :); 

hv = mean(Y_VALS([3,5,7,11], :));
wv = mean(Y_VALS([1,2,4,6,8], :));

% plot(x_vals, hv, 'r')
% hold on 
% plot(x_vals, wv, 'k')

%%

nWT = 5; 
nHET = 4;

mean_WT = wv; 
mean_HET = hv; 

x = x_vals;

semWT = std(speed_WT)/sqrt(nWT); 
y1 = mean_WT+semWT;
y2 = mean_WT-semWT;
     
semHET = std(speed_HET)/sqrt(nHET); 
y3 = mean_HET+semHET;
y4 = mean_HET-semHET;

figure
plot(x, y1, 'w')
hold on
plot(x, y2, 'w')
patch([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x_vals, (mean_WT)', 'k', 'LineWidth', 1.3)

plot(x, y3, 'w')
hold on 
plot(x, y4, 'w')
patch([x fliplr(x)], [y3 fliplr(y4)],  col, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(x_vals, (mean_HET)', 'Color', col, 'LineWidth', 1.3)


[Qpre, p, sm, varcov] = fit_logistic(x_vals,mean_WT);
%  p = [1.9898, 77.6541, 5.0272];
% sm = [0.0034, 0.2228, 0.0918];

[Qpre, p, sm, varcov] = fit_logistic(x_vals,mean_HET);
%  p = [2.0367, 69.6833, 4.0547];
% sm = [6.44-04, 0.0344, 0.0095];



x0 = p(1); 
Qinf = p(2); 
k = p(3); 
x = 0:0.1:4;

for j = 1:numel(x)
    val = x(j);
    y(j) = (Qinf/(1+exp(-k*(val-x0))));
end 

figure; plot(x,y, 'k')
hold on 

plot(x,y, 'r')









% To normalise divde all EC values by p(1) - since p(1) is the midpoint and
% if you divide it by itself you will get 1!! :) 













