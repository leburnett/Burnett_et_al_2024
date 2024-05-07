% Kmeans

% Multiloom
% spikes = all_av_hist(:, 1:3515);
% 
% ncl = 10; 
% 
% k = kmeans(spikes, ncl);
% all_hist(:,3521) = k; 
% 
% figure
% for i = 1:ncl
%     
%     rows = find(k ==i);
%     subplot(ncl, 1, i)
%     plot(mean(spikes(rows, :)), 'Color', [0.3 0.3 0.3]);
%    
% end 
% 
% allwt = numel(find(all_hist(:, 3519)==1));
% allhet = numel(find(all_hist(:, 3519)==0));
% 
% figure
% for i = 1:ncl
%     
%     rows = find(k ==i);
%     
%     nrows = numel(rows);
%     
%     wtvals = find(all_hist(rows, 3519)==1);
%     hetvals = find(all_hist(rows, 3519)==0);
%     
%     subplot(ncl, 1, i)
%     b = bar([numel(wtvals)/allwt, numel(hetvals)/allhet], 'Horizontal', 'on', 'FaceColor', 'Flat');
%     b.CData(2,:) = [0 0 0];
%     b.CData(1,:) = [1 1 1];
%     box off
%     
% end 



% 5 looms

spikes = all_av_hist(:, 1:350);

ncl = 10; 

k = kmeans(spikes, ncl);
all_av_hist(:,356) = k; 

figure
for i = 1:ncl
    
    rows = find(k ==i);
    subplot(ncl, 1, i)
    plot(mean(spikes(rows, :)), 'Color', [0.3 0.3 0.3]);
   
end 

allwt = numel(find(all_av_hist(:, 355)==1));
allhet = numel(find(all_av_hist(:, 355)==0));

% PER GENOTYPE

% figure
% for i = 1:ncl
%     
%     rows = find(k ==i);
%     
%     nrows = numel(rows);
%     
%     wtvals = find(all_av_hist(rows, 355)==1);
%     hetvals = find(all_av_hist(rows, 355)==0);
%     
%     subplot(ncl, 1, i)
%     b = bar([numel(wtvals)/allwt, numel(hetvals)/allhet], 'Horizontal', 'on', 'FaceColor', 'Flat');
%     b.CData(2,:) = [0 0 0];
%     b.CData(1,:) = [1 1 1];
%     box off
%     
% end 

% % % % % % % % % % % % % % % % % % % % % % % % % 

% % PER ANIMAL - Ptchd1
% 
% % het_animals = ["7614", "7476", "7790", "7269","1970", "2833", "1385", "1386", "1388", "1394", "2709"]; 
% 
%   all_animals  = unique(all_av_hist(:, 352));
%   n_animals = numel(all_animals);
%       na1 = numel(find(all_av_hist(:, 352)==all_animals(1)));
%         na2 = numel(find(all_av_hist(:, 352)==all_animals(2)));
%             na3 = numel(find(all_av_hist(:, 352)==all_animals(3)));
%                 na4 = numel(find(all_av_hist(:, 352)==all_animals(4)));
%                     na5 = numel(find(all_av_hist(:, 352)==all_animals(5)));
%                         na6 = numel(find(all_av_hist(:, 352)==all_animals(6)));
%                             na7 = numel(find(all_av_hist(:, 352)==all_animals(7)));
%                                 na8 = numel(find(all_av_hist(:, 352)==all_animals(8)));
% 
% figure
% for i = 1:ncl
%     
%     rows = find(k ==i);
%     nrows = numel(rows);
%     
%     a1 = find(all_av_hist(rows, 352)==all_animals(1));
%     a2 = find(all_av_hist(rows, 352)==all_animals(2));
%     a3 = find(all_av_hist(rows, 352)==all_animals(3));
%     a4 = find(all_av_hist(rows, 352)==all_animals(4));
%     a5 = find(all_av_hist(rows, 352)==all_animals(5));
%     a6 = find(all_av_hist(rows, 352)==all_animals(6));
%     a7 = find(all_av_hist(rows, 352)==all_animals(7));
%     a8 = find(all_av_hist(rows, 352)==all_animals(8));
%    
%     
%     subplot(ncl, 1, i)
%     b = bar([numel(a1)/na1, numel(a2)/na2, numel(a3)/na3, numel(a4)/na4, numel(a5)/na5,numel(a6)/na6, numel(a7)/na7,numel(a8)/na8], 'Horizontal', 'on', 'FaceColor', 'Flat');
%     b.CData(1,:) = [0 0 0];
%     b.CData(2,:) = [0 0 0];
%     b.CData(3,:) = [1 1 1];
%     b.CData(4,:) = [1 1 1];
%     b.CData(5,:) = [1 1 1];
%     b.CData(6,:) = [0 0 0];
%     b.CData(7,:) = [0 0 0];
%     b.CData(8,:) = [1 1 1];
%     box off
%     
% end 


% PER ANIMAL - Setd5

% het_animals = ["7614", "7476", "7790", "7269","1970", "2833", "1385", "1386", "1388", "1394", "2709"]; 

  all_animals  = unique(all_av_hist(:, 352));
  n_animals = numel(all_animals);
      na1 = numel(find(all_av_hist(:, 352)==all_animals(1)));
        na2 = numel(find(all_av_hist(:, 352)==all_animals(2)));
            na3 = numel(find(all_av_hist(:, 352)==all_animals(3)));
                na4 = numel(find(all_av_hist(:, 352)==all_animals(4)));
                    na5 = numel(find(all_av_hist(:, 352)==all_animals(5)));
                        na6 = numel(find(all_av_hist(:, 352)==all_animals(6)));

figure
for i = 1:ncl
    
    rows = find(k ==i);
    nrows = numel(rows);
    
    a1 = find(all_av_hist(rows, 352)==all_animals(1));
    a2 = find(all_av_hist(rows, 352)==all_animals(2));
    a3 = find(all_av_hist(rows, 352)==all_animals(3));
    a4 = find(all_av_hist(rows, 352)==all_animals(4));
    a5 = find(all_av_hist(rows, 352)==all_animals(5));
    a6 = find(all_av_hist(rows, 352)==all_animals(6));
   
    subplot(ncl, 1, i)
    b = bar([numel(a1)/na1, numel(a2)/na2, numel(a3)/na3, numel(a4)/na4, numel(a5)/na5,numel(a6)/na6], 'Horizontal', 'on', 'FaceColor', 'Flat');
    b.CData(1,:) = [0 0 0];
    b.CData(2,:) = [1 1 1];
    b.CData(3,:) = [1 1 1];
    b.CData(4,:) = [0 0 0];
    b.CData(5,:) = [1 1 1];
    b.CData(6,:) = [0 0 0];
    box off
    
end 















%%

% all_animals  = unique(all_av_hist(:, 352));
% n_animals = numel(all_animals);
% 
% for i = 1:n_animals
%   
%     ani = all_animals(i);
%     all_ani = find(all_av_hist(:, 352)==ani);
%     geno = all_av_hist(all_ani(1), 355);
%     if geno == 1
%         col = 'k';
%     else 
%         col = 'r';
%     end 
%     
%     subplot(n_animals,1, i)
%     plot(mean(all_av_hist(all_ani, 1:350)), 'Color', col); 
%     title(ani)
%     
% end 
% 
% 
% all_ani = find(all_av_hist(:, 352)==1386 | all_av_hist(:, 352)==1390);
% all_ani = find(all_av_hist(:, 352)==1970 | all_av_hist(:, 352)==1971 | all_av_hist(:, 352)==7475);
% all_av_hist(all_ani, :) = [];
% 
% save('210901_Spiking_Looms_Setd5_HC_N6.mat', 'all_av_hist');
% 
% 
