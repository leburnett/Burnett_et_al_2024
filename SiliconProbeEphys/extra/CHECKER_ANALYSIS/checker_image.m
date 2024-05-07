function Cgen = checker_image(ncol,nrow, distr, fraction)
switch distr
    case 'uniform'        
        Cgen=uint8(rand(ncol, nrow)*255);        
    case 'normal'        
%         Cgen=uint8(randn(ncol, nrow)*255+fraction);  
        Cgen=randn(ncol, nrow);  
        minC=min(min(Cgen));
        maxC=max(max(Cgen));
        Cgen=uint8((Cgen-minC)/(maxC-minC)*255+(fraction-127.5));
    case 'sparse'
        indrange=ncol*nrow;
        samples=[];
        for i=1:2*fraction
            resample=1;
            while resample
                ind=max(1,min(round(rand*indrange),indrange));
                if ~ismember(ind,samples)
                    samples=cat(1,samples,ind);
                    resample=0;
                end
            end                    
        end
        Cgen=ones(ncol, nrow, 'uint8')*127;
        for i=1:fraction
               [I_row, I_col] = ind2sub(size(Cgen),samples(i));
                Cgen(I_row, I_col)=255;
                [I_row, I_col] = ind2sub(size(Cgen),samples(end-i+1));
                Cgen(I_row, I_col)=0;
        end
    otherwise % black and white with the given fraction        
        Cgen=rand(ncol, nrow);
        Cgen(Cgen>fraction)=255;
        Cgen(Cgen<=fraction)=0;
end