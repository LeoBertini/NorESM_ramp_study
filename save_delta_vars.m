

varlist={'ph';'o2';'templvl'}  ;

for varIDX= 1 %length(varlist)
    tic
    varname = varlist{varIDX};
    sprintf("Doing delta for %s",varname)

    [d.delta1_481, d.year1_field, d.year481_field] = do_delta(varname,'complete'); %calling function
    %[d.delta_1_140, d.year1_field, d.year140_field] = do_delta(varname,'year1',1,'year2',140);

    delta = cell(481,1);
    %extract one layer for year 1
    for year = 1:481
        sprintf('%s year %d', varlist{varIDX}, year)
        tic
        field = NaN(size(d.delta1_481,1),size(d.delta1_481,2),size(d.delta1_481,3));
        for line = 1:size(d.delta1_481,1)
            for col = 1:size(d.delta1_481,2)
                for k = 1: size(d.delta1_481,3)
                    if ~isnan(d.delta1_481{line,col,k}(1))
                        field(line,col,k) = d.delta1_481{line,col,k}(year);
                    else
                        field(line,col,k) = NaN;
                    end
                end
            end
        end
        delta(year,1) = {field};
        toc
    end
    
     addpath  /Volumes/LaCie_Leonardo/NorESM
        
%     figure
%     subplot(2,1,1)
%     pcolor(field(:,:,5)); shading flat
%     title ('Delta 1')
%     colorbar
%     colormap(redblue(20))
%     caxis([-1 1])
% 
%     subplot(2,1,2)
%     pcolor(field(:,:,20)); shading flat
%     title ('Delta 2')
%     colorbar
%     colormap(redblue(20))
%     caxis([-1 1])
    fprintf('Saving data')
    tic
    save(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Binded_var_files/NorESM_deltas_%s.mat',varlist{varIDX}),'delta','-v7.3');
    toc
    fprintf('Data saved')
    
end

