

clear all
close all

%extracting  values for the surface layers for each year..

parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
parea=parea';

addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered
addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
addpath /Volumes/LaCie_Leonardo/NorESM

domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected_lat_bands.mat');
domain_mask=domain_mask.domain_mask_lat_bands;
pcolor(domain_mask)

%Code for latitudinal bands in mask
% tropical = 9
% subtropical = 10
% subpolar = 11

numerator=nan(384/2,320); %arrays of NaNs where values are appended to when there's a value to be added
denominator=nan(384/2,320);


varlist={'templvl','ph', 'o2', 'AOU', 'omegac'};
varIDX=5;

ocn_series=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_1.mat',varlist{varIDX}));
ocn_series=ocn_series.new_series;

lat_code = {"tropical","subtropical","subpolar"};
loop_idx = 1;

for lat_band = [9 10 11]   
    
    sprintf("Latitudinal Band is %s", lat_code{loop_idx})

    for year=1:481 
       for l=384/2:384 
        for c=1:320

            if  domain_mask(l,c)==lat_band && length(ocn_series{l,c})==481

                numerator(l,c)=ocn_series{l,c}(year)*parea(l,c);
                denominator(l,c)=parea(l,c);    
            end
        end

       end

       %storing things in a structure

       if lat_band == 9
       seriesSURF.tropical(year)=nansum(numerator(:))/nansum(denominator(:));

       elseif lat_band == 10
              seriesSURF.subtropical(year)=nansum(numerator(:))/nansum(denominator(:));

       elseif lat_band == 11
              seriesSURF.subpolar(year)=nansum(numerator(:))/nansum(denominator(:));
       end 

        numerator=nan(384/2,320); %clear contents 
        denominator=nan(384/2,320); 
    end
    loop_idx = loop_idx + 1;
end

%%

%doing mean for the whole water column and depths

numerator2ocn=nan(384/2,320,70);
denominator2ocn=nan(384/2,320,70);  
depthbds=ncread('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/N1850AERCNOC_f19_g16_CTRL_02.micom.hbgcm.1000-12.nc','depth_bnds');

for k=1:length(depthbds)
thickness(k)=depthbds(2,k)-depthbds(1,k);
end

loop_idx = 1;

for lat_band = [9 10 11]   

tic
    for year=1:481
        for k=1:70
        sprintf("Variable is = %s, Lat band code = %s ; Year = %d ; depth level = %d", varlist{varIDX}, lat_code{loop_idx}, year, k)
        
            ocn_series=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},k));
            ocn_series=ocn_series.new_series;

            for l=384/2:384
                for c=1:320

                    if  domain_mask(l,c) == lat_band && length(ocn_series{l,c}) == 481
                        numerator2ocn(l,c,k) = ocn_series{l,c}(year)*parea(l,c)*thickness(k);
                        denominator2ocn(l,c,k)=parea(l,c)*thickness(k);
                    end

                end
            end
            
            %volume weighted time series evolution per depth (fig 2 plot)
            if lat_band == 9
                ocn_volseries_depth.tropical(k,year)=nansum(numerator2ocn(:,:,k),'all')/nansum(denominator2ocn(:,:,k),'all');
            elseif lat_band == 10
                ocn_volseries_depth.subtropical(k,year)=nansum(numerator2ocn(:,:,k),'all')/nansum(denominator2ocn(:,:,k),'all');
            elseif lat_band == 11
                ocn_volseries_depth.subpolar(k,year)=nansum(numerator2ocn(:,:,k),'all')/nansum(denominator2ocn(:,:,k),'all');
            end
            
        end

        %volume weighted timeseries (single line - fig 1 plot)
        if lat_band == 9
            ocn_volseries.tropical(year)=nansum(numerator2ocn(:,:,:),'all')/nansum(denominator2ocn(:,:,:),'all');
        elseif lat_band == 10
            ocn_volseries.subtropical(year)=nansum(numerator2ocn(:,:,:),'all')/nansum(denominator2ocn(:,:,:),'all');
        elseif lat_band == 11
            ocn_volseries.subpolar(year)=nansum(numerator2ocn(:,:,:),'all')/nansum(denominator2ocn(:,:,:),'all');
        end

        
    end
    
    numerator2ocn=nan(384/2,320,70);
    denominator2ocn=nan(384/2,320,70); 
    loop_idx = loop_idx +1;
end
toc



save(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_%s_lat_bands',varlist{varIDX}),...
    'seriesSURF','ocn_volseries','ocn_volseries_depth' )
