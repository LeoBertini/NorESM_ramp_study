
clear all
addpath /Volumes/LaCie_Leonardo/NorESM/

%% Loading Ocean Mask
domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected_lat_bands.mat');
mask_natl=domain_mask.domain_mask_lat_bands;

figure
pcolor(mask_natl)
shading flat

depth_lvls=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth_lvls=depth_lvls.depths;

%%
varlist={'templvl';'ph';'o2';'AOU';'omegac'};
for varIDX=1:length(varlist)
tic
year=1;

var_field=cell(size(mask_natl,1) ,size(mask_natl,2),size(depth_lvls,1));


for depth=1:size(depth_lvls,1)


ocn_series=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},depth));
ocn_series=ocn_series.new_series;



for line = 1:size(ocn_series,1)
    for col = 1:size(ocn_series,2)
        
        if ~isnan(ocn_series{line,col}(year)) && mask_natl(line,col) == 9 || mask_natl(line,col) == 10 || mask_natl(line,col) == 11
           
            var_field{line,col,depth} = round(ocn_series{line,col}(:),2);
            
        
        else
            var_field{line,col,depth} = NaN;
        end
        
    end
end

end


%% saving binded file 
save(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Binded_var_files/NorESM_filtered_%s.mat',varlist{varIDX}),'var_field','-v7.3');
toc
%takes ~30 min per variable
end


%% unpacking to see if it works
% depth = 50;
% extract=NaN(384,320);
% for line = 1:384
%     for col = 1:320
%         extract(line,col) = var_field{line,col,depth}(1);
%     end
% end
% 
% figure
% pcolor(extract)
% 
% figure
% plot(var_field{300,100,1})
% hold on
% plot(var_field{300,100,10})
% hold on
% plot(var_field{300,100,15})