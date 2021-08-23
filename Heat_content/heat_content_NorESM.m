
% heat capcitiy 3850 J/(kg C)
Heat_Capacity = 3850;
% heat content (Joules)  = mass (volume* mean density = 1024 kg m3) * sw specific heat capacity * Temperature
% volume = grid point area * thickness

%import temperature NorESM series and unpack
varname = 'templvl';
ocn_mask=load('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected_lat_bands.mat');
ocn_mask=ocn_mask.domain_mask_lat_bands;

parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
area=parea';

addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered
addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
addpath /Volumes/LaCie_Leonardo/NorESM

depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';
depthbds=ncread('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/N1850AERCNOC_f19_g16_CTRL_02.micom.hbgcm.1000-12.nc','depth_bnds');

%% preparing mask %% working on a mask
grid_filename='/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc';
plat=ncread(grid_filename, 'plat'); 
plon=ncread(grid_filename, 'plon'); 

plat1=plat(200:end,:);
plat2=plat(1:199,:);
new_plat=[plat1 ; plat2]';

plon1=plon(200:end,:);
plon2=plon(1:199,:);
new_plon=[plat1 ; plat2]';

%#####reshapping MAP GRID so atlantic is in the middle... 
domain_mask=load('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat');
domain_mask = domain_mask.domain_mask';
domain_mask(domain_mask==8)=2;
pcolor(domain_mask)
    
%% Heat content time series

tropN_series=NaN(481,1);
tropS_series=NaN(481,1);
subtropN_series=NaN(481,1);
subtropS_series=NaN(481,1);


for time=1:481
    tic
    dummy = NaN(384, 320, 70);
    sprintf('model year = %d', time)
    
    for k = 1:70
        load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_templvl_k_%d',k))

        for l=1:384
                for c=1:320 
                
                if length(new_series{l,c})>1
                    dummy(l,c,k) = new_series{l,c}(time);
                else 
                    dummy(l,c,k) = new_series{l,c}(1);
                end
                
                end
        end
    end
    
 %% gridcell volumes = gridpoint area *thickness    
if time == 1    %then get the volume of gridcells
   
    for k=1:size(depthbds,2)
    thickness(k)=depthbds(2,k)-depthbds(1,k);
    end

    %calculating grid point volumes and doing volume mask
    vol = NaN(size(area,1),size(area,2),size(thickness,2));

    for k = 1:size(thickness,2)
        for line= 1:size(area,1)
            for col = 1:size(area,2)

                if sum(~isnan(dummy(line,col,:))) ~=0

                vol(line,col,k) = area(line,col)*thickness(k); %[m3]
                else
                vol(line,col,k) = NaN;
                end
            end

        end
    end
end

%% do heat content
Heat_content= NaN(384,320,70);
% heat content (Joules)  = mass (volume* mean density = 1024 kg m3) * sw specific heat capacity * Temperature
for l=1:384
    for c=1:320
        for k=1:70
            if ~isnan(dummy(l,c,k))
            
                Heat_content(l,c,k) = vol(l,c,k)*1024*Heat_Capacity*dummy(l,c,k);
            else
                Heat_content(l,c,k) = NaN;
            end
        end
    end
end


%% lat bands

tropN = NaN(384,320,70); %this is 0-30N
subtropN = NaN(384,320,70); %this is 30-60N
tropS = NaN(384,320,70); %this is 0-30S
subtropS = NaN(384,320,70);%this is 30-60S


for l=1:384
    for c =1:320
        for k = 1:70
            
        if domain_mask(l,c) ==2 && new_plat(l,c) >=0 && new_plat(l,c) <=30 %this is 0-30N
            tropN(l,c,k) = Heat_content(l,c,k);
            
        elseif domain_mask(l,c) ==2 && new_plat(l,c) >=30 && new_plat(l,c) <=60 %this is 30-60N
            subtropN(l,c,k) = Heat_content(l,c,k);
        
        elseif domain_mask(l,c) ==2 && new_plat(l,c) <=0 && new_plat(l,c) >=-30 %this is 0-30S   
            tropS(l,c,k) = Heat_content(l,c,k);
        
        elseif domain_mask(l,c) ==2 && new_plat(l,c) <=-30 && new_plat(l,c) >=-60 %this is 30-60S
            subtropS(l,c,k) = Heat_content(l,c,k);
            
        end
        
        end
        
    end
end

tropN_series(time,1) = nansum(tropN,'all');
subtropN_series(time,1) = nansum(subtropN,'all');
tropS_series(time,1) = nansum(tropS,'all');
subtropS_series(time,1) = nansum(subtropS,'all');
toc
end

save('/Volumes/LaCie_Leonardo/NorESM/Heat_content/NAtl_heat_content_series.mat',...
       'tropN_series', 'subtropN_series', 'tropS_series', 'subtropS_series')


%% figure - plot panel A 4 lines

load ('/Volumes/LaCie_Leonardo/NorESM/Heat_content/NAtl_heat_content_series.mat')

tropN_series = tropN_series/10e22;
tropS_series = tropS_series/10e22;
subtropN_series = subtropN_series/10e22;
subtropS_series = subtropS_series/10e22;

pchange_tropN = ((tropN_series(end)-tropN_series(1))/tropN_series(1))*100;
pchange_tropS = ((tropS_series(end)-tropS_series(1))/tropS_series(1))*100;
pchange_subtropN = ((subtropN_series(end)-subtropN_series(1))/subtropN_series(1))*100;
pchange_subtropS = ((subtropS_series(end)-subtropS_series(1))/subtropS_series(1))*100;

delta_tropN = tropN_series(end)-tropN_series(1);
detla_tropS = tropS_series(end)-tropS_series(1);
delta_subtropN = subtropN_series(end)-subtropN_series(1);
delta_subtropS = subtropS_series(end)*6-subtropS_series(1)*6;

figure
plot(1:481,tropN_series, 'LineWidth',2)
hold on
plot(1:481,subtropN_series,'LineWidth',2)
hold on
plot(1:481,tropS_series,'LineWidth',2)
hold on
plot(1:481,subtropS_series*6,'LineWidth',2)
hold on
l2=xline(140,'--','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1, 'FontSize',12);
l3=xline(280,'--','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1, 'FontSize',12);
hold on
lgd = legend('NAtl (0-30\circN)','NAtl (30-60\circN)','SAtl (0-30\circS)','SAtl (30-60\circS)');
lgd.Location = 'southeast';
lgd.FontSize = 10;
hold on
ylabel({'Heat Content (10^2^2 Joules)';""})
xlabel('Model year')
ylim ([10 35])
grid on
grid minor
set(gca,'FontSize',14) % Creates an axes and sets its FontSize to 14
set(gcf,'color','w')
print('/Volumes/LaCie_Leonardo/NorESM/Heat_content/heat_content_evolution.png','-dpng','-r300')

    % a) North altantic between 0-30N and 30-60N
    
    % b) South atlantic between 0-30S and 30-60S