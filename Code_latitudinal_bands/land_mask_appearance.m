
clear all
close all

grid_filename='/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc';
addpath(genpath('/Volumes/LaCie_Leonardo/NorESM'))
atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');


%% ADJUSTING LANDMASK APPEARANCE
pmask=ncread(grid_filename, 'pmask');
pmask=double(pmask); %converting information to double to pass NaN conditional
pmask(pmask==0)=NaN; %replacing continents (where there is 0, with NaN) to appear as black on the plots
pmask(pmask==1)=0;


plat=ncread(grid_filename, 'plat'); plat=plat';
plon=ncread(grid_filename, 'plon'); plon=plon';


figure;
pcolor (pmask');
shading flat
set(gca, 'color', 'k');
colormap white
colorbar

%#####reshapping MAP GRID so atlantic is in the middle... 
pmask_bit1=pmask(200:end,:); %position 200 is where Pacific is
pmask_bit2=pmask(1:199,:);
new_pmask=[pmask_bit1 ; pmask_bit2];

figure;
pcolor (new_pmask');
shading flat
set(gca, 'color', 'k');
colormap white
colorbar

%% Subdividing mask in different oceans

%% Loading Raw Ocean Mask
domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/noresm_mertraoceans_gx1v6.dat');

%fixing map appearance
%#####reshapping MAP GRID so atlantic is in the middle... 
domainmask_bit1=domain_mask(200:end,:); %position 200 is where Pacific is
domainmask_bit2=domain_mask(1:199,:);
domain_mask=[domainmask_bit1 ; domainmask_bit2];
figure
pcolor(domain_mask')
shading flat
%% excluding mediterranean from Atlantic, dead and black seas)

%assingning code for black and dead seas
for l=150:250
    for c=280:319      
if domain_mask(l,c)==1
    domain_mask(l,c)=5;
end
    end
end
    

%%assining code for mediterranean
for l=153:190
    for c=275:309      
if domain_mask(l,c)==2
    domain_mask(l,c)=6;
end
    end
end    
   

%assinging code to the baltic
for l=169:189
    for c=328:350      
if domain_mask(l,c)==2
    domain_mask(l,c)=7;
end
    end
end    
   
%%code to southern ocean 
for l=97:180
    for c=1:106      
if domain_mask(l,c)==1
    domain_mask(l,c)=8;
end
    end
end    


domain_mask(154:156,309)=2; %fixing some overlaps
domain_mask(169:173,349:351)=2;

figure
ax=pcolor(domain_mask')
shading flat
grid on
h=colorbar
colormap(parula(9))
caxis([0 9])
h.Ticks=.5:1:8;
h.TickLabels={'Land';'Southern & Arctic Ocn';'Atlantic';'Pacific';'Idian';'Black and Dead seas';'Mediteranean';'Baltic';'SO'};


%% Assigning code to latitudinal bands (Subpolar - 40-65; Subtropical 20-40; Tropical 0-20

domain_mask_lat_bands = domain_mask';


%solving lat limitis
%(x-191)*90 = y*191 , where x is the lat limit on old grid and y is lat from 0 - 90 deg
% x =(y*191 + 191*90)/90

tropical_lim = round((20*191 + 191*90)/90) + 24;
subtropical_lim = round((40*191 + 191*90)/90) + 26;
subpolar_lim = round((65*191 + 191*90)/90) + 22;

%assinging code to latitudinal bands

for l=189:tropical_lim %going from equator to 20N
    for c=30:250      
if domain_mask_lat_bands(l,c)==2
    domain_mask_lat_bands(l,c)=9;
end
    end
end 

for l=tropical_lim:subtropical_lim %going from 20N to 40N
    for c=30:250      
if domain_mask_lat_bands(l,c)==2
    domain_mask_lat_bands(l,c)=10;
end
    end
end 

for l=subtropical_lim:subpolar_lim %going from 40N to 65N
    for c=30:250      
if domain_mask_lat_bands(l,c)==2
    domain_mask_lat_bands(l,c)=11;
end
    end
end 


figure
ax=pcolor(domain_mask_lat_bands)
shading flat
grid on
h=colorbar
colormap(parula(11))
caxis([0 11])
h.Ticks=.5:1:13;
h.TickLabels={'Land';'Southern & Arctic Ocn';'Atlantic';'Pacific';'Idian';'Black and Dead seas';'Mediteranean';'Baltic';...
    'Southern O.'; 'Tropical'; 'Suptropical'; 'Subpolar'};




%% Ddelimiting box on map
% 
% lines_subPOL=[380 345]; %this is North Atl subPolar (Norwegian Sea etc)
% lines_NortAt=[344 320]; %this is North Atl (UK, Canada, Ireland)
% lines_subTro=[319 280]; %this is North Atl subPolar (Gulf of Maine, Sargasso Sea, Portugal, Morroco)
% 
% columns_band=[0 320];
% 
% domain_mask=domain_mask';
% 
% figure
% pcolor(domain_mask)
% %drawing box on map
% hold on
% subpolarbox=polyshape([columns_band(1) columns_band(1) columns_band(2) columns_band(2)], [lines_subPOL(2) lines_subPOL(1) lines_subPOL(1) lines_subPOL(2)]);
% NorAtlbox=polyshape([columns_band(1) columns_band(1) columns_band(2) columns_band(2)], [lines_NortAt(2) lines_NortAt(1) lines_NortAt(1) lines_NortAt(2)]);
% subtropbox=polyshape([columns_band(1) columns_band(1) columns_band(2) columns_band(2)], [lines_subTro(2) lines_subTro(1) lines_subTro(1) lines_subTro(2)]);
% plot(subpolarbox)
% hold on
% plot(NorAtlbox)
% hold on
% plot(subtropbox)
% hold off
% shading flat


save('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat','domain_mask')
save('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected_lat_bands.mat','domain_mask_lat_bands')

%%%%% ======= DOING GLOBE PROJECTION =========%%

field1=domain_mask_lat_bands; 
%unflipping
bit1=field1(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=field1(:,1:121);   %bit2 (1:199,:,:);
field2=cat(2,bit1,bit2);


long_new=shift_atl(plon);
lat_new=shift_atl(plat);
mask_globe=shift_atl(field2);

% for l=1:384
%     for c=1:320
%         
%         if mask_globe(l,c)==2 && (lat_new(l,c)>=0 && lat_new(l,c)<=65)
%             mask_globe(l,c)=9;
%         end
%         
%         if mask_globe(l,c)==0 
%             mask_globe(l,c)=NaN;
%         end
%         
%         if field2(l,c)==0 
%         field2(l,c)=NaN;
%         end
%     end
% end



figure
m_proj('Robinson','lon',[-100 15], 'lat', [0 65]) 
map2=m_pcolor(shift_atl(plon),shift_atl(plat),mask_globe)%,mask_globe)
%m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field2))%,mask_globe)
%map2.FaceAlpha=0.6;

map2.FaceLighting='gouraud';

        
caxis([0 11])
colormap(parula(12));        
h=colorbar
h.Ticks=1.5:1:9.5;
h.TickLabels={'Southern & Arctic';'Atlantic';'Pacific';'Idian';'Black & Dead seas';'Med.';'Baltic';'North Atlantic'};
%h.Location='southoutside';

m_coast('patch',[.77 .77 .77],'edgecolor','k');
m_grid('linest','-','yaxis','middle','xaxis','middle','xtick',-180:30:180,'ytick',[0 20 40 65],'backgroundcolor','k', ...
    'xticklabels',[], 'yticklabels',[0 20 40 65],'linewidth',5);
%m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);

%patch(.55*[-1 1 1 -1],.25*[-1 -1 1 1]-.55,'w'); 
hold on
%text(0,-.55,'M\_Map','fontsize',25,'color','b',...
%    'verticalalignment','middle','horizontalalignment','center');


set(gcf,'color','w')
supersizeme(1.5)
print('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/NATL_NORESM_grid_lat_bands.png','-dpng','-r300')

