%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SECTIONS PEAK VALUE AS COLOURED AND TPeak as contour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PRODUCES SURFACE AND SECTION PLOTS FOR THE RESULTS OF GLOBAL MAX or MIN
%ANALYSES

% Author: Leonardo Bertini
% Modified: 29th April 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CODE

depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';

depth_interp=min(depth):200:max(depth);

parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
plon=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plon') ;
plat=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plat') ;

plat=plat';
plon=plon';
plon=plon;

addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered
addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
addpath /Volumes/LaCie_Leonardo/NorESM

atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

mode = 'peak'; %can be set to either 'peak, delta_end or delta_rup

n = 4;

if n==1
    levels = [500 750 1000];
elseif n==2
    levels = [1250 1500 1750];
elseif n==3
    levels = [2000 2500 3000];
elseif n==4
    levels = [1000 1500 2500];
end

varlist={'ph';'o2';'omegac';'detoc';'templvl';'AOU'}  ;
for varIDX=5

if strcmp(varlist{varIDX},'ph')==1
    varname='pH';
    varunit='';
    axmin=7.3;
    axmax=8.3;
    colorscheme_val=parula(20);
    isolinessct=[20 50 70 100 140 180 240 280 340 400];
    isolinessfc=[20 50 100 150 180 200 240 300 380];

    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'o2')==1
    varname='DO';
    varunit='(\mumol O_2 kg^-^1)';
    axmin=-50;
    axmax=50;
    colorscheme_val=parula(20);
    colorscheme_time=parula(20);
    axmin_time=100;
    axmax_time=400;
    isolinessct=[20 50 70 100 140 180 240 280 340 400];
    isolinessfc=[80 150 240 380];
%120 320
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'omegac')==1
    varname='\Omega_C';
    varunit='';
    axmin=0;
    axmax=2;
    colorscheme_val=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/colorscheme_redblue_leo_100shades.mat');%flip(redblue,2);
    colorscheme_val=colorscheme_val.colorscheme;
    isolinessct=20:40:280;
    isolinessfc=[80 150 240 380];

    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;

elseif strcmp(varlist{varIDX},'detoc')==1
    varname='POC flux';
    varunit='(mgC m^-^2 day^-^1)';
    colorscheme_val=parula(20);
    axmin=0;
    axmax=0.5;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    isolinessct=[50 100 150 200 250];
    isolinessfc=[80 150 240 380];

    
elseif strcmp(varlist{varIDX},'templvl')==1
    varname='Temp';
    varunit='(\circC)';
    colorscheme_val=parula(12);
    axmin=6;
    axmax=12;
    colorscheme_time=parula(15);
    axmin_time=100;
    axmax_time=400;
    isolinessct=20:40:380;
    isolinessfc=20:40:380;

    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=10;

elseif strcmp(varlist{varIDX},'AOU')==1
    varname='AOU';
    varunit='(\mumol O2 kg^-^1)';
    colorscheme_val=parula(10);
    axmin=0;
    axmax=200;
   isolinessct=[60 140 200 280 380 ];
   isolinessfc=[20 60 100 140 200 380];

    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
 
end

 if strcmp(mode,'peak')
    fields=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}));
    fieldZERO = fields.PeakVALUE;
    %year1 = load year 1 data to do delta tmax
    fieldTpeak=fields.TpeakMAX; %Tmax values
    fieldTpeak1=permute(fieldTpeak,[3,1,2]);
    bit1=fieldTpeak1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
    bit2=fieldTpeak1(:,:,1:121);   %bit2 (1:199,:,:);
    fieldTpeak2=cat(3,bit1,bit2);
    
    figname = sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/TracerTmax_panels/Delta_%s_at_tmax_%d.png',varlist{varIDX},n);
    ttl_text = sprintf('  ) %s_t_m_a_x_ at',varname);

    delta_tpeak= NaN(384,320,70);
    for line = 1:384
        for col =1:320
            for k =1:70
                if ~isnan(fieldTpeak2(k,line,col))
                    delta_tpeak(line,col,k) = fieldZERO(line,col,k) - fieldTpeak2(k,line,col);
                else
                    delta_tpeak(line,col,k) = NaN;
                end
            end
        end
    end
    
    fieldZERO = fields.PeakVALUE;
    

 


 elseif strcmp(mode,'delta_end') || strcmp(mode,'delta_rup')
     data = load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Binded_var_files/NorESM_deltas_%s.mat',varlist{varIDX}));
   
     if strcmp(mode,'delta_rup')
        fieldZERO = data.delta_1_140;
        ttl_text = sprintf('  ) \\Delta%s at',varname);
        figname = sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/TracerTmax_panels/Delta_%s_RUP_%d.png',varlist{varIDX},n);
     elseif strcmp(mode,'delta_end')
        fieldZERO = data.delta_1_481;
        ttl_text = sprintf('  ) \\Delta%s at',varname);
        figname = sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/TracerTmax_panels/Delta_%s_EndExt_%d.png',varlist{varIDX},n);
     end
     
     if strcmp(varlist{varIDX},'o2')
    %if O2 convert to umol O2 kg-1 and do Delta based on Year 1
    fieldZERO=((fieldZERO)/1024)*1000000;
     end
 end

        
  
%     data = load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Binded_var_files/NorESM_deltas_%s.mat',varlist{varIDX}));
%     
%     
%     year1_o2=year1_o2.year1_o2;
%     fieldZERO=((fields.PeakVALUE)/1024)*1000000;
% 
field1=permute(fieldZERO,[3,1,2]);

%field1(field1<0)=NaN;
%unflipping
bit1=field1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=field1(:,:,1:121);   %bit2 (1:199,:,:);
field2=cat(3,bit1,bit2);

%%

k=[find(depth==levels(1)) find(depth==levels(2)) find(depth==levels(3))];  
hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[10.0894   31.3972   16.0161   27.1286]); %
set(hFig,'visible','on')
m1=0.01; p1=.002; s1=0.005; pb=.02; 

% surface 1
subaxis(3,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
get(gca,'Position');
set(gca,'Position',[0.03    0.7067+0.03    0.4875    0.2533]);

field3=squeeze(field2(k(1),:,:));

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
end

colormap(parula(12)) % [177/255 47/255 245/255]
caxis([axmin axmax])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on

% % if want to add Tcontour plot
% fieldTpeak3=squeeze(fieldTpeak2(k(1),:,:));
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
%     fieldTpeak3=smooth2d(round(fieldTpeak3));
% end
% [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3));
% 
% c.LevelList=isolinessfc; %levels
% c.LineColor=[.9 .9 .9];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',8) 

%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('%s %i m',ttl_text,round(depth(k(1)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])



%  surface 2
subaxis(3,2,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]);
set(gca,'Position',[0.03    0.7067-0.25    0.4875    0.2533]);

field3=squeeze(field2(k(2),:,:));

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
end

colormap(parula(12)) % [177/255 47/255 245/255]
caxis([axmin axmax])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on


%Grid box
%m_grid('backgroundcolor',[.2 .2 .2],'box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('%s %i m',ttl_text,round(depth(k(2)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])


% surface 3
subaxis(3,2,5,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
set(gca,'Position',[0.03    0.7067-0.25-0.28    0.4875    0.2533]);

field3=squeeze(field2(k(3),:,:));

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
end

colormap(parula(12)) % [177/255 47/255 245/255]
caxis([axmin axmax])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on

%Grid box
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
ttl=title(sprintf('%s %i m',ttl_text,round(depth(k(3)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])

% colorbar square 1
h=colorbar;
h.Location='southoutside';
h.Position=[0.06 0.2700-0.12 0.4493 0.0177]
caxis([axmin axmax])
h.Label.String=sprintf('%s_t_m_a_x %s',varname,varunit);
set(h,'fontsize',12);

%Print fig
set(gcf,'color','w')
hFig.InvertHardcopy='off';
print(figname,'-dpng','-r300')



%% New Fig
fieldTpeak=fields.TpeakMAX; %Tmax values
fieldTpeak1=permute(fieldTpeak,[3,1,2]);
bit1=fieldTpeak1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldTpeak1(:,:,1:121);   %bit2 (1:199,:,:);
fieldTpeak2=cat(3,bit1,bit2);


hFig = figure(2); clf;
set(hFig, 'Units','centimeters','Position',[10.0894   31.3972   16.0161   27.1286]); %
set(hFig,'visible','on')
m1=0.01; p1=.002; s1=0.005; pb=.02; 
axmax_time = 400
% surface 4

subaxis(3,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
get(gca,'Position');
set(gca,'Position',[0.03    0.7067+0.03    0.4875    0.2533]);




fieldTpeak3=squeeze(fieldTpeak2(k(1),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTpeak3=smooth2d(round(fieldTpeak3));
end

if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3)); shading interp
end

colormap(colorscheme_time) % [177/255 47/255 245/255]
caxis([axmin_time axmax_time])

m_coast('patch',[.77 .77 .77],'edgecolor','k');

hold on

%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('  ) tmax %s at %i m',varname,round(depth(k(1)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])


% surface 5
subaxis(3,2,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
get(gca,'Position');
set(gca,'Position',[0.03    0.7067-0.25    0.4875    0.2533]);

fieldTpeak3=squeeze(fieldTpeak2(k(2),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTpeak3=smooth2d(round(fieldTpeak3));
end

if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3)); shading interp
end

colormap(colorscheme_time) % [177/255 47/255 245/255]
caxis([axmin_time axmax_time])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on


%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('  ) tmax %s at %i m',varname,round(depth(k(2)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])

% surface 6
subaxis(3,2,5,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
get(gca,'Position');
set(gca,'Position',[0.03    0.7067-0.25-0.28    0.4875    0.2533]);

fieldTpeak3=squeeze(fieldTpeak2(k(3),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTpeak3=smooth2d(round(fieldTpeak3));
end

if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3)); shading interp
end

colormap(colorscheme_time) % [177/255 47/255 245/255]
caxis([axmin_time axmax_time])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on


%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('  ) tmax %s at %i m',varname,round(depth(k(3)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])

% colorbar square 2
h=colorbar;
h.Location='southoutside';
h.Position=[0.06 0.2700-0.12 0.4493 0.0177]
caxis([axmin_time axmax_time])
h.Label.String=sprintf('tmax %s (%s)',varname,'yr');
set(h,'fontsize',12);

%Print fig
set(gcf,'color','w')
hFig.InvertHardcopy='off';

print(sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/TracerTmax_panels/tmax_%s_%d.png',varlist{varIDX},n),'-dpng','-r300')

%% FIg Initial DO layers

fieldTpeak=year1_o2; %Tmax values
fieldTpeak1=permute(fieldTpeak,[3,1,2]);
bit1=fieldTpeak1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldTpeak1(:,:,1:121);   %bit2 (1:199,:,:);
fieldTpeak2=cat(3,bit1,bit2);


hFig = figure(3); clf;
set(hFig, 'Units','centimeters','Position',[10.0894   31.3972   16.0161   27.1286]); %
set(hFig,'visible','on')
m1=0.01; p1=.002; s1=0.005; pb=.02; 

% surface 7

subaxis(3,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
get(gca,'Position');
set(gca,'Position',[0.03    0.7067+0.03    0.4875    0.2533]);

fieldTpeak3=squeeze(fieldTpeak2(k(1),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTpeak3=smooth2d(round(fieldTpeak3));
end

if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3)); shading interp
end

colormap(parula(20)) % [177/255 47/255 245/255]
caxis([0 400])

m_coast('patch',[.77 .77 .77],'edgecolor','k');

hold on

%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('  ) PI DO at %i m', round(depth(k(1)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])


% surface 8
subaxis(3,2,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
get(gca,'Position');
set(gca,'Position',[0.03    0.7067-0.25    0.4875    0.2533]);

fieldTpeak3=squeeze(fieldTpeak2(k(2),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTpeak3=smooth2d(round(fieldTpeak3));
end

if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3)); shading interp
end

colormap(parula(20)) % [177/255 47/255 245/255]
caxis([0 400])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on


%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('  ) PI DO at %i m', round(depth(k(2)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])

% surface 9
subaxis(3,2,5,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
get(gca,'Position');
set(gca,'Position',[0.03    0.7067-0.25-0.28    0.4875    0.2533]);

fieldTpeak3=squeeze(fieldTpeak2(k(3),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTpeak3=smooth2d(round(fieldTpeak3));
end

if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3)); shading interp
end

colormap(parula(20)) % [177/255 47/255 245/255]
caxis([0 400])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on


%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(sprintf('  ) PI DO at %i m', round(depth(k(3)))),'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])

% colorbar square 2
h=colorbar;
h.Location='southoutside';
h.Position=[0.06 0.2700-0.12 0.4493 0.0177]
caxis([0 400])
h.Label.String=sprintf('DO %s', varunit);
set(h,'fontsize',12);

%Print fig
set(gcf,'color','w')
hFig.InvertHardcopy='off';


print(sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/TracerTmax_panels/PI_Value_%s_%d.png',varlist{varIDX},n),'-dpng','-r300')


%% section

% subplot(3,2,[2,4,6])
% get(gca,'Position');
% % % %option2
% % % subplot(3,2,[2,4,6])
% % % get(gca,'Position')
% % % set(gca,'Position',[0.5 0.09760 0.4073 0.7914])
% 
% % bottomflag=-1111;
% % PI_mean_1(PI_mean_1==bottomflag)=NaN; %this makes the points representing the bottom 'transparent'
% % 
% % %changing all the negative flags to avoid the concentric rings when interpolating
% % PI_mean_1(PI_mean_1<0)=600;
% 
% %getting section data
% [sect_fld_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(field2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% 
% for smoothnumber=1:smoothnumbsect %smoothing 3 times
%     sect_fld_original=smooth2d(sect_fld_original);
% end
% 
% %section plot
% if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% pcolor(sect_lat,-sect_dep,sect_fld_original*12.0107*5*1000);
% else
% pcolor(sect_lat,-sect_dep,sect_fld_original);
% end
% shading interp
% hold on
% 
% %Tpeak
% [sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTpeak2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% 
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times
%     sect_ctr_original=smooth2d(sect_ctr_original);
% end
% 
% [M,c]=contour(sect_lat,-sect_dep,sect_ctr_original);
% 
% c.LevelList=isolinessct; %levels
% c.LineColor=[.9 .9 .9];%'grey';
% c.LineStyle='-';
% c.LineWidth=2;
% clabel(M,c,'FontSize',10) 
% hold on
% title(['d) North Atlantic section of Tracer_t_m_a_x ' sprintf('%s', varname)],'fontsize',14)
% 
% %contourlines
% % [sect_fld_original,C,h]=makesect_atl_noresm_fixed_z_COUNTOUR_Tmin(field2,depth,plat,atl_ind);
% % h.LevelList=[250 350 450 480];
% % h.Color=[.7 .7 .7];
% % clabel(C,h,'Color','white','FontSize',8);
% % clabel(C,h,'LabelSpacing',300)
% % h.LineStyle='-';
% 
% 
% colormap(colorscheme) % [177/255 47/255 245/255]
% caxis([axmin axmax])
%  
% h=colorbar;
% caxis([axmin axmax])
% h.Label.String=sprintf('Tracer_t_m_a_x %s %s',varname,varunit);
% h.Location='southoutside';
% set(h,'fontsize',12);
% 
% xlabel('Latitude','fontsize',12)
% ylabel('Depth (m)','fontsize',12)
% 
% if varIDX==4 %if detoc then change axis to end at 2000m deep
%     axis([-10 65 -1800 0])
% else
%     axis([-10 65 -6000 0])
% end
% set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on','xtick',-10:10:60,'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'  });
% set(gca,'Position',[0.40    0.18    0.55    0.77]);

%saving image


close all
end











% k=[find(depth==500) find(depth==1000) find(depth==2000)];
% 
% fig = figure('Units','centimeters','Position',[0,0,21,29.7]); %that's a4 size
% 
% %% surface 1
% subplot(2,3,1)
% m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 
% 
% field3=squeeze(field2(k(1),:,:));
% 
% for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
%     field3=smooth2d(field3);
% end
% 
% 
% if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
% else
% m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
% end
% 
% colormap(colorscheme) % [177/255 47/255 245/255]
% caxis([axmin axmax])
% 
% m_coast('patch',[.7 .7 .7],'edgecolor','k');
% hold on
% 
% %Tcontour
% fieldTpeak3=squeeze(fieldTpeak2(k(1),:,:));
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
%     fieldTpeak3=smooth2d(round(fieldTpeak3));
% end
% [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3))
% 
% c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',10) 
% 
% %Grid box
% m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
% hold on
% m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
% title(sprintf('Depth= %.2f m',depth(k(1))))
% 
% %% surface 2
% subplot(2,3,2)
% m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 
% 
% field3=squeeze(field2(k(2),:,:));
% 
% for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
%     field3=smooth2d(field3);
% end
% 
% 
% if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
% else
% m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
% end
% 
% colormap(colorscheme) % [177/255 47/255 245/255]
% caxis([axmin axmax])
% 
% m_coast('patch',[.7 .7 .7],'edgecolor','k');
% hold on
% 
% %Tcontour
% fieldTpeak3=squeeze(fieldTpeak2(k(2),:,:));
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
%     fieldTpeak3=smooth2d(round(fieldTpeak3));
% end
% [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3))
% 
% c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',10) 
% 
% %Grid box
% m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
% hold on
% m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
% title(sprintf('Depth= %.2f m',depth(k(2))))
% 
% %% surface 3
% 
% subplot(2,3,3)
% m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 
% 
% field3=squeeze(field2(k(3),:,:));
% 
% for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
%     field3=smooth2d(field3);
% end
% 
% if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
% else
% m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
% end
% 
% colormap(colorscheme) % [177/255 47/255 245/255]
% caxis([axmin axmax])
% 
% m_coast('patch',[.7 .7 .7],'edgecolor','k');
% hold on
% 
% %Tcontour
% fieldTpeak3=squeeze(fieldTpeak2(k(3),:,:));
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
%     fieldTpeak3=smooth2d(round(fieldTpeak3));
% end
% [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTpeak3))
% 
% c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',10) 
% 
% %Grid box
% m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
% hold on
% m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
% title(sprintf('Depth= %.2f m',depth(k(3))))
% 
% %% section
% subplot(2,3,4:6)
% get(gca,'Position')
% set(gca,'Position',[0.1300    0.08    0.7750    0.5])
% 
% bottomflag=-1111;
% 
% field2(field2==bottomflag)=NaN; %this makes the points representing the bottom 'transparent'
% 
% %changing all the negative flags to avoid the concentric rings when interpolating
% field2(field2<0)=600;
% 
% %getting section data
% [sect_fld_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(field2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% 
% for smoothnumber=1:smoothnumbsect %smoothing 3 times
%     sect_fld_original=smooth2d(sect_fld_original);
% end
% 
% %section plot
% if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% pcolor(sect_lat,-sect_dep,sect_fld_original*12.0107*5*1000);
% else
% pcolor(sect_lat,-sect_dep,sect_fld_original);
% end
% shading interp
% hold on
% 
% %Tcontour
% [sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTpeak2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% 
% for smoothnumber=1:smoothnumbcontour %smoothing 
%     sect_ctr_original=smooth2d(sect_ctr_original);
% end
% 
% [M,c]=contour(sect_lat,-sect_dep,sect_ctr_original)
% 
% c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',10) 
% hold on
% 
% %contourlines
% % [sect_fld_original,C,h]=makesect_atl_noresm_fixed_z_COUNTOUR_Tmin(field2,depth,plat,atl_ind);
% % h.LevelList=[250 350 450 480];
% % h.Color=[.7 .7 .7];
% % clabel(C,h,'Color','white','FontSize',8);
% % clabel(C,h,'LabelSpacing',300)
% % h.LineStyle='-';
% 
% 
% colormap(colorscheme) % [177/255 47/255 245/255]
% caxis([axmin axmax])
%  
% h=colorbar;
% caxis([axmin axmax])
% h.Location='northoutside';
% h.Label.String=sprintf('%s Peak Value %s',varname,varunit);
% h.Location='northoutside';
% 
% xlabel('Latitude (\circ)')
% ylabel('Depth (m)')
% if varIDX==4 %if detoc then change axis to end at 2000m deep
%     axis([-10 65 -2000 0])
% else
%     axis([-10 65 -6000 0])
% end
% set(gca, 'color','k' );





% 
% 
% 
% 
% 
% 
% %% if it's pH
% if varIDX==1
%     
%     %pcolor
%     fig=figure('Units','centimeters','Position',[0,0,29.7,21]); %that's a4 size
%     [sect_fld_original]=makesect_atl_noresm_fixed_z(field2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
%     
%     %contourlines of percentage change
%     statfields=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/%s_stats_updated_v1D_novertinterp',varlist{varIDX}));
%     percentage_change_field=statfields.Percentage_change_hidrogen_ion;
%     
%     percentage_change_field=permute(percentage_change_field,[3,1,2]);
%     
%     hold on
%     %unflipping
%     bit1=percentage_change_field(:,:,120:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
%     bit2=percentage_change_field(:,:,1:119);   %bit2 (1:199,:,:);
%     percentage_change_field=cat(3,bit1,bit2);
%     
%     [countourlinefield,C,h]=makesect_atl_noresm_fixed_z_COUNTOUR_Tmin(percentage_change_field,depth,plat,atl_ind);
%     
%     %h.LevelList=7.5:0.10:8.2;
%     
%     shading interp
%     colormap(parula(20)) % [177/255 47/255 245/255]
%     
%     h=colorbar;
%     h.Location='northoutside';
%     h.Label.String=sprintf('%s Minimum',varlist{varIDX});
%     
%     %caxis([7.5 8.0])
%     %caxis([0 0.5])
%     %h.Ticks=7.5:0.05:8.0;
%     %h.TickLabels(1)={'\leq.5'};
%     caxis([7.0 8.0])
%     h.TickLabels(end)={'\geq8.0'};
%     h.TickLabels(1)={'\leq7.0'};
%     
%     xlabel('Latitude')
%     ylabel('Depth (m)')
%     axis([-10 75 -6000 0])
%     set(gca, 'color', 'k');
%     %     hold on
%     %      %yline(-depth(k),'--','LineWidth',1)
%     %      print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Figures/Section_North_Atlantic_%s_%s.png',timevarname,varlist{varIDX}),'-dpng','-r300')
%     %
%     
% else
%     %% if it is omega
%     
%     %pcolor
%     fig=figure('Units','centimeters','Position',[0,0,29.7,21]); %that's a4 size
%     [sect_fld_original]=makesect_atl_noresm_fixed_z(field2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
%     hold on
%     %single contour line of 1
%     [sect_fld_original,C,h]=makesect_atl_noresm_fixed_z_COUNTOUR_Tmin(field2,depth,plat,atl_ind)
%     h.LevelList=1;
%     h.LineColor='red';
%     h.LineWidth=1;
%     h.LineStyle='--';
%     clabel(C,h,'color','red')
%     
%     hold on
%     %%pre industrial mean field contour lines
%     PI_field=load('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_omegac_first_5yr_mean.mat');
%     PI_field=PI_field.PI_mean;
%     
%     PI_field=permute(PI_field,[3,2,1]);
%     
%     %unflipping
%     %      bit1=PI_field(:,:,120:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
%     %      bit2=PI_field(:,:,1:119);   %bit2 (1:199,:,:);
%     %      PI_field=cat(3,bit1,bit2);
%     
%     [countourlinefield,C,h]=makesect_atl_noresm_fixed_z_COUNTOUR_Tmin(PI_field,depth,plat,atl_ind);
%     clabel(C,h,'color','white')
%     h.LineColor='k';
%     h.LevelList=0.5:0.25:2;
%     
%     
%     shading interp
%     colormap(parula(10)) % [177/255 47/255 245/255]
%     
%     h=colorbar;
%     h.Location='northoutside';
%     h.Label.String=sprintf('%s Minimum',varname);
%     
%     %caxis([7.5 8.0])
%     %caxis([0 0.5])
%     caxis([0.5 1.5])
%     %h.Ticks=7.5:0.05:8.0;
%     h.TickLabels(end)={'\geq1.5'};
%     h.TickLabels(1)={'\leq0.5'};
%     %h.TickLabels(1)={'\leq.5'};
%     
%     xlabel('Latitude')
%     ylabel('Depth (m)')
%     axis([-10 75 -6000 0])
%     set(gca, 'color', 'k');
%     
%     lgnd= legend('','Chemical Equilibrium','\Omega_C Pre-Industrial mean');
%     set(lgnd,'color','white');
%     lgnd.Location='southeast';
%     
%     % hold on
%     %yline(-depth(k),'--','LineWidth',1)
%     %print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Figures/Section_North_Atlantic_%s_%s.png',timevarname,varlist{varIDX}),'-dpng','-r300')
%     
%     
%     
%     
%     % close all
%     %
%     %
%     %
%     %     drawnow
%     %         frame = getframe(fig);
%     %         im{k} = frame2im(frame);
%     %         hold all
%     %
%     % pause( 0.5 )
%     % end
%     % close;
%     %
%     % filename = 'Section_North_Atlantic_pH_Minimum.gif'; % Specify the output file name
%     %
%     %  for k = 1 : (length(depth)-1)
%     %     [A,map] = rgb2ind(im{k},256);
%     %     if k == 1
%     %         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     %     else
%     %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
%     %     end
%     %  end
%     %
%     %  close all
%     
%     
%     %section
%     