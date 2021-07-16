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

addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered
addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
addpath /Volumes/LaCie_Leonardo/NorESM

atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');
varlist={'ph';'o2';'templvl'}  ;

n = 1;

    if n ==1
        levels = [500 750 1000];
    elseif n==2
        levels = [1250 1500 1750];
    elseif n==3
        levels = [2000 2500 3000];
    end



for varIDX = [1]

if strcmp(varlist{varIDX},'ph')==1
    varname='pH';
    varunit='';
    axmin=-.8;
    axmax=.8;
    colorscheme_val=redblue(32);
    isolinessct=[20 50 70 100 140 180 240 280 340 400];
    isolinessfc=[20 50 100 150 180 200 240 300 380];
    bartick_vals = [axmin -.5 -.2 0 .2 .5 axmax];
    bartick_labels = {string(axmin);'-0.5';'-0.2';'0';'0.2';'0.5';string(axmax)};

    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'o2')==1
    varname='DO';
    varunit='(\mumol O_2 kg^-^1)';
    axmin=-50;
    axmax=50;
    colorscheme_val=redblue(50);
    isolinessct=[20 50 70 100 140 180 240 280 340 400];
    isolinessfc=[80 150 240 380];
    bartick_vals = [axmin -30 -10 0 10 30 axmax];
    bartick_labels = {string(axmin);'-30';'-10';'0';'10';'30';string(axmax)};

    
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    

elseif strcmp(varlist{varIDX},'templvl')==1
    varname='T';
    varunit='(\circC)';
    colorscheme_val=redblue(24)
    colorscheme_val=colorscheme_val(9:end,:);
    axmin=-2;
    axmax=10;
    isolinessct=20:40:380;
    isolinessfc=20:40:380;
    bartick_vals = [axmin 0 2 4 6 8 axmax];
    bartick_labels = {string(axmin);'0';'2';'4';'6';'8';string(axmax)};

    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=10;
    
end

    folder_name = sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/Deltas/%s',varlist{varIDX});
    if not(isfolder(folder_name))
        mkdir(folder_name)
    end
    
    %load delta data using matfile function
    m = matfile(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Binded_var_files/NorESM_deltas_%s.mat',varlist{varIDX}));
    strucvar = who(m);
    delta_struc = m.(strucvar{1})(:,1);
    
    
    for time = 1:481
    
    figname = sprintf('%s/Delta_%s_at_year_%d_depth_group_%d.png',folder_name,varlist{varIDX},time, n);

    %extract data from structure at specific year
    fieldZERO = delta_struc.year{time,1};


     if strcmp(varlist{varIDX},'o2')
        %if O2 convert to umol O2 kg-1 and do Delta based on Year 1
        fieldZERO=((fieldZERO)/1024)*1000000;
     end

     ttl_text = sprintf('  ) \\Delta%s at',varname);

    field1=permute(fieldZERO,[3,1,2]);

    %field1(field1<0)=NaN;
    %unflipping
    bit1=field1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
    bit2=field1(:,:,1:121);   %bit2 (1:199,:,:);
    field2=cat(3,bit1,bit2);

    %% Ploting figure
    k=[find(depth==levels(1)) find(depth==levels(2)) find(depth==levels(3))];  
    hFig = figure(1); clf;
    set(hFig, 'Units','centimeters','Position',[10.0894   31.3972   16.0161   27.1286]); %
    set(hFig,'Visible','on')
    m1=0.01; p1=.002; s1=0.005; pb=.02; 

    % surface 1
    subaxis(3,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
    m_proj('Lambert','lon',[-100 15], 'lat', [0 65]) ;
    get(gca,'Position');
    set(gca,'Position',[0.03+0.01    0.7067+0.03    0.4875    0.2533]);

    field3=squeeze(field2(k(1),:,:));

    for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
        field3=smooth2d(field3);
    end


    if varIDX==4 %if POC then calculate data in mgC m-2 day-1
    m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
    else
    m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
    end

    colormap(colorscheme_val) % [177/255 47/255 245/255]
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
    set(gca,'Position',[0.03+0.01    0.7067-0.25    0.4875    0.2533]);

    field3=squeeze(field2(k(2),:,:));

    for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
        field3=smooth2d(field3);
    end


    if varIDX==4 %if POC then calculate data in mgC m-2 day-1
    m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
    else
    m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
    end

    colormap(colorscheme_val) % [177/255 47/255 245/255]
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
    set(gca,'Position',[0.03+0.01    0.7067-0.25-0.28    0.4875    0.2533]);

    field3=squeeze(field2(k(3),:,:));

    for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
        field3=smooth2d(field3);
    end


    if varIDX==4 %if POC then calculate data in mgC m-2 day-1
    m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
    else
    m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
    end

    colormap(colorscheme_val) % [177/255 47/255 245/255]
    caxis([axmin axmax])

    m_coast('patch',[.77 .77 .77],'edgecolor','k');
    hold on

    %Grid box
    m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
    set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
    hold on
    m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
    ttl=title(sprintf('%s %i m',ttl_text,round(depth(k(2)))),'fontsize',14);
    set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])

    % colorbar square 1
    h=colorbar;
    h.Location='southoutside';
    h.Position=[0.06 0.2700-0.12 0.4493 0.0177]
    h.Ticks = bartick_vals
    h.TickLabels = bartick_labels
    caxis([axmin axmax])
    h.Label.String=sprintf('\\Delta%s %s',varname,varunit);
    set(h,'fontsize',12);
    
    hold on
    anot = annotation('textbox', [0.2237    0.1716    0.1310    0.0219], 'String', sprintf('Year = %03d',time))
    anot.BackgroundColor = 'k'
    anot.Color = 'w'
    anot.FontWeight = 'bold'
    anot.FontSize = 9;
    
    %Print fig
    set(gcf,'color','w')
    hFig.InvertHardcopy='off';
    %save fig  
    print(figname,'-dpng','-r300')
    end
end