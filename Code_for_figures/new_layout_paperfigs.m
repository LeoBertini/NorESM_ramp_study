%% LOADING INITIAL STUFF

depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';

depth_interp=min(depth):100:max(depth);

%parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
plon=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plon') ;
plat=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plat') ;

plat=plat';
plon=plon';
plon=plon;

addpath(genpath('/Volumes/LaCie_Leonardo/NorESM'))
addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
addpath /Volumes/LaCie_Leonardo/NorESM


domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat');
domain_mask=domain_mask.domain_mask;

atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

varlist={'ph';'o2';'omegac';'templvl';'AOU'};

varIDX=5

  %% PI_mean

    PI_mean=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_first_5yr_mean.mat',varlist{varIDX}));
   
    %if O2 convert to umol O2 kg-1
    if varIDX==2
    PI_mean=((PI_mean.PI_mean)/1024)*1000000;
    else
    PI_mean=PI_mean.PI_mean;
    end
    
    PI_mean=permute(PI_mean,[2,1,3]);
    
    %unflipping domain mask
    bit1=domain_mask(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
    bit2=domain_mask(:,1:121);   %bit2 (1:199,:,:);
    domain_mask_1=cat(2,bit1,bit2);
    
    for l=1:384
    for c=1:320
        for k=1:70
            if domain_mask_1(l,c)~=2 
                PI_mean(l,c,k)=NaN;
            end
        end
    end
    end

%pcolor(PI_mean(:,:,1));shading flat
PI_mean_1=permute(PI_mean,[3,1,2]);   

%% Importing TOD
timevarname='TOD';
fieldTOD=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}),timevarname);
fieldTOD=fieldTOD.TOD;

fieldTOD1=permute(fieldTOD{1,2},[3,1,2]); 
%unflipping
bit1=fieldTOD1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldTOD1(:,:,1:121);   %bit2 (1:199,:,:);
fieldTOD2=cat(3,bit1,bit2);

%pcolor(squeeze(fieldTOD2(1,:,:)));shading flat

%flag of no departure
for k=1:70
    for l=1:384
        for c=1:320
            if domain_mask_1(l,c)==2 && ~isnan(PI_mean_1(k,l,c)) && isnan(fieldTOD2(k,l,c))
                fieldTOD2(k,l,c)=-10; %flag of no departure
            end
        end
    end
end 

%% IMPORTING VERTICAL AVERAGES OF TOD
load(sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/TOD_vertprof_%s',varlist{varIDX}))

%% Importing Trevocery dataset
%Possible flags included in the Trecovery dataset
badflagsnegslope=-7777;
badflagtooearly=-4444;
bottomflag=-1111;

%IMPORTING DATA
timevarname='Trecovery';
fields=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}),timevarname);
fieldZERO=fields.Trecovery{1,2};%2 is the one based on 2sd envelope

fieldZERO(fieldZERO==badflagsnegslope)=120; %this makes all different flags the same value %Definite no recovery
fieldZERO(fieldZERO==badflagtooearly)=120; %this makes all different flags the same value %Gets

fieldTrec1=permute(fieldZERO,[3,1,2]); 
%unflipping
bit1=fieldTrec1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldTrec1(:,:,1:121);   %bit2 (1:199,:,:);
fieldTrec2=cat(3,bit1,bit2);

fieldTrec2(fieldTrec2<0 & fieldTrec2~=bottomflag)=600; %this makes all negative flags which are not bottom positive to indicate very late or no recovery


%% IMPORTING VERTICAL AVERAGES OF Trec
load(sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/Trec_vertprof_%s',varlist{varIDX}))

%% IMPORTING PCHANGE for time windows AND DOING VERTICAL AVERAGES

Pchange_data=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D_vertinterp.mat',varlist{varIDX}));

%calculating mean pchange per period at depth
for numb=1:3

         dummy1=eval(sprintf('Pchange_data.Percentage_change_mean_pair%d',numb));
         bit1=dummy1(:,122:end,:); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
         bit2=dummy1(:,1:121,:);   %bit2 (1:199,:,:);
         dummy2=cat(2,bit1,bit2);
          
         dummy3=permute(dummy2,[3,1,2]);
         
         Pchange_flipped(numb)={dummy3};
  
            for i=1:size(dummy3,1)  %calculating the mean for each depth layer for a set range of latitudes
            
            tgtlayer=dummy3(i,:,:);
            tgtlayer=tgtlayer(:);
            indxlat=find(plat<=65 & plat>=0);

            Pchange_sd_(numb,i)={nanstd(tgtlayer(~isnan(tgtlayer(indxlat))))};
            Pchange_vertprof(numb,i)={nanmean(tgtlayer(~isnan(tgtlayer(indxlat))))};

            end      
end


%% FIGURE PARAMETERS

answer=[];
axminTrec=120;
axmaxTrec=720;

%%%Trec special colormap
    pinky=pink(20);
    use=pinky(5:10,:);use=flip(use);
    dummycolor=parula(13);
    colorscheme2=[dummycolor(2:10,:);use];


    if strcmp(varlist{varIDX},'ph')==1
        %call message input
        opts.Interpreter = 'tex';
        answer = inputdlg({'Type 1 for pH or Type 2 for [H^+]:'},...
            'pH mode Plots',[1 50],{'1'},opts);
        
        if cell2mat(answer)=='2'
            varname='[H^+]';
            %calculating mean pchange per period at depth
            for numb=1:3
                
                dummy1=eval(sprintf('Pchange_data.Percentage_change_hidrogen_ion_pair%d',numb));
                bit1=dummy1(:,122:end,:); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
                bit2=dummy1(:,1:121,:);   %bit2 (1:199,:,:);
                dummy2=cat(2,bit1,bit2);
                
                dummy3=permute(dummy2,[3,1,2]);
                
                Pchange_flipped(numb)={dummy3};
                
                for i=1:size(dummy3,1)  %calculating the mean for each depth layer for a set range of latitudes
                    
                    tgtlayer=dummy3(i,:,:);
                    tgtlayer=tgtlayer(:);
                    indxlat=find(plat<=65 & plat>=0);
                    
                    Pchange_sd_(numb,i)={nanstd(tgtlayer(~isnan(tgtlayer(indxlat))))};
                    Pchange_vertprof(numb,i)={nanmean(tgtlayer(~isnan(tgtlayer(indxlat))))};
                    
                end
            end
            
            varunit='';
            axminpchange=-20;
            axmaxpchange=+20;
            colorscheme3=redblue(50);
            contourlevels=280:20:480;
            smoothnumbsurface=3;
            smoothnumbsect=3;
            tickjump=5;
            axminred=-10;
            axmaxred=90;
            
        else
            
            varname='pH';
            
            for numb=1:3
                
                dummy1=eval(sprintf('Pchange_data.Percentage_change_mean_pair%d',numb));
                bit1=dummy1(:,122:end,:); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
                bit2=dummy1(:,1:121,:);   %bit2 (1:199,:,:);
                dummy2=cat(2,bit1,bit2);
                
                dummy3=permute(dummy2,[3,1,2]);
                
                Pchange_flipped(numb)={dummy3};
                
                for i=1:size(dummy3,1)  %calculating the mean for each depth layer for a set range of latitudes
                    
                    tgtlayer=dummy3(i,:,:);
                    tgtlayer=tgtlayer(:);
                    indxlat=find(plat<=65 & plat>=0);
                    
                    Pchange_sd_(numb,i)={nanstd(tgtlayer(~isnan(tgtlayer(indxlat))))};
                    Pchange_vertprof(numb,i)={nanmean(tgtlayer(~isnan(tgtlayer(indxlat))))};
                end
            end
            
            varunit='';
            axminTOD=0;
            axmaxTOD=100;
            
            colorscheme=parula(20);
            contourlevels_layers=[20 50 70 100];
            contourlevels_section=[10 20 40 80 100];
            smoothnumbsurface=7;
            smoothnumbsect=1;
            smoothnumbcontour=7;
            isolinessct=[200 400 480 650 750 1000];
            labelpos=[2.5   50.000  0];
            
            axmaxred=50;
            axminred=0;
            axmaxblue=1400;
            axminblue=200;
            axminpchange=-3;
            axmaxpchange=0;
            axminpchangesect=-2.5;
            axmaxpchangesect=2.5;
            tickjumpTOD=20;
            tickjumpTrec=80;
            tickjumppchange=.5;
            colorscheme3=redblue(50);

        end
        
        
    elseif strcmp(varlist{varIDX},'o2')==1
        varname='DO';
        varunit='(\mumol O_2 kg^-^1)';
        axminTOD=0;
        axmaxTOD=120;
        colorscheme=parula(30);
        contourlevels_layers=[40 80 140];
        contourlevels_section=[20 40 80 120];
        smoothnumbsurface=7;
        smoothnumbsect=1;
        smoothnumbcontour=7;
        labelpos=[2.5   50.000  0];
        
        axmaxred=120;
        tickjumpTOD=20;
        
       axminpchange=-10;
       axmaxpchange=30;
       axminpchangesect=-20;
       axmaxpchangesect=20;
       tickjumppchange=5;
       colorscheme3=redblue(40);
       isolinessct=[200 480 750 1000];

        
    elseif strcmp(varlist{varIDX},'omegac')==1
        varname='\Omega_C';
        varunit='';
        axminTOD=0;
        axmaxTOD=100;
        colorscheme=parula(20)%flip(redblue,2);
        contourlevels_layers=[20 50 80 100];
        contourlevels_section=[5 20 40 80 100];
        smoothnumbsurface=7;
        smoothnumbsect=7;
        smoothnumbcontour=7;
        labelpos=[2.5 50.0000  0];

        axminred=0;
        axmaxred=80;
        axminblue=200;
        axmaxblue=1600;
        
        tickjumpTOD=20;
        tickjumpTrec=80;
        
       axminpchange=-40;
       axmaxpchange=10;
       axminpchangesect=-20;
       axmaxpchangesect=20;
                   
       tickjumppchange=5;
       colorscheme3=redblue(40);
       isolinessct=[200 320 480 750 1000];
        
        
    elseif strcmp(varlist{varIDX},'templvl')==1
        varname='T';
        varunit='(\circC)';
        colorscheme=parula(20);
        axminTOD=0;
        axmaxTOD=100;
        contourlevels_layers=[20 40 60 100];
        contourlevels_section=[20 40 60 80];
        smoothnumbsurface=12;
        smoothnumbsect=1;
        smoothnumbcontour=15;
        labelpos=[2.5   50  0];
        
        axminred=0;
        axmaxred=100;
        axminblue=200;
        axmaxblue=1600;
        
        tickjumpTOD=20;
        tickjumpTrec=80;
        
       axminpchange=-10;
       axmaxpchange=90;
       axminpchangesect=-40;
       axmaxpchangesect=40;
                   
       tickjumppchange=10;
       colorscheme3=redblue(40);
       isolinessct=[200 480 600 1000];
        
    elseif strcmp(varlist{varIDX},'AOU')==1
        varname='AOU';
        varunit='(\mumol O2 kg^-^1)';
        colorscheme=parula(10);
        axminTOD=0;
        axmaxTOD=200;
        contourlevels_layers=[20 40 120];
        contourlevels_section=[20 40 80 120];
        
        smoothnumbsurface=7;
        smoothnumbsect=1;
        smoothnumbcontour=7;
        labelpos=[2.3937  100.0001         0];
         
        axminred=0;
        axmaxred=140;
        axminblue=200;
        axmaxblue=1400;
        
        tickjumpTOD=20;
        tickjumpTrec=80;
        
       axminpchange=-30;
       axmaxpchange=70;
       axminpchangesect=-40;
       axmaxpchangesect=40;
                   
       tickjumppchange=10;
       colorscheme3=redblue(40);
       isolinessct=[200 340 480 1000];
        
    end


%% Figure 3X3 (3 lines and 3 columns)

hFig1 = figure(1); clf;
set(hFig1, 'Units','pixels','Position',[0,0,800,700]); %
set(hFig1,'visible','on')

%PANEL Vertical data
panel1=subplot(8,6,[1,2,7,8,13,14,19,20])
%filled area
    x1=smooth(cell2mat(TOD_vertprof(1,:)));
    x3=smooth(cell2mat(TOD_vertprof(3,:)));
    X = [x1(1:65).', fliplr(x3(1:65).')];
    Y = [-depth(1:65), fliplr(-depth(1:65))];
    shaded1=fill(smooth(X), smooth(Y), [255/255 204/255 204/255]);
    shaded1.FaceAlpha=0.4;
    shaded1.LineStyle='none';
    
    hold on
    x2=smooth(cell2mat(TOD_vertprof(2,:)));
    p2=line(smooth(x2(1:65)),-depth(1:65),'linewidth',2,'color','r','LineStyle','-'); %TOD mean
    
    grid on
    ax1 = gca; % current axes
    ax1.XColor = 'r';
    ax1.YColor = 'r';
    ax1.LineWidth=2;
    ax1.XLim=[axminred axmaxred];
    ax1.YLim=[-6000 0];
    ax1.GridColor=[.5 .5 .5];
    ax1.YMinorTick='on';
    ax1.MinorGridColor=[.85 .85 .85];
    ax1.MinorGridLineStyle='-';
    ax1.MinorGridAlpha=0.3;
    ax1.GridAlpha=0.4;
    ax1.FontName='Helvetica Neue';
    ax1.FontSize=12;
    ax1.YLabel.String='Depth (km)';
    ax1.XLabel.String={'ToD (yr)'};
    ax1.XTick=linspace(axminred,axmaxred,5);
    ax1.YTick=-6000:1000:-500;
    
    ax1.YTickLabel={'-6  '; '-5  ';'-4  ';'-3  ';'-2  ';'-1  '};
    
    
    grid(gca,'minor')
   
    
    ax1_pos = ax1.Position; % position of first axes
    
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');
    hold on
    
    
    x10=smooth(cell2mat(Trec_vertprof(1,:)));
    x30=smooth(cell2mat(Trec_vertprof(3,:)));

    %filled area
    X = [x10(1:65).', fliplr(x30(1:65).')];
    Y = [-depth(1:65), fliplr(-depth(1:65))];
    shaded2=fill(smooth(X), smooth(Y), [102/255 178/255 225/255])
    shaded2.FaceAlpha=0.4;
    shaded2.LineStyle='none';


    hold on
    x20=smooth(cell2mat(Trec_vertprof(2,:)));
    p20=line(smooth(x20(1:65)),-depth(1:65),'linewidth',2,'color','b','LineStyle','-'); %TREC
    grid on
    
    ax2 = gca; % current axes
    ax2.XColor = 'b';
    ax2.YColor = 'b';
    ax2.LineWidth=2;
    ax2.XLim=[axminblue axmaxblue];
    ax2.YLim=[-6000 0];
    ax2.GridColor=[.5 .5 .5];
    ax2.YMinorTick='on';
    ax2.MinorGridColor=[.85 .85 .85];
    ax2.MinorGridLineStyle='-';
    ax2.MinorGridAlpha=0.3;
    ax2.GridAlpha=0.4;
    ax2.FontName='Helvetica Neue';
    ax2.FontSize=12;
    ax2.YLabel.String='';
    ax2.XLabel.String={'Trec (yr)'};
    ax2.XTick=linspace(axminblue,axmaxblue,5);
    ax2.YTick=-6000:1000:-500;
    ax2.YTickLabel={'-6  '; '-5  ';'-4  ';'-3  ';'-2  ';'-1  '};
    grid(gca,'minor')
    ax2.XMinorGrid='off';

    
%     lgd=legend([p2 p20],'ToD','Trec')
%     lgd.Location='southeast';
%     lgd.Box='on';
%     lgd.Color=[1 1 1];
%     
%     
%     uistack(lgd,'top')

%title(sprintf('d) %s ToD vertical profile NAtl',varname),'fontsize',14)


%% PANEL TOD section
panel2=subplot(8,6,[4,5,6,10,11,12])

%getting section data
[sect_fld_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTOD2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:7 %smoothing 7 times
    sect_fld_original=smooth2d(sect_fld_original);
end

%section plot
pcolor(sect_lat,-sect_dep,sect_fld_original);
shading interp
hold on

%contourlines
[sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTOD2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:7 %smoothing 7 times
    sect_ctr_original=smooth2d(sect_ctr_original);
end

[M,c]=contour(sect_lat,-sect_dep,sect_ctr_original);

c.LevelList=contourlevels_section; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
c.LineWidth=1.5;
clabel(M,c,'FontSize',10) 

colormap(colorscheme)
caxis([axminTOD axmaxTOD])

h=colorbar;
h.Location='eastoutside';
h.Label.String=sprintf('%s ToD (yr)',varname);
h.Label.Position=labelpos;
h.Ticks=[0:tickjumpTOD:axmaxTOD];
%h.TickLabels(end)={['\geq' num2str(axmax)]};
h.Label.FontSize=12;

xlabel('Latitude (\circ)')
ylabel('Depth (km)')
title(['c) ' sprintf('%s ToD ',varname)],'fontsize',14)

axis([-10 65 -6000 0])


set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in', ...
    'TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on','xtick',-10:10:60, ... 
    'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'},...
        'ytick',-6000:1000:0,'ytickLabel',{'-6  ';'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '});


%========INSET=================
insetaxis=axes('Position',[.84 .59+0.08 .06 .18]); hold on;
m_proj('hammer-aitoff','lon',[-70 -5],'lat',[-10 75]);%ATL
lonx=-63:-5;latx2=25:65;field=NaN(length(latx2),length(lonx));field(:,-45-lonx(1):-30-lonx(1))=1;
m_plot(plon(atl_ind),plat(atl_ind),'-r','linewidth',2);
m_coast('patch',[.7 .7 .7]); m_grid('ytick',0:15:60,'yaxislocation','right','xtick',-90:30:0,'yticklabel',[],'xticklabel',[],'fontsize',7);
%========INSET=================

%%
panel3=subplot(8,6,[22,23,24,28,29,30])
fieldTrec2(fieldTrec2==bottomflag)=NaN; %this makes the points representing the bottom 'transparent'

%changing all the negative flags to avoid the concentric rings when interpolating
fieldTrec2(fieldTrec2<0)=600;

%getting section data
[sect_fld_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTrec2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:7 %smoothing 7 times
    sect_fld_original=smooth2d(sect_fld_original);
end

%section plot
pcolor(sect_lat,-sect_dep,sect_fld_original);
shading interp
hold on

%contourlines
[sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTrec2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:7 %smoothing 7 times
    sect_ctr_original=smooth2d(sect_ctr_original);
end

[M,c]=contour(sect_lat,-sect_dep,sect_ctr_original);

c.LevelList=isolinessct; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
c.LineWidth=1.5;
clabel(M,c,'FontSize',10) 

colormap(gca,colorscheme2)
caxis([axminTrec axmaxTrec])

h=colorbar;
h.Location='eastoutside';
h.Label.String=sprintf('%s Trec (yr)',varname);
%h.Label.Position=labelpos;
h.Ticks=[120 160:tickjumpTrec:axmaxTrec];
h.TickLabels(1)={''};
h.Label.FontSize=12;

xlabel('Latitude (\circ)')
ylabel('Depth (km)')
title(['d) ' sprintf('%s Trec ',varname)],'fontsize',14)

axis([-10 65 -6000 0])


set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in', ...
    'TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on','xtick',-10:10:60, ... 
    'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'},...
        'ytick',-6000:1000:0,'ytickLabel',{'-6  ';'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '});

%%
clear x* y* X Y
panel4=subplot(8,6,[31,32,37,38,43,44])

%filled area
    x4=smooth(cell2mat(Pchange_vertprof(3,:)));
    x5=smooth(cell2mat(Pchange_vertprof(3,:)));
    X = [x4(:).', fliplr(x5(:).')];
    Y = [-depth_interp, fliplr(-depth_interp)];
    shaded1=fill(X, Y, [255/255 204/255 204/255]);
    shaded1.FaceAlpha=0.4;
    shaded1.LineStyle='none';
    
 
    
    x1=smooth(cell2mat(Pchange_vertprof(1,:)));
    p1=line(x1(1:56),-depth_interp(1:56),'linewidth',2.5,'color',[0 158/255 115/255],'LineStyle',':','MarkerSize',1); %TOD mean
    hold on
    
    x2=smooth(cell2mat(Pchange_vertprof(2,:)));
    p2=line(x2(1:56),-depth_interp(1:56),'linewidth',3,'color',[230/255 159/255 0],'LineStyle','-.'); %TOD mean
    hold on
    
    x3=smooth(cell2mat(Pchange_vertprof(3,:)));
    p3=line(x3(1:56),-depth_interp(1:56),'linewidth',3,'color',[0 114/255 178/255],'LineStyle','-'); %TOD mean
    
    zeroline=xline(0);
    zeroline.LineWidth=3;
    
    lgd=legend([p1 p2 p3],{'250-280';'350-380';'450-480'})
    title(lgd,{'Simulation','period (yr-yr)'},'Fontsize',10)
    lgd.Location='southeast';
    
    grid on
    ax1 = gca; % current axes
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    ax1.LineWidth=2;
    ax1.XLim=[axminpchange axmaxpchange];
    ax1.YLim=[-6000 0];
    ax1.GridColor=[.5 .5 .5];
    ax1.YMinorTick='on';
    ax1.MinorGridColor=[.85 .85 .85];
    ax1.MinorGridLineStyle='-';
    ax1.MinorGridAlpha=0.3;
    ax1.GridAlpha=0.4;
    ax1.FontName='Helvetica Neue';
    ax1.FontSize=12;
    ax1.YLabel.String='Depth (km)';
    ax1.XLabel.String='Percentage change (%)';
    ax1.XTick=linspace(axminpchange,axmaxpchange,5);
    ax1.YTick=-6000:1000:-500;
    ax1.YTickLabel={'-6  '; '-5  ';'-4  ';'-3  ';'-2  ';'-1  '};

    grid minor

title(sprintf('%s_%%_c_h_a_n_g_e',varname),'fontsize',12)

%% pchange section
panel5=subplot(8,6,[40,41,42,46,47,48])

field2=Pchange_flipped{1,3};
[sectfield,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(field2,depth_interp,plat,atl_ind); %this gets the orginial section field
dummyfield=sectfield;

for smoothnumber=1:smoothnumbsect %smoothing 7 times at the specific depth level
    dummyfield=smooth2d(dummyfield);
end

pcolor(sect_lat,-sect_dep,sectfield)
shading interp
hold on

%% stippling mask for section


pairs={'pair1';'pair2';'pair3'};
pairIDX=3; %window 450-480

Pflags=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D_vertinterp.mat',varlist{varIDX})...
    ,sprintf('Pflags_%s',pairs{pairIDX}));
Pflags1=eval(sprintf('Pflags.Pflags_%s',pairs{pairIDX}));

%flipping and breaking atlantic
Ptest_ind=permute(Pflags1,[3,1,2]);
bit1=Ptest_ind(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=Ptest_ind(:,:,1:121);   %bit2 (1:199,:,:);
Ptest_ind=cat(3,bit1,bit2);


% 1st make section field
[sect_lat2,sect_depth,section_stpl_field]=makesect_atl_noresm_fixed_z_stippl(Ptest_ind,depth_interp,plat,atl_ind);
mask2=section_stpl_field==100; %this is the stippling mask where significant differences were labelled with flag=100
mask2(1,:)=0; mask2(2,:)=0;
stipple(sect_lat,-sect_dep,mask2,'color','k','density',500,'markersize',5);
hold on

colormap(gca,colorscheme3)
caxis([axminpchangesect axmaxpchangesect])

h=colorbar;
h.Location='eastoutside';
h.Label.String=sprintf('%s_%%_c_h_a_n_g_e %s',varname,'(%)');
h.Ticks=[axminpchangesect:tickjumppchange:axmaxpchangesect];
h.Label.FontSize=12;


xlabel('Latitude (\circ)')
ylabel('Depth (km)')
title(['e) ' sprintf('%s_%%_c_h_a_n_g_e ',varname)],'fontsize',14)

axis([-10 65 -6000 0])



set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in', ...
    'TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on','xtick',-10:10:60, ... 
    'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'},...
        'ytick',-6000:1000:0,'ytickLabel',{'-6  ';'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '});

increment=.035; %%using to adjust postions
decrement=0.04;
wideincre=0.17;

panel1.Position=[0.08    0.5323+decrement    0.20    0.3927-decrement]
ax2.Position=panel1.Position;

panel4.Position(1)=panel1.Position(1); 
panel4.Position(3)=panel1.Position(3);
panel4.Position(4)=panel1.Position(4);

panel2.Position=[0.40                   0.7434-increment          0.3356+wideincre    0.1816+increment];
panel3.Position=[panel2.Position(1)     0.4267-increment+0.02          0.3356+wideincre    0.1816+increment]; 
panel5.Position=[panel2.Position(1)     0.1100          0.3356+wideincre    0.1816+increment];

panel2.XLabel.String=[];
panel3.XLabel.String=[];

set(gcf,'color','w')
hFig1.InvertHardcopy='off';

a=annotation('textbox','String','a)','fontsize',14,'Edgecolor','none','FontName','Helvetica Neue','FontWeight','bold');
a.Position=[0.0220     0.9206    0.0239    0.0348]

b=annotation('textbox','String','b)','fontsize',14,'Edgecolor','none','FontName','Helvetica Neue','FontWeight','bold');
b.Position=[0.0220     0.9206-0.46    0.0239    0.0348]


anotfase=annotation('textbox','String','@ End of Extension (years 450-480) ','fontsize',12,'Edgecolor','none','FontName','Helvetica Neue',...
    'FontWeight','bold','Color','w','HorizontalAlignment','center');
anotfase.Position=[0.7500    0.1271    0.1525    0.0348]


%print(sprintf('%s1.png',varname),'-dpng','-r300')
print(sprintf('%s1.pdf',varname),'-dpdf','-r300', '-bestfit')
print(sprintf('%s1.svg',varname),'-dsvg','-r300')

% print('Omegac1.pdf','-dpdf','-r300')
% print('Omegac1.svg','-dsvg','-r300')

