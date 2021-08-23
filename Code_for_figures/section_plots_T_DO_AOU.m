clear all;

addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered
addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
addpath /Volumes/LaCie_Leonardo/NorESM/Animations/
addpath /Volumes/LaCie_Leonardo/NorESM

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


atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

GPATH='/Volumes/LaCie_Leonardo/NorESM/Animations/Panel_figs/section_data/';

varlist={'templvl';'DO'; 'AOU'};
ncol=10;


%%%%%GETTING SECTION DATA FOR EACH TIME STEP    

load([GPATH 'sectionfields_t_o2_AOU_timestep_1.mat']); 
ini_T=sect_fld_original_TEMP; ini_DO=(sect_fld_original_o2/1024)*1000000; ini_AOU = sect_fld_original_AOU;
for smoothnumber=1:7
  ini_T=smooth2d(ini_T);
  ini_DO=smooth2d(ini_DO);
  ini_AOU=smooth2d(ini_AOU);
end


tic
for varIDX=1

    folder_name = sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/Sections/%s',varlist{varIDX});
    
    if not(isfolder(folder_name))
        mkdir(folder_name)
    end

    if strcmp(varlist{varIDX},'templvl')==1

        varname='Temp';
        varunit='(\circC)';
        axmin=0;
        axmax=15;
        colorscheme_val=parula(15);
        bartick_vals = axmin:5:axmax;
        bartick_labels_vals = {'0'; '5'; '10'; ['\geq' num2str(max(caxis))]};


        colorscheme=redblue(40);
        colorscheme_delta = colorscheme;
        %colorscheme_delta=[colorscheme(4:8,:);colorscheme(15,:); colorscheme(17,:);colorscheme(20:30,:)];
        axmin_delta=-4;
        axmax_delta=4; 
        bartick_vals_delta = [axmin_delta -3 -2 -1 0 1 2 3 axmax_delta];
        bartick_labels_delta = {string(axmin_delta);'-3';'-2';'-1';'0';'1';'2';'3';string(axmax_delta)};

        smoothnumbsect=3;


    elseif strcmp(varlist{varIDX},'DO')==1

        varname='DO';
        varunit='(\mumol O_2 kg^-^1)';
        axmin = 100;
        axmax = 340;
        bartick_vals=axmin:60:axmax;
        bartick_labels_vals={'100' ;  '160' ;  '220'  ; '280'  ; '340'; ['\geq' num2str(max(caxis))]};
        colorscheme_val=parula(12);

        axmin_delta=-50;
        axmax_delta=50;
        colorscheme_delta=redblue(50);
        bartick_vals_delta = [axmin_delta -30 -10 0 10 30 axmax_delta];
        bartick_labels_delta = {string(axmin_delta);'-30';'-10';'0';'10';'30';string(axmax_delta)};

        smoothnumbsurface=3;
        smoothnumbsect=1;
        smoothnumbcontour=7;

    elseif strcmp(varlist{varIDX},'AOU')==1

 
        varname='AOU';
        varunit='(\mumol O_2 kg^-^1)';
        axmin = 0;
        axmax = 200;
        bartick_vals=axmin:40:axmax;
        bartick_labels_vals={'0' ;  '40' ;  '80'  ; '120'  ; '160'; ['\geq' num2str(max(caxis))]};
        colorscheme_val=parula(10);

        axmin_delta=-50;
        axmax_delta=50;
        colorscheme_delta=redblue(50);
        bartick_vals_delta = [axmin_delta -30 -10 0 10 30 axmax_delta];
        bartick_labels_delta = {string(axmin_delta);'-30';'-10';'0';'10';'30';string(axmax_delta)};

        smoothnumbsurface=3;
        smoothnumbsect=1;
        smoothnumbcontour=7;

    end

    
    

    for tstep=1:480
        fprintf("Year %d var %s \n ",tstep,varlist{varIDX} )
        load([GPATH sprintf('sectionfields_t_o2_AOU_timestep_%d.mat',tstep)]);
        

        %% Doing section series
        fig=figure('visible','off');
        axis tight manual


        if strcmp(varlist{varIDX},'templvl')==1
            sect_fld_original = sect_fld_original_TEMP;
        elseif strcmp(varlist{varIDX},'DO')==1
            sect_fld_original = (sect_fld_original_o2/1024)*1000000;
        elseif strcmp(varlist{varIDX},'AOU')==1
            sect_fld_original = sect_fld_original_AOU;
        end



        %
        for smoothnumber=1:smoothnumbsect %smoothing 3 times
            sect_fld_original=smooth2d(sect_fld_original);
        end
        %

        %%section plot
        pcolor(sect_lat_TEMP,-sect_dep_TEMP,sect_fld_original);shading interp


        set(gca,'FontSize',18,'LineWidth',2);
        ylabel('Depth (Km)')
        axis([-10 65 -6000 0])
        set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in',...
            'TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
            'ytick',-6000:1000:0,'yticklabel',{'-6' '-5' '-4' '-3' '-2' '-1' '0'  },...
            'xtick',-10:10:60,'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'  });

        grid on


        colormap(colorscheme_val)
        caxis([axmin axmax])
        h=colorbar;
        h.Location='eastoutside';
        h.Ticks = bartick_vals;
        h.Label.String=sprintf('%s %s',varname,varunit);
        set(h,'fontsize',12);

        H=annotation('textbox','String', sprintf('Year %d', tstep));
        H.BackgroundColor='w';
        H.Position=[0.6834    0.1411    0.1320    0.0588];
        H.HorizontalAlignment='center';
        H.FontSize=12;
        H.VerticalAlignment='middle';


         %========INSET=================
        a=axes('Position',[.70 .205 .1 .18]); hold on;
        m_proj('hammer-aitoff','lon',[-70 -5],'lat',[-10 75]);%ATL
        lonx=-63:-5;latx2=25:65;field=NaN(length(latx2),length(lonx));field(:,-45-lonx(1):-30-lonx(1))=1;
        m_plot(plon(atl_ind),plat(atl_ind),'-r','linewidth',2);
        m_coast('patch',[.7 .7 .7]); m_grid('ytick',0:15:60,'yaxislocation','right','xtick',-90:30:0,'yticklabel',[],'xticklabel',[],'fontsize',7);
        %========INSET=================

        set(gcf,'color','w')
        figname = sprintf('%s/Section_%s_RawValue_at_year_%d.png',folder_name,varlist{varIDX},tstep);

        %Print fig
        set(gcf,'color','w')
        set(gcf, 'InvertHardCopy', 'off');
        %save fig  
        print(figname,'-dpng','-r300')



%         %% Doing section delta
%         load([GPATH sprintf('sectionfields_t_o2_AOU_timestep_%d.mat',tstep)]);
% 
% 
%         if strcmp(varlist{varIDX},'templvl')==1
%             delta = sect_fld_original_TEMP-ini_T;
%         elseif strcmp(varlist{varIDX},'DO')==1
%             delta =(sect_fld_original_o2/1024 *1000000)-ini_DO;
%         elseif strcmp(varlist{varIDX},'AOU')==1
%             delta = sect_fld_original_AOU-ini_AOU;
%         end
% 
% 
%         fig=figure('visible','off');
%         axis tight manual
%         set(gcf,'nextplot','replacechildren');
% 
%          %
%         for smoothnumber=1:smoothnumbsect %smoothing 3 times
%             sect_fld_original=smooth2d(delta);
%         end
%         %
% 
%         %%section plot
%         pcolor(sect_lat_TEMP,-sect_dep_TEMP,sect_fld_original);shading interp
% 
% 
%         set(gca,'FontSize',18,'LineWidth',2);
%         ylabel('Depth (Km)')
% 
%         axis([-10 65 -6000 0])
%         set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in',...
%             'TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
%             'ytick',-6000:1000:0,'yticklabel',{'-6' '-5' '-4' '-3' '-2' '-1' '0'  },...
%             'xtick',-10:10:60,'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'  });
% 
%         grid on
% 
%         colormap(colorscheme_delta)
%         caxis([axmin_delta axmax_delta])
%         h=colorbar;
%         h.Location='eastoutside';
%         h.Ticks = bartick_vals_delta;
%         h.TickLabels = bartick_labels_delta;
%         h.Label.String=sprintf('\\Delta%s %s',varname,varunit);
%         set(h,'fontsize',12);
% 
%         H=annotation('textbox','String', sprintf('Year %d', tstep));
%         H.BackgroundColor='w';
%         H.Position=[0.6834    0.1411    0.1320    0.0588];
%         H.HorizontalAlignment='center';
%         H.FontSize=12;
%         H.VerticalAlignment='middle';
% 
%         %========INSET=================
%         a=axes('Position',[.70 .205 .1 .18]); hold on;
%         m_proj('hammer-aitoff','lon',[-70 -5],'lat',[-10 75]);%ATL
%         lonx=-63:-5;latx2=25:65;field=NaN(length(latx2),length(lonx));field(:,-45-lonx(1):-30-lonx(1))=1;
%         m_plot(plon(atl_ind),plat(atl_ind),'-r','linewidth',2);
%         m_coast('patch',[.7 .7 .7]); m_grid('ytick',0:15:60,'yaxislocation','right','xtick',-90:30:0,'yticklabel',[],'xticklabel',[],'fontsize',7);
%         %========INSET=================
% 
%         figname = sprintf('%s/Section_%s_Delta_at_year_%d.png',folder_name,varlist{varIDX},tstep);
% 
%         %Print fig
%         set(gcf,'color','w')
%         set(gcf, 'InvertHardCopy', 'off');
%         %save fig  
%         print(figname,'-dpng','-r300')
% 

    end

end
toc

%% Doing video animation

fprintf('Doing Animations \n')
tic
%do animation    
varlist={'templvl';'DO'; 'AOU'};
  for varIDX = 1%:length(varlist)
      workingDir = sprintf('/Volumes/LaCie_Leonardo/NorESM/PAPER_FIGS/Sections/%s',varlist{varIDX}); %where fig files are saved
      
      
      %%Normal Section Video
      outdir = '/Volumes/LaCie_Leonardo/NorESM/Animations/Panel_figs/Sections/';
            if not(isfolder(outdir))
            mkdir(outdir)
            end
      
      outputVideo1 = VideoWriter(['/Volumes/LaCie_Leonardo/NorESM/Animations/Panel_figs/Sections/' sprintf('Section_RawValue_%s.avi',varlist{varIDX})]); %name of the output file
      outputVideo1.FrameRate = 15; %number of frames per seconds 480/15 approx 30 sec
      open(outputVideo1)
      

      for time = 1:480
       img = imread([workingDir sprintf('/Section_%s_RawValue_at_year_%d.png',varlist{varIDX},time)]);
       writeVideo(outputVideo1,img)
      end
      close(outputVideo1)
      
      
      %% Delta section video
      outputVideo2 = VideoWriter(['/Volumes/LaCie_Leonardo/NorESM/Animations/Panel_figs/Sections/' sprintf('Section_Delta_%s.avi',varlist{varIDX})]); %name of the output file
      outputVideo2.FrameRate = 15; %number of frames per seconds 480/15 approx 30 sec
      open(outputVideo2)
      
      for time = 1:480
       img = imread([workingDir sprintf('/Section_%s_Delta_at_year_%d.png',varlist{varIDX},time)]);
       writeVideo(outputVideo2,img)
      end
      close(outputVideo2)
      
      
  end
fprintf('Finished')
toc

    
    