

%%
%=====================================================================================================================================================================
%=====================================================================================================================================================================
%=====================================================================================================================================================================

depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';

model = '';
mode = {"Tropical","Subtropical", "Subpolar"};
variables = {"templvl", "pH","o2","AOU","omegac"};


for i = 1:length(variables)
dummy=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/Testing_new_vol_avg/area_and_vol_NorESM_%s_lat_bands.mat', variables{i}));

%dummy=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_%s_lat_bands.mat', variables{i}));
data(1,i) = {dummy.vol_weighted_avg.tropical};
data(2,i) = {dummy.vol_weighted_avg.subtropical};
data(3,i) = {dummy.vol_weighted_avg.subpolar};

end


new_avg = [data{1,1} + data{1,2} + data{1,3}]./3;


for m = 1:length(mode)
    
if mode{m} == "Tropical"
    
    var_temperature = data{1,1};
    var_pH = data{1,2};
    var_o2 = data{1,3};
    var_AOU = data{1,4};
    var_omegaC = data{1,5};
    
    letters = {"a";"b";"c";"d";"e"};

elseif mode{m} == "Subtropical"
    var_temperature = data{2,1};
    var_pH = data{2,2};
    var_o2 = data{2,3};
    var_AOU = data{2,4};
    var_omegaC = data{2,5};
    
    letters = {"f"; "g"; "h"; "i";"j"};
    
elseif mode{m} == "Subpolar"
    
    fill_gap_temp = data{3,1};
    var_temperature = [fill_gap_temp(1:62,:); repmat(fill_gap_temp(62,:),70-62,1)];
    
    var_pH = data{3,2};
    var_o2 = data{3,3};
    
    fill_gap_AOU = data{3,4};
    var_AOU = [fill_gap_AOU(1:62,:); repmat(fill_gap_AOU(62,:),70-62,1)];
   
    var_omegaC = data{3,5};
    
    letters = {"k";"l";"m";"n";"o"};
    
end
    
    
fig=figure('Visible', 'off', 'color', 'white','units','centimeters', 'Position', [   -2.8331   28.5750   15.2047   32.1028])

%Temperature section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
panel1=subplot(5,1,1)
set(gca, 'Units', 'centimeters')
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_templvl_lat_bands.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
pcolor(1:481,-depth,var_temperature);shading interp

caxis([0 15])
h=colorbar;
h.Label.String='Temp (\circC)';
h.FontSize=12;
h.Ticks=[0 5 10 15];
h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
colormap(gca,parula(15))
xlim([0 480])
ylim([-5000 0])
%xlabel('Time (yr)')
ylabel('Depth (Km)')
title(sprintf('%s) %s Evolution of NAtl mean Temp',letters{1},model),'FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0,'ytickLabel',{'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '})% ...
    %'Position',[    0.1501    0.85    0.6364    0.1137]);
hold on
l1=xline(140,'--','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);
l2=xline(280,'--','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);

%%
%pH section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
panel2=subplot(5,1,2)
set(gca, 'Units', 'centimeters')
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_ph.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
pcolor(1:481,-depth,var_pH);shading interp

caxis([7.2 8.2])
h=colorbar;
h.Label.String='pH';
h.FontSize=12;
h.Ticks=[7.2 7.7 8.2];

h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
colormap(gca,parula(10));
xlim([0 480])
ylim([-5000 0])

%xlabel('Time (yr)')
ylabel('Depth (Km)')
title(sprintf('%s) %s Evolution of NAtl mean pH',letters{2},model),'FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0,'ytickLabel',{'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '})% ...
   % 'Position',[    0.1501    0.68    0.6364    0.1137])
hold on
l1=xline(140,'--','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);
l2=xline(280,'--','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);

%%
%O2 section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
panel3=subplot(5,1,3)
set(gca, 'Units', 'centimeters')
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_o2.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
pcolor(1:481,-depth,(var_o2/1024)*1000000);shading interp

caxis([100 340])
h=colorbar;
h.Label.String='DO (\mumol O_2 kg^-^1)';
h.FontSize=12;

h.Location='eastoutside';
h.Ticks=100:60:340;
h.TickLabels(end)={['\geq' num2str(max(caxis))]};
colormap(gca,parula(12));
xlim([0 480])
ylim([-5000 0])

%xlabel('Time (yr)')
ylabel('Depth (Km)')
title(sprintf('%s) %s Evolution of NAtl mean DO',letters{3},model),'FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0,'ytickLabel',{'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '})% ...
    %'Position',[    0.1501    0.50    0.6364    0.1137])
hold on
l1=xline(140,'--','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);
l2=xline(280,'--','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);

%%
%AOU section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
panel4=subplot(5,1,4)
set(gca, 'Units', 'centimeters')
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_AOU.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
pcolor(1:481,-depth,var_AOU);shading interp

caxis([0 200])
h=colorbar;
h.Ticks=0:40:200
h.Label.String='AOU (\mumol O_2 kg^-^1)';
h.FontSize=12;
h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
colormap(gca,parula(10))
xlim([0 480])
ylim([-5000 0])

%xlabel('Time (yr)')
ylabel('Depth (Km)')
title(sprintf('%s) %s Evolution of NAtl mean AOU',letters{4},model),'FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0,'ytickLabel',{'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '})% ...
%set(gca,'Position',[0.45,0.11,0.42,0.3412]);
hold on 
l1=xline(140,'--','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);
l2=xline(280,'--','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);


%OMEGA
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
panel5=subplot(5,1,5)
set(gca, 'Units', 'centimeters')
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_omegac.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
var_omegaC= var_omegaC.*0.95;
pcolor(1:481,-depth,var_omegaC);shading interp

caxis([0 2])
colorscheme=flip(redblue);
colorscheme2=[repmat(colorscheme(1,:),6,1); colorscheme(1:50,:);colorscheme(72:end,:)];
colormap(gca,colorscheme2)


% colorscheme=flip(redblue);
% colorscheme2 = [repmat(colorscheme(1,:),38,1);...
%                 colorscheme(1:35,:);...
%                 colorscheme(75:120,:);...
%                 repmat(colorscheme(120,:),25,1)];
          
colormap(gca,colorscheme2)
%colormap(gca,flip(redblue(80),2))

h=colorbar;
h.Label.String='\Omega_c';
h.Ticks=0:0.5:2;
h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
h.FontSize=12;

xlim([0 480])
ylim([-5000 0])
xlabel('Time (yr)')
ylabel('Depth (Km)') 
title(sprintf('%s) %s Evolution of NAtl mean \\Omega_c',letters{5},model),'FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0,'ytickLabel',{'-5  ';'-4  ';'-3  ';'-2  ';'-1  ';' 0  '})% ...
%set(gca,'Position',[0.45,0.11,0.42,0.3412]);
hold on
l1=xline(140,'--','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);
l2=xline(280,'--','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold','Linewidth',1);

set(gcf,'color','w')

print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/fig2_NAtl_%s_%s.png',model, mode{m}),'-dpng','-r300')
print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/fig2_NAtl_%s_%s.svg',model, mode{m}),'-dsvg','-r300')

end
%=====================================================================================================================================================================
%=====================================================================================================================================================================
%=====================================================================================================================================================================

% % % %% VIDEO by screenshots
% % %    figure
% % %    axis tight manual 
% % %    set(gca,'nextplot','replacechildren');
% % %    v = VideoWriter('RUP.avi');
% % %    open(v)
% % %    for i = 1:140
% % %        % Example of plot
% % %        set(gcf,'Position',[10 40 1200 900],'color','w')
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        rampup=plot(1:i,co2annual(1:i),'r','linewidth',10);
% % %        yline(co2start,'--', 'Pre-industrial baseline 284.7 ppm','LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% % %        if i==140
% % %        yline(max(co2annual),'--',['CO_2 peak ' sprintf('%.2f ppm',max(co2annual))],'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% % %        end
% % %        lgd=legend([rampup],'Ramp-up')
% % %        set(lgd,'Position',[0.195,0.80,0.1133,0.0576],'FontSize',20)
% % %        xlabel('Time (yr)')
% % %        ylabel('Atmospheric CO_2 (ppm)')
% % %        grid on
% % %        title('a) Atmospheric CO_2 forcing','FontSize',20)
% % %        set(gca,'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
% % %            'xtick',0:50:500)
% % %        % Store the frame and save picture file
% % %        print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RUP/img_%d',i),'-dpng','-r300')
% % %        M=getframe(gcf);% leaving gcf out crops the frame in the movie.
% % %        writeVideo(v,M);
% % %    end
% % %    
% % %   rampdownseries=co2annual(140:280);
% % %   set(gca,'nextplot','replacechildren');
% % %    v = VideoWriter('RDOWN.avi');
% % %    open(v)  
% % %    for j=1:length(rampdownseries)
% % %        % Example of plot
% % %        set(gcf,'Position',[10 40 1200 900],'color','w')
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        rampdown=plot(140:140+j-1,rampdownseries(1:j),'b','linewidth',4);
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        yline(co2start,'--', 'Pre-industrial baseline 284.7 ppm','LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% % %        yline(max(co2annual),'--',['CO_2 peak ' sprintf('%.2f ppm',max(co2annual))],'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% % %        lgd=legend([rampdown],'Ramp-down');
% % %        set(lgd,'Position',[0.40,0.80,0.1133,0.0576],'FontSize',14)
% % %        xlabel('Time (yr)')
% % %        ylabel('Atmospheric CO_2 (ppm)')
% % %        grid on
% % %        title('a) Atmospheric CO_2 forcing','FontSize',14)
% % %        set(gca,'FontName','Helvetica','fontsize',14,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
% % %            'xtick',0:50:500)
% % %        set(gcf,'color','w')
% % %        % Store the frame
% % %        K=getframe(gcf); % leaving gcf out crops the frame in the movie.
% % %        writeVideo(v,K);
% % %    end
% % % 
% % %    close(v)
% % %  
% % %    
% % %  extension=co2annual(280:end);
% % %   set(gca,'nextplot','replacechildren');
% % %    v = VideoWriter('Extention.avi');
% % %    open(v)  
% % %    for j=1:length(extension)
% % %        % Example of plot
% % %        set(gcf,'Position',[10 40 1200 900],'color','w')
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        rampdown=plot(280:280+j-1,extension(1:j),'c','linewidth',4);
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %         yline(co2start,'--', 'Pre-industrial baseline 284.7 ppm','LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% % %        yline(max(co2annual),'--',['CO_2 peak ' sprintf('%.2f ppm',max(co2annual))],'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% % %        lgd=legend([rampdown],'Extention');
% % %        set(lgd,'Position',[0.60,0.80,0.1133,0.0576],'FontSize',14)
% % %        xlabel('Time (yr)')
% % %        ylabel('Atmospheric CO_2 (ppm)')
% % %        grid on
% % %        title('a) Atmospheric CO_2 forcing','FontSize',14)
% % %        set(gca,'FontName','Helvetica','fontsize',14,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
% % %            'xtick',0:50:500)
% % %        set(gcf,'color','w')
% % %        % Store the frame
% % %        K=getframe(gcf); % leaving gcf out crops the frame in the movie.
% % %        writeVideo(v,K);
% % %    end
% % % 
% % %    close(v)
% % % 
% % % %%
% % % %% VIDEO making by image series
% % % %1
% % %    fig=figure;
% % %    axis tight manual 
% % %    set(fig, 'Visible', 'off');
% % %    set(gca,'nextplot','replacechildren');
% % %    for i = 1:140
% % %        % Example of plot 
% % %        set(gcf,'Position',[10 40 1200 900],'color','w')
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        rampup=plot(1:i,co2annual(1:i),'r','linewidth',10);
% % %        yline(co2start,'--', 'Pre-industrial baseline 284.7 ppm','LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold');
% % %        if i==140
% % %        yline(max(co2annual),'--',['CO_2 peak ' sprintf('%.2f ppm',max(co2annual))],'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold');
% % %        end
% % %        lgd=legend([rampup],'Ramp-up');
% % %        set(lgd,'Position',[0.195,0.80,0.1133,0.0576],'FontSize',20)
% % %        xlabel('Time (yr)')
% % %        ylabel('Atmospheric CO_2 (ppm)')
% % %        grid on
% % %        title('a) Atmospheric CO_2 forcing','FontSize',20)
% % %        set(gca,'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
% % %            'xtick',0:50:500)
% % %        % Store the frame and save picture file
% % %        print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RUP/img_%d',i),'-dpng','-r300')
% % %    end
% % %    
% % %    %creating video from image sequence
% % %   workingDir = '/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RUP/';
% % %   outputVideo = VideoWriter('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RUP/RUP_sequence.avi');
% % %   outputVideo.FrameRate = 30;
% % %   open(outputVideo)
% % %   
% % %   for ii = 1:140
% % %    img = imread(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RUP/img_%d.png',ii));
% % %    writeVideo(outputVideo,img)
% % %   end
% % % 
% % %   close(outputVideo)
% % %   
% % %   %2
% % %    fig=figure;
% % %    axis tight manual 
% % %    set(fig, 'Visible', 'off');
% % %    set(gca,'nextplot','replacechildren');
% % %   for j=1:length(rampdownseries)
% % %        % Example of plot
% % %        set(gcf,'Position',[10 40 1200 900],'color','w')
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        rampdown=plot(140:140+j-1,rampdownseries(1:j),'b','linewidth',10);
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        yline(co2start,'--', 'Pre-industrial baseline 284.7 ppm','LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold');
% % %        yline(max(co2annual),'--',['CO_2 peak ' sprintf('%.2f ppm',max(co2annual))],'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold');
% % %        lgd=legend([rampdown],'Ramp-down');
% % %        set(lgd,'Position',[0.40,0.80,0.1133,0.0576],'FontSize',20)
% % %        xlabel('Time (yr)')
% % %        ylabel('Atmospheric CO_2 (ppm)')
% % %        grid on
% % %        title('a) Atmospheric CO_2 forcing','FontSize',20)
% % %        set(gca,'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
% % %            'xtick',0:50:500)
% % %        set(gcf,'color','w')
% % %        % Store the frame
% % %        print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RDOWN/img_%d',140+j-1),'-dpng','-r300')
% % %    end
% % %   
% % %   %creating video from image sequence
% % %   workingDir = '/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RDOWN/';
% % %   outputVideo = VideoWriter('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RDOWN/RDOWN_sequence.avi');
% % %   outputVideo.FrameRate = 30;
% % %   open(outputVideo)
% % %   
% % %   for ii = 140:280
% % %    img = imread(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/RDOWN/img_%d.png',ii));
% % %    writeVideo(outputVideo,img)
% % %   end
% % % 
% % %   close(outputVideo)
% % %   
% % %   
% % %   %3
% % %    fig=figure;
% % %    axis tight manual 
% % %    set(fig, 'Visible', 'off');
% % %    set(gca,'nextplot','replacechildren');
% % %    
% % %     extension=co2annual(280:end);
% % %   set(gca,'nextplot','replacechildren');
% % %    v = VideoWriter('Extention.avi');
% % %    open(v)  
% % %    for j=1:length(extension)
% % %        % Example of plot
% % %        set(gcf,'Position',[10 40 1200 900],'color','w')
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %        rampdown=plot(280:280+j-1,extension(1:j),'c','linewidth',10);
% % %        xlim([0 500])
% % %        ylim([200 1350])
% % %         yline(co2start,'--', 'Pre-industrial baseline 284.7 ppm','LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold');
% % %        yline(max(co2annual),'--',['CO_2 peak ' sprintf('%.2f ppm',max(co2annual))],'LabelVerticalAlignment','bottom','FontSize',12,'FontWeight','bold');
% % %        lgd=legend([rampdown],'Extention');
% % %        set(lgd,'Position',[0.60,0.80,0.1133,0.0576],'FontSize',20)
% % %        xlabel('Time (yr)')
% % %        ylabel('Atmospheric CO_2 (ppm)')
% % %        grid on
% % %        title('a) Atmospheric CO_2 forcing','FontSize',20)
% % %        set(gca,'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
% % %            'xtick',0:50:500)
% % % 
% % %        % Store the frame
% % %        print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/EXTENSION/img_%d',280+j-1),'-dpng','-r300')
% % % 
% % %    end
% % % 
% % % %creating video from image sequence
% % %   workingDir = '/Volumes/LaCie_Leonardo/NorESM/Initial_figs/EXTENSION/';
% % %   outputVideo = VideoWriter('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/EXTENSION/EXT_sequence.avi');
% % %   outputVideo.FrameRate = 30;
% % %   open(outputVideo)
% % %   
% % %   for ii = 280:480
% % %    img = imread(sprintf('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/EXTENSION/img_%d.png',ii));
% % %    writeVideo(outputVideo,img)
% % %   end
% % % 
% % %   close(outputVideo)
  