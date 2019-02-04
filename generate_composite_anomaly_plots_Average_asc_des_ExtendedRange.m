clear all
close all
clc
fclose all;

addpath /data02/transfer/Chris/mfile_library/
addpath /data02/transfer/Chris/mfile_library/m_map/
addpath /data01/pisco/ucsb/shared-files/oceanography/mfewings/matlab_bits/ % path to c2h.m
addpath ../NARR_analysis/

set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultaxesfontweight','bold')
set(0,'defaultaxeslinewidth',1)

zooooom = 0; % 1 = this switch no longer used

% take average of grandMeans confidence for ascending and descending
gnd_dir = '/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/mat/';
g_asc = load([gnd_dir,'Qscat_JPL_ascending_grandMean.mat']); % this file created in generate_data_mean_plots.m 
g_des = load([gnd_dir,'Qscat_JPL_descending_grandMean.mat']);
gU95 = sqrt(g_asc.gndmean.rot_Uconf95.^2 + g_des.gndmean.rot_Uconf95.^2)./2;
gV95 = sqrt(g_asc.gndmean.rot_Vconf95.^2 + g_des.gndmean.rot_Vconf95.^2)./2;
ny = g_asc.ny;
nx = g_asc.nx;

wdw = 6:10; % days (enter a day here -10 to 5 and also edit 'sequence' below to list files in order)

for ww = 1:length(wdw)

sequence = sprintf('%02d',wdw(ww)+11); % used in fname to line the sequence up in directories

% SINCE WE'VE SAVED OUT THESE FILES WE CAN PERFORM THE AVERAGE FOR CONFU V
% HERE AS WELL AS AVERAGING DU AND DV
indir_asc = '/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/anomaly_component_composite/JPL/ne_pacific/ascending/mat/';
indir_des = '/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/anomaly_component_composite/JPL/ne_pacific/descending/mat/';
in_asc = load([indir_asc,'Qscat_JPL_windStress_anomaly',sequence,'_window_',num2str(wdw(ww)),'_ascending_Data.mat']);
in_des = load([indir_des,'Qscat_JPL_windStress_anomaly',sequence,'_window_',num2str(wdw(ww)),'_descending_Data.mat']);
dU = (in_asc.dU + in_des.dU)./2; % rotated anomaly
dV = (in_asc.dV + in_des.dV)./2;

acu = in_asc.confU(:,1)-in_asc.dU; % reconstruct just the error estimate here
acv = in_asc.confV(:,1)-in_asc.dV;
dcu = in_des.confU(:,1)-in_des.dU; 
dcv = in_des.confV(:,1)-in_des.dV;

confU(:,1) = dU + (sqrt(acu.^2 + dcu.^2)./2); % rotated anomaly
confU(:,2) = dU - (sqrt(acu.^2 + dcu.^2)./2);
confV(:,1) = dV + (sqrt(acv.^2 + dcv.^2)./2);
confV(:,2) = dV - (sqrt(acv.^2 + dcv.^2)./2);

slon = in_asc.qlon;
slat = in_asc.qlat;

lowDClon = [in_asc.lowDClon;in_des.lowDClon];
lowDClat = [in_asc.lowDClat;in_des.lowDClat];

Unotsig   = ~((confU(:,1)<-gU95) | (confU(:,2)>gU95));
Vnotsig   = ~((confV(:,1)<-gV95) | (confV(:,2)>gV95));

Unotsig(isnan(confU(:,1))) = 0;
Vnotsig(isnan(confV(:,1))) = 0;      

%keyboard % data file saved out for melanie and all plotting function commented below
matdir = ['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/anomaly_component_composite/JPL/ne_pacific/average/mat/'];
qlon = slon(:);
qlat = slat(:);
save([matdir,'Qscat_JPL_windStress_anomaly',sequence,'_window_',num2str(wdw(ww)),'_average_Data'],'qlon','qlat','dU','dV','Unotsig','Vnotsig','confU','confV','lowDClon','lowDClat','nx','ny')
clear qlon qlat data_mean_*


%% we will generate 2 figures
%  one for dU and another for dV

figure(1)  % dU
set(gcf,'units','normalized','position',[0 0 1 1],'color','w','PaperPosition',[0 0 8.25 10],'renderer','painters')
  
cmap = c2h(100); % cold (blue) to red (hot) white in the middle
cmap = brighten(cmap,-.5);
colormap(cmap) 

    subplot('position',[.1 .1 .55 .8]) 
    
        if zooooom
            %m_proj('lambert','long',[220 250],'lat',[25 51]);
        else
            m_proj('lambert','long',[228 242],'lat',[30 45]);
        end
        
        hold on
                
        m_plot(lowDClon,lowDClat,'.','color',[.9 .9 .9]) % plot regions of low data coverage
        
        m_usercoast('/home/pisco/mfewings/nasaUpwellRelax/data/N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
        
        if zooooom
            %m_grid('box','on','tickdir','out','linestyle','none', ...
               %'xtick',[210 250],'ytick',[30 51]);
        else
            m_grid('box','on','tickdir','out','linestyle','none', ...
               'xtick',[228 242],'ytick',[30 45]);
        end
            
        % set up contour variables
        X = reshape(slon,ny,nx);
        X = X+.05; % centers pcolor boxes
        Y = reshape(slat,ny,nx);
        Z = reshape(dU,ny,nx); % include an offset to center on c2h colormap 
    
        pp = m_pcolor(X,Y,Z);
        set(pp,'linestyle','none')
        caxis([-.2 .2])        
        
        [~,ch] = m_contour(X,Y,Z,[-0.05,-0.05]); % add contour lines for +- 0.05 Pa
        set(ch,'linewidth',2,'linecolor','b');
        
        [~,ch] = m_contour(X,Y,Z,[0.05,0.05]);
        set(ch,'linewidth',2,'linecolor','r');
        
        m_plot(slon(Unotsig),slat(Unotsig),'+','markersize',4	,'color',[.65 .65 .65])    
 
        title(num2str(wdw(ww)),'fontsize',20,'fontweight','bold','verticalalignment','bottom')
        
        fname = ['qwkskt_JPL_windStress_anomaly',sequence,'_window_',num2str(wdw(ww)),'_AscDesAverage_dU'];       
        
        set(gca,'layer','top')
        
%         xtk = linspace(-.2,.2,26); 
%         xtklbl = {'-0.20';'';'';'';'';'-0.12';'';'';'';'';'-0.04';'';'';'';'';'0.04';'';'';'';'';'0.12';'';'';'' ;'';'0.20'};
%         colorbar_cmap(xtk,xtklbl, [.75 .2 .015 .6],'wind stress U anomaly (Pa)',cmap)
        
        if zooooom
            %fname = [fname,'_zoom'];
        end
        
    drawnow
    saveas(gcf,['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/anomaly_component_composite/JPL/ne_pacific/average/U/png/',fname],'png')
    print('-depsc2',['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/anomaly_component_composite/JPL/ne_pacific/average/U/eps/',fname])    
    

close all

    
figure(1)  % dV
set(gcf,'units','normalized','position',[0 0 1 1],'color','w','PaperPosition',[0 0 8.25 10],'renderer','painters')
  
cmap = c2h(100); % cold (blue) to red (hot) white in the middle
cmap = brighten(cmap,-.5);
colormap(cmap) 

    subplot('position',[.1 .1 .55 .8]) 
    
        if zooooom
            %m_proj('lambert','long',[220 250],'lat',[25 51]);
        else
            m_proj('lambert','long',[228 242],'lat',[30 45]);
        end
        
        hold on
        
        m_plot(lowDClon,lowDClat,'.','color',[.9 .9 .9])
        
        m_usercoast('/home/pisco/mfewings/nasaUpwellRelax/data/N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
        
        if zooooom
            %m_grid('box','on','tickdir','out','linestyle','none', ...
               %'xtick',[210 250],'ytick',[30 51]);
        else
            m_grid('box','on','tickdir','out','linestyle','none', ...
               'xtick',[228 242],'ytick',[30 45]);
        end         
    
        % set up contour variables
        X = reshape(slon,ny,nx);
        X = X+.05; % centers pcolor boxes
        Y = reshape(slat,ny,nx);
        Z = reshape(dV,ny,nx); % include an offset to center on c2h colormap 
    
        pp = m_pcolor(X,Y,Z);
        set(pp,'linestyle','none')
        caxis([-.2 .2])        
        
        [~,ch] = m_contour(X,Y,Z,[-0.05,-0.05]); % add contour lines for +- 0.05 Pa
        set(ch,'linewidth',2,'linecolor','b');
        
        [~,ch] = m_contour(X,Y,Z,[0.05,0.05]);
        set(ch,'linewidth',2,'linecolor','r');
    
        m_plot(slon(Vnotsig),slat(Vnotsig),'+','markersize',4,'color',[.65 .65 .65])    
        
        title(num2str(wdw(ww)),'fontsize',20,'fontweight','bold','verticalalignment','bottom')
        
        fname = ['qwkskt_JPL_windStress_anomaly',sequence,'_window_',num2str(wdw(ww)),'_AscDesAverage_dV'];       
    
        set(gca,'layer','top')
        
%         xtk = linspace(-.2,.2,26); 
%         colorbar_cmap(xtk,xtklbl, [.75 .2 .015 .6],'wind stress V anomaly (Pa)',cmap)

        if zooooom
            %fname = [fname,'_zoom'];
        end
        
    drawnow
    saveas(gcf,['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/anomaly_component_composite/JPL/ne_pacific/average/V/png/',fname],'png')
    print('-depsc2',['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/anomaly_component_composite/JPL/ne_pacific/average/V/eps/',fname])    
        

close all
   
clear Unotsig Vnotsig f aU aV dU dV in_asc* in_des*

end % of ww wdw loop


    
