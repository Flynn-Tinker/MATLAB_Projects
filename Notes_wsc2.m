%% Notes on Chris's Step 2
% Wind Stress Curl
%
% Averages the grandmeans into image
%
% Saving everything in May folder
%
% Last edited 03.21.16
%
% Edited for June values
%
clear all
close all
fclose all;
clc

%addpath /data/data02/transfer/Chris/mfile_library/
%addpath /data/data02/transfer/Chris/mfile_library/m_map/
%gnd_dir = '/data/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/mat/';
%fname = 'qwkskt_JPL_windStress_1999-2009_asc_desc_average_w_vectors'; % the output figure name

%%%gnd_dir = '/Volumes/Lacie_Primary/kayla/oaflux_research/quikscat/May/';
gnd_dir = '/Volumes/Lacie_Primary/kayla/oaflux_research/quikscat/wsc/';

wsmax = 0.2; % max wind stress (N/m2) for color scaling

fname = 'norotate_qwkskt_JPL_windStress_2002-2009_average_w_vectors';
fname2 = 'norotate_qwkskt_JPL_windStressCurl_2002-2009_average_w_vectors';

%{
zooooom = 0; % 1 = make axes a bit smaler

if zooooom
    fname = [fname,'_zoom'];
end

ndbc(1).abb = '46011';
ndbc(1).ll = [-120.992, 35];
ndbc(2).abb = '46023';
ndbc(2).ll = [-120.967, 34.714];
ndbc(3).abb = '46054';
ndbc(3).ll = [-120.462, 34.274];
ndbc(4).abb = '46062';
ndbc(4).ll = [-121.01, 35.101];
%}

A = load([gnd_dir,'Qscat_JPL_ascending_grandMean_wsc.mat']); % load up the ascending and descending results and average them
D = load([gnd_dir,'Qscat_JPL_descending_grandMean_wsc.mat']); % end with Chris for May

%%

% note: grand mean is now in units of wind stress
% Finding the Average, though I thought you couldn't do it this way?
wsU = (A.gndmean.U + D.gndmean.U)./2;
wsV = (A.gndmean.V + D.gndmean.V)./2;
%%% Are ascending and descending wsc different...yea they should be,
%%% different swaths
wsC = (A.gndmean.S + D.gndmean.S)./2;

% A and D the same
lon = A.gndmean.lon;
lat = A.gndmean.lat;

nx = A.nx;
ny = A.ny;

lowDClon = lon(isnan(A.gndmean.U+D.gndmean.U));
lowDClat = lat(isnan(A.gndmean.U+D.gndmean.U));

clear gndmean


figure(1) % this will show data coverage
set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultaxesfontweight','bold')
set(0,'defaultaxeslinewidth',1)


set(gcf,'units','normalized','position',[0 0 1 1],'color','w','PaperPosition',[0 0 8.25 10],'renderer','painters')
cmap = colormap(jet(100));

windst = double(sqrt(wsU.^2+wsV.^2));

subplot('position',[.1 .1 .55 .7])

%{
        if zooooom
            m_proj('lambert','long',[228 242],'lat',[30 45]);
        else
            m_proj('lambert','long',[210 250],'lat',[25 50]);
        end
%}

m_proj('lambert','long',[210 250],'lat',[25 50]);
m_plot(lowDClon,lowDClat,'.','color',[.9 .9 .9]); 
% [.68 .75 .68] plot regions of low data coverage
m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.85.*[1 1 1],'edgecolor','k')

%{
        if zooooom
            m_grid('box','on','tickdir','out','linestyle','none', ...
               'xtick',[228 242],'ytick',[30 45]);
        else
            %m_grid('box','on','tickdir','out','linestyle','none', ...
               %'xtick',[220 250],'ytick',[30 60]);
            m_grid('box','on','tickdir','out','linestyle','none', ...
               'xtick',[220 250],'ytick',[25 50]);
        end
%}

m_grid('box','on','tickdir','out','linestyle','none', ...
    'xtick',[220 250],'ytick',[25 50]);

hold on

X = reshape(lon,ny,nx);
Y = reshape(lat,ny,nx);
Z = reshape(windst,ny,nx);

pp = m_pcolor(X,Y,Z);
set(pp,'linestyle','none')
caxis([0 .2])

% set all wind vectors to unit speed for wdir display
% if zooooom;step=15;else step=20;end

step = 20;
sub = 1:step:length(lon);
sub = sub(1:step:length(sub));

junk = atan2(wsV,wsU);
dU = cos(junk);
dV = sin(junk);

m_quiver(lon(sub),lat(sub),dU(sub),dV(sub),1,'color','k')

%{
for ii = 1:length(ndbc)
    m_plot(360+ndbc(ii).ll(1),ndbc(ii).ll(2),'ko','markerfacecolor','w','markersize',3)
    
end
%}

set(gca,'layer','top')

xtk = 0:.02:.2;
tklbl = {'0';'';'0.04';'';'0.08';'';'0.12';'';'0.16';'';'0.2'};
colorbar_h(xtk,tklbl, [.2 .87 .35 .02],'wind stress (Pa)',7)

saveas(gcf,['/Volumes/Lacie_Primary/kayla/oaflux_research/figures/quikscat/',fname],'png')


%     saveas(gcf,['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/png/',fname],'png')
%     saveas(gcf,['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/fig/',fname],'fig')
%     print('-depsc2',['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/eps/',fname])

close all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot grand mean average asc and des Wind Stress Curl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1) % this will show data coverage
set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultaxesfontweight','bold')
set(0,'defaultaxeslinewidth',1)


set(gcf,'units','normalized','position',[0 0 1 1],'color','w','PaperPosition',[0 0 8.25 10],'renderer','painters')
%cmap = colormap(jet(100));
cmap = c2h(100); % cold (blue) to red (hot) white in the middle
%cmap = brighten(cmap,-.5);
cmap(1,:) = [0 0 0]; % force offscale to black and magenta
cmap(100,:) = [1 0 1];
colormap(cmap)


windst = double(sqrt(wsU.^2+wsV.^2));

subplot('position',[.1 .1 .55 .7])

%{
        if zooooom
            m_proj('lambert','long',[228 242],'lat',[30 45]);
        else
            m_proj('lambert','long',[210 250],'lat',[25 50]);
        end
%}

m_proj('lambert','long',[210 250],'lat',[25 50]);
m_plot(lowDClon,lowDClat,'.','color',[.9 .9 .9]); 
% [.68 .75 .68] plot regions of low data coverage
m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.85.*[1 1 1],'edgecolor','k')

%{
        if zooooom
            m_grid('box','on','tickdir','out','linestyle','none', ...
               'xtick',[228 242],'ytick',[30 45]);
        else
            %m_grid('box','on','tickdir','out','linestyle','none', ...
               %'xtick',[220 250],'ytick',[30 60]);
            m_grid('box','on','tickdir','out','linestyle','none', ...
               'xtick',[220 250],'ytick',[25 50]);
        end
%}

m_grid('box','on','tickdir','out','linestyle','none', ...
    'xtick',[220 250],'ytick',[25 50]);

hold on

X = reshape(lon,ny,nx);
Y = reshape(lat,ny,nx);
Z = reshape(wsC,ny,nx);

pp = m_pcolor(X,Y,Z);
set(pp,'linestyle','none')
caxis([-.0000005 .0000005])

% set all wind vectors to unit speed for wdir display
% if zooooom;step=15;else step=20;end

step = 20;
sub = 1:step:length(lon);
sub = sub(1:step:length(sub));

junk = atan2(wsV,wsU);
dU = cos(junk);
dV = sin(junk);

m_quiver(lon(sub),lat(sub),dU(sub),dV(sub),1,'color','k')

%{
for ii = 1:length(ndbc)
    m_plot(360+ndbc(ii).ll(1),ndbc(ii).ll(2),'ko','markerfacecolor','w','markersize',3)
    
end
%}

set(gca,'layer','top')

%{
xtk = 0:.02:.2;
tklbl = {'0';'';'0.04';'';'0.08';'';'0.12';'';'0.16';'';'0.2'};
colorbar_h(xtk,tklbl, [.2 .87 .35 .02],'wind stress (Pa)',7)
%}

B = colorbar;
set(B, 'Position', [.9 .3 .02 .5]) % left,bottom,width,height
xlabel(B, '10^-7  N/m^(-3)') %'m/s') %'^oC') %'W / m^2') 'g / Kg') 'Pa'


saveas(gcf,['/Volumes/Lacie_Primary/kayla/oaflux_research/figures/quikscat/',fname2],'png')

