%% Monthly Climatology for Argo Floats
%
% Last edited 3.15.16
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a map of the California Current System with all of the mean
% values (MLD and Temperautre) using density and temperature algorithms.
% We look at the std for the MLD and Temperature as well to estimate
% What the composite mfile will look like.
%
% INPUT:
% monthlyclim.mat from http://mixedlayer.ucsd.edu
%
% Current Set up:
% Load mat file
% Keep JJA summer months in the corresponding lat/lon
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

load('monthlyclim.mat')

% remove values that we will not be using here:
%
% density algorithm maximum MLD
% density algorithm median MLD
% density threshold maximum MLD
% density threshold median MLD
% temperature algorithm median MLD
% temperature threshold median MLD

% mean mixed layer salinity and potential density calculated with the density algorithm mean MLD
% mean mixed layer salinity and potential density calculated with the density threhsold mean MLD

clear mld_da_max mld_da_median mld_dt_max mld_dt_median mld_ta_median mld_tt_median
clear mlpd_da mlpd_dt mls_dt mls_da

% lat and lon
lat = latm(40:71, 115:141); % 25 to 50 deg
%lat = lat';
%lat = flipud(lat);
lon = lonm(40:71, 115:141); % -140 to -110 deg
lon = lon+360;
%lon = lon';

clear lonm latm

% density algorithm mean MLD
mld_da_mean = mld_da_mean(6:8,40:71,115:141);

% density algorithm MLD standard deviation
mld_da_std = mld_da_std(6:8,40:71,115:141);

% density threshold mean MLD
mld_dt_mean = mld_dt_mean(6:8,40:71,115:141);

% density threshold MLD standard deviation
mld_dt_std = mld_dt_std(6:8,40:71,115:141);

% temperature algorithm mean MLD
mld_ta_mean = mld_ta_mean(6:8,40:71,115:141);

% temperature algorithm sdt MLD
mld_ta_std = mld_ta_std(6:8,40:71,115:141);

% temperature threshold mean MLD
mld_tt_mean = mld_tt_mean(6:8,40:71,115:141);

% temperature threshold sdt MLD
mld_tt_std = mld_tt_std(6:8,40:71,115:141);

% mean mixed layer temperature calculated with the density algorithm mean MLD
temp_da = mlt_da(6:8,40:71,115:141);
clear mlt_da

% mean mixed layer temperature calculated with the density threhsold mean MLD
temp_dt = mlt_dt(6:8,40:71,115:141);
clear mlt_dt

% num??
num = num(6:8,40:71,115:141);

% lat 37:50 lon :
da_mean_N = nanmean(mld_da_mean(:,:,14:27));
da_mean_N = reshape(da_mean_N,32,14);
da_mean_N = nanmean(da_mean_N(:));

% lat 25:37 lon:
da_mean_S = nanmean(mld_da_mean(:,:,1:13));
da_mean_S = reshape(da_mean_S,32,13);
da_mean_S = nanmean(da_mean_S(:));

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot each set of data on a map.
% 3 subplots (June, July, and August)
%
% Input a varibale from above into var
% Edit Colorbar label
%

for mm = 1:3 % 3 months
    
    subplot(1,3,mm)
    
    %cmap = c2h(100);
    %colormap(cmap)
    
    m_proj('lambert','long',[220 250],'lat',[25 50]);
    hold on
    
    m_grid('box','on','tickdir','out','linestyle','none', ...
        'xtick',[220 250],'ytick',[25 50], ...
        'xticklabel','','yticklabel','');
    
    var_name = 'mld_da_mean'; % mld temp dt
    var = mld_da_mean(mm,:,:);
    var = reshape(var,32,27);
    %var = var';
    %var = flipud(var);
    
    
    pp = m_pcolor(lon,lat,var);
    set(pp,'linestyle','none')
    caxis([5 25])
    %caxis([0 20])
    %caxis([12 26])
    %colorbar
    
    if mm == 1
        title('June')
    elseif mm == 2
        title('July')
    else
        title('August')
    end
    
    
    m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
    [pbx,pby] = political_boundaries;
    m_plot(pbx,pby,'-','color','k')
    
    clear var
end % of ww window loop


outdir = '/Volumes/Lacie_Primary/kayla/oaflux_research/figures/argo/';

packfig(1,3)

suptitle(var_name);
B = colorbar;
set(B, 'Position', [.91 .3 .03 .4]) % left,bottom,width,height
xlabel(B, 'm') %'^oC' m

saveas(gcf,[outdir,var_name,'_argo_climo_re_JJA.png'],'png')

%%
%%%%%%%%%%%%%%%%%%%
%
% Plot summer mean values to find appropriate h
%
%%%%%%%%%%%%%%%%%%%


m_proj('lambert','long',[220 250],'lat',[25 50]);
hold on

m_grid('box','on','tickdir','out','linestyle','none', ...
    'xtick',[220 250],'ytick',[25 50], ...
    'xticklabel','','yticklabel','');

var_name = 'mld_da_mean'; % mld temp dt
var = nanmean(mld_da_mean(:,:,:),1);
var = reshape(var,32,27);
%var = var';
%var = flipud(var);


pp = m_pcolor(lon,lat,var);
set(pp,'linestyle','none')
caxis([10 45])
%caxis([0 20])
%caxis([12 26])
%colorbar

title('Summer')


m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
[pbx,pby] = political_boundaries;
m_plot(pbx,pby,'-','color','k')




outdir = '/Volumes/Lacie_Primary/kayla/oaflux_research/figures/argo/';


suptitle(var_name);
B = colorbar;
set(B, 'Position', [.91 .3 .03 .4]) % left,bottom,width,height
xlabel(B, 'm') %'^oC' m

saveas(gcf,[outdir,var_name,'_argo_climo_re_summer.png'],'png')

%%

clear all
close all
clc

matdir = '/Volumes/Lacie_Primary/kayla/oaflux_research/argo/';
load([matdir,'argo_data_JJA']);

load('/Users/kayla/Documents/MATLAB/matfiles/others_matfiles/chris/auto/Auto_etimes_1982_now_10m.mat'); % began using this file when we changed to 10m winds and decided to use only auto etimes
arrive = etimes;
clear etimes

%{
%%% Fit events using 500mb heights from day -1 to +3
load('positive_anomaly_flag_500mbht_1982_2010.mat')
keep = pos_flag(:,2)==1;
arrive = pos_flag(keep,1);
clear keep pos_flag
%}

arrive = arrive(arrive>datenum(2002,6,6) & arrive<datenum(2009,8,25)); % limit this range to available SST
junk = datevec(arrive);
keep = junk(:,2)>5 & junk(:,2)<9; % limit to June-Aug here
arrive = arrive(keep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find Argo composite climatology 2002-2009 JJA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(220:250,25:50); % output grid

count = 1;

% Loop through each arrival
for aa = 1:size(arrive,1)
    
    e_idx = find(profiledate>=(arrive(aa,1)-6) & profiledate<=(arrive(aa,1)+5));
    %
    % e indx
    %
    for ff = 1:length(e_idx)
        % profile latitude
        e_profilelat(:,count) = profilelat(e_idx(ff));
        
        % profile longitude
        e_profilelon(:,count) = profilelon(e_idx(ff));
        
        % profile date
        e_profiledate(:,count) = profiledate(e_idx(ff));
        
        % float number (WMO number)
        e_floatnumber(:,count) = floatnumber(e_idx(ff));
        
        % profile number
        e_profilenumber(:,count) = profilenumber(e_idx(ff));
        
        % status (delayed mode corresponds to 1)
        e_status(:,count) = status(e_idx(ff));
        
        % density algorithm MLD
        e_da_mld(:,count) = da_mld(e_idx(ff));
        
        % density threshold MLD
        e_dt_mld(:,count) = dt_mld(e_idx(ff));
        
        % temperature algorithm MLD
        e_ta_mld(:,count) = ta_mld(e_idx(ff));
        
        % temperature threshold MLD
        e_tt_mld(:,count) = tt_mld(e_idx(ff));
        
        % mixed layer temperature calculated with the density algorithm MLD
        e_temp_da(:,count) = temp_da(e_idx(ff));
        
        % mixed layer temperature calculated with the density threhsold MLD
        e_temp_dt(:,count) = temp_dt(e_idx(ff));
        
        count = count+1;
    end
    clear e_idx
    
end

%%
m_proj('lambert','long',[220 250],'lat',[25 50]);
hold on

m_grid('box','on','tickdir','out','linestyle','none', ...
    'xtick',[220 250],'ytick',[25 50], ...
    'xticklabel','','yticklabel','');

%pp = m_pcolor(profilelon,profilelat);
%set(pp,'linestyle','none')
m_plot(e_profilelon,e_profilelat,'.')

m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
[pbx,pby] = political_boundaries;
m_plot(pbx,pby,'-','color','k')

%outdir = '/Volumes/Lacie_Primary/kayla/oaflux_research/figures/argo/';
%saveas(gcf,[outdir,'argofloat_map.png'],'png')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Will not be using some terms here, but will determine JJA climo
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Removed:
%
%   e_floatnumber e_profilenumber e_status
%   e_dt_mld e_tt_mld e_ta_mld e_temp_dt
%
%   Separate into JJA data
%
%

e_string = datestr(e_profiledate,'mmm.dd,yyyy');

for mm = 1:3 % 3 months
    
    subplot(1,3,mm)
    
    %cmap = c2h(100);
    %colormap(cmap)
    
    m_proj('lambert','long',[220 250],'lat',[25 50]);
    hold on
    
    m_grid('box','on','tickdir','out','linestyle','none', ...
        'xtick',[220 250],'ytick',[25 50], ...
        'xticklabel','','yticklabel','');
    
    var_name = 'composite da temp'; % mld
    
    if mm == 1
        title('June')
        jun_idx = strmatch('Jun',e_string);
        F = scatteredInterpolant(e_profilelon(jun_idx)',e_profilelat(jun_idx)',...
            e_da_mld(jun_idx)');
        jun_mld = F(X,Y);
        clear F
        
        F = scatteredInterpolant(e_profilelon(jun_idx)',e_profilelat(jun_idx)',...
            e_temp_da(jun_idx)');
        jun_temp = F(X,Y);
        clear F
        pp = m_pcolor(X,Y,jun_temp);
    elseif mm == 2
        title('July')
        jul_idx = strmatch('Jul',e_string);
        F = scatteredInterpolant(e_profilelon(jul_idx)',e_profilelat(jul_idx)',...
            e_da_mld(jul_idx)');
        jul_mld = F(X,Y);
        clear F
        
        F = scatteredInterpolant(e_profilelon(jul_idx)',e_profilelat(jul_idx)',...
            e_temp_da(jul_idx)');
        jul_temp = F(X,Y);
        clear F
        pp = m_pcolor(X,Y,jul_temp);
    else
        title('August')
        aug_idx = strmatch('Aug',e_string);
        F = scatteredInterpolant(e_profilelon(aug_idx)',e_profilelat(aug_idx)',...
            e_da_mld(aug_idx)');
        aug_mld = F(X,Y);
        clear F
        
        F = scatteredInterpolant(e_profilelon(aug_idx)',e_profilelat(aug_idx)',...
            e_temp_da(aug_idx)');
        aug_temp = F(X,Y);
        clear F
        pp = m_pcolor(X,Y,aug_temp);
    end
    
    set(pp,'linestyle','none')
    %caxis([0 50])
    %caxis([0 20])
    caxis([12 26])
    %colorbar
    
    m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
    [pbx,pby] = political_boundaries;
    m_plot(pbx,pby,'-','color','k')
    
    clear var
end % of ww window loop


outdir = '/Volumes/Lacie_Primary/kayla/oaflux_research/figures/argo/';

packfig(1,3)

suptitle(var_name);
B = colorbar;
set(B, 'Position', [.91 .3 .03 .4]) % left,bottom,width,height
xlabel(B, '^oC') %'^oC' m

saveas(gcf,[outdir,var_name,'_argo_climo_JJA.png'],'png')



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Yearly climatology for each lat/lon. Mean of each lat/lon for JJA
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%








