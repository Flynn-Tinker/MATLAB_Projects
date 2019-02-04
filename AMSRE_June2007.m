%
% Last edited 09.18.15
% Main computer
%
% Output AMSR-E SST from 25-50 deg N, 110-140 deg W.

clear all
close all
clc

% How many years is your data set?
% EDIT: years, Filename
%
years = (2007);
Lengthyears = length(years);

% LIS MAP
fs = 12;
% PP = ([0 0 18 10]);

%{
m_proj('mercator','latitude',[25 50], ...
    'longitude',[-140 -110])
m_gshhs_f('save','/Users/kayla425/Documents/MATLAB/WestCoastMap.mat');
%}


load WestCoastMap.mat



% Creates Year cell array to easily find event files
for yy = 1%:Lengthyears
    
    Year{yy} = num2str(years(yy));
    FileName{yy} = dir('/Users/kayla425/Desktop/West_Coast_AMSRE');
    % Using dir leaves top 3 as (. .. and DS_Store) , so we need to delete those
    FileName{1,yy}(1:3,:) = [];
    
    % Gives all of the available files for that year
    
    % Will try one single file initially:
    
    numbfiles(yy) = length(FileName{1,yy});
    numbplots(yy) = floor(numbfiles(yy)/20); %removes some files
    
    
    m = 5; %row
    n = 4; %column
    aa = 1;
    %}
    
    %stop_try = numbplots(yy) - 8;
    % I need to read all of the files, but only plot 20 at a time
    for ff = 1:numbplots(yy)
        clf
        for gg = 1:20;
            clear sst lat lon
            Subsetted{gg} = ['/Users/kayla425/Desktop/West_Coast_AMSRE/',FileName{1,yy}(gg,1).name];
            % for sinlge plot Users/kayla425/Documents/MATLAB/West_Coast_AMSRE/
            [sst,lon,lat,time,dt,bias,sigma,rjct,conf,prox] = ...
                readL2Pcore(Subsetted{gg});
            %aa = aa+1;
            
            % Want SST in Degree C
            sst = double(sst)-273.15;
            sst(sst<-1.5) = nan; % Delete cloud top
            lat = double(lat);
            lon = double(lon);
            Datenumb = datenum(1981,1,1,0,0,double(time));
            
            
            t = subplot(m,n,gg);
            
            % THIS is for if the packfig IS NOT used
            if (gg == 1 || gg == 6 || gg == 11 || gg == 16)
                t = subplot(m,n,gg);
                set(gca,'XAxisLocation','top')
            else
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
            end
            
            if gg == 1
                ylabel('Latitude','Fontsize',16)
                xlabel('Longitude','Fontsize',16)
            end
            
            
            %{
            % Land %
            Land = prox;
            Land (Land > .9) = [];
            Land_percent = (length(Land)/(length(prox)))*100;
            
            % Non-Land %
            Non_Land = prox;
            Non_Land (Non_Land < .9) = [];
            Non_Land_percent = (length(Non_Land)/(length(prox)))*100;
            
            % Prox_5 %
            Prox_5 = prox;
            Prox_5(Prox_5 <4.9) = [];
            Prox_5_percent(gg) = (length(Prox_5)/(length(prox)))*100;
            %}
            
            hold on
            % Curently set to the highest = 5
            %sst(double(prox)<5) = NaN; % PROX MAX IS 4??
            
            % I dont want to delete images yet
            %{
        [a,e] = size(lon);
        if (Prox_5_percent(ff) < 40 || isnan(Prox_5_percent(ff)))...
            sst = nan(a,e);
        end
            %}
            
            
            % HOLD ON FOR THIS
            
            sst(~(lat >=25 & lat <= 50 & lon >= -140 & lon <= -110)) = nan;
            prox(~(lat >=25 & lat <= 50 & lon >= -140 & lon <= -110)) = nan;
            %if ~isempty(sst_x) 
            % removed gg
            min_temp(gg) = nanmin(sst(:));
            max_temp(gg) = nanmax(sst(:));
            lat_x = lat;
            lon_x = lon;
            lat_x(~(lat >=25 & lat <= 50 & lon >= -140 & lon <= -110)) = nan;
            lon_x(~(lat >=25 & lat <= 50 & lon >= -140 & lon <= -110)) = nan;
            pcolor(lon_x,lat_x,sst)
            % warning off
            shading flat
            
            
            %{
            %else 
            pcolor(lon,lat,sst)
            % warning off
            shading flat
            %end
            %}
          
            % For Strings
        %{
    Event_dn = cat(1,DNumb{1,yy});
    % Remove Nans which were added where times didn't
    % lay within +/-5 days of the Event
    nanval_dn = find(isnan(Event_dn)); % Find where nans are located
    FileName{1,yy}(nanval_dn) = []; % remove non Event times from FileName
    Event_dn(isnan(Event_dn)) = []; % remove nans from dn
    % Create a datestring for plotting
    lengDN = length(Event_dn);
    Event_dstr = datestr(Event_dn,'mmmddyyyy HH:MM:SS');
    Event_dstr = cellstr(Event_dstr);
        %}
            Subset_dstr = datestr(Datenumb,'mmmddyyyy HH:MM:SS');
            Subset_dstr = cellstr(Subset_dstr);
            text(0.25,0.98,Subset_dstr,'Fontsize',5,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            
            set(gca,'LineWidth',1,'FontSize',10,'Tick','out')
            
            ind = find(isnan(ncst(:,1)));
            for k = 1:length(ind)-1
                fill(ncst(ind(k)+1:ind(k+1)-1,1),ncst(ind(k)+1:ind(k+1)-1,2),.7*[1 1 1],'LineStyle','none')
            end
            
            % Edited axis
            axis([-140 -110 25 50]);
            
            % Land line
            plot(ncst(:,1),ncst(:,2),'k','LineWidth',0.5);
            % fix map to correct for deformation
            % home/set path/ add subfolder
            c_basescale;
            
            
            % For colorbar
            % Confused?
            %
            %{
            if ~isempty(sst) 
            min_temp(gg) = nanmin(sst(:));
            max_temp(gg) = nanmax(sst(:));
            end
            %}
            
            
            %}
            %colorbar
            
        end
        
        min_temp = nanmean(min_temp(:));
        max_temp = nanmean(max_temp(:));
        if ~isnan(min_temp) && ~isnan(max_temp) && ...
                    (min_temp < max_temp)
                caxis([min_temp max_temp]);
        end
        
        %end
        %FileDatenum{yy} = Datenumb;
        %
        
        
        % For percent coverage
        %{
    if max(lengDN) < .1
        continue
    end
        %}
        
        %caxis([min(mean_sst_rect_out) max(mean_sst_rect_in)]);
        B = colorbar;
        set(B, 'Position', [.92 .11 .04 .8])
        xlabel(B, '^oC')
        % x-placement; y placement; width;length
        set(gcf,'PaperPosition',[0 0 11 8.5])
        % orient landscape
        % packboth(m,n);
        
        %{
        h = packfig(m,n);
        for ii = 1:length(h);
            axes(h(ii));
            if ii == 1
                set(gca,'XAxisLocation','top')
                set(gca,'xticklabel', [-74 -73 -72 -71],'Fontsize',7)
                xlabel('Longitude','Fontsize',16)
            elseif ii > 1
                set(gca,'xtick',[]);
                set(gca,'ytick',[]);
            end
        end
        %}
        
        
        
        % Delete empty figures
        %{
    subplotMax = m*n;
    subplotsTODelte = subplotMax - max(lengDN);
    if max(lengDN) < n*4;
        set(gcf,'xtick',[])
        xtick([]);
    end
    ff = [];
    delete(subplot(m,n,(max(lengDN)+1):28));
    delete(subplot(m,n,29:35));
    set(gca,'XTick',-122:-119);
        %}
        
        suptitle(Subset_dstr);
        plot_num = num2str(ff);
        
        print('-dpng', ['/Users/kayla425/Desktop/West_Coast_AMSRE_figures/plot2007',plot_num]) ;
        % '-depsc2'
        %clear gg Subset_dstr
    end
end












