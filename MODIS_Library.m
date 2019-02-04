% Last edited 09.17.15
%
% map showing location of Santa Barbara map
% Maps one event +/- 5 days
%
clear all
close all
clc



% Choose year of data to load and then start and end of event date

% right now I have to do a year of data at a time because of my lable
% variable


%load Events_dn_09to13_0529
load Events_dn_0629


for zz = 13%:length(Event)
    Event(zz)
    load('DeltaT_06to13_0625.mat','numbfiles','FileName','Year','DNumb','mean_sst_rect_out');
    
    % Why 09to13?
    % DeltaT_09to13_0529.mat
    
    % MAB MAP
    fs = 12;
    % PP = ([0 0 18 10]);
    
    %m_proj('mercator','latitude',[34 35.4], ...
    %'longitude',[-121.2 -119.6])
    %m_gshhs_f('save','/Users/kayla/Documents/MATLAB/CaliMap.mat');
    
    load CaliMap.mat
    
    clear sst lat lon
    
    figure(1)
    clf % clears figure 1
    m = 5; % rows
    n = 7; %colm
    
    % Select DN that match event time
    % How does it know which SST to use?
    % I can use the new DN for DeltaT, but not sure about SST?
    
    
    
    for yy = 1:length(Year)
        
        %{
    FileName = dir(['/Volumes/LaCie/MODIS_A/',Year{yy}]);
    
    % Using dir leaves top 3 as (. .. and DS_Store) , so we need to delete those
    FileName(1:3,:) = [];
    % Gives all of the available files for that year
        %}
        
        for cc = 1:numbfiles(yy)
            if (DNumb{1,yy}(1,cc) >= starting_dn(zz)) && (DNumb{1,yy}(1,cc) <= ending_dn(zz))
                %starting_dn <= DNumb{1,1}(1,ff) <= ending_dn
                % DNumb{1,1}(1,ff) = DNumb{1,1}(1,ff);
                continue
            else
                DNumb{1,yy}(1,cc) = NaN;
            end
        end
        
        % Create a new variable name
        New_dn{1,yy} = DNumb{1,yy};
        
        nanval{1,yy} = find(isnan(New_dn{1,yy}));
        FileName{1,yy}(nanval{1,yy},:) = [];
        New_dn{yy}(isnan(New_dn{1,yy})) = [];
        
        lengDN(yy) = length(New_dn{yy});
        New_dstr{yy} = datestr(New_dn{1,yy},'mmmddyyyy HH:MM:SS');
        New_dstr{yy} = cellstr(New_dstr{1,yy});
        
        
        
        %numbfiles = length(FileName);
        
        for ff = 1:lengDN(yy);
            Subsetted{ff} = ['/Volumes/LaCie/kayla/MODIS_A/',Year{yy},'/',FileName{1,yy}(ff,1).name];
            [sst,lon,lat,time,dt,bias,sigma,rjct,conf,prox] = ...
                readL2Pcore(Subsetted{ff});
            
            
            % Want SST in Degree C
            sst = double(sst)-273.15;
            lat = double(lat);
            lon = double(lon);
            
            
            
            t = subplot(m,n,ff);
            
            if (ff == 1 | ff == 8 | ff == 15 | ff == 22)
                t = subplot(m,n,ff);
            else
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
            end
            
            if ff == 15
                ylabel('Latitude','Fontsize',16)
            end
            
            if ff == 18
                xlabel('Longitude','Fontsize', 16)
            end
            
            hold on
            
            
            
            
            % Curently set to the highest = 5
            sst(double(prox)<5) = NaN;
            pcolor(lon,lat,sst)
            % warning off
            shading flat
            
            % Five{1,tt} = num2str(Five{1,tt});
            text(0.25,0.98,New_dstr{1,yy}(ff,1),'Fontsize',5,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            % text(0.8,0.87,Numb{tt},'Fontsize',5.5,'Units', 'Normalized', 'VerticalAlignment', 'Top')
            % text(0.8,0.8,Plot{tt},'Units', 'Normalized', 'VerticalAlignment', 'Top')
            
            
            
            
            
            set(gca,'LineWidth',1,'FontSize',10,'Tick','out')
            
            ind = find(isnan(ncst(:,1)));
            for k = 1:length(ind)-1
                fill(ncst(ind(k)+1:ind(k+1)-1,1),ncst(ind(k)+1:ind(k+1)-1,2),.7*[1 1 1],'LineStyle','none')
            end
            
            % text(.5,.5,'December','horizontalAlignment','center');
            
            
            % Edited axis
            axis([-121.2 -119.6 34 35.4]);
            
            % color scale range
            % Make seasonal variation definition
            %{
            LabNew_dstr = New_dstr{1,yy};
            
            if (LabNew_dstr{zz,1}(1:3) == 'Nov')
                caxis([10 17]);
            elseif (LabNew_dstr{zz,1}(1:3) == 'Dec')
                caxis([10 17]);
            elseif (LabNew_dstr{zz,1}(1:3) == 'Jan')
                caxis([10 17]);
            elseif (LabNew_dstr{zz,1}(1:3) == 'Feb')
                caxis([10 17]);
            elseif (LabNew_dstr{zz,1}(1:3) == 'Mar')
                caxis([10 17]);
            elseif (LabNew_dstr{zz,1}(1:3) == 'Apr')
                caxis([10 17]);
            elseif (LabNew_dstr{zz,1}(1:3) == 'May')
                caxis([10 17]);
            else
                caxis([13 20]);
            end
            %}
            
            plot(ncst(:,1),ncst(:,2),'k','LineWidth',0.5);
            % fix map to correct for deformation
            % home/set path/ add subfolder
            c_basescale;
            % Draw rectangles on map
            Rect_out = rectangle('Position',[-120.75,34.55,.1,.35]);
            Rect_in = rectangle('Position',[-120.3,34.37,.8,.08]);
            
            % dbstop if error ???
            
        end
        
        %clear LabNew_dstr
    end
    
    %{
    if isempty(ff)==1
            continue
    end
    %}
    
        %caxis([min(mean_sst_rect_out) max(mean_sst_rect_in)]);
    B = colorbar;
    set(B, 'Position', [.92 .11 .04 .8])
    xlabel(B, '^oC')
    % x-placement; y placement; width;length
    orient landscape
    % packboth(m,n);
    
    
    h = packfig(m,n);
    for ii = 1:length(h);
       axes(h(ii));
       if ii == 1
           set(gca,'XAxisLocation','top')
           set(gca,'xticklabel', [-121 -120.5 -120 -119.5],'Fontsize',7)
           xlabel('Longitude','Fontsize',16)
       elseif ii > 1
           set(gca,'xtick',[]);
           set(gca,'ytick',[]);
       end
    end
    
    suptitle(Event{zz});
    
    % print('-depsc2', ['/Volumes/LaCie/kayla/Figures/Library_0903/',Event{zz}]) ;
    
end

