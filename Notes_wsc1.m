%% Going through Chris Step 1 To better understand
% This is for no rotate and WSC
%
% Last edited 03.22.16
%
% Saving everything in May folder
%
% Editing for JJA 02-09
%
clear all
close all
fclose all;
clc

%{
addpath /data/data02/transfer/Chris/mfile_library/
addpath /data/data02/transfer/Chris/mfile_library/m_map/
indir = '/data/data02/transfer/Chris/raw_data/QuikSCAT/v3/';
gnd_dir = '/data/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/mat/';
%}
indir = '/Volumes/Lacie_Secondary/data/quikscat/';
%%%gnd_dir = '/Volumes/Lacie_Primary/kayla/oaflux_research/quikscat/May/';
gnd_dir = '/Volumes/Lacie_Primary/kayla/oaflux_research/quikscat/wsc/';
wsmax = 0.2; % max wind stress (N/m2) for color scaling
%%% wind stress curl max? +/- 5 *10^-7 N/m^3
Data_coverage_switch = 1;  % enable to remove locations with less than min_n finite points
min_n = 10;  % minimum number of points required to average (typically 10)
a_or_d = 'ascending';

%{
zooooom = 0; % 1 = make axes a bit smaler
if zooooom
    %fname = [fname,'_zoom'];
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

reprocess_switch = 1;
if reprocess_switch
    
    %year = num2cell(1999:2009);
    year = num2cell(2002:2009);
    for yy = 1:length(year);
        qlist{yy} = fuf([indir,num2str(year{yy}),'/','qs_l2b_*_v3_*.nc'],0,'normal'); % this takes a while...
        %list = [list, lis(:)];
    end
    
    list = [qlist{1} ; qlist{2} ; qlist{3}; qlist{4}; qlist{5}; qlist{6}; qlist{7}; qlist{8}];% qlist{9}; qlist{10}; qlist{11}];
    
    %list = fuf([indir,'qs_l2b_*_v3_*.nc'],0,'normal'); % this takes a while...
    
    % determine how many layers to allocate for sU and sV
    junk = char(list);
    junk = str2num(junk(:,21:22));
    %junk = find(junk>4 & junk<9);
    junk = find(junk>5 & junk<9);
    nfiles = round(length(junk)*0.2);
    % ~80% of files are skipped due to lon/lat or time ranges
    clear junk
    
    % add preallocations and gridding parameters
    [X,Y] = meshgrid(220:.1:250,25:.1:50);
    sU = zeros((size(X,1)*size(X,2)),nfiles,'single');
    sV = zeros((size(X,1)*size(X,2)),nfiles,'single');
    stime = zeros(1,nfiles,'single');
    nrawpts = zeros(1,nfiles,'single');
    clear nfiles year
    
    count = 0; % will have to impose a counter rather than using ii to accommodate skipped files (memory issue)
    
    for ii = 1:length(list)
        
        junk = str2num(char(list{ii}(21:22))); % a quick date check here
        junk2 = str2num(char(list{ii}(17:20))); % year check
        if junk<6 || junk>8 && junk2>2001 % junk<5 || junk>8
            clear junk
            continue
        end
        
        year = num2str(list{ii}(17:20));
        [~, W] = read_nc_file_struct([indir,year,'/',list{ii}]); % read a daily file
        
        %[~, W] = read_nc_file_struct([indir,list{ii}]); % read a daily file
        
        % inspect and see if we have data within our date and location ranges
        lonidx = W.lon(:)>=220 & W.lon(:)<=250; % locate data within our longitude range
        latidx = W.lat(:)>=25 & W.lat(:)<=50;  % and data within our lat range
        flagidx = W.flags(:)<16384; % guarantee 4 beams
        
        mtime = datenum(1999,1,1,0,0,double(W.time));
        mtime = repmat(mtime.',size(W.lat,1),1); % same size it for use below
        mvec = datevec(mtime(:));
        
        switch a_or_d
            case 'ascending'
                %tidx = (mvec(:,2)>4 & mvec(:,2)<9) & mvec(:,4)<8;
                tidx = (mvec(:,2)>5 & mvec(:,2)<9) & mvec(:,4)<8;
            case 'descending'
                %tidx = (mvec(:,2)>4 & mvec(:,2)<9) & mvec(:,4)>8;
                tidx = (mvec(:,2)>5 & mvec(:,2)<9) & mvec(:,4)>8;
        end
        
        keep = (lonidx+latidx+tidx+flagidx)==4;  % all have to be true (1) to keep the point
        clear lon* lat* mvec tidx flagidx
        
        if sum(keep(:))<=2  % no data within our time/location range (turns out that triscatteredinterp errors out if only 1 point available.)
            clear keep W
            continue
        else
            count = count+1;
            disp(list{ii})
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Calculate wind stress for keep swaths
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        wspd = W.retrieved_wind_speed(:);
        wspd(~keep) = nan;
        zvec = 10.*ones(length(wspd),1);
        ts = 15.*ones(length(wspd),1);
        
        [~,~,~,~,tau] = mf_dragNC35(zvec,wspd,ts);
        
        [wU,wV] = SpeedDir2UV(tau.', W.retrieved_wind_direction(:));
        clear tau
        
        % now interp these data onto our 0.1 degree grid
        F = scatteredInterpolant(double(W.lon(:)),double(W.lat(:)),-wU);
        sU(:,count) = single(F(X(:),Y(:)));
        clear F
        F = scatteredInterpolant(double(W.lon(:)),double(W.lat(:)),-wV);
        sV(:,count) = single(F(X(:),Y(:)));
        clear F
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %
        % Calculate terms of wind stress curl; need lat gridded first
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %a = diff(X); % .1
        rad = pi/180;
        %dlat = mean(a);
        dlat = 0.1;
        clear a
        
        % Why 111,176
        deltay = dlat*111176; % all values are the same
        
        nx = 251;
        ny = 301;
        
        clear dd ee
        
        new_Y = Y-360;
        
        for dd=1:nx
            for ee=1:ny
                long(dd, ee)=new_Y(dd, ee)*111176*cos(X(dd, ee)*rad);
                %long(i,j)=lon(j)*6378137*rad*cos(lat(i)*rad);
                % [m] earth radious in meters= 6,378,137.0 m.. from wikipedia.
            end % endfor
        end % endfor
        
        clear dd ee
        
        Tx = sU(:,count);
        Tx = reshape(Tx,nx,ny);
        Ty = sV(:,count);
        Ty = reshape(Ty,nx,ny);
        
        
        % Centeral difference method in x and y
        for dd=2:nx-1
            for ee=2:ny-1
                
                curlZ(dd, ee)=(Ty(dd, ee+1)-Ty(dd, ee-1))./...
                    (2*(long(dd, ee+1)-long(dd, ee-1))) - ...
                    (Tx(dd+1, ee)-Tx(dd-1, ee))./(2*deltay) ;
                
            end % endfor
        end % endfor
        
        clear dd ee
        % Make wsc same size
        
        %re_wsc = reshape(curlZ,250,300);
        re_wsc = [NaN(250,1),curlZ];
        re_wsc = cat(1,NaN(1,301),re_wsc);
        
        curl(:,count) = re_wsc(:);
        
        % Make wsc same size
        %{
        [a,b] = size(sT);
        wsc = [];
        re_wsc = [];
        for ss = 1:b
            re_wsc = sT(:,ss);
            re_wsc = reshape(re_wsc,250,300);

            re_wsc = [NaN(250,1),re_wsc];
            re_wsc = cat(1,NaN(1,301),re_wsc);
            wsc = [wsc, re_wsc(:)];
            clear re_wsc
        end
        clear sT re_wsc
        sT = wsc;
        %}
        
        clear Tx Ty
        
        
        % took nanmean, to edit shape
        stime(count) = nanmean(mtime(keep));
        nrawpts(count) = sum(keep(:));
        % a talley of the number of raw points that went into a gridded field
        
        clear W wU wV keep mtime
        
        %     if count>=50 % uncomment for debugging memory issues
        %         break
        %     end
        
    end % of ii interp loop
    
    
    % % % % - one time fix of spike on 20060505
    % % % % save([gnd_dir,'Look4Wildpts_ascending.mat'],'sU','sV','stime','-v7.3')
    % % % load([gnd_dir,'Look4Wildpts_ascending.mat'])
    % % % count = 2944;
    % Not sure if this did what I wanted
    % or do I need to use keyboard?
    % DID NOT need, gndmeans matched Chris without this:
    %{
    junku = squeeze(sU(:,:,1767));
    junkv = squeeze(sV(:,:,1767));
    nbad = junku>0.3;
    junku(nbad) = nan;
    junkv(nbad) = nan;
    sU(:,:,1767) = junku;
    sV(:,:,1767) = junkv;
    clear junk* nbad
    % % % % keyboard  % can fix spike in sU and sV layer 1767 here
    %}
    %
    stime = stime(1:count); % get rid of any unused preallocated layers
    sU = sU(:,1:count); % matrix
    sV = sV(:,1:count);
    curl = curl(:,1:count);
    nrawpts = nrawpts(1:count);
    
    lon = X(:);
    lat = Y(:);
    % columns
    nx = size(X,2);
    % rows
    ny = size(X,1);
    nr = length(lon);
    clear X Y
    
    
    % due to memory constraints we're going to have to loop through each
    % grid point and perform the calculations sequentially
    
    wsU = nan(nr,1); % preallocate output (doubles fine here)
    wsV = nan(nr,1);
    wsc = nan(nr,1);
    
    %dU = nan(nr,1);
    %dV = nan(nr,1);
    stdU = nan(nr,1);
    stdV = nan(nr,1);
    stdS = nan(nr,1);
    %rot_stdU = nan(nr,1);
    %rot_stdV = nan(nr,1);
    
    for yy = 1:nr % rows
        
        % apply data coverage switch
        npts = sum(isfinite(sU(yy,:)),2);
        if npts < min_n;
            sU(yy,:) = nan;
            sV(yy,:) = nan;
            curl(yy,:) = nan;
            clear npts
            continue
        end
        
        % calculate mean wind stress for each row
        wsU(yy) = nanmean(sU(yy,:));
        wsV(yy) = nanmean(sV(yy,:));
        wsc(yy) = nanmean(curl(yy,:));
        
        % get standard deviation for each point
        stdU(yy) = nanstd(sU(yy,:));
        stdV(yy) = nanstd(sV(yy,:));
        stdS(yy) = nanstd(curl(yy,:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Remove rotation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        
        
        % we rotate U and V to line up the vector with
        % the mean and calculate a
        % rotated standard deviation:
        
        %{
        ra = atan2(wsV(yy),wsU(yy));
        junkU = (sU(yy,:).*cos(ra)) + (sV(yy,:).*sin(ra));
        junkV = (sV(yy,:).*cos(ra)) - (sU(yy,:).*sin(ra));
        rot_stdU(yy) = nanstd(junkU);
        rot_stdV(yy) = nanstd(junkV);
        dU(yy) = nanmean(junkU);
        dV(yy) = nanmean(junkV);
        clear junk*
        %}
        
    end % of yy loop
    
    
    % calculate df based on Melanie's 11/4 email
    df_vec = getdof(stime,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Save Values
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    gndmean.U = wsU(:);  % note: grand mean is now in units of wind stress
    gndmean.V = wsV(:);
    gndmean.S = wsc(:);
    %gndmean.dU = dU; %nanmean(dU,3); gndmean.dU = gndmean.dU(:);
    %gndmean.dV = dV; %nanmean(dV,3); gndmean.dV = gndmean.dV(:);
    %gndmean.rstdU = rot_stdU; % added these in to generate std plots (3/2013)
    %gndmean.rstdV = rot_stdV;
    % Can the same df_vec be used for curl?
    gndmean.Uconf95 = stdU.*tinv(.975,df_vec)./sqrt(df_vec);
    gndmean.Vconf95 = stdV.*tinv(.975,df_vec)./sqrt(df_vec);
    gndmean.Sconf95 = stdS.*tinv(.975,df_vec)./sqrt(df_vec);
    %gndmean.rot_Uconf95 = rot_stdU.*tinv(.975,df_vec)./sqrt(df_vec);
    %gndmean.rot_Vconf95 = rot_stdV.*tinv(.975,df_vec)./sqrt(df_vec);
    gndmean.lon = lon;
    gndmean.lat = lat;
    gndmean.nrawpts = nrawpts;
    gndmean.stime = stime;
    
    save([gnd_dir,'Qscat_JPL_',a_or_d,'_grandMean_wsc.mat'],'gndmean','nx','ny')
    
end % of reprocess_switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load previous calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([gnd_dir,'Qscat_JPL_',a_or_d,'_grandMean_wsc.mat'])
%
wsU = gndmean.U;  % note: grand mean is now in units of wind stress
wsV = gndmean.V;
wsS = gndmean.S;

lon = gndmean.lon;
lat = gndmean.lat;

% Is this ok??
lowDClon = lon(isnan(gndmean.U));
lowDClat = lat(isnan(gndmean.U));

clear gndmean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1)
set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultaxesfontweight','bold')
set(0,'defaultaxeslinewidth',1)
set(gcf,'units','normalized','PaperPosition',[0 0 8.25 10],'color','w','renderer','painters')
cmap = colormap(jet(100));

windst = sqrt(wsU.^2+wsV.^2);

subplot('position',[.1 .1 .55 .7])
%{
if zooooom
    m_proj('lambert','long',[228 242],'lat',[30 45]);
else
    m_proj('lambert','long',[220 250],'lat',[25 50]);
end
%}
m_proj('lambert','long',[220 250],'lat',[25 50]);

m_plot(lowDClon,lowDClat,'.','color',[.9 .9 .9])
% [.68 .75 .68] plot regions of low data coverage

%{
if zooooom
    m_grid('box','on','tickdir','out','linestyle','none', ...
        'xtick',[228 242],'ytick',[30 45]);
else
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
%{
if zooooom
    [X,Y] = meshgrid(linspace(229,241,8),linspace(31,44,10));
else
    [X,Y] = meshgrid(linspace(210,250,19),linspace(25,50,12));
end
%}

[X,Y] = meshgrid(linspace(210,250,19),linspace(25,50,12));

junk = atan2(wsV,wsU);
dU = cos(junk);
dV = sin(junk);

ZdU = griddata(lon,lat,dU,X,Y);
ZdV = griddata(lon,lat,dV,X,Y);
scl = 0.75;

m_quiver(X(:),Y(:),ZdU(:).*scl,ZdV(:).*scl,0,'color','k')

%{
for ii = 1:length(ndbc)
    m_plot(360+ndbc(ii).ll(1),ndbc(ii).ll(2),'ko','markerfacecolor','w','markersize',3)
    
end
%}

m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.85.*[1 1 1],'edgecolor','k')


xtk = 0:.02:.2;
tklbl = {'0';'';'0.04';'';'0.08';'';'0.12';'';'0.16';'';'0.2'};
colorbar_h(xtk,tklbl, [.2 .87 .35 .02],'wind stress (Pa)',7)

fname = ['norotate_qwkskt_JPL_windStress_2002-2009_',a_or_d,'_average_w_vectors'];
saveas(gcf,['/Volumes/Lacie_Primary/kayla/oaflux_research/figures/quikscat/',fname],'png')


% saveas(gcf,['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/png/',fname],'png')
% saveas(gcf,['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/fig/',fname],'fig')
% print('-depsc2',['/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/quickscat/all_year_average/ne_pacific/JPL/eps/',fname])



%%%%%%%%%%%%%%%%%%%%%%
%
% Plot wind stress curl with wind stress vectors
%
%%%%%%%%%%%%%%%%%%%%%%

close all

figure(1)
set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultaxesfontweight','bold')
set(0,'defaultaxeslinewidth',1)
set(gcf,'units','normalized','PaperPosition',[0 0 8.25 10],'color','w','renderer','painters')
%cmap = colormap(jet(100));
cmap = c2h(100); % cold (blue) to red (hot) white in the middle
%cmap = brighten(cmap,-.5);
cmap(1,:) = [0 0 0]; % force offscale to black and magenta
cmap(100,:) = [1 0 1];
colormap(cmap)


windst = sqrt(wsU.^2+wsV.^2);

subplot('position',[.1 .1 .55 .7])
%{
if zooooom
    m_proj('lambert','long',[228 242],'lat',[30 45]);
else
    m_proj('lambert','long',[220 250],'lat',[25 50]);
end
%}
m_proj('lambert','long',[220 250],'lat',[25 50]);

m_plot(lowDClon,lowDClat,'.','color',[.9 .9 .9])
% [.68 .75 .68] plot regions of low data coverage

%{
if zooooom
    m_grid('box','on','tickdir','out','linestyle','none', ...
        'xtick',[228 242],'ytick',[30 45]);
else
    m_grid('box','on','tickdir','out','linestyle','none', ...
        'xtick',[220 250],'ytick',[25 50]);
end
%}
m_grid('box','on','tickdir','out','linestyle','none', ...
    'xtick',[220 250],'ytick',[25 50]);

hold on

X = reshape(lon,ny,nx);
Y = reshape(lat,ny,nx);
Z = reshape(wsS,ny,nx);

pp = m_pcolor(X,Y,Z);
set(pp,'linestyle','none')
caxis([-.0000005 .0000005])

% set all wind vectors to unit speed for wdir display
%{
if zooooom
    [X,Y] = meshgrid(linspace(229,241,8),linspace(31,44,10));
else
    [X,Y] = meshgrid(linspace(210,250,19),linspace(25,50,12));
end
%}

[X,Y] = meshgrid(linspace(210,250,19),linspace(25,50,12));

junk = atan2(wsV,wsU);
dU = cos(junk);
dV = sin(junk);

ZdU = griddata(lon,lat,dU,X,Y);
ZdV = griddata(lon,lat,dV,X,Y);
scl = 0.75;

m_quiver(X(:),Y(:),ZdU(:).*scl,ZdV(:).*scl,0,'color','k')

%{
for ii = 1:length(ndbc)
    m_plot(360+ndbc(ii).ll(1),ndbc(ii).ll(2),'ko','markerfacecolor','w','markersize',3)
    
end
%}

m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.85.*[1 1 1],'edgecolor','k')

%{
xtk = 0:.02:.2;
tklbl = {'0';'';'0.04';'';'0.08';'';'0.12';'';'0.16';'';'0.2'};
colorbar_h(xtk,tklbl, [.2 .87 .35 .02],'wind stress (Pa)',7)
%}
%{
xtk = 0:.02:.2;
tklbl = {'0';'';'0.04';'';'0.08';'';'0.12';'';'0.16';'';'0.2'};
colorbar_h(xtk,tklbl, [.2 .87 .35 .02],'wind stress (Pa)',7)
%}
B = colorbar;
set(B, 'Position', [.9 .3 .02 .5]) % left,bottom,width,height
xlabel(B, '10^-7  N/m^(-3)') %'m/s') %'^oC') %'W / m^2') 'g / Kg') 'Pa'


fname = ['wsc_norotate_qwkskt_JPL_windStress_2002-2009_',a_or_d,'_average_w_vectors'];
saveas(gcf,['/Volumes/Lacie_Primary/kayla/oaflux_research/figures/quikscat/',fname],'png')







