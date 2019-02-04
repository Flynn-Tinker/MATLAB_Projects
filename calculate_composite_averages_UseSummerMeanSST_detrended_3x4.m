clear all
close all
fclose all;
clc

addpath /data/data02/transfer/Chris/mfile_library/
addpath /data/data02/transfer/Chris/mfile_library/m_map/
addpath /data/data02/transfer/Chris/mfile_library/jlab/

%a_or_d = 'ascending'; % let's skip this switch for now

nr = 251; % gridded data rows ( found by loading a file)
nc = 301; % gridded data columns

outdir = [get_path('data01_sbc'),'internal/research/Collaborative_Research/upwelling_relaxation/May_August/AMSRE_SST/composite_anomaly_UseSummerMean/'];

%load('May_Aug_allData_regression_results.mat')
load('June_Aug_allData_regression_results_QF4.mat')


%%% Melton index
% load('/data02/users/gots/AutoPoleFlow/Auto_etimes_1982_now_10m.mat'); % began using this file when we changed to 10m winds and decided to use only auto etimes
% arrive = etimes;
% clear etimes

%%% Fit events using 500mb heights from day -1 to +3
load('positive_anomaly_flag_500mbht_1982_2010.mat')
keep = pos_flag(:,2)==1;
arrive = pos_flag(keep,1);
clear keep pos_flag




datadir = '/data/data02/transfer/Chris/raw_data/MODIS_SST/gridded_mat_qualityflag_4/';
list = fuf([datadir,'AMSRE-REMSS-L2P-amsr_l2b_v05_*.mat'],0,'normal');
junk = char(list);
% can use ftime to grab a subset of files to process
ftime = datenum(str2num(junk(:,36:39)),str2num(junk(:,40:41)),str2num(junk(:,42:43)),str2num(junk(:,45:46)),str2num(junk(:,47:48)),str2num(junk(:,49:50)));
clear junk


arrive = arrive(arrive>datenum(2002,6,11) & arrive<datenum(2009,9,20)); % limit this range to available SST
junk = datevec(arrive);
%keep = junk(:,2)>4 & junk(:,2)<9; % limit to May-Aug here
keep = junk(:,2)>5 & junk(:,2)<9; % limit to June-Aug here
arrive = arrive(keep);

wdw = -10:10; % time windows

CSST = nan(nr*nc,length(wdw)); % this will hold the averaged composite data
Cstd = CSST;
dcSST = CSST; % will use for data coverage

figure(1) % this will show data coverage
    set(gcf,'units','normalized','position',[.1 .1 .8 .8])
    set(gcf,'PaperPosition',[0 0 10 8.25],'color','w','renderer','painters') 


for ww = 5:16 %10:21  %1:length(wdw)
    
    disp(['processing window ',num2str(wdw(ww))])
    
    csst = nan(nr*nc,length(arrive)*5);  % we'll load up each applicable file and then take the average ( should not be more than 5 frames per day)
    count = 1;

    sequence = sprintf('%02d',wdw(ww)+11); % used in fname to line the sequence up in directories

    for aa = 1:size(arrive,1) % loop through each arrival    
    
        fidx = find(ftime>=(arrive(aa,1)+wdw(ww)-0.5) & ftime<(arrive(aa,1)+(wdw(ww)+1))); 
%         %fidx0 = find(ftime>=(arrive(aa,1)+wdw(11)) & ftime<(arrive(aa,1)+(wdw(11)+1))); % use day zero
%         fidx0 = find(ftime>=(arrive(aa,1)+wdw(6)-0.5) & ftime<(arrive(aa,1)+(wdw(6)+1))); % use day -5
%     
%         if isempty(fidx) || isempty(fidx0)
%             clear fidx
%             continue
%         end  
%         
%         % loop through fidx0 to generate the value used to calc the anomaly
%         num = 1;
%         dzro = [];
%         for ff = 1:length(fidx0)            
%             junk = datevec(ftime(fidx0(ff)));
%             if junk(4)>=3 & junk(4)<14  % limit to night-time gmt (20:00 to 7:00 PDT)
%                 
%                 load([datadir,list{fidx0(ff)}]) 
%                 
%                 % detrend SST(:) here
%                 yearday = ftime(fidx0(ff)) - datenum(junk(1),1,1);
%                 yval = (yearday.*gndmean.M) + gndmean.B;
%                 dzro(:,num) = SST(:)-yval(:);
%                 clear yval yearday
%                 
%                 %dzro(:,num) = SST(:);
%                 
%                 num=num+1;
%             end
%             clear SST junk 
%         end % of ff fidx0 loop
%         day0 = nanmean(dzro,2);
%         clear dzro num fidx0
        
        
        if true %~isempty(day0)   
        for ff = 1:length(fidx)
            
            junk = datevec(ftime(fidx(ff)));
            if junk(4)>=3 & junk(4)<14  % limit to night-time gmt (20:00 to 7:00 PDT)
                
                load([datadir,list{fidx(ff)}])  
                
                % detrend SST here
                yearday = ftime(fidx(ff)) - datenum(junk(1),1,1);
                yval = (yearday.*gndmean.M) + gndmean.B;
                %csst(:,count) = (SST(:)-yval(:)) - day0;
                csst(:,count) = (SST(:)-yval(:));% - gndmean.Tbar(:);
                
                
                
                
                %csst(:,count) = SST(:)-gndmean.SSTbar(:);
                
                
                count = count+1;
            end
            clear SST junk      

        end % of ff fidx loop
        clear fidx day0
        end

     end % of aa arrive loop
     
     CSST(:,ww) = nanmean(csst,2);
     Cstd(:,ww) = nanstd(csst,0,2);
     dcSST(:,ww) = sum(isfinite(csst),2);
     frames(ww,1) = count;
     
        % % calculate 95% confidence limits for these data
    conf(:,1) = CSST(:,ww) + (Cstd(:,ww).*tinv(.975,dcSST(:,ww))./sqrt(dcSST(:,ww))); 
    conf(:,2) = CSST(:,ww) - (Cstd(:,ww).*tinv(.975,dcSST(:,ww))./sqrt(dcSST(:,ww)));    
    
    %SSTnotsig  = ~((conf(:,1)<-gndmean.Tconf95(:)) | (conf(:,2)>gndmean.Tconf95(:)));
    SSTnotsig  = ~((conf(:,1)<0) | (conf(:,2)>0));
    SSTnotsig(isnan(conf(:,1))) = 0;
       
    
    %subplot(3,4,(ww-9))
    subplot(3,4,(ww-4))

    cmap = c2h(100);
    cmap = brighten(cmap,0);
    colormap(cmap)
     
%     m_proj('lambert','long',[228 242],'lat',[30 45]);        
%     hold on
%     
%     m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
%     m_grid('box','on','tickdir','out','linestyle','none', ...
%            'xtick',[228 242],'ytick',[30 45], ...
%            'xticklabel','','yticklabel','');                 


    m_proj('lambert','long',[220 250],'lat',[25 50]);        
    hold on
    
    m_usercoast('N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
    m_grid('box','on','tickdir','out','linestyle','none', ...
           'xtick',[220 250],'ytick',[25 50], ...
           'xticklabel','','yticklabel','');      
                
    [pbx,pby] = political_boundaries;                
    m_plot(pbx,pby,'-','color','k')          
    
    X = lon; 
    Y = lat;
    Z = reshape(CSST(:,ww),nr,nc);
    
    pp = m_pcolor(X,Y,Z);
    set(pp,'linestyle','none')
    caxis([-1 1]) 
    
%     %     % plot the non-sig points
%     notsig = [lon(:),lat(:),SSTnotsig];
%     notsig = notsig(SSTnotsig,:); % keep only not sig lon/lat
%     m_plot(notsig(:,1),notsig(:,2),'+','markersize',.1,'color',[.6 .6 .6])
    
    
    %title(num2str(wdw(ww)))
    
%     %fname = ['AMSRE_SST_composite_anomaly_UseDayZero_',sequence,'_window_',num2str(wdw(ww))];          
%     fname = ['AMSRE_SST_composite_anomaly_UseDayminus5_',sequence,'_window_',num2str(wdw(ww)),'_detrended_fit'];          
%     
%     saveas(gcf,[outdir,fname],'png')
%         
%     close all
     
end % of ww window loop

packboth(3,4)

saveas(gcf,[outdir,'SST_anomaly_composite_summerMean_JJA-6+5_extendedRange_pcolor.png'],'png')


return

% figure(1)  
% set(gcf,'PaperPosition',[0 0 8.25 10],'color','w','renderer','painters') 
% 
%         cmap = c2h(100); % cold (blue) to red (hot) white in the middle
%         %cmap = brighten(cmap,-.5);
%         colormap(cmap) 
%         
%         cb = colorbar;
%         set(cb,'yaxislocation','left','ylim',[1 99],'ytick',linspace(1,99,21),'yticklabel',{'-1';'';'';'';'';'-0.5';'';'';'';'';'0';'';'';'';'';'0.5';'';'';'';'';'1'})
%         
%         junk = get(cb,'position');
%         set(cb,'position',[junk(1:2) .025 junk(4)])
%         
%         
%         set(gca,'visible','off')
%         
%         print('-depsc2',[outdir,'SST_colorbar'])

