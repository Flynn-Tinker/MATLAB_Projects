clear all
close all
fclose all;
clc

addpath /data02/transfer/Chris/mfile_library/
addpath /data02/transfer/Chris/mfile_library/m_map/

%a_or_d = 'ascending'; % let's skip this switch for now

nr = 251; % gridded data rows ( found by loading a file)
nc = 301; % gridded data columns

load('/data02/users/gots/AutoPoleFlow/Auto_etimes_1982_now_10m.mat'); % began using this file when we changed to 10m winds and decided to use only auto etimes
arrive = etimes;
clear etimes

datadir = '/data02/pisco/modis_sst/ne_pacific/gridded_mat/';
list = fuf([datadir,'AMSRE-REMSS-L2P-amsr_l2b_v05_*.mat'],0,'normal');
junk = char(list);
% can use ftime to grab a subset of files to process
ftime = datenum(str2num(junk(:,36:39)),str2num(junk(:,40:41)),str2num(junk(:,42:43)),str2num(junk(:,45:46)),str2num(junk(:,47:48)),str2num(junk(:,49:50)));
clear junk

outdir = '/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/AMSRE_SST/composite/png/';

arrive = arrive(arrive>datenum(2002,6,11) & arrive<datenum(2009,9,20)); % limit this range to available SST
junk = datevec(arrive);
keep = junk(:,2)>4 & junk(:,2)<9; % limit to May-Aug here
arrive = arrive(keep);

wdw = -10:10; % time windows

CSST = nan(nr*nc,length(wdw)); % this will hold the averaged composite data
dcSST = CSST; % will use for data coverage

for ww = 1:length(wdw)
    
    disp(['processing window ',num2str(wdw(ww))])
    
    csst = nan(nr*nc,length(arrive)*5);  % we'll load up each applicable file and then take the average ( should not be more than 5 frames per day)
    count = 1;

    sequence = sprintf('%02d',wdw(ww)+11); % used in fname to line the sequence up in directories

    for aa = 1:size(arrive,1) % loop through each arrival    
    
        fidx = find(ftime>=(arrive(aa,1)+wdw(ww)) & ftime<(arrive(aa,1)+(wdw(ww)+1))); 
    
        if isempty(fidx) 
            clear fidx
            continue
        end    
            
        for ff = 1:length(fidx)
            
            junk = datevec(ftime(fidx(ff)));
            if junk(4)>=3 & junk(4)<14  % limit to night-time gmt (20:00 to 7:00 PDT)
                
                load([datadir,list{fidx(ff)}])            
                csst(:,count) = SST(:);
                count = count+1;
            end
            clear SST junk            

        end % of ff fidx loop
        clear fidx
        

     end % of aa arrive loop
     
     CSST(:,ww) = nanmean(csst,2);
     dcSST(:,ww) = sum(isfinite(csst),2);
     
         
     % I guess we can plot this stuff up?     
    figure(1) % this will show data coverage
    set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultaxesfontweight','bold')
    set(0,'defaultaxeslinewidth',1)
    set(gcf,'PaperPosition',[0 0 8.25 10],'color','w','renderer','painters') 

    cmap = colormap(jet(100));
     
    m_proj('lambert','long',[228 242],'lat',[30 45]);        
    hold on
    
    m_usercoast('/home/pisco/mfewings/nasaUpwellRelax/data/N_Amer_coast_NASA_l.mat','patch',0.9.*[1 1 1],'edgecolor','k')
    m_grid('box','on','tickdir','out','linestyle','none', ...
           'xtick',[228 242],'ytick',[30 45]);                 
                
    [pbx,pby] = political_boundaries;                
    m_plot(pbx,pby,'-','color','k')          
    
    X = lon; 
    Y = lat;
    Z = reshape(CSST(:,ww),nr,nc);
    
    pp = m_pcolor(X,Y,Z);
    set(pp,'linestyle','none')
    caxis([12 20]) 
    
    title(num2str(wdw(ww)))
    
    fname = ['AMSRE_SST_composite_',sequence,'_window_',num2str(wdw(ww))];          
    
    saveas(gcf,[outdir,fname],'png')
    
    close all
     
end % of ww window loop




