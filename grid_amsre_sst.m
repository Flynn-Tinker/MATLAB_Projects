clear all
close all
clc
fclose all;

addpath /data/data02/transfer/Chris/mfile_library/

indir = '/data/data02/transfer/Chris/raw_data/MODIS_SST/raw/';
outdir = '/data/data02/transfer/Chris/raw_data/MODIS_SST/gridded_mat_qualityflag_4/';

list = fuf([indir,'AMSRE-REMSS-L2P-amsr_l2b_v05_*.nc'],0,'normal');

[X,Y] = meshgrid(220:.1:250,25:.1:50); % output grid

for ii = 1:length(list)
    
    disp([num2str(ii),' ',list{ii}]) 

    [sst,lon,lat,time,dt,bias,sigma,rjct,conf,prox] = readL2Pcore([indir,list{ii}]);
    dn = datenum(1981,1,1,0,0,double(time));
    
    sst = double(sst)-273.15;
    lat = double(lat);
    lon = double(lon)+360;
           
    sst(double(prox)<4) = NaN; % quality flag
    
    nbad = lon<-32000;  % missing value is -32768
    sst(nbad) = NaN; 
    lat(nbad) = NaN; 
    lon(nbad) = NaN; 
        
    % there is a problem along the swath edges where triscattered interp is
    % connecting the dots through arcs in the swath.  by forcing a nan
    % along all sides of the available data this is eliminated.  can't
    % think of a more elegant way of doing it right now.
    [nr,nc] = size(sst);
    for cc = 1:nc
        nidx = find(isfinite(sst(:,cc)),1,'first');
        sst(nidx,cc) = NaN;
        nidx = find(isfinite(sst(:,cc)),1,'last');
        sst(nidx,cc) = NaN;
    end
    clear cc nc
    for rr = 1:nr
        nidx = find(isfinite(sst(rr,:)),1,'first');
        sst(rr,nidx) = NaN;
        nidx = find(isfinite(sst(rr,:)),1,'last');
        sst(rr,nidx) = NaN;
    end   
    clear rr nr nidx
    
    lon = lon(:);
    lat = lat(:);
    sst = sst(:);
    
    keep = (isfinite(lon)+isfinite(lat))==2;
    lon = lon(keep);
    lat = lat(keep);
    sst = sst(keep);
    
    
    % now interp these data onto our 0.1 degree grid
    F = TriScatteredInterp(double(lon),double(lat),sst);
    SST = single(F(X,Y));
    clear F
    
    if sum(isfinite(SST(:)))<2 % need at the very least two finite data points
        disp('skipped')
        clear sst lat lon dn mtime
        continue        
    end
    
    lon = X;
    lat = Y;
    mtime = dn;
    
    % save out with a more precise timestamp
    save([outdir,list{ii}(1:end-11),datestr(dn,30),'.mat'],'lon','lat','mtime','SST')
                
    clear sst lat lon dn mtime            
     
end % of ii list loop


 

