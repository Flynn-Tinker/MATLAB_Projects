clear all
close all
fclose all;
clc

addpath /data/data02/transfer/Chris/mfile_library/
addpath /data/data02/transfer/Chris/mfile_library/m_map/

%a_or_d = 'ascending'; % let's skip this switch for now

nr = 251; % gridded data rows ( found by loading a file)
nc = 301; % gridded data columns
count = 1;

datadir = '/home/pisco/gots/data/MODIS_SST/gridded_mat_qualityflag_4/';
list = fuf([datadir,'AMSRE-REMSS-L2P-amsr_l2b_v05_*.mat'],0,'normal');
junk = char(list);
% can use ftime to grab a subset of files to process
ftime = datenum(str2num(junk(:,36:39)),str2num(junk(:,40:41)),str2num(junk(:,42:43)),str2num(junk(:,45:46)),str2num(junk(:,47:48)),str2num(junk(:,49:50)));
clear junk

%outdir = '/data01/sbclter/internal/research/Collaborative_Research/upwelling_relaxation/May_August/AMSRE_SST/composite_anomaly/png/';

% preallocate space
y = nan(nr,nc,length(ftime),'single'); % this will hold the averaged composite data
x = nan(length(ftime),1);

gndmean.SSTbar = nan(nr,nc,'single');
gndmean.M = nan(nr,nc,'single');
gndmean.B = nan(nr,nc,'single');
gndmean.Tbar = nan(nr,nc,'single'); % NOTE: all stats performed on residuals
gndmean.Tstd = nan(nr,nc,'single');
gndmean.Tn = nan(nr,nc,'single');
gndmean.Tconf95 = nan(nr,nc,'single');
gndmean.r2 = nan(nr,nc,'single');

% load up all the required files
for ii = 1:length(ftime)
    
    junk = datevec(ftime(ii));
    %if (junk(4)>=3 & junk(4)<14) & (junk(2)>4 & junk(2)<9)  % limit to night-time gmt (20:00 to 7:00 PDT) and May/Aug
    if (junk(4)>=3 & junk(4)<14) & (junk(2)>5 & junk(2)<9)  % limit to night-time gmt (20:00 to 7:00 PDT) and May/Aug

        load([datadir,list{ii}])
        y(:,:,count) = SST;
        x(count,1) = mtime;
        count = count+1;
        
    end
    clear SST mtime junk
    

end % of ii list loop
clear ii
        

% now we're going to loop through each grid location and calculate M and B
for rr = 1:nr
    for cc = 1:nc
        
        Y = squeeze(y(rr,cc,:));
        X = x(isfinite(Y));
        Y = Y(isfinite(Y));
        
        if length(Y)>10 % let's set a 10 point minimum
            
            junk = datevec(X);
            junk = junk(:,1);
            yyyy = datenum(junk,ones(length(junk),1),ones(length(junk),1),zeros(length(junk),1),zeros(length(junk),1),zeros(length(junk),1));
            X = X-yyyy; % zero based yearday
            clear yyyy junk
            
            P = polyfit(X,Y,1); % get the equation for the line
            gndmean.M(rr,cc) = P(1);
            gndmean.B(rr,cc) = P(2);
            
            rY = Y - ((X.*gndmean.M(rr,cc)) + gndmean.B(rr,cc)); % calculate the residuals
            
            gndmean.SSTbar(rr,cc) = nanmean(Y);
            gndmean.Tbar(rr,cc) = nanmean(rY);  % these should be close to zero
            gndmean.Tstd(rr,cc) = nanstd(rY);
            gndmean.Tn(rr,cc) = sum(isfinite(rY));
            
            df_vec = getdof(X,2);
            gndmean.Tconf95(rr,cc) = gndmean.Tstd(rr,cc).*tinv(.975,df_vec)./sqrt(df_vec);
            
            R = corrcoef(X,Y);
            gndmean.r2(rr,cc) = R(1,2)^2;            
          
        end % of if branch
       
    clear X Y P df_vec rY R
    end % of cc column loop
    
end % of rr row loop

gndmean.lon = lon;
gndmean.lat = lat;


%save('May_Aug_allData_regression_results.mat','gndmean')
save('June_Aug_allData_regression_results_QF4.mat','gndmean')








