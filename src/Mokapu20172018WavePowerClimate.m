% script that pulls and plots spectral data for a cdip buoy
% from the CDIP Thredds server: see http://thredds.cdip.ucsd.edu

% The majority of historical cdip buoy data is stored in a single
% archive netcdf file (where XXX is the 3 digit cdip buoy number):
%
% http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/XXXp1/XXXp1_historic.nc
%
% but the most recent few months are found in the buoy realtime "rt" file
%
%  http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/XXXp1_rt.nc


close all
clearvars

% specify dataset and start/end dates of interest

% KANEOHE 198, WETS 225, 098 MAKAPU
% Get any 433 historic data for 2018 (1/1/2018 - 12/31/2018) 
dset = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/098p1/098p1_historic.nc';
% Display variable names in the netcdf file in the command window
ncdisp(dset,'/','min')
% just get May 2017 - Apr 2018 data
dstart = datenum(2017,5,1,0,0,0);
dend = datenum(2018,4,30,23,59,59);

% load time variable, find time indices needed
wtime = ncread(dset,'waveTime');
wtime = double(wtime)/86400 + datenum(1970,1,1);
dindices = find(wtime >= dstart & wtime <= dend);
dlength = size(dindices,1);
nstart_1d = dindices(1);
ncount_1d = dlength;
ptime = wtime(dindices(1):dindices(ncount_1d));  
  
% load qc flags; records good where qcflag = 1 
qcflag = ncread(dset,'waveFlagPrimary',nstart_1d,ncount_1d);
mask_1d = qcflag==1;
  
% load frequency and frequency bandwidths
freq = ncread(dset,'waveFrequency');
fbw = ncread(dset,'waveBandwidth');
  
% load wave energy densities
nstart_2d = [1, dindices(1)];
ncount_2d = [size(freq,1), dlength];
a0 = ncread(dset,'waveEnergyDensity',nstart_2d,ncount_2d);
mask_2d = repmat(mask_1d,1,size(freq,1))';
a0(mask_2d~=1) = NaN;
  
% get Hs from a0 freq spectra
hs=4*sqrt(fbw'*a0)';
  
% load 2d var mean direction (dmf)
dmf = ncread(dset,'waveMeanDirection',nstart_2d,ncount_2d);
mask_2d = repmat(mask_1d,1,size(freq,1))';
dmf(mask_2d~=1) = NaN;
% convert back to a1, b1 fourier coeffs for directional moment
% integrations
a1=cosd(dmf);b1=sind(dmf);
  
%
% get wave energy (aka energy flux, aka power) direction, de (weighted by a0(f)/f)
%
a1e=((1./(freq)')*(a1.*(a0.*fbw)))./(fbw'*a0);% a1 weigthed by a0(f)/f
b1e=((1./(freq)')*(b1.*(a0.*fbw)))./(fbw'*a0);% b1 weigthed by a0(f)/f
de=atan2d(b1e,a1e);de=de';% energy direction in degrees
  
% get wave energy period (aka energy flux, aka power), te (weighted by a0(f)/f)
te=((1./(freq)')*(a0.*fbw))./(fbw'*a0);te=te';
te(mask_1d~=1) = NaN;  % mask data times that did not pass qc

% make month and year vectors
months=str2num(datestr(ptime,'mm'));
years=str2num(datestr(ptime,'yyyy'));

% wave power in kW/m of wavelength
pw=0.5*(hs.^2).*te; 



%% -----------------------------------------------

%------------------------------------------
%  make plot of monthly mean wave power 
%   vs direction for individual months
%   from May 2017 to Apr 2018
%-----------------------------------------

figure('position',[244    64   681   741]);
  
  um=unique(months,'stable');
  for m=1:12
  subplot(4,3,m);
  p1=plot(de,pw,'k.');
  if( m == 7)
      ylabel(['{ }                  ',...
          '        Wave Power (kW/m)'],...
          'fontsize',14,'fontweight','demi');end
  if(m == 11);xlabel('Wave Power Direction','fontsize',14,'fontweight','demi');end
  
  set(p1,'color',[.7 .7 .7])
  set(gca,'xlim',[-45 122.5],'xgrid','on');
  i=find(months == um(m));hold on;plot(de(i),pw(i),'k.','markersize',10)
  title(datestr(ptime(i(1)),'mmm yyyy'));
  %i=find(months == m & te > 11 );
  %if(~isempty(i));hold on;plot(dp(i),hs(i),'b.','markersize',10);end
  demon(m)=atan2d(sum(pw(i).*sind(de(i)))/sum(pw(i)),sum(pw(i).*cosd(de(i)))/sum(pw(i)));
  pwmon(m)=nansum(pw(i));
  mpwmon(m)=nanmean(pw(i));
  minpwmon(m)=nanmin(pw(i));
  maxpwmon(m)=nanmax(pw(i));
  
  
  pmax=max(pw);
  %plot([0 0],[0 1.1*pmax],'k-');
   %plot([90 90],[0 1.1*pmax],'k-');
   if( um(m) > 4 && um(m) < 11 )
   plot([demon(m) demon(m)],[0 1.1*pmax],'g-','linewidth',2);
   title(datestr(ptime(i(end)),'mmm yyyy'),'color','g',...
   'fontsize',12,'fontweight','bold');
   else
   plot([demon(m) demon(m)],[0 1.1*pmax],'r-','linewidth',2);
   title(datestr(ptime(i(end)),'mmm yyyy'),'color','r',...
   'fontsize',12,'fontweight','bold');
   end
   
   plot([49.2788 49.2788],[0 1.1*pmax],'k--','linewidth',2); %long-term mean dir
  
  set(gca,'ylim',[0 1.1*pmax]);
  xticks([-45 0 45 90]);xticklabels({'NW','N','NE','E'});
  set(gca,'fontweight','bold');
  end
  
  pfile='MokapPowerMay2017Apr2018.png';
fprintf(1,' Making image %s\n',pfile);
print(gcf,'-dpng','-r300','-loose',pfile);
  
 
  
