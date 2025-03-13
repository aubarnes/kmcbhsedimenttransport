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
dstart = datenum(2000,1,1,0,0,0);
dend = datenum(2019,6,30,23,59,59);

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


%------------------

% find unique years
uy=unique(years);

% step through years in data set
ym=0;
for y=uy(1):uy(end)
% step through months in years
  for m=1:12
      ym=ym+1;
      nbr=48*(datenum(y,m+1,1,0,0,0)-datenum(y,m,1,0,0,0));
      i=find(years == y & months == m);
      if(~isempty(i))
      n=length(hs(i) >= 0); % number of good  data records in month
      fprintf('%6i %6i %8i %8i %5.3f\n',[y m n nbr n/nbr]);
      if(n/nbr > 0.95)
      demon(ym)=atan2d(sum(pw(i).*sind(de(i)))/sum(pw(i)),sum(pw(i).*cosd(de(i)))/sum(pw(i)));
      pwrtime(ym)=ptime(i(1));
      pwmon(ym)=nansum(pw(i));
      mpwmon(ym)=nanmean(pw(i));
      minpwmon(ym)=nanmin(pw(i));
      maxpwmon(ym)=nanmax(pw(i));
      else
        pwrtime(ym)=NaN;
        demon(ym)=NaN;
        pwmon(ym)=NaN;
        mpwmon(ym)=NaN;
        minpwmon(ym)=NaN;
        maxpwmon(ym)=NaN;
      end
      else
        fprintf('No Data %6i %4i\n',[y m]);
        pwrtime(ym)=NaN;
        demon(ym)=NaN;
        pwmon(ym)=NaN;
        mpwmon(ym)=NaN;
        minpwmon(ym)=NaN;
        maxpwmon(ym)=NaN;
      end
      
% skip months with less than 95% data
   end
end

% energy direction for entire data set
det=atan2d(nansum(pwmon.*sind(demon))/nansum(pwmon),nansum(pwmon.*cosd(demon))/nansum(pwmon));
mp=nanmean(mpwmon);
  

%------------------------------------------
%  make plot of monthly mean wave power 
%   and mean direction
%-----------------------------------------

pwrmonths=str2num(datestr(pwrtime,'mm'));

figure('position',[2006 453 592 575]);
subplot(2,1,2);hold on;
plot(pwrtime,demon,'k-');plot(pwrtime,demon,'k.','markersize',15);
plot([ptime(1) ptime(end)],[det det],'k--');
datetick('x','yyyy','keeplimits');
set(gca,'xlim',[ptime(1) ptime(end)],'xgrid','on','box','on');

yVal = ylim;
dy=0.01*(yVal(2)-yVal(1));
y1 = [yVal(1)+dy,yVal(2)-dy,yVal(2)-dy,yVal(1)+dy];

for y=uy(1):uy(end)
x1 = datenum(datetime(y,10,1));
x2 = datenum(datetime(y+1,3,31));
patch([x1 x1 x2 x2],y1,[.9 .9 .9],'linestyle','none');
ax=gca;
ax.Children = ax.Children([end,1:end-1]);
end

smr=find(pwrmonths > 4 & pwrmonths < 11);
wnt=find(pwrmonths < 5 | pwrmonths > 10);
plot(pwrtime,demon,'k-');
plot(pwrtime(smr),demon(smr),'g.','markersize',15);
plot(pwrtime(wnt),demon(wnt),'r.','markersize',15);
plot([ptime(1) ptime(end)],[det det],'k--');

set(gca,'ylim',yVal);
xlabel('Date');
ylabel('Wave Power Direction')
set(gca,'fontsize',12,'fontweight','demi');
yticks([0 22.5 45 67.5 90]);yticklabels({'N','NNE','NE','ENE','E'});

%------------------

subplot(2,1,1);hold on;
plot(pwrtime,mpwmon,'k-');plot(pwrtime,mpwmon,'k.','markersize',15);
plot([ptime(1) ptime(end)],[nanmean(pw) nanmean(pw)],'k--');
datetick('x','yyyy','keeplimits');
set(gca,'xlim',[ptime(1) ptime(end)],'xgrid','on','box','on');

yVal = ylim;
dy=0.01*(yVal(2)-yVal(1));
y1 = [yVal(1)+dy,yVal(2)-dy,yVal(2)-dy,yVal(1)+dy];


for y=uy(1):uy(end)
x1 = datenum(datetime(y,10,1));
x2 = datenum(datetime(y+1,3,31));
patch([x1 x1 x2 x2],y1,[.9 .9 .9],'linestyle','none');
ax=gca;
ax.Children = ax.Children([end,1:end-1]);
end
plot(pwrtime,mpwmon,'k-');
p1=plot(pwrtime(smr),mpwmon(smr),'g.','markersize',15);
p2=plot(pwrtime(wnt),mpwmon(wnt),'r.','markersize',15);
plot([ptime(1) ptime(end)],[nanmean(pw) nanmean(pw)],'k--');

set(gca,'ylim',yVal);
xlabel('Date');
ylabel('Wave Power (kW/m)');
set(gca,'fontsize',12,'fontweight','demi');
legend([p1 p2],'May-Oct','Nov-Apr');

pfile='MokapuMonthlyPower.png';
fprintf(1,' Making image %s\n',pfile);
print(gcf,'-dpng','-r300','-loose',pfile);


%-----------------------------
% make power vs direction plot
%------------------------------
k=0;
  for i=-180:5:175
      k=k+1;
      j=find(de >= i & de < i+5);
      tpw(k)=sum(pw(j))/(sum(pw)*5);
      j=find(de >= i & de < i+5 & months > 4 & months < 11);
      tpwsum(k)=sum(pw(j))/(sum(pw)*5);tpwwin(k)=tpw(k)-tpwsum(k);
  end
  a1p=cosd(-180:5:175);b1p=sind(-180:5:175);
  
  deyr=atan2d(sum(pw.*sind(de))/sum(pw),sum(pw.*cosd(de))/sum(pw));
  j=find(months > 3 & months < 10);
  desum=atan2d(sum(pw(j).*sind(de(j)))/sum(pw(j)),...
      sum(pw(j).*cosd(de(j)))/sum(pw(j)));
  j=find(months < 4 | months > 9);
  dewin=atan2d(sum(pw(j).*sind(de(j)))/sum(pw(j)),...
      sum(pw(j).*cosd(de(j)))/sum(pw(j)));
  
  
  figure;
  p1=plot(-180:5:175,tpw,'k-','linewidth',2);hold on;
  plot(-180:5:175,tpw,'k.','markersize',15);
  p2=plot(-180:5:175,tpwwin,'r-','linewidth',2);
  plot(-180:5:175,tpwwin,'r.','markersize',15);
  p3=plot(-180:5:175,tpwsum,'g-','linewidth',2);
  plot(-180:5:175,tpwsum,'g.','markersize',15);
  
 
  set(gca,'xlim',[-45 122.5]);
  set(gca,'ylim',[0 1.1*max(tpw)]);
  xticks([-45 0 45 90]);xticklabels({'NW','N','NE','E'});
  set(gca,'xgrid','on');
 
  ylabel('Wave Power Density (kW/m-deg)');
  xlabel('Wave Power Direction');
  plot([deyr deyr],[0 1.1*max(tpw)],'k--','linewidth',2);
  plot([desum desum],[0 1.1*max(tpw)],'g--','linewidth',2);
  plot([dewin dewin],[0 1.1*max(tpw)],'r--','linewidth',2);
  set(gca,'fontsize',12,'fontweight','demi');
  legend([p1 p2 p3],'Total','Winter','Summer');
  
%pfile='MokapuMay17Apr18Power.png';
pfile='MokapPowervsDir.png';
fprintf(1,' Making image %s\n',pfile);
print(gcf,'-dpng','-r300','-loose',pfile);

%% -----------------------------------------------

%------------------------------------------
%  make plot of monthly mean wave power 
%   vs direction for individual months
%   from Aug 2018 to Jul 2019
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
  
 
  
