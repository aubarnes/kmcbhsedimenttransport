% Example code to 
%
% - load VosData for HI in VosData.mat 
% - step through the North Beach coastsat transects and
%       - calculate the annual mean shoreline positions
%       - save the desired range of years for all transects (ComparisonYears)
%       - remove the transect means to make annual position anomalies
%          (AnnualAnoms)
%       - Make various plots of AnnualAnoms

VosReach=16; % north beach reach number
VosTransects=0:23;  % north beach transect numbers
ComparisonYears=2000:2024; % Oct-Sep beach years to include

load VosDataHI.mat

Nt=0; % transect counter
AnnualAnoms=[]; % annual mean shoreline position array

for tn=VosTransects % transect loop

    % find the transect in VosData
 idx=find(strcmp({VosData.VosName},...
     ['usa_HI_' num2str(VosReach,'%4.4i') '-' num2str(tn,'%4.4i')]));

   % reduce the transect time series to annual means
 [BeachYears,AnnualMeans,Sparse]=...
     GetBeachYearSeasonallyWeightedAnnualMeans(VosData(idx).VosDatetimes,VosData(idx).VosX);
 
 % find and save the annual mean BeachYears that match the desired ComparisonYears
 [Lia,Lib]=ismember(BeachYears,ComparisonYears);
 Nt=Nt+1;
 AnnualAnoms(Nt,1:numel(ComparisonYears))=NaN; % initialize as NaNs
 % insert results for the comparison years and turn annual mean positions
 % into position anomalies by subtracting their mean
 AnnualAnoms(Nt,Lib(Lib > 0))=AnnualMeans(Lia)-mean(AnnualMeans(Lia),'omitnan'); 
 
 % turn into anomalies
 %AnnualAnoms(Nt,:)=AnnualAnoms(Nt,:)-mean(AnnualAnoms(Nt,:),'omitnan');

end

% make smoothed running mean version of anomaly matrix 
% (3 transects in alongshore, 3 years in time) 
SmoothedAnom=movmean(movmean(AnnualAnoms,3,2,'omitnan'),3,1);

% year to year change
Change=diff(AnnualAnoms,1,2);
SmoothChange=movmean(Change,3,1,'omitnan'); % only smooth in alongshore, not time

% set end of runway transects 6,7,8 (indices 7-9, ) to NaNs for both
%   versions of the anomaly matrix
AnnualAnoms(7:9,:)=NaN;
SmoothedAnoms(7:9,:)=NaN;
Change(7:9,:)=NaN;
SmoothChange(7:9,:)=NaN;

% Make a couple figures

figure('position',[ 10        188        1096         609]);
imagesc(ComparisonYears,0:numel(VosTransects)-1,AnnualAnoms,'AlphaData',~isnan(AnnualAnoms));
polarmap;colormap(flipud(colormap));
cb=colorbar;cb.Label.String='Anomaly (m)';
set(gca,'color','k','fontsize',16,'xtick',2000:2024)
xlabel('Oct-Sep Year');
ylabel('CoastSat Transect Number');
title('North Beach Annual Anomalies (Transect 0 = Pyramid; Black = insufficient data)')

%

figure('position',[ 100        148        1096         609]);

imagesc(ComparisonYears,0:numel(VosTransects)-1,SmoothedAnom,'AlphaData',~isnan(AnnualAnoms));
polarmap;colormap(flipud(colormap));
cb=colorbar;cb.Label.String='Anomaly (m)';
set(gca,'color','k','fontsize',16,'xtick',2000:2024)
xlabel('Oct-Sep Year');
ylabel('CoastSat Transect Number');
title({'North Beach Annual Anomalies (Transect 0 = Pyramid; Black = insufficient data)',...
    '3 year and 3 transect Smoothing'})

%
figure('position',[ 200        108        1096         609]);
plot(ComparisonYears,ComparisonYears*0,'k-')
hold on;
p(1)=plot(ComparisonYears,mean(AnnualAnoms,'omitnan'),'k.-','linewidth',2,...
    'DisplayName','All Transects (Transects 0-5;9-24)');

p(2)=plot(ComparisonYears,mean(AnnualAnoms(1:12,:),'omitnan'),'r.-','linewidth',2,...
     'DisplayName','Pyramid Rock (Transects 0-5;9-11)');
p(3)=plot(ComparisonYears,mean(AnnualAnoms(13:24,:),'omitnan'),'g.-','linewidth',2,...
    'DisplayName','North Beach (Transects 12-24)');
grid on;
set(gca,'fontsize',16,'xtick',2000:2024);
xlabel('Oct-Sep Year');
ylabel('Anomaly (m)');
title('Alongshore Averaged Annual Anomalies');
legend(p,'location','northwest')

figure('position',[ 289     80        1096         609]);

imagesc(ComparisonYears(2:end),0:numel(VosTransects)-1,SmoothChange,'AlphaData',~isnan(SmoothChange));
polarmap;colormap(flipud(colormap));
cb=colorbar;cb.Label.String='Anomaly Change (m/yr)';
set(gca,'color','k','fontsize',16,'xtick',2000:2024)
xlabel('Oct-Sep Year');
ylabel('CoastSat Transect Number');
title('Annual Anomaly Change, with 3 transect alongshore smoothing')

%%
figure('position',[ 343    52     1096     609]);
plot(ComparisonYears,ComparisonYears*0,'k-')
hold on;
p(1)=plot(ComparisonYears(2:end),mean(Change,'omitnan'),'k.-','linewidth',2,...
    'DisplayName','All Transects (Transects 0-5;9-24)');

p(2)=plot(ComparisonYears(2:end),mean(Change(1:12,:),'omitnan'),'r.-','linewidth',2,...
     'DisplayName','Pyramid Rock (Transects 0-5;9-11)');
p(3)=plot(ComparisonYears(2:end),mean(Change(13:24,:),'omitnan'),'g.-','linewidth',2,...
    'DisplayName','North Beach (Transects 12-24)');
grid on;
set(gca,'fontsize',16,'xtick',2000:2024);
xlabel('Oct-Sep Year');
ylabel('Anomaly Change (m/yr)');
title('Alongshore Averaged Annual Anomaly Change');
legend(p,'location','northwest')
