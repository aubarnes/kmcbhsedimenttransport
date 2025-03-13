
% Saves Online CoastSat shoreline time series for HI CoastSat transects
%
% Makes VosDataHI.mat containing the struct array VosData
%
% Other files needed: 
%  VosHI.mat (HI coastsat transect info)
%
% Saved VosData structural array fields
%
% VosData.VosName : (Online CoastSat Transect ID that is closest to the Mop) 
% VosData.Orientation : (Transect orientation based on its UTM end point coords)
% VosData.VosDatetimes : (vector of CoastSat data dates as matlab datetimes)
% VosData.VosX : (vector of shoreline xshore location, m, CoastSat transect coordinates)  
% 
%  Notes:
%
%      The VosData.VosDatetimes and VosData.VosX vectors are unchanged valid
%      time series info from the CoastSat website. Downloaded data for a transect
%      can include dates where there is no valid data ('None' in data
%      field for specific dates).  These are NOT included in VosData.VosDatetimes
%      and VosData.VosX

%%
clearvars
load VosHI.mat % load alongshore CoastSat transect info

% loop through CoastSat transects
nmt=0;
for nm=1:size(VosHI,2)

 fprintf('%i %s :',nm,VosHI(nm).Name);

 nmt=nmt+1;
 VosData(nmt).VosName=VosHI(nm).Name;
 VosData(nmt).Orientation=VosHI(nm).OrientationUTM;
 
 % get CoastSat transect shoreline time series from website
 %url=['http://coastsat.wrl.unsw.edu.au/time-series/' VosHI(nm).Name '/'];
 url=['http://coastsat.space/time-series/' VosHI(nm).Name '/'];
 s=webread(url);

 % split returned data be newlines
 ss=strsplit(s,'\n');

% add any valid data to the struct array time series
k=0; % data counter
for n=1:size(ss,2)
    if ~isempty(ss{n})
       st=strsplit(regexprep(ss{n},',',' '),' ');
       k=k+1;
       VosData(nmt).VosDatetimes(k)=datetime([st{1} ' ' st{2}]);
       VosData(nmt).VosX(k)=str2double(st{3});
    end
end

% check for and remove any time series shoreline position NaNs
idx=find(isnan([VosData(nmt).VosX]));
if numel(idx) > 0
    %fprintf(' *** Removing %i CoastSat shoreline "None"s\n',...
    %    numel(idx));
    VosData(nmt).VosDatetimes(idx)=[];
    VosData(nmt).VosX(idx)=[];
end

fprintf('CoastSat N= %i\n',numel([VosData(nmt).VosX]));
 
% end Mop-CoastSat Transect loop
end

fprintf('Saving Struct array VosData to VosDataHI.mat\n')
save VosDataHI.mat VosData