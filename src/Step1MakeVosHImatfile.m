%%  Makes VosHI.mat  
%%  Reads CoastSat Transect info in CoastSat_transect_layer.geojson 
%%  Reorders the mainland transect info from Ensenada to the OR border 
%%    going from S to N.
%%  Places info in struct array VosHI
%
%  - needs the file CoastSat_transect_layer.geojson to be in path.
%  - includes the function deg2utm.m by Rafael Palacios
%
%  - skips transects on islands
%  - reorders some transect reaches that go N->S to S -> N
%  - adds UTM coords of transect end point.
%  - calculates new transect orientation based on UTM coords
%  - saves struct array of ordered Ensenda MX-> CA/OR border transect info 
%     in VosHI in VosHI.mat
%
%  VosHI.mat is used by BuildVosDataFinal.m to download an store the
%   corresponding online CoastSat transect shoreline position time series.
%

clearvars

% read in coastsat transect meta data and parse info to get HA transects

s=fileread('CoastSat_transect_layer.geojson');
ss=strsplit(s,'TransectId'); % split into individual transect character strings

n=0; % valid transect counter

%  loop through transect strings and parse out transect details
for nn=2:size(ss,2)
    st=strsplit(ss{nn});
    VosName=regexprep(st{2}, {'"',','}, '');
    
    if contains(VosName,'usa_HI')
             
    ReachNum=round(str2num(VosName(8:11)));
    TranNum=round(str2num(VosName(13:16)));
    
    % skip reach 50 transect 52 which can't be downloaded
    if ReachNum == 50 && TranNum == 52
        % skip
    else
    
    % save parsed transect info the the struct array VosHI
    %if ~isempty(ReachNum) % valid reach number means valid transect info
        n=n+1; % transect counter
    VosHI(n).Name=regexprep(st{2}, {'"',','}, '');
    VosHI(n).ReachNum=ReachNum;
    VosHI(n).TranNum=TranNum;
    VosHI(n).Orientation=str2double(regexprep(st{6}, {'"',','}, ''));
    VosHI(n).Slope=str2double(regexprep(st{8}, {'"',','}, ''));
    VosHI(n).Trend=str2double(regexprep(st{10}, {'"','},'}, ''));
    VosHI(n).BackLon=str2double(regexprep(st{15}, {'[[',','}, ''));
    VosHI(n).BackLat=str2double(regexprep(st{16}, {']',','}, ''));
    VosHI(n).OffLon=str2double(regexprep(st{17}, {'[',','}, ''));
    VosHI(n).OffLat=str2double(regexprep(st{18}, {']]','}},'}, ''));

    % UTM coords of transect back beach point
    [XutmV,YutmV,UTMzone]=deg2utm(VosHI(n).BackLat,VosHI(n).BackLon);
    VosHI(n).BackYutm=YutmV; 
    VosHI(n).BackXutm=XutmV;
     
    % UTM coords of transect offshore point
    [XutmV,YutmV,UTMzone]=deg2utm(VosHI(n).OffLat,VosHI(n).OffLon);
    VosHI(n).OffYutm=YutmV;  
    VosHI(n).OffXutm=XutmV;
     
     % save UTM zone and new UTM orientation
     VosHI(n).UTMzone=UTMzone;
     orient=90-atan2d(VosHI(n).OffYutm-VosHI(n).BackYutm,...
                      VosHI(n).OffXutm-VosHI(n).BackXutm);
     if orient < 0;orient=orient+360;end
     VosHI(n).OrientationUTM=orient;
     VosHI(n).BackMop=NaN; % default is no nearest Mop
       
     end % end valid transect if then

    end % end HI transect if then
   
% end transect info loop
end

% save struct array in mat file
save VosHI.mat VosHI

%% ------------------------------------------------------------------------
function  [x,y,utmzone] = deg2utm(Lat,Lon)
% -------------------------------------------------------------------------
% [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Description: Function to convert lat/lon vectors into UTM coordinates (WGS84).
% Some code has been extracted from UTM.m function by Gabriel Ruiz Martinez.
%
% Inputs:
%    Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
%    Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84
%
% Outputs:
%    x, y , utmzone.   See example
%
% Example 1:
%    Lat=[40.3154333; 46.283900; 37.577833; 28.645650; 38.855550; 25.061783];
%    Lon=[-3.4857166; 7.8012333; -119.95525; -17.759533; -94.7990166; 121.640266];
%    [x,y,utmzone] = deg2utm(Lat,Lon);
%    fprintf('%7.0f ',x)
%       458731  407653  239027  230253  343898  362850
%    fprintf('%7.0f ',y)
%      4462881 5126290 4163083 3171843 4302285 2772478
%    utmzone =
%       30 T
%       32 T
%       11 S
%       28 R
%       15 S
%       51 R
%
% Example 2: If you have Lat/Lon coordinates in Degrees, Minutes and Seconds
%    LatDMS=[40 18 55.56; 46 17 2.04];
%    LonDMS=[-3 29  8.58;  7 48 4.44];
%    Lat=dms2deg(mat2dms(LatDMS)); %convert into degrees
%    Lon=dms2deg(mat2dms(LonDMS)); %convert into degrees
%    [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Author: 
%   Rafael Palacios
%   Universidad Pontificia Comillas
%   Madrid, Spain
% Version: Apr/06, Jun/06, Aug/06, Aug/06
% Aug/06: fixed a problem (found by Rodolphe Dewarrat) related to southern 
%    hemisphere coordinates. 
% Aug/06: corrected m-Lint warnings
%-------------------------------------------------------------------------
% Argument checking
%
error(nargchk(2, 2, nargin));  %2 arguments required
n1=length(Lat);
n2=length(Lon);
if (n1~=n2)
   error('Lat and Lon vectors should have the same length');
end
% Memory pre-allocation
%
x=zeros(n1,1);
y=zeros(n1,1);
utmzone(n1,:)='60 X';
% Main Loop
%
for i=1:n1
   la=Lat(i);
   lo=Lon(i);
   sa = 6378137.000000 ; sb = 6356752.314245;
         
   %e = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sa;
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
   %alpha = ( sa - sb ) / sa;             %f
   %ablandamiento = 1 / alpha;   % 1/f
   lat = la * ( pi / 180 );
   lon = lo * ( pi / 180 );
   Huso = fix( ( lo / 6 ) + 31);
   S = ( ( Huso * 6 ) - 183 );
   deltaS = lon - ( S * ( pi / 180 ) );
   if (la<-72), Letra='C';
   elseif (la<-64), Letra='D';
   elseif (la<-56), Letra='E';
   elseif (la<-48), Letra='F';
   elseif (la<-40), Letra='G';
   elseif (la<-32), Letra='H';
   elseif (la<-24), Letra='J';
   elseif (la<-16), Letra='K';
   elseif (la<-8), Letra='L';
   elseif (la<0), Letra='M';
   elseif (la<8), Letra='N';
   elseif (la<16), Letra='P';
   elseif (la<24), Letra='Q';
   elseif (la<32), Letra='R';
   elseif (la<40), Letra='S';
   elseif (la<48), Letra='T';
   elseif (la<56), Letra='U';
   elseif (la<64), Letra='V';
   elseif (la<72), Letra='W';
   else Letra='X';
   end
   a = cos(lat) * sin(deltaS);
   epsilon = 0.5 * log( ( 1 +  a) / ( 1 - a ) );
   nu = atan( tan(lat) / cos(deltaS) ) - lat;
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   ta = ( e2cuadrada / 2 ) * epsilon ^ 2 * ( cos(lat) ) ^ 2;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   xx = epsilon * v * ( 1 + ( ta / 3 ) ) + 500000;
   yy = nu * v * ( 1 + ta ) + Bm;
   if (yy<0)
       yy=9999999+yy;
   end
   x(i)=xx;
   y(i)=yy;
   utmzone(i,:)=sprintf('%02d %c',Huso,Letra);
end
end