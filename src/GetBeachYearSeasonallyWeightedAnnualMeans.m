function [BeachYears,AnnualMeans,Sparse]=GetBeachYearSeasonallyWeightedAnnualMeans(xdatetime,xpos)

vflag=0; % verbose flag; 1 = print details to screen ; 0 = quiet 

% Calculates annual mean beach widths or transect positions for 
% Oct 1- Sep 30 beach years from an input time series of values

%  Input
%   xdatetime =  timeseries dates as matlab datetimes (not datenums)
%   xpos = shoreline positions
%   
%     (eg. BeachYears=1985:2021)
%   beach years without sufficient data for a mean estimates will 
%   have NaN values.

%  Output
%    BeachYears = vector of beach years in time series
%    AnnualMeans = annual mean beach widths/position vector matching
%    Sparse = vector of sparse data flags 0 = no problems; 1 = sparse 
%             data algorithm was invoked.
%             

% figure out the range of beach years covered by the data (even if the
%   first and/or last year might not have enough coverage to make an
%   estimate.


% years with data
BeachYears=unique(year(xdatetime+calmonths(3)));
% use a continuous vector of years
BeachYears=BeachYears(1):BeachYears(end);

AnnualMeans=NaN(size(BeachYears));
AvgMonthCount=zeros(size(BeachYears));
Sparse=zeros(size(BeachYears));

% reduce to individual year-month means (ymmeans) 
    %xdatetime=datetime(sdates,'convertfrom','datenum');
    y=year(xdatetime);mn=month(xdatetime);
    uy=unique(y); %unique years
    ymmean(1:numel(uy(1):uy(end)),1:12)=NaN; % initialize matrix as NaNs
    ymcount(1:numel(uy(1):uy(end)),1:12)=0; % initialize matrix as zeros

    yrs=uy(1):uy(end); %loop through unique years
    for ny=uy
        for m=unique(mn(y == ny)) % loop through unique months in year
            ymmean(ny-uy(1)+1,m)=mean(xpos(y == ny & mn == m),'omitnan');
            idx=find(y == ny & mn == m);
            if ~isempty(idx)
               ymcount(ny-uy(1)+1,m)=numel(idx);
            end
        end
    end
   
    % use year-month means to calculate beach year annual means
    for y=yrs(2:end)
        iy=find(yrs == y); % find this year
        ipy=find(yrs == y-1); % find previous year
        iny=find(yrs == y+1); % find next year
        q1=NaN;q2=NaN;q3=NaN;q4=NaN;wmean=NaN;smean=NaN;
        if ~isempty(ipy)
            if vflag == 1 
            fprintf('%i %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n',...
                y,[ymmean(ipy,10:12) ymmean(iy,1:10)]);
            end
            q1=mean(ymmean(ipy,10:12),'omitnan'); % Oct-Dec prev year
            q2=mean(ymmean(iy,1:3),'omitnan'); % Jan-Mar year
            q3=mean(ymmean(iy,4:6),'omitnan'); % Apr-Jun year
            q4=mean(ymmean(iy,7:9),'omitnan'); % Jul-Sep year

%% ---------------------------------------------
% if all 4 quarters do not have data, use sparse data 
% balancing algorithm
            if sum(isnan([q1 q2 q3 q4])) > 0 % check for missing quarters
                idx=find(BeachYears == y);
                Sparse(idx)=1;
                if vflag == 1 
                fprintf('%i %s %4.1f %4.1f %4.1f %4.1f\n',y,'Making Sparse Data Adjustment of these 4 Qs: ',q1,q2,q3,q4)
                end
            
            % Goal of sparse data estimate is to have one (min) width
            %   in winter q1 or q2, and one (max) width in summer q3 or q4, to 
            %   minimize any two-season bias

            % sparse data sep-oct boundary adjustment 
            % 1. if q4 has no data but Oct has data, use Oct mean for q4
            %    and set q1 to mean of Nov Dec of previous year
            %     
            if isnan(q4) & ~isnan(ymmean(iy,10))
                q4=ymmean(iy,10);
                % take Oct of beginning of the beach year out of play
                q1=mean(ymmean(ipy,11:12),'omitnan');
               
                if vflag == 1 
                fprintf('%i %s  %3i %3i %3i %3i\n',y,'1. filling q4 using next Oct data',q1,q2,q3,q4)
                end
                
                
            end

            % 2. if q2 has no data but Apr has data use Apr mean for q2
            %     and make q3=mean(May Jun).
            if isnan(q2) & ~isnan(ymmean(iy,4))
                q2=ymmean(iy,4);
                % recalc q3 without april
                q3=mean(ymmean(iy,5:6),'omitnan');
               
                if vflag == 1 
                fprintf('%i %s %3i %3i %3i %3i\n',y,'2. adjusting q2 to use Apr data',q1,q2,q3,q4)
                end
                
    
            end

            % keep larger of q3 and q4 for summer
                if ~isnan(q3) & q4 >= q3
                    q3=NaN;
                elseif q3 > q4
                    q4=NaN;
                end

                % keep smaller of q2 and q1 for winter
                if ~isnan(q1) & q2 <= q1
                    q1=NaN;
                elseif q1 < q2
                    q2=NaN;
                end

            % if only have q1 and q4
            wmean=min([q1 q2]); % winter min mean
            smean=max([q3 q4]); % summer max mean
            else
            % normal two-season means when all q's have data    
            wmean=mean([q1 q2],'omitnan'); % winter mean
            smean=mean([q3 q4],'omitnan'); % summer mean
            end
        end
        if vflag == 1 
            fprintf('%i Q %4.1f %4.1f %4.1f %4.1f\n',y,q1,q2,q3,q4)
            fprintf('%i S %4.1f %4.1f\n',y,wmean,smean)
        end
         
        % make an annual mean if have both winter and summer values    
     
            ymean=mean([wmean smean]);
            % add to single vectors of data for all mops in reach
            idx=find(BeachYears == y);
            if ~isempty(idx)
                AnnualMeans(idx)=ymean;
            end
        if vflag == 1 
        fprintf('%i Y %4.1f \n',y,ymean)
        end
    end
    
end