%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate glacier melt balance over UAV survey area of Fountain Glacier

% Read in the AWS data
% Headers are ignored by the matlab function xlsread
clear all       % clears matlab memory
data=xlsread('Bylot_AWSdata_nohead.xlsx');      % AWS data from 2015-2016
ndata=size(data);
ndat=ndata(1,1);

filepath = 'G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ebal_data'; %path to data directory

albedofile1 = [filepath,'\j21albLLadj.tif'];
albedofile2 = [filepath,'\j23albLLadj.tif'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start procssing the AWS data
year=data(:,1);
dday=round(data(:,2),2);
Tn=data(:,5);
T=data(:,6);
Tx=data(:,7);

% convert SW data from MJ/m2 to W/m2
mj2w=1e6/3600;

% model coefficients
MF = 0.0002;
RF = 0.0026;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Process the AWS data for the 1-hour energy balance

% constants for the energy balance model - these need to be set properly
Tmelt=0;    % C

% miscellaneous QC or diagnostic thresholds
pdd_threshold=0.;
frost_threshold=0.;
SW_threshold=10;      % W/m2
max_albedo=0.9;

% set up solar radiation file location
location = [filepath,'\July2016_all_update'];
list = dir(location);
filenames = {list.name}; % first two cells blank

%% calculate the 1-hour and daily energy balance at all points
Nt=24;      % integrate 1-hour data; Nt = time steps per day
% Just do calculations for July 2016 (dday 183-215)
jul1=9014;
aug1=9758;
np=31*2;       % number of periods, month of July
ndat_period=12;    % 12 data points per period
for m=1:np
    startdat=1+(m-1)*ndat_period;   % start of 12-hr period
    enddat=m*ndat_period;           % end of 12-hr period
    startdatd=startdat+jul1-1;      % index, hourly dataset
    enddatd=enddat+jul1-1;        % index, hourly dataset
    day(m)=floor(dday(startdatd));   % current day
    
    % create arrays for cumulative temperature and radiation
    period_absrad = zeros(3294,7792);
    period_temp = zeros(3294,7792);
    period_tres = zeros(3294,7792);
    
    % run hourly loop
    for n=startdatd:enddatd
        k=n-jul1+1; % index for radiation
        hour = 1+n-startdatd;
        
        % hourly radiation balance
        if (n<9547)==1 %for days up to July 23
            ice_albedo=geotiffread(albedofile1);
        else
            ice_albedo=geotiffread(albedofile2);
        end
        
        [dSWin,R]=geotiffread([location '\' filenames{k+2}]); % hourly incoming radiation
        SWnet=dSWin.*(1-ice_albedo);
    
        % calculate temperature distribution
        dT = distribute_temp(SWnet(1885,911),SWnet,T(n));
        dT(dT<Tmelt)=0; % if temp below 0, no PDD accumulated
        
        % calculate residual temperature at AWS
        tf = 9.7721;
        Tres_aws = T(n)-(SWnet(1885,911)*tf);
        Tres_aws(Tres_aws<0)=0;
        
        % cumulative variables for 12-hr period
        period_absrad = period_absrad+SWnet;
        period_temp = period_temp+dT;
        period_tres = period_tres+Tres_aws;
    end
    % calculate melt in m
    period_melt = MF*period_tres+ RF*period_absrad;
    
    % save melt grid
    period(m) = 2-mod(m,2); % 1/2 for am/pm
    newname = ['G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ETI\DistributedMelt\melt_' num2str(day(m)) '_' num2str(period(m))]
    geotiffwrite(newname,period_melt,R)
      
    % save AWS location
    aT(m)=period_temp(1885,911);
    aTres(m)=period_tres(1885,911);
    albedo(m)=ice_albedo(1885,911);
    aSWnet(m)=period_absrad(1885,911);
    amelt(m)=period_melt(1885,911);

    waitbar(m/np)
end

Bylot_results = double([day' period' aT' aTres' albedo' aSWnet' amelt']);

save 'Distributed_ETI_July2016.dat' Bylot_results -ascii

%% calculate total melt for July 21-24 (noon-noon)
location = 'G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ETI\DistributedMelt\';
list = dir(location);
filenames = {list.name}; % first two cells blank
sday = 203;
eday = 206;

total_melt = zeros(3294,7792);

nr = length(filenames);
for i=3:nr
    strings = split(filenames{i},{'_','.'});
    day = str2num(strings{2});
    per = str2num(strings{3});
    if day==sday && per==2
        disp([day per])
        [melt,R] = geotiffread([location filenames{i}]);
    elseif day==eday && per==1
        disp([day per])
        [melt,R] = geotiffread([location filenames{i}]);
    elseif day>sday && day<eday
        disp([day per])
        [melt,R] = geotiffread([location filenames{i}]);
    else
        melt = zeros(3294,7792);
        disp([day per])
    end
    
    total_melt = total_melt + melt;
    
    n_amelt(i-2)=melt(1885,911);
end

total_melt = total_melt*-1;

imagesc(total_melt)
contour(total_melt)
contourf(total_melt)
pcolor(total_melt)
surf(total_melt, 'edgecolor', 'none'); view(2);

total_melt(total_melt<-100) = -9999;
geotiffwrite('G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ETI\j2124_ETI.tif',total_melt,R);

