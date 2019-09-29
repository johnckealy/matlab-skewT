function Plot_skewT(sonde,S_num)
%function to plot a skewT logP diagram from sounding data. The input is a
%structure which must contain the fields temperature, pressure and
%dewpointT (in that spelling). The second input specifies the sonde's number, or index if the
%structure "sonde" is multidimensional. Optional additional inputs for 
%wind barbs are wind_u and wind_v. 
%The code uses the lift parcel method, each point being ascended
%through the atmosphere. Stability indices (CAPE, KI, LI) are also
%included. Precipitable water is another optional input if desired.  


%J. Kealy 2011

if nargin<2
    S_num=1;
end

Tsonde=sonde(S_num).temperature;
Psonde=sonde(S_num).pressure;
Tdsonde=sonde(S_num).dewpointT;
if isfield(sonde,'wind_u')==1 & isfield(sonde,'wind_v')==1;
wind_u=-sonde(S_num).wind_u; %an error with the sonde structs used reverses the 
wind_v=-sonde(S_num).wind_v; %sign of the wind vectors. Remove '-' for other uses. 
end

%-------start function here-----------

%set up the figure
figure
set(gca,'yscale','log','ydir','reverse')
hold on
set(gca,'fontweight','bold')
set(gca,'ytick',100:100:1000)
%set(gca,'xgrid','on')
%set(gca,'ygrid','on')
axis([-48 56 100 1050]);
W=[0.0001,0.000138,0.00017,0.00021,0.00031,0.00046,0.0007,0.0009,0.00119,...
    0.00165,0.00225,0.00355,0.00475,0.00676,0.01,0.1455,0.0204,0.05];
s=1;

%plot the background: adiabats, temperature lines, w-lines.
for T=-100:5:200  %background T
    [T_dry,~,WinT] = lift_parcel(T,NaN,1000,'dry',W(s)); %find dry adiabat
    [T_moist,p,~] = lift_parcel(T,NaN,1000,'moist',W(s)); %find moist adiabat
    T_moist_skew=T_moist-(30.*log(0.001.*p)); %skew the moist adiabat
    T_dry_skew=T_dry-(30.*log(0.001.*p));   %skew the dry adiabat
    T_skew=T-(30.*log(0.001.*p));  % skew the T
    T_moist_skew(T_moist_skew>50)=NaN; 
    T_dry_skew(T_dry_skew>50)=NaN; 
    T_skew(T_skew<-50)=NaN; 
    T_moist_skew(T_moist_skew<-50)=NaN; 
    T_dry_skew(T_dry_skew<-50)=NaN; 
    T_skew(T_skew>50)=NaN;
    WinT(WinT<-50)=NaN; 
    WinT(WinT>50)=NaN;
    plot(T_moist_skew(1:56:end),p(1:56:end),'r');
    plot(WinT,p,'g--');
    plot(T_skew,p,'k');   %plot them up
    plot(T_dry_skew,p,'b');
    if s==max(size(W))
        continue;
    else
        s=s+1;
    end
end


%plot up the sounding data
Tsonde_skew=Tsonde-(30.*log(0.001.*Psonde)); % skew the T
Tdsonde_skew=Tdsonde-(30.*log(0.001.*Psonde)); % skew the Td
plot(Tsonde_skew,Psonde,'k','linewidth',2)
plot(Tdsonde_skew,Psonde,'k--','linewidth',2)


%plot the convective parcel ascent
MLT = mean(Tsonde(1:15));
MLTd = mean(Tdsonde(1:15));
MLP = mean(Psonde(1:15));
[Tparcel,p,~] = lift_parcel(MLT,MLTd,MLP,'parcel',W);
Tparcel_skew=Tparcel-(30.*log(0.001.*p)); % skew the parcel T
plot(Tparcel_skew,p,'b','linewidth',2)


%plot up the wind barbs
pbarb=100:25:1030;
if exist('wind_u','var')==1 & exist('wind_v','var')==1;
for i=1:max(size(pbarb))
    [f1 indbarb]=min(abs(Psonde-pbarb(i))); 
    %The minus on wind_v is to compensate for the reversed y axis
    plot_windbarb(wind_u(indbarb),-wind_v(indbarb),Psonde(indbarb),53,0.3,'ms');
end
end


%Calculate the CAPE
CAPE = 0;
for P=1:max(size(p))-1
    [f1 indPsonde]=min(abs(Psonde-p(P))); %for CAPE
    Tenv = Tsonde(indPsonde);
    Tpar = Tparcel(P);
    Tenv_K = Tenv +273.15;
    Tpar_K = Tpar +273.15;
    R=287;
    AlphaEnv=R*Tenv_K/p(P);
    AlphaPar=R*Tpar_K/p(P);
    if (Tpar_K > Tenv_K) 
        C = (AlphaPar - AlphaEnv).*(p(P)-p(P+1));
        CAPE = CAPE + C;
    end
end

%Calculate the LI
[f1 indPs500]=min(abs(Psonde-500));  %for the LI
Tenv500 = Tsonde(indPs500)+273;
[f1 indp500]=min(abs(p-500));  %for the LI
Tpar500 = Tparcel(indp500)+273;
LI = Tenv500 - Tpar500;

%find the K-index
[f1 indPenv850]=min(abs(Psonde-850));
T850e = Tsonde(indPenv850)+273;
[f1 indPenv700]=min(abs(Psonde-700));
T700e = Tsonde(indPenv700)+273;
Td850e = Tdsonde(indPenv850)+273;
Td700e = Tdsonde(indPenv700)+273;

KI = T850e - Tenv500 + Td850e - (T700e - Td700e);
KI = KI-273.15;

%add precipitable water
if isfield(sonde,'ph2o')==1
    PW=sonde(S_num).ph2o;
end

%add time and date to title if available
if isfield(sonde,'fdoy')
    T=datestr(sonde(S_num).fdoy(1)-0.4166666);
    T=[T(1:6),' ',T(13:17),' (HST)'];
end
    
    
%titles and labels
titlelin2=['\rm{Sonde No.',num2str(S_num),' CAPE=',num2str(CAPE),' LI=',num2str(LI),' KI=',num2str(KI),' PW=',num2str(PW),'}'];
title(['\bf{SkewT-LogP Chart }' '\rm' T  char(10) titlelin2]);
xlabel('Temperature (C)')
ylabel('Pressure (hPa)')


