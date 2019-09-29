function L=getLatentHeat(Tp)
%Enter the parcel temperature in celcius to find the appropriote latent
%heat value


%           LATENT HEAT LOOK UP TABLE
if (Tp<-40)
    L=2603;
elseif (Tp<-30) & (Tp>=-40)
    L=2575;
elseif (Tp<-20) & (Tp>=-30)
    L=2549;
elseif (Tp<-0) & (Tp>=-20)
    L=2525;
elseif (Tp<5) & (Tp>=0)
    L=2501;
elseif (Tp<10) & (Tp>=5)
    L=2489;
elseif (Tp<15) & (Tp>=10)
    L=2477;
elseif (Tp<20) & (Tp>=15)
    L=2466;
elseif (Tp<25) & (Tp>=20)
    L=2453;
elseif (Tp<30) & (Tp>=25)
    L=2442;
elseif (Tp<35) & (Tp>=30)
    L=2430;
elseif (Tp<40) & (Tp>=35)
    L=2418;
elseif (Tp>=40)
    L=2406;
end

L=L*1000;