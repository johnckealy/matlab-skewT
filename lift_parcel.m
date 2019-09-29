function [Tparcel_C,p,W_inTcords] = lift_parcel(T_surf,Td_surf,P_surf,ascent_type,W)
% Takes in the surface T, Td and Ps and lifts the parcel adiabatically
% through the standard atmosphere. 
% inputs are in degrees celcius, hPa and ascent_type is a string equal to
% either 'dry', 'moist', or 'parcel' for the lapse rate to be followed.
%P_surf is in mb, T values are in degree celcius.


T_surf=T_surf + 273.15;
Tsd=Td_surf + 273.15;
p=P_surf:-1:100;
stop=0;
stop2=0;

%find the 1000mb pressure index
[f1 indp1000]=min(abs(p-1000));


%Tp(1)=T_surf;   %initial parcel temperature
for P=1:max(size(p))-1;   %integrate through standard atmosphere to 16km
    
    %find the w profile
    es = W*p(P)/(W + 0.622);
    W_inTcords(P) = (243.5*log(es/6.1121) )/ (log(es/6.1121) - 17.67);
        
    if P==1
        Tp(P)=T_surf;
    end
    
    %find the equivalent potential temperature. Needed for the MA nudging.     
    if P==indp1000
        es_surf=6.1121*exp(17.67*(Tp(P)-273)/((Tp(P)-273)+243.5));
        L=getLatentHeat(Tp(P));
        Ws_surf=0.622*(es_surf/(1000-es_surf));
        Te1000=Tp(P) + (L/1004)*Ws_surf;
    end
    %find the equivalent potental temperature of the dew point
    if Tp(P)<=Tsd & stop==0
        es_Td=6.1121*exp(17.67*(Tp(P)-273)/((Tp(P)-273)+243.5));
        L=getLatentHeat(Tp(P));
        Ws_Td=0.622*(es_Td/(p(P)-es_Td));
        Te_dew=Tp(P) + (L/1004)*Ws_Td;
        Td_theta = Te_dew*((1000/p(P)).^(2/7));
        stop=1;
    end

    T_theta1(P) = T_surf*((1000/p(P)).^(-2/7));
    T_theta2(P) = T_surf*((1000/p(P+1)).^(-2/7));
    GammaD=T_theta1(P) - T_theta2(P);
    
    %--------find Ws-----------------------
    Tp(P)=Tp(P)-273.15;

    es=6.1121*exp(17.67*Tp(P)/(Tp(P)+243.5));
    %es=es*100;
    
    L=getLatentHeat(Tp(P));
    
    Tp(P)=Tp(P)+273.15;
    
    Ws=0.622*(es/(p(P)-es));
    
    %nudge the moist adiabat
    if strcmp(ascent_type,'moist')==1 
        Te_parcel = Tp(P) + (L/1004)*Ws;
        Te = Te1000*((1000/p(P)).^(-2/7));
        
        if Te_parcel>Te
        Te_diff = Te_parcel - Te;
        Tp(P) = Tp(P) - Te_diff;
        end
    end
    
    if (strcmp(ascent_type,'parcel')==1 & Tp(P)<=Tsd)
        Te_parcel = Tp(P) + (L/1004)*Ws;
        if stop2==0
        Thetae_dew = Te_parcel*((1000/p(P)).^(2/7));
        stop2=1;
        end
        Te = Thetae_dew*((1000/p(P)).^(-2/7));
        
        if Te_parcel>Te
            Te_diff = Te_parcel - Te;
            Tp(P) = Tp(P) - Te_diff;
        end
    end
    
    
    
    
    % This is the moist adiabatic lapse rate
    GammaS=GammaD.*((1+(L.*Ws)/(287.*Tp(P))))/((1+(0.622*Ws*(L.^2))/(287*1004*(Tp(P).^2))));
    
    
    %--------------------------------------
    
    if strcmp(ascent_type,'moist')==1 | (strcmp(ascent_type,'parcel')==1 & Tp(P)<=Tsd)
        Gamma=GammaS;   %Initiate moist ascent
    else
        Gamma=GammaD;
    end
    
    Tp(P+1)=Tp(P)-Gamma;
    Tparcel(P)=Tp(P);
    Tparcel(P+1)=NaN;
   
end

W_inTcords(max(size(p)))=NaN;
Tparcel_C=Tparcel-273.15;


