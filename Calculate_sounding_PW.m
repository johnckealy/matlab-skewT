%calculate PW from a sounding for a given height
clear

load km1103_sonde

H = sonde(9).height;
T = sonde(9).temperature + 273.15;
Td = sonde(9).dewpointT;

for layer = 1:10:10000
    
    
Td_layer = Td(layer:layer+10);
meanTd = nanmean(Td_layer);

T_layer = T(layer:layer+10);
meanT = nanmean(T_layer);

H1 = H(layer);
H2 = H(layer+10);
meanH = (H1+H2)/2;
layer_thickness = H2 - H1;


e = 6.112 * exp ( (17.67 * meanTd) / (meanTd + 243.5) ); 

if H(layer) < 1000  
    vapor_density = 14;
else
    vapor_density =  (e / (461.5*meanT))*10^5;
end
    
pw(layer) = (vapor_density)*(layer_thickness)*(1/1000000)*(100/1);
pw_mm(layer) = 10.*pw(layer);


end

PW = nansum(pw_mm)