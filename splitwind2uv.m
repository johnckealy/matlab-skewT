function [u,v]=splitwind2uv(WD,WS)
% Splits a given wind speed and direction into constituent wind components,
% intended for use with plot_windbarb.m. 
% Inputs ==>  WD : Wind direction in degrees. Clockwise from 0degress=North
%             WS : Wind speed, units are only important when using in 
%             conjunction with plot_windbarb.m

% J. Kealy 2011


%find the u and v components
if WD>0 & WD<=90
    theta=WD;
    phi=90-WD;
    u=WS*sind(theta);
    v=WS*sind(phi);
elseif WD>90 & WD<=180
    theta=WD-90;
    phi=180-WD;
    v=-WS*sind(theta);
    u=WS*sind(phi);
elseif WD>180 & WD<=270
    theta=WD-180;
    phi=270-WD;
    u=-WS*sind(theta);
    v=-WS*sind(phi);
elseif WD>270 & WD<=360
    theta=WD-270;
    phi=360-WD;
    v=WS*sind(theta);
    u=-WS*sind(phi);
else
    disp('Error')
end

end