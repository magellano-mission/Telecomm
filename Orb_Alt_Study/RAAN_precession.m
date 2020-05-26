clc; clearvars;

Re = 3390;              %[km] Mars surface radius
J2 = 1960e-6;           %[-]
inc = deg2rad(-15);       %[rad]
mu = 42828;             %[km3/s2]

a = linspace(4000,8000,100);
e = linspace(0,0.5,100);


%add in some check to see if the periapsis is below the surface
%a = 6000;
%e = 0.5;

T = 2.*pi.*sqrt(a.^3./mu);
n = 2*pi./T;

wp = zeros(length(a),length(e));
rp = zeros(length(a),length(e));

for i = 1:length(a)
    for j = 1:length(e)
        rp(i,j) = a(i)*(1-e(j));
        if rp(i,j) > Re + 150
            wp(i,j) = -3/2*(Re^2/((a(i)*(1-e(j)^2))^2))*J2*n(i)*cos(inc);
            wp(i,j) = rad2deg(wp(i,j))*86400;
        else
            rp(i,j) = NaN;
            wp(i,j) = NaN;
        end
    end
end

subplot(1,2,1)
levels1 = -5:0.2:0;
contour(e,a,wp,levels1,'ShowText','on')
xlabel('Eccentricity')
ylabel('Semi-major Axis [km]')
title('Daily RAAN Shift [degrees]')
grid on

subplot(1,2,2)
levels2 = 250:250:5000;
contour(e,a,rp-Re,levels2,'ShowText','on')
xlabel('Eccentricity')
ylabel('Semi-major Axis [km]')
title('Altitude of Periapsis [km]')
grid on



        