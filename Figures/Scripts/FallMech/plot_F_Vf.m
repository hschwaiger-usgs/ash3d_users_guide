clear all;

  global grav;
  global rho_air;
  global eta_air;
  global lam_air;

grav        =    9.81;        % m/s
rho_ash     = 2000.0;         % kg/m^3
d_0         =   50.0;         % diameter of reference particle (in microns)
F_0         =    1.0;
% Get atmospheric conditions
height =     20000.0; %m
[Pres Temp rho_air eta_air lam_air] = Func_IntStdAtmos(height);

Vf_0  = Func_FV_ash_WH_slip(rho_ash,d_0,F_0);

Vol0 = (4.0/3.0)*pi*(0.5*d_0)^3;
fac = 0.5;  % The minimum size of B as a fraction (<1.0) of d_0

% The ellipsoid is characterize by the long, intermediate, and short diameters given
% by A,B,C respectively
%  B and C will be independent variables and A will be determined either by assuming
%  a constant da or a constant Vol0
n   = 101;
B   = linspace(d_0*fac,d_0,n);
% C is over the same space as B, but masked to only those values <= B
for i = 1:n
 for j = 1:n
  if(B(j)>B(i))
    C(i,j) = NaN;
  else
    C(i,j) = B(j);
  end
 end
end


%  Assuming a constant da (Eq 6 of WH) of d_0
%   da = (1/3) * (A + B + C)
for i = 1:n
 for j = 1:n
  A_a(i,j)   = 3.0*d_0 - B(i) - C(i,j);
   % Here is the volume of the ellipse with the specified da
  Vol_a(i,j) = (4.0/3.0)*pi*(0.125*A_a(i,j)*B(i)*C(i,j));
   % And here is the F_a resulting from the assumption that da is constant
  F_a(i,j) = (B(i) + C(i,j))/(2.0*A_a(i,j));
 end
end


% Assuming a constant Vol of Vol0
%  Volume of ellipse = (4.0/3.0) pi (ABC)/8   : the 8 is since we are using diameter instead of radius
for i = 1:n
 for j = 1:n
  A_v(i,j) = (6.0*Vol0)/(pi*B(i)*C(i,j));
   % Here is the da of ellipse with the specified Vol0
  da_v(i,j) = (A_v(i,j) + B(i) + C(i,j))/3.0;
   % And here is the F_v resulting from the assumption that Vol0 is constant
  F_v(i,j) = (B(i) + C(i,j))/(2.0*A_v(i,j));
 end
end

% Finally, calculate the fall velocities
for i = 1:n
 for j = 1:n
   Vf_a(i,j) = Func_FV_ash_WH_slip(rho_ash,d_0,F_a(i,j));
   Vf_v(i,j) = Func_FV_ash_WH_slip(rho_ash,da_v(i,j),F_v(i,j));
 end
end



F_int = linspace(0,1.0,21);
V_int = linspace(0,1.0,21);
figure;
[c, h] = contour(B,B,F_a',F_int,'LineWidth',4);
clabel(c, h, 0.25:0.25:1.0, 'fontsize', 12);
axis([d_0*fac d_0 d_0*fac d_0])
axis equal;
title('F at constant da')
xlabel('B')
ylabel('C')
print -dpng F_Const_da.png

figure;
[c, h] = contour(B,B,Vf_a'/Vf_0,V_int,'LineWidth',4);
clabel(c, h, 0.4:0.2:1.0, 'fontsize', 12);
axis([d_0*fac d_0 d_0*fac d_0])
axis equal;
title('Vf/Vf_0 at constant da')
xlabel('B')
ylabel('C')
print -dpng Vfnorm_Const_da.png


figure;
[c, h] = contour(B,B,F_v',F_int,'LineWidth',4);
clabel(c, h, 0.25:0.25:1.0, 'fontsize', 12);
axis([d_0*fac d_0 d_0*fac d_0])
axis equal;
title('F at constant Vol')
xlabel('B')
ylabel('C')
print -dpng F_Const_Vol.png

figure;
[c, h] = contour(B,B,Vf_v'/Vf_0,V_int,'LineWidth',4);
clabel(c, h, 0.4:0.2:1.0, 'fontsize', 12);
axis([d_0*fac d_0 d_0*fac d_0])
axis equal;
title('Vf/Vf_0 at constant Vol')
xlabel('B')
ylabel('C')
print -dpng Vfnorm_Const_Vol.png


