function [vset,Re] = Func_FV_ash_WH_slip(rho_ash,diam_part,F)

  global grav;
  global rho_air;
  global eta_air;
  global lam_air;

      % Calculates fall velocity (Wilson and Huang)
      % Given:
      %  particle density  : rho_ash : kg/m3
      %  particle diameter : d1      : um
      %  shape factor      : F       : dimensionless
      % Returns:
      %  particle fall vel : vest : m/s
      %  Reynolds number   : Re   : dimensionless


  %Cc = 1.0+(2.0*6.5e-8./dp).*(1.257+0.4.*exp(-1.1.*dp/(2.0*6.5e-8))); %Eq 8.36 of Seinfeld
  %vt = dp.^2*rhoa*g.*Cc/(18.0*mu_a);

  d1_mks = diam_part*1.0e-6;
  Kn = 2.0*lam_air/d1_mks;
  DahnekeFac = 1.0;
  Kna = Kn/DahnekeFac;
  Cslip = 1.0+ Kna * (1.257 + 0.4*exp(-1.1/Kna));

  vold = 1.0;                               % assume an initial settling velocity of 1 m/s

  Re = rho_air*vold*d1_mks/eta_air;               % Eq. 4  of Wilson79

  Ffac1 = F^(-0.828);
  Ffac2 = sqrt(1.07-F);

  % Get initial Cd from Wilson79, Eq 12
  Cd = (24./Re)*Ffac1 + Ffac2;

  vset = sqrt((4.0*rho_ash*d1_mks*grav)/(3.0*rho_air*Cd));   % Eq. 15 of WoodsBursik91

  while ((abs(vold-vset)/vset)>0.001)
    vold = vset;
    Re = rho_air*vold*d1_mks/eta_air;
    Cd = (24.0/Re)*Ffac1 + Ffac2;
    Cd = Cd/Cslip;
    vset = sqrt((4.0*rho_ash*d1_mks*grav)/(3.0*rho_air*Cd));
  end

  % If no slip is used, then use the explicit solution
  d3=d1_mks*d1_mks*d1_mks;
  vset2 = (2.44948974278318*sqrt(d3*Ffac2*grav*rho_air*rho_ash + ...
               54.0*eta_air*eta_air*Ffac1*Ffac1)-18.0*eta_air*Ffac1)/(d1_mks*Ffac2*rho_air)/3.0;
  vset = vset2;
end
