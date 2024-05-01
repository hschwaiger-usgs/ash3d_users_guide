function [vset,Re] = Func_FV_ash_Gans_slip(rho_ash,diam_part,F,G)

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
  %Cd = (24./Re)*Ffac1 + Ffac2;
  beta=2*F/(1+G);
  gamma=beta*G;
  p = 1.6075;
  phi_sphere = (beta*gamma)^(2.0/3.0) * ((beta^p + gamma^p + (beta*gamma)^p)/3.0)^(-1.0/p);
  K1=3.0/(1.0+2.0*phi_sphere^(-0.5));
  K2=1.84148*(-log10(phi_sphere))^0.5743;
  K2 = 10.0^K2;
  Cd1 = 1.0 + 0.1118 * (Re*K1*K2)^0.6567;
  Cd1 = Cd1 * 24.0 /(K1*Re);
  Cd2 = 0.4305*K2/(1.0+3305.0/(Re*K1*K2));
  Cd = Cd1+Cd2;

  vset = sqrt((4.0*rho_ash*d1_mks*grav)/(3.0*rho_air*Cd));   % Eq. 15 of WoodsBursik91

  while ((abs(vold-vset)/vset)>0.001)
    vold = vset;
    Re = rho_air*vold*d1_mks/eta_air;
    Cd1 = 1.0 + 0.1118 * (Re*K1*K2)^0.6567;
    Cd1 = Cd1 * 24.0 /(K1*Re);
    Cd2 = 0.4305*K2/(1.0+3305.0/(Re*K1*K2));
    Cd = Cd1+Cd2;
    %Cd = Cd/Cslip;
    vset = sqrt((4.0*rho_ash*d1_mks*grav)/(3.0*rho_air*Cd));
  end

end
