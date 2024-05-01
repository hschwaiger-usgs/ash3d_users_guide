function [Pres Temp rho_air eta_air lam_air] = Func_IntStdAtmos(H)

      % Given a height in m, returns:
      %  Pres    : Pa
      %  Temp    : K
      %  rho_air : kg/m3
      %  eta_air : Pa s
      %  lam_air : m

      eta0      = 1.8325e-5;    %! Ref visc (Pa s)
      suthcons  = 120.0;        %! Sutherland Constant (K)
      suthtref  = 296.16;       %! Sutherland Ref temperature (K)
      temper0   = 300.0;        %! Surf temp (K)
      pres0     = 101325.0;     %! Surf pres (Pa)
      skinz     = 7000.0;       %! pres skin depth (m)
      lpsr      = -0.007;       %! lapse rate (K/m)

      R_GAS_IDEAL  = 8.3144621;    %! Ideal gas constant (J /(kg K))
      R_GAS_DRYAIR = 286.98;       %! Specific gas constant of R=286.98 J /(kg K)
      MB_DRY_AIR   = 0.028966;     %! Molecular weight of dry air in kg/mol
      BoltzK       = 1.380658e-23; %! Boltzmann's constant kg m2 s-2 K-1 molec-1

      height_in_km = H * 1.0e-3;  

      % http://en.wikipedia.org/wiki/International_Standard_Atmosphere
      H_bm = [ 0.0  11.019 20.063 32.162 47.350 51.413 71.802 86.00];
      l_bm = [-6.5   0.0    1.0    2.8    0.0   -2.8   -2      0.0  ];
      t_bm = [15.0 -56.5  -56.5  -44.5   -2.5   -2.5  -58.5  -86.2];

      for i = 8:-1:2
       if height_in_km<=H_bm(i)
        interval = i;
       end
      end
      frac   = (height_in_km-H_bm(interval-1))/(H_bm(interval)-H_bm(interval-1));
      Temp_C = t_bm(interval-1)+ frac*(t_bm(interval)-t_bm(interval-1));
      Temp = Temp_C + 273.0;

        % Eq 1.8 of Wallace and Hobbs
      Pres = pres0 * exp(-H/skinz);
      %Temp = temper0 + lpsr*H;
        % Get the density (kg/m^3) of dry air via the ideal gas law and
        % Specific gas constant of R=286.98 J /(kg K)
      rho_air = Pres/(R_GAS_DRYAIR*Temp);

        % Get the dynamic viscosity (kg/(m s)) of air via
        % Sutherland's equation (Jacobson05 p. 102 Eq 4.54))
      eta_air = eta0*((suthcons+suthtref)/(Temp+suthcons))*(Temp/suthtref)^1.5;

        % Get mean-free-path of dry air
        % Eq. 9.6 of Seinfeld and Pandis
      lam_air = (2.0*eta_air)/(Pres * sqrt(8.0*MB_DRY_AIR/(pi*R_GAS_IDEAL*Temp)));
end 
