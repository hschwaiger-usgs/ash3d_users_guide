clear all;

  global grav;
  global rho_air;
  global eta_air;
  global lam_air;

grav        =    9.81;        % m/s
rho_ash     = 2000.0;         % kg/m^3
d_0         =  500.0;
%d_0         =   64.0;         % diameter of reference particle (in microns)
% Get atmospheric conditions
height =     5000.0; %m
[Pres Temp rho_air eta_air lam_air] = Func_IntStdAtmos(height);

f = linspace(0.01,1,100);
g = linspace(0.01,1,100);
beta  = 2*f./(1+g);
gamma = 2*f.*g./(1+g);

n   = 100;
for i = 1:n
 Vf_wh(i) = Func_FV_ash_WH_slip(rho_ash,d_0,f(i));
 for j = 1:n
  Beta(i,j) =2*f(i)/(1+g(j));
  Gamma(i,j)=Beta(i,j)*g(j);
  Vf_gs(i,j) = Func_FV_ash_Gans_slip(rho_ash,d_0,f(i),g(j));
  %Vf_gs2(i,j)= Func_FV_ash_Gans_slip(rho_ash,d_0,f_2(i),g_2(j));
 end
end

[F,G] = meshgrid(f,g);
figure;
subplot(1,2,1),contour(g,f,Vf_gs)
%subplot(1,2,1),surf(G,F,Vf_gs);
shading interp
axis square
xlabel('G')
ylabel('F')
subplot(1,2,2),plot(Vf_wh,f,'b-',Vf_gs(:,100),f,'r-',Vf_gs(:,50),f,'r--',Vf_gs(:,10),f,'r:')
axis square
legend('W/H','Ganser (G=1.0)','Ganser (G=0.5)','Ganser (G=0.1)')
xlabel('V_s')
ylabel('F')

print('FigFV','-depsc')

