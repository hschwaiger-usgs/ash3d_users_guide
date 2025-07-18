clear all;

% Constants
Ric = 0.25;
vk  = 0.4;
nz  = 101;

% Assume a friction velocity
ustr = 1.0;

% Assume a ABL thickness
H  = 2000.0;
h  = 0.1*H;
z  = linspace(0,1.25*H,nz);

% Assume a Ri profile
Rimin=-0.10;
Risrf= 0.05;
Rimax= 0.10;
zs=0;
for i=1:nz
 if z(i)<h
  zs = i;
  % Linear in surface layer since Ri = z/L
  Ri(i) = Rimin + (Risrf-Rimin)*z(i)/h;
 elseif z(i)<H
  Ri(i) = Risrf + (Rimax-Risrf)*(z(i)-h)/(H-h);
 else
  Ri(i) = Rimax;
 end
end
% Get MO length
for i=1:nz
 if Ri(i)<0
  L(i) = z(i)/Ri(i);
 else
  L(i) = (z(i)/Ri(i))*(1-5*Ri(i));
 end
end
zeta = z./L;

% Loop over z and calculate Kv
for i=1:nz
 % Scaling terms for M-L
 sc1(i) = 1.0;
 sc2(i) = 1.0;
 sc3(i) = exp(-2*z(i)/H);
 if z(i)<=H
  sc4(i) = (1 - z(i)/H)^2;
  sc5(i) = (1 - z(i)/H)^1;
 else
  sc4(i) = 0;
  sc5(i) = 0;
 end
 % phi values
 if zeta(i)<0 % unstable
  p1=1.0;a1=4;c1=16.0;phisrf_BD(i) = p1*(1-c1*zeta(i))^(-1.0/a1);
  p1=1.0;a1=3;c1=15.0;phisrf_Cl(i) = p1*(1-c1*zeta(i))^(-1.0/a1);
  p1=1.0;a1=4;c1=16.0;phiml_BA(i)  = p1*(1-c1*zeta(i))^(-1.0/a1);
  p1=1.0;a1=3;c1= 7.0;phiml_TM(i)  = p1*(1-c1*zeta(i))^(-1.0/a1);
  p1=1.0;a1=2;c1=13.0;phiml_Ul(i)  = p1*(1-c1*zeta(i))^(-1.0/a1);
 else % stable
  p1=1.0 ;b1=5.0;phisrf_BD(i) = p1 + b1 * zeta(i);
  p1=1.0 ;b1=5.0;phisrf_Cl(i) = p1 + b1 * zeta(i);
  p1=0.74;b1=4.7;phiml_BA(i)  = p1 + b1 * zeta(i);
  p1=1.0 ;b1=5.0;phiml_TM(i)  = p1 + b1 * zeta(i);
  p1=1.0 ;b1=9.2;phiml_Ul(i)  = p1 + b1 * zeta(i);
 end
 % Diffusivity values
 Kv1(i) = ustr*vk*z(i)/phisrf_BD(i) * sc1(i);
 Kv2(i) = ustr*vk*z(i)/phisrf_Cl(i) * sc2(i);
 Kv3(i) = ustr*vk*z(i)/phiml_BA(i)  * sc3(i);
 Kv4(i) = ustr*vk*z(i)/phiml_TM(i)  * sc4(i);
 Kv5(i) = ustr*vk*z(i)/phiml_Ul(i)  * sc5(i);
end

figure;
subplot(1,3,1),plot(Ri,z,'b-',"linewidth",4)
axis([1.5*Rimin 1.5*Rimax 0 1.25*H])
xlabel('Ri')
ylabel('Height (m)')
subplot(1,3,2),plot(Kv3(1:nz),z(1:nz),'b-',"linewidth",4);hold on;
subplot(1,3,2),plot(Kv4(1:nz),z(1:nz),'m-',"linewidth",4);
subplot(1,3,2),plot(Kv5(1:nz),z(1:nz),'k-',"linewidth",4);
subplot(1,3,2),plot(Kv1(1:zs),z(1:zs),'ro',"linewidth",2);
subplot(1,3,2),plot(Kv2(1:zs),z(1:zs),'g+',"linewidth",2);
hold off;
axis([0 150 0 1.25*H])
xlabel('K_v (m^2/s)')


subplot(1,3,3),plot(sc3(1:nz),z(1:nz),'b-',"linewidth",4);hold on;
subplot(1,3,3),plot(sc4(1:nz),z(1:nz),'m-',"linewidth",4);
subplot(1,3,3),plot(sc5(1:nz),z(1:nz),'k-',"linewidth",4);
subplot(1,3,3),plot(sc1(1:zs),z(1:zs),'ro',"linewidth",2);
subplot(1,3,3),plot(sc2(1:zs),z(1:zs),'g+',"linewidth",2);
legend('BA_m','TM_,','Ul_m','BD_s','Cl_s')
axis([0 1.1 0 1.25*H])
hold off;
xlabel('Mixed-Layer scale factor')

print -dpng -color "-S1000,500" Kv_ML.png

