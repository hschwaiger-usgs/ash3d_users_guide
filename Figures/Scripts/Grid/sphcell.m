clear all;

R1 = 1000.0;
R2 = 1200.0;

d2r = pi/180;

ns = 101;

lam1 = 40.0; lam1 = lam1 * d2r;
lam2 = 60.0; lam2 = lam2 * d2r;
phi1 =  0.0; phi1 = phi1 * d2r;
phi2 = 15.0; phi2 = phi2 * d2r;

R   = linspace(R1,R2,ns);
lam = linspace(lam1,lam2,ns);
phi = linspace(phi1,phi2,ns);

% Bottom face of cell
b1x = R1*cos(phi1)*cos(lam);
b1y = R1*sin(phi1)*cos(lam);
b1z = R1*sin(lam);
b2x = R1*cos(phi)*cos(lam1);
b2y = R1*sin(phi)*cos(lam1);
b2z = R1*sin(lam1)*ones(1,ns);
b3x = R1*cos(phi2)*cos(lam);
b3y = R1*sin(phi2)*cos(lam);
b3z = R1*sin(lam);
b4x = R1*cos(phi)*cos(lam2);
b4y = R1*sin(phi)*cos(lam2);
b4z = R1*sin(lam2)*ones(1,ns);

% Upper face of cell
u1x = R2*cos(phi1)*cos(lam);
u1y = R2*sin(phi1)*cos(lam);
u1z = R2*sin(lam);
u2x = R2*cos(phi)*cos(lam1);
u2y = R2*sin(phi)*cos(lam1);
u2z = R2*sin(lam1)*ones(1,ns);
u3x = R2*cos(phi2)*cos(lam);
u3y = R2*sin(phi2)*cos(lam);
u3z = R2*sin(lam);
u4x = R2*cos(phi)*cos(lam2);
u4y = R2*sin(phi)*cos(lam2);
u4z = R2*sin(lam2)*ones(1,ns);

% z (radial) legs of cell
r1x = R*cos(phi1)*cos(lam1);
r1y = R*sin(phi1)*cos(lam1);
r1z = R*sin(lam1);
r2x = R*cos(phi1)*cos(lam2);
r2y = R*sin(phi1)*cos(lam2);
r2z = R*sin(lam2);
r3x = R*cos(phi2)*cos(lam1);
r3y = R*sin(phi2)*cos(lam1);
r3z = R*sin(lam1);
r4x = R*cos(phi2)*cos(lam2);
r4y = R*sin(phi2)*cos(lam2);
r4z = R*sin(lam2);


figure;
plot3(b1x,b1y,b1z,'k-',"linewidth",5);
hold on;
plot3(b2x,b2y,b2z,'k-',"linewidth",5);
plot3(b3x,b3y,b3z,'k-',"linewidth",5);
plot3(b4x,b4y,b4z,'k-',"linewidth",5);
plot3(u1x,u1y,u1z,'k-',"linewidth",5);
plot3(u2x,u2y,u2z,'k-',"linewidth",5);
plot3(u3x,u3y,u3z,'k-',"linewidth",5);
plot3(u4x,u4y,u4z,'k-',"linewidth",5);
plot3(r1x,r1y,r1z,'k-',"linewidth",5);
plot3(r2x,r2y,r2z,'k-',"linewidth",5);
plot3(r3x,r3y,r3z,'k-',"linewidth",5);
plot3(r4x,r4y,r4z,'k-',"linewidth",5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
lam1 = 20.0; lam1 = lam1 * d2r;
lam2 = 80.0; lam2 = lam2 * d2r;
% Draw lines of longitude
for iphi1 = -15:15:30;
  phi = iphi1 * d2r;
  lam = linspace(lam1,lam2,ns);

  b1x = R1*cos(phi)*cos(lam);
  b1y = R1*sin(phi)*cos(lam);
  b1z = R1*sin(lam);
  plot3(b1x,b1y,b1z,'k-',"linewidth",1);
end

% Draw z (radial) extensions of longitude lines
for iphi1 = -15:15:30;
  phi = iphi1 * d2r;
  R3 = R1 + 2*(R2-R1);
  R   = linspace(R1,R3,ns);
  lam = lam2;

  b1x = R*cos(phi)*cos(lam);
  b1y = R*sin(phi)*cos(lam);
  b1z = R*sin(lam);
  plot3(b1x,b1y,b1z,'k-',"linewidth",1);
end

% Draw lines of latidute
phi1 =-15.0; phi1 = phi1 * d2r;
phi2 = 30.0; phi2 = phi2 * d2r;
for ilam1 = 20:20:80
  lam1 = ilam1 * d2r;
  phi = linspace(phi1,phi2,ns);
  b1x = R1*cos(phi)*cos(lam1);
  b1y = R1*sin(phi)*cos(lam1);
  b1z = R1*sin(lam1)*ones(1,ns);
  plot3(b1x,b1y,b1z,'k-',"linewidth",1);
end

delR = R2-R1;
r=[R1 R1+delR R1+2*delR];
% Draw latitude lines on z-plane
for i = 1:3;
  phi = linspace(phi1,phi2,ns);
  R = r(i);
  lam = lam2;
  b1x = R*cos(phi)*cos(lam);
  b1y = R*sin(phi)*cos(lam);
  b1z = R*sin(lam)*ones(1,ns);
  plot3(b1x,b1y,b1z,'k-',"linewidth",1);
end

% Set up to draw hatch area
lam1 = 40.0; lam1 = lam1 * d2r;
lam2 = 60.0; lam2 = lam2 * d2r;
phi1 =  0.0; phi1 = phi1 * d2r;
phi2 = 15.0; phi2 = phi2 * d2r;
% Longitude hatches
hdens = 10;
R   = linspace(R1,R2,hdens);
lam = linspace(lam1,lam2,ns);
phi = linspace(phi1,phi2,ns);
for i = 1:hdens
  b1x = R(i)*cos(phi1)*cos(lam);
  b1y = R(i)*sin(phi1)*cos(lam);
  b1z = R(i)*sin(lam);
  plot3(b1x,b1y,b1z,'k-',"linewidth",1);
end
rdens = 30;
% Radial hatches
R   = linspace(R1,R2,ns);
lam = linspace(lam1,lam2,rdens);
phi = phi1;
for i = 1:rdens
  b1x = R*cos(phi)*cos(lam(i));
  b1y = R*sin(phi)*cos(lam(i));
  b1z = R*sin(lam(i));
  plot3(b1x,b1y,b1z,'k-',"linewidth",1);
end

%Annotations
fon = "Times";
fsize = 20;
tphi =  -7.5 * d2r; tlam = 10.0 * d2r; tr   = R1;
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'i-1',"fontname",fon,"fontsize",fsize)
tphi =  7.5 * d2r; tlam = 10.0 * d2r; tr   = R1;
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'i',"fontname",fon,"fontsize",fsize)
tphi =  22.5 * d2r; tlam = 10.0 * d2r; tr   = R1;
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'i+1',"fontname",fon,"fontsize",fsize)

tphi =  -22 * d2r; tlam = 30.0 * d2r; tr   = R1;
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'j-1',"fontname",fon,"fontsize",fsize)
tphi =  -23 * d2r; tlam = 50.0 * d2r; tr   = R1;
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'j',"fontname",fon,"fontsize",fsize)
tphi =  -30 * d2r; tlam = 70.0 * d2r; tr   = R1;
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'j+1',"fontname",fon,"fontsize",fsize)

tphi =  40 * d2r; tlam = 80.0 * d2r; tr   = R1 + 0.5*(R2-R1);
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'k',"fontname",fon,"fontsize",fsize)
tphi =  40 * d2r; tlam = 80.0 * d2r; tr   = R1 + 1.5*(R2-R1);
tx = tr*cos(tphi)*cos(tlam); ty = tr*sin(tphi)*cos(tlam); tz = tr*sin(tlam);
text(tx,ty,tz,'k+1',"fontname",fon,"fontsize",fsize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis equal;
axis off;
view(75,20);
hold off;
print('Spherical_cell.eps','-color','-depsc')

