clear all;

kap = 0.4;
lam1 = 30;
lam2 = 50;
lam3 = 100;
lam4 = 150;

z = linspace(0,5000,101);

num = kap*z;
den1 = 1+kap*z/lam1;
den2 = 1+kap*z/lam2;
den3 = 1+kap*z/lam3;
den4 = 1+kap*z/lam4;

el1 = num./den1;
el2 = num./den2;
el3 = num./den3;
el4 = num./den4;

plot(z,el1,'r-',"linewidth",4);hold on;
plot(z,el2,'g-',"linewidth",4);
plot(z,el3,'b-',"linewidth",4);
plot(z,el4,'k-',"linewidth",4);
hold off;
axis([0 5000 0 150]);
xlabel('Height (m)',"fontweight","bold","fontsize",20)
ylabel('l_m (m)',"fontweight","bold","fontsize",20)
set(gca,"fontsize",20)
h=get (gcf, "currentaxes");
set(h,"fontweight","bold","linewidth",2)

h2=legend('\lambda=30m','\lambda=50m','\lambda=100m','\lambda=150m','location', 'northwest');
set(h2,"fontweight","bold","fontsize",20)
set(h2,"Box","off")

print -dpng -color "-S700,500" Kz_MixLen.png

