clear all;

dat1 = load('GS_Spurr_McGim.dat');
dat2 = load('GS_Webserver_dep.dat');
dat3 = load('GS_Ganser.dat');
dat4 = load('GS_Gaussian.dat');

mfr1=dat1(:,1);
dim1=dat1(:,2);
den1=dat1(:,3);
shp1=dat1(:,4);
phi1=dat1(:,5);

mfr2=dat2(:,1);
dim2=dat2(:,2);
den2=dat2(:,3);
shp2=dat2(:,4);
phi2=dat2(:,5);
nag2=find(den2>600.0);
agg2=find(den2<=600.0);

mfr3=dat3(:,1);
dim3=dat3(:,2);
den3=dat3(:,3);
shp3=dat3(:,4);
phi3=dat3(:,5);

mfr4=dat4(:,1);
dim4=dat4(:,2);
den4=dat4(:,3);
shp4=dat4(:,4);
phi4=dat4(:,5);


%%%%%%%%%%%%%%%%%%%%
figure; hold on;

lt=2;
plot(phi1,mfr1,'-go','LineWidth',lt)
plot(phi2(nag2),mfr2(nag2),'-ro','LineWidth',lt)
plot(phi3,mfr3,'-bo','LineWidth',lt)
plot(phi4,mfr4,'-ko','LineWidth',lt)

plot(phi2(agg2),mfr2(agg2),'-rx','LineWidth',lt,'MarkerSize',15)
plot(phi3(4),mfr3(4),'-bx','LineWidth',lt,'MarkerSize',15)
legend('Spurr','WS dep','Ganser','Gaussian');

hold off;
xlabel('\phi')
ylabel('Mass Fraction')

%saveas(gcf,'FigInputGSD.epsc')
print('FigInputGSD','-depsc')
