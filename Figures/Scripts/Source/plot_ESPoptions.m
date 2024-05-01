clear all;

Suz02raw=load('Suz_02.dat');
Suz04raw=load('Suz_04.dat');
Suz12raw=load('Suz_12.dat');
pointraw=load('point.dat');
lineraw=load('line.dat');
Umbrelraw=load('Umbrel.dat');

z = Suz02raw(:,1);
Sz02=Suz02raw(:,2);
Sz04=Suz04raw(:,2);
Sz12=Suz12raw(:,2);
Pnt=pointraw(:,2);
Lin=lineraw(:,2);
Um1=Umbrelraw(:,2);
Um2=Umbrelraw(:,3)*9;


plot(Sz04,z,'b-',Sz12,z,'r-',Um1,z,'g-',Um2,z,'g*',Lin,z,'m-')
