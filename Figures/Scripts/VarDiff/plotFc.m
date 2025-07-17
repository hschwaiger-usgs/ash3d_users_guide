clear all;

% Note: for unstable conditions, a non-local scheme is probably better
vk  = 0.4;
Ric = 0.25;
lm  = 100;

n1=100;
n2=51;
Ri=linspace(-1,1,n1);
z0 = 0.1;
z = [1 5 10 50 100 500 1000 5000];

% Louis (1979)
b=9.4;
for i=1:n1
 for j=1:8
  if(Ri(i)<0)
   % Unstable atmosphere
   a = vk/log(z(j)/z0);          % Eq. 13
   c = 7.4*a*a*b*sqrt(z(j)/z0);  % Eq. 20
   Fc1(i,j) = 1 - b*Ri(i)/(1+c*sqrt(abs(Ri(i)))); % Eq. 14
  else
   % Stable atmosphere
   Fc1(i,j) = (1 + 0.5*b*Ri(i))^-2.0; % Eq. 15
  end
 end
end

% Stull (1988) p.209
% Jacobson (2005) p.251
for i=1:n1
  if(Ri(i)<=0)
    % Unstable atmosphere
    Fc2(i) = (1-18*Ri(i))^(0.5); % Stull Tab. 6.4
  elseif(Ri(i)<=Ric)
    % Weakly unstable atmosphere 
    Fc2(i) = (Ric-Ri(i))/Ric; % Jac Eq. 8.70 or Stull Tab. 6.4
  else
    % Stable atmosphere
    Fc2(i)=0.0;
  end
end

% Betts et al 1996 (also Hong 2006)
for i=1:n1
  if(Ri(i)<0)
   % Unstable atmosphere
   Fc3(i) = 1-8*Ri(i)/(1+1.746*sqrt(-Ri(i))); % Eq. A20b
  else
   % Stable atmosphere
   Fc3(i) = (1 + 5*Ri(i))^(-2);  % Eq. A18
  end
end

% Collins et al (2004)
for i=1:n1
  if(Ri(i)<0)
   % Unstable atmosphere
   Fc4(i) = (1-18*Ri(i))^(0.5); % Eq. 4.464
  else
   % Stable atmosphere
   Fc4(i) = 1/(1 + 10*Ri(i)+80*Ri(i)*Ri(i));  % Eq. 4.465
  end
end

iz = 4;
z(iz)
plot(Ri(:),Fc1(:,iz),'r-',"linewidth",4);hold on;
plot(Ri(:),Fc2(:),'g-',"linewidth",4);
plot(Ri(:),Fc3(:),'b:',"linewidth",4);
plot(Ri(:),Fc4(:),'k:',"linewidth",4);
hold off;
xlabel('Ri',"fontweight","bold","fontsize",20)
ylabel('F',"fontweight","bold","fontsize",20)
set(gca,"fontsize",20)
h=get (gcf, "currentaxes");
set(h,"fontweight","bold","linewidth",2)

h2=legend('Louis (1979)','Stull (1988)','Betts et al (1996)','Collins et al (2004)','location', 'north');
set(h2,"fontweight","bold","fontsize",20)
set(h2,"Box","off")

print -dpng -color "-S700,500" Kz_Fc.png

