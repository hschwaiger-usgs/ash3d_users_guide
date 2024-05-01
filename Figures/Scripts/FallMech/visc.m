clear all;

temp=linspace(250,330,101);
eta=1.8325*0.00001*(416.16./(temp+120)).*(temp./296.16).^1.5;

plot(temp,eta)
