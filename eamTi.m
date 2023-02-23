function [frho,F,phi,nr,nrho,dr,drho]=eamTi(a0,dNN,nNN)
 
nr=1001;
nrho=1001;
 
re=dNN(1);
 
%--------------------------------------------------------------------------
rhostart=0.0 ; rhostop =3.35989;
rho=linspace(rhostart,rhostop,nrho);
 
rstart=0.0 ; rstop =6.1953;
r=linspace(rstart,rstop,nr);
 
dr=(r(end)-r(1))/(nr-1);
drho=(rho(end)-rho(1))/(nrho-1);
%--------------------------------------------------------------------------
 
 
% Effective embedding function for H
 
b=[20.44510 -2.897859 52.89785 0.412562];
 
for i=1:length(rho)
  F(i)=b(1)*rho(i)+b(2)+1.0d0/(b(3)*rho(i)+b(4));
end
 
 
%------------------------------------------------
% Normalized electron density for Ti
%------------------------------------------------
 
 
alpha=1.0;                              %alpha=1.211 et beta=2.32 d apres Foiles 1987
beta=3.255;  
 
%alpha=1.21;
%beta=2.32;
for i=1:length(r)
  frho(i)=alpha*exp(-beta*r(i));
end
 
for i=1:length(rho)
  F(i)=F(i)-rho(i)*15.0;
end
 
%------------------------------------------------
% Effective pair potentiel for Ti-Ti
%------------------------------------------------
 
p=3.83d0;
 
for i=2:length(r)
  
  if (r(i) < 2.0)
    Z=(1.0d0-r(i)/2.0d0)^p;
  else
    Z=0.0d0;
  end
  
  phi(i)=14.4*(Z*Z)/r(i);           %kc=14.4 eV.Ang.e^-2
  
  phi(i)=phi(i)+2.0*frho(i)*15.0;
  
end
 
phi(1)=phi(2);
 
figure
plot(r,frho);
title('rho(r)')
figure
plot(rho,F)
title('F(rho)')
axis([0.0 0.1 -7.5 10])
figure
plot(r/10,phi)
title('Phi(r)')
axis([0.15 0.55 -0.3 0.1])