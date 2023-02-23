clear all; close all; clc; 
 
sc=0 ; fcc=1 ; bcc=2 ; hcp=3;
 
 
%Parametre de maille Cai  : Al=4.05 / Ag=4.09 / Au=4.08 / Cu=3.615 / Ni=3.52 / Pd=3.89 / H=2.65
%Parametre de maille Ruda : H=2.65 / Ti=2.95
 
a0=2.95;
%Reseau : sc / fcc / bcc / hcp
lattice=hcp;   
 
rcut=6.9
 
%--------------------------------------------------------------------------
%Shell multiplicity and NN distance
%d'apres "Crystallography and Surface Structure: An Introduction for Surface
%Scientits and nanoscientists, Klaus Hermann, (p. 54)
 
% SC : Shell multiplicity and NN distance
if (lattice == sc)
  disp('Lattice SC');
  nNN(1)=6  ; dNN(1)=sqrt(1);
  nNN(2)=12 ; dNN(2)=sqrt(2);
  nNN(3)=8  ; dNN(3)=sqrt(3);
  nNN(4)=6  ; dNN(4)=sqrt(4);
  nNN(5)=24 ; dNN(5)=sqrt(5);
  nNN(6)=24 ; dNN(6)=sqrt(6);
end
 
% FCC : Shell multiplicity and NN distance
if (lattice == fcc)
  disp('Lattice FCC');
  nNN(1)=12 ; dNN(1)=sqrt(1/2);
  nNN(2)=6  ; dNN(2)=sqrt(2/2); 
  nNN(3)=24 ; dNN(3)=sqrt(3/2);
  nNN(4)=12 ; dNN(4)=sqrt(4/2);
  nNN(5)=24 ; dNN(5)=sqrt(5/2);
  nNN(6)=8  ; dNN(6)=sqrt(6/2);
end
 
% BCC : Shell multiplicity and NN distance
if (lattice == bcc)
  disp('Lattice BCC');
  nNN(1)=8  ; dNN(1)=sqrt(3/4);
  nNN(2)=6  ; dNN(2)=1;
  nNN(3)=12 ; dNN(3)=sqrt(2);
  nNN(4)=24 ; dNN(4)=sqrt(11/4);
  nNN(5)=8  ; dNN(5)=sqrt(3);
  nNN(6)=6  ; dNN(6)=2;
end
 
% Perfect HCP i.e. with c/a=sqrt(8/3) : Shell multiplicity and NN distance
if (lattice == hcp)
  disp('Lattice HCP');
  nNN(1)=12 ; dNN(1)=1;
  nNN(2)=6  ; dNN(2)=sqrt(2);
  nNN(3)=2  ; dNN(3)=sqrt(8/3);
  nNN(4)=18 ; dNN(4)=sqrt(3);
  nNN(5)=12 ; dNN(5)=sqrt(11/3);
  nNN(6)=6  ; dNN(6)=2;
end
 
dNN=a0*dNN;
 
 
[frho,fF,fphi,nr,nrho,dr,drho]=eamTi(a0,dNN,nNN);

%--------------------------------------------------------------------------
 
disp(['rcut = ' num2str(rcut) '  -  d6NN = ' num2str(dNN(6)) ]);;
 
E=0.0;
rhoTOT=0.0;
 
for i=1:6
  
  if (dNN(i) > rcut), break, end
    
  iweight = floor(dNN(i)/dr);
  weight  = dNN(i)/dr - iweight;
 
  rho=(1.0d0-weight)*frho(iweight+1)+weight*frho(iweight+2);
  rhoTOT=rhoTOT+nNN(i)*rho;  
  
  phi=(1.0d0-weight)*fphi(iweight+1)+weight*fphi(iweight+2);
  E=E+0.5*nNN(i)*phi;
  
  if (i == 1), rho_e=rho; end
  
end
 
iweight = floor(rhoTOT/drho);
weight  = rhoTOT/drho - iweight;
 
F=(1.0d0-weight)*fF(iweight+1)+weight*fF(iweight+2);
 
Ecoh=E+F;
 
disp(['a0 = ' num2str(a0) '  -  Ecohesion = ' num2str(Ecoh)]);
 
r=0:dr:dr*(nr-1);
rho=0:drho:drho*(nrho-1);
 
%re=dNN(1);
%
%figure
%plot(r,frho)
%figure
%plot(rho/rho_e,fF)
%axis([0.0 4.5 -5 4])
%figure
%plot(r/re,fphi)
%axis([0.5 3.0 -0.5 0.1])