% Script to plotEAM functions

close all, clear all, clc

eam=readEAM('Ti1.eam.fs','setfl')

embed=reshape(eam.embed',1,[]);
elecden=reshape(eam.elecden',1,[]);
pair=reshape(eam.pair',1,[]);

rho=0:eam.drho:(eam.nrho-1)*eam.drho;
r=0:eam.dr:(eam.nr-1)*eam.dr; 

pair=pair./r;

figure(1)
plot(r,pair'); axis([2 7 -1 5]); title('Pair')

figure(2)
plot(r,elecden); axis([2 7 -1 5]); title('Elecden')

figure(3)
plot(rho,embed); axis([0 60 -10 10]); title('Embed')

