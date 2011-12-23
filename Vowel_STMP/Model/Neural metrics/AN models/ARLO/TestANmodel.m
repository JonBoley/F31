% File: TestANmodel.m
%
%
%

clear

[noise,fs,N]=wavread('romnoise.wav');

tdres=1/fs;

No=10; %dB SPL
W=50e3;

OAL=No+10*log10(W);
RMS=sqrt(mean(noise.^2));

CurrSPL=20*log10(RMS/20e-6);
LEVfact=10^((OAL-CurrSPL)/20);
noiseNEW=noise*LEVfact;

RMSnew=sqrt(mean(noiseNEW.^2));
NewSPL=20*log10(RMSnew/20e-6);

cf=1000;
spont=50;
model=1;
species=9;
ifspike=1;

sout = an_arlo([tdres,cf,spont,model,species,ifspike],noiseNEW');

nrep=50;
[sptime,nspikes] = sgmodel([tdres, nrep],sout);