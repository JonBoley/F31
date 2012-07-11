function [p,f] = oct3bank2(Pref,refdB,x,Fs);
% OCT3BANK Simple one-third-octave filter bank.
%    OCT3BANK(X) plots one-third-octave power spectra of signal vector X.
%    Implementation based on ANSI S1.11-1986 Order-3 filters.
%    Sampling frequency = Fs . Restricted one-third-octave-band
%    range (from 100 Hz to 10000 Hz). RMS power is computed in each band
%    and expressed in dB with 0 as reference level.
%
%    [P,F] = OCT3BANK(X) returns two length-21 row-vectors with
%    the RMS power (in dB) in P and the corresponding preferred labeling
%    frequencies (ANSI S1.6-1984) in F.
%
%    See also OCT3DSGN, OCT3SPEC, OCTDSGN, OCTSPEC.

% Author: Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
%         couvreur@thor.fpms.ac.be
% Last modification: Aug. 23, 1997, 10:30pm.

% References:
%    [1] ANSI S1.1-1986 (ASA 65-1986): Specifications for
%        Octave-Band and Fractional-Octave-Band Analog and
%        Digital Filters, 1993.
%    [2] S. J. Orfanidis, Introduction to Signal Processing,
%        Prentice Hall, Englewood Cliffs, 1996.
%
% Last Revised: 03/29/08

if nargin < 4,
    origindir = cd;
    [filename, pathname] = uigetfile('*.wav', 'Select file to analyze');
    if isequal(filename,0) | isequal(pathname,0)
        cd(origindir);
        close
        return
    end;
    [x,Fs] = wavread(strcat(pathname,char(filename)));
end;

N = 3; 					% Order of analysis filters.
f = [ 100 125 160, 200 250 315, 400 500 630, 800 1000 1250, ...
    1600 2000 2500, 3150 4000 5000, 6300 8000 10000]; % Preferred labeling freq.

i = find(f<(0.88.*(Fs./2)));
F = f(1:length(i));
nF = length(F);

ff = (1000).*((2^(1/3)).^[-10:nF-11]); 	% Exact center freq.
P = zeros(1,nF);
m = length(x);

% Design filters and compute RMS powers in 1/3-oct. bands
% 10000 Hz band to 1600 Hz band, direct implementation of filters.
for i = nF:-1:13
    [B,A] = oct3dsgn(ff(i),Fs,N);
    y = filter(B,A,x);
    P(i) = sum(y.^2)/m;
end
% 1250 Hz to 100 Hz, multirate filter implementation (see [2]).
[Bu,Au] = oct3dsgn(ff(15),Fs,N); 	% Upper 1/3-oct. band in last octave.
[Bc,Ac] = oct3dsgn(ff(14),Fs,N); 	% Center 1/3-oct. band in last octave.
[Bl,Al] = oct3dsgn(ff(13),Fs,N); 	% Lower 1/3-oct. band in last octave.
for j = 3:-1:0
    x = decimate(x,2);
    m = length(x);
    y = filter(Bu,Au,x);
    P(j*3+3) = sum(y.^2)/m;
    y = filter(Bc,Ac,x);
    P(j*3+2) = sum(y.^2)/m;
    y = filter(Bl,Al,x);
    P(j*3+1) = sum(y.^2)/m;
end

% Convert to decibels.
if nargin < 2,
    Pref = min(P); 				% Reference level for dB scale.
    refdB = 0;
end;

if (Pref == -99) & (refdB == -99)
    Pref = min(P); 				% Reference level for dB scale.
    refdB = 0;
end;  
Pref = Pref.^2;
idx = (P>0);
P(idx) = refdB+10*log10(P(idx)/Pref);
P(~idx) = NaN*ones(sum(~idx),1);

% Generate the plot
if (nargout == 0)
    bar(P);
    ax = axis;
    axis([0 nF+1 ax(3) ax(4)])
    set(gca,'XTick',[2:3:nF]); 		% Label frequency axis on octaves.
    set(gca,'XTickLabel',F(2:3:nF));  % MATLAB 4.1c
    xlabel('Frequency band [Hz]'); ylabel('Power [dB]');
    if nargin < 4,
    title(strcat('One-third-octave spectrum for -',filename));
    else
    title('One-third-octave spectrum');
    end;
    % Set up output parameters
elseif (nargout == 1)
    p = P;
elseif (nargout == 2)
    p = P;
    f = F;
end