

% ofdm decoder
clear all
close all
clc


FILTERSPEC = '.wav';
TITLE = 'Pick a Recorded OFDM IQ wav file';
FILE = 'D:\gnu_radio_work';

[FILENAME, PATHNAME, FILTERINDEX] = uigetfile(FILTERSPEC, TITLE, FILE);
xstring = [PATHNAME FILENAME];

[x,fswav] = audioread(xstring);
Tswave = 1/fswav;

% tx bandwidth
fs = 100e3;  % final sample rate to convert IQ file to 

Nguard = 512;   % time between pulse and preamble fft data

Nchar = 256;  % make power of 2
bits = Nchar * 8;
N2 = bits/2;
Nfft = bits + 1024;  % bits bins plus 2 512 guard bins 


PreN = 8192;
Pulseseconds = PreN / fs;

sigtime = Nguard + 2*Nfft + 2*Nfft + 2*Nguard;





x = x(:,1) + 1i*x(:,2);
x = x.';

pw = Pulseseconds * fswav;
pw = round(pw);
pwvec = ones(1,pw);
xdet = fftfilt(pwvec,abs(x));
figure(1)
plot(abs(xdet))
title('abs of tone pulse');
[maxval,maxindex] = max(abs(xdet));

a = maxindex - (pw - 1);
b = maxindex;

tonesig = x(a:b);
figure(2)
plot(real(tonesig))
hold on
plot(imag(tonesig))
hold off

% compute freq of tonesig
td1 = [0 tonesig(1:end-1)]; % delay by one sample 
td1 = conj(td1);
z = tonesig .* td1;
z = angle(z);
z = (1/(2*pi)) * (z / Tswave);
z = z(56:end-56);
foffset = mean(z);

% toneNfft = round(4*fswav);
% tonesig = fft(tonesig,toneNfft);
% tonesig = fftshift(tonesig);
% 
% fvec = [0:toneNfft - 1];  % 0 to tonefft
% fvec = fvec / (toneNfft );  % 0 to 1
% fvec = fvec  *  fswav;  % 0 to fswave
% fvec = fvec - fswav/2 ;  % -fswave/2 to fswave/2
% 
% 
% figure(22)
% plot(fvec/1e3, abs(tonesig));
% title('fft of tone signal');
% xlabel('frequency khz');
% [maxtone,toneindex] = max(abs(tonesig));
% toneindex = toneindex - 1;
% foffset = toneindex - toneNfft/2;
% foffset = foffset * (fswav /toneNfft);


n = length(x);
t = [1:n]/fswav;
x = x .* exp(1i*2*pi*-foffset*t);   % freq shift sig to zero hz




figure(3)
plot(real(x))
hold on
plot(imag(x))
hold off


n = gcd(fs,fswav);
p = fs / n;
q = fswav / n;
x = resample(x,p,q);   % resample file to expected sample rate

hpw = ones(1,PreN);
xdet = filter(hpw,1,abs(x));
[maxv, maxi] = max(abs(xdet));
xstart = maxi ;
xend = xstart + sigtime ;


figure(32)
xplot = abs(x) / max(abs(x));
xplot = 20 * log10(xplot);
yplot = length(xdet);
yplot = zeros(1,yplot);
yplot(maxi) = -40;
plot(xplot)
hold on
plot(yplot)
hold on





x = x(xstart:xend);
a = Nguard  + round(Nfft/2) + 1 ;
b = a + Nfft   - 1 ;
xp = x(a:b);

a1 =  Nguard + 2*Nfft + round(Nfft/2) + 1 ;
b2 = a1 + Nfft  - 1 ;
xd = x(a1:b2);







Xp = fft(xp);
Xp = fftshift(Xp);   % pilot tones  amp = 1, phase = 0

figure(6)
plot(20*log10(abs(Xp)))
title('preamble fft');
hold on

Xd = fft(xd);
Xd = fftshift(Xd);       % data tones,  bit1 = 180 deg,  bit0 = 0 deg

plot(20*log10(abs(Xd)))
title('data fft');
hold off






det = angle(Xp ./ Xd);          % ratio to get delta phase diff
det = abs(det);

thres = pi/2;

% det = ( det >= thres) .* 1   +  (det < thres) .* 0;
det = det(513:end);
det = det(1:bits);   % bits for char  

figure(8)
plot(det)


det(det <= thres) = 1;
det(det > thres) = 0;


    figure(105)
    plot(det)
    title('decoded bits vector')



det = reshape(det,8,[]);
det = det';


charmes = zeros(1,Nchar);
w = [7:-1:0];
w = 2.^w;


for k = 1:Nchar
    x = det(k,:);
    x = x .* w;
    x = sum(x);
    charmes(k) = x;
end



    clc

   xstr =  char(charmes)

   msgbox(xstr,'replace')









