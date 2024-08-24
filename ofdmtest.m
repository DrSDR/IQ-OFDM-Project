

clear all
close all


% pulse at start for sync and freq offset correction
N = 8192;  %number of samples for sync pulse
pream = ones(1,N);
fs = 100e3;  % sample rate of ofdm signal , bw = 100khz 

% gap between sync pulse and xp data vector 
guardN = 512;

% two seconds of pause inserted at start and end of iq file 
td = 2;
tdsamples = round(td * fs);
tdvec = zeros(1,tdsamples);


% number of char for text message 
Nchar = 256;  % make power of 2,  number of char for text
bits = Nchar * 8;
Nfft = bits + 1024;  % bits plus guard freq bins, gbins 512 + 512 



% make pilot vector with random phases 
xp = 23 * randn(1,Nfft);



% pilot vector with 512 guard on each side
xp = exp(1i*xp);  % mag = 1, phases = random 
xp(1:512) = 0;  % first 512 bins guard bins = 0 
xp(end-511:end) = 0;  % last 512 bins guard bins = 0


% get message from txt file
fid = fopen('C:\SDR_Work\textfile.txt');
mtxt = fread(fid);  % vector of bytes,  each char = a number, one char = one byte or 8 bits 
fclose(fid);
mtxt = mtxt';
if length(mtxt) >= Nchar
    mtxt = mtxt(1:Nchar);
else
    z = length(mtxt);
    z = Nchar - z;
    mtxt = [mtxt zeros(1,z)];
end
mbits = dec2bin(mtxt,8);   %  Nchar x 8 matrix 
mbits = reshape(mbits',1,[]);  % 1 x nbits  vector 
mbits = mbits(1:bits);  % 1010101010101 as string vector  
txbits = zeros(1,bits);  % transmit vector or freq bins  
% make BPSK data vector 
for k = 1:bits
    if mbits(k) == '1'
        txbits(k) = 1;  % bit one ,  zero phase:  exp(j * 0) = 1
    else
        txbits(k) = -1;  % bit zero,  180 phase:  exp(j*pi) = -1 
    end
end

% make data vector bpsk,  1 = 1,  0 = -1
xd = txbits;
xd = [zeros(1,512)   xd    zeros(1,512) ];  %512 guard on each side





% relate data vector to pilot vector
xd = xp .* xd;  % e^a * e^b =  e^(a + b)





% pilot time vector
xp = ifftshift(xp);
xp = ifft(xp);
xp = xp / max(abs(xp));

% data time vector
xd = ifftshift(xd);
xd = ifft(xd);
xd = xd / max(abs(xd));

% put the whole movie together
xz = zeros(1,guardN);
xt = [tdvec pream xz xp xp xd xd tdvec];

% ensure 90 % bandwidth, reduce out of band signals
hlpf = fir1(64,0.9);
xt = filter(hlpf,1,xt);
xt = xt / max(abs(xt));  % 1xN

datafile = [ real(xt).'  imag(xt).' ];  %  Nx2


% write iq file as wav file , use file in gnuradio flow blocks: 
 audiowrite('c:\SDR_Work\ofdmtest100khzIQ256Char.wav',datafile,fs,'BitsPerSample',16);


% % make a sim iq rx file
% 
% xt = resample(xt,4,1);
% N = length(xt);
% nvec = randn(1,N) + 1i*randn(1,N);
% 
% xt = 10^(-10/20) * xt + 10^(-40/20) * nvec ; % create 20db SNR 
% 
% t = [1:N]/(4 * fs);
% xt = xt .* exp(1i*2*pi*1235.45*t);  % apply freq offset  
% fs = 4*fs;
% datafile = [real(xt).' imag(xt).']; 
% 
% audiowrite('c:\SDR_Work\ofdmSimRX100khz256Char.wav',datafile,fs,'BitsPerSample',16);
% 




