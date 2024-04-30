%%signal1:
lamda = 3/24;
fs = 60e6;
Ts = 1/fs;
CIT = 0.5;
load("sursignal_010.mat");
load("symbolsignal_010.mat")

n_total = length(sursignal);
P = 10;
V = zeros(n_total,P);
for p = 1:P
    M = [zeros(1,p-1),symbolsignal.'];
    V(:,p) = M(1:n_total);
end
W = (V'*V)\(V'*sursignal);
U = sursignal-V*W;

signal_number = 4;
fd = -80:2:80;
FFT = zeros(7,n_total);
symbolsignal_2 = symbolsignal;
sportsignal_2 = U;
symbolsignal_2 = symbolsignal_2.';
sportsignal_2 = sportsignal_2.';
nn = length(symbolsignal_2);
n = 1:length(symbolsignal_2);

doppler_time = zeros(signal_number,length(fd));
% for n_tau = 0:6
%     for n_fd = 1:length(fd)
%         ref = [zeros(1,n_tau),symbolsignal_2];
%         ref = ref(1:nn);
%         COR(n_tau+1,n_fd) = sum(conj(ref).*sportsignal_2.*exp(-1j*2*pi*n*Ts*fd(n_fd)));
%     end
% end

for n_tau = 0:6
    ref = [zeros(1,n_tau),symbolsignal_2];
    ref = ref(1:nn);
    FFT(n_tau+1,:) = abs(fftshift(fft(conj(ref).*sportsignal_2)));
end
COR=FFT(:,(-40+n_total/2+1:40+n_total/2+1));

[maxValue, maxIndex] = max(COR(1,:));
% [row, col] = ind2sub(size(COR), maxIndex);
doppler_time(1,:) = COR(1,:)/maxValue;
surf(lamda*fd,3*10^8*Ts*(0:6),abs(COR),'CDataMapping','scaled');
% mesh(Ts*(1:delta_n:n_total),lamda*fd,abs(doppler_time.'),'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');

% pcolor(Ts*((1:delta_n:n_total)+delta_n),lamda*fd,abs(doppler_time.')),colorbar
% yticks(lamda*fd);
%% signal2
lamda = 3/24;
fs = 60e6;
Ts = 1/fs;
CIT = 0.5;
n_total = length(sursignal1);
fd = -80:2:80;
COR = zeros(7,length(fd));
FFT = zeros(7,n_total);
symbolsignal_2 = symblesignal2;
sportsignal_2 = sursignal2;
symbolsignal_2 = symbolsignal_2.';
sportsignal_2 = sportsignal_2.';
nn = length(symbolsignal_2);
n = 1:length(symbolsignal_2);
for n_tau = 0:6
    ref = [zeros(1,n_tau),symbolsignal_2];
    ref = ref(1:nn);
    FFT(n_tau+1,:) = abs(fftshift(fft(conj(ref).*sportsignal_2)));
end
COR=FFT(:,(-40+n_total/2+1:40+n_total/2+1));
[maxValue, maxIndex] = max(COR(1,:));
% [row, col] = ind2sub(size(COR), maxIndex);
doppler_time(2,:) = COR(1,:)/maxValue;
% surf(lamda*fd,3*10^8*Ts*(0:6),abs(COR),'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');
%% signal3
lamda = 3/24;
fs = 60e6;
Ts = 1/fs;
CIT = 0.5;
n_total = length(sursignal1);
fd = -80:2:80;
FFT = zeros(7,n_total);
symbolsignal_2 = symblesignal3;
sportsignal_2 = sursignal3;
symbolsignal_2 = symbolsignal_2.';
sportsignal_2 = sportsignal_2.';
nn = length(symbolsignal_2);
n = 1:length(symbolsignal_2);

for n_tau = 0:6
    ref = [zeros(1,n_tau),symbolsignal_2];
    ref = ref(1:nn);
    FFT(n_tau+1,:) = abs(fftshift(fft(conj(ref).*sportsignal_2)));
end
COR=FFT(:,(-40+n_total/2+1:40+n_total/2+1));
[maxValue, maxIndex] = max(COR(1,:));
% [row, col] = ind2sub(size(COR), maxIndex);
doppler_time(3,:) = COR(1,:)/maxValue;
% surf(lamda*fd,3*10^8*Ts*(0:6),abs(COR),'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');
%% signal4
lamda = 3/24;
fs = 60e6;
Ts = 1/fs;
CIT = 0.5;
n_total = length(sursignal1);
fd = -80:2:80;
FFT = zeros(7,n_total);
symbolsignal_2 = symblesignal4;
sportsignal_2 = sursignal4;
symbolsignal_2 = symbolsignal_2.';
sportsignal_2 = sportsignal_2.';
nn = length(symbolsignal_2);
n = 1:length(symbolsignal_2);
for n_tau = 0:6
    ref = [zeros(1,n_tau),symbolsignal_2];
    ref = ref(1:nn);
    FFT(n_tau+1,:) = abs(fftshift(fft(conj(ref).*sportsignal_2)));
end
COR=FFT(:,(-40+n_total/2+1:40+n_total/2+1));
[maxValue, maxIndex] = max(COR(1,:));
% [row, col] = ind2sub(size(COR), maxIndex);
doppler_time(4,:) = COR(1,:)/maxValue;
% surf(lamda*fd,3*10^8*Ts*(0:6),abs(COR),'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');
%% dopplor-time
surf(CIT*(1:signal_number),lamda*fd,abs(doppler_time.'),'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');

