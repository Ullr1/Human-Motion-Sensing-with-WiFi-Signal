clc;clear;
% % Y = zeros(4,4);
% % k = 0;
% % x = 1:4;
% % y = 1:4;
% % for i = 1:4
% %     for j = 1:4
% %         Y(i,j)=k;
% %         k=k+1;
% %     end
% %     k=0;
% % end
% % figure;
% % mesh(y,x,Y,'CDataMapping','scaled');
% % colorbar
% % set(gca,'YDir','normal');
% 
% % f = fft(symbolsignal);
% % n = 1:length(sportsignal);
% % figure(1)
% % plot(n,sportsignal2)
% % plot(n,10*log10(abs(f)));
% % n = 1:length(sportsignal);
% % f = fft(sportsignal);
% % figure(2)
% % plot(n,10*log10(abs(f)));
% % n = 1:length(sportsignal2);
% % f = fft(sportsignal2);
% % figure(3)
% % plot(n,10*log10(abs(f)));
% lamda = 3/24;
% fs = 60e6;
% Ts = 1/fs;
% CIT = 0.5;
% % delta_n = fs*CIT;
% n_total = length(sursignal);
% fd = -120:0.1:120;
% % floor(n_total/delta_n)
% % doppler_time = zeros(n_total,length(fd));
% % for index = 1:delta_n:n_total-delta_n
%     COR = zeros(7,length(fd));
%     symbolsignal_2 = sursignal1;
%     sportsignal_2 = sursignal2;
%     sportsignal2_2 = sportsignal2(index:index+delta_n);
%     symbolsignal_2 = symbolsignal_2.';
%     sportsignal_2 = sportsignal_2.';
%     sportsignal2_2 = sportsignal2_2.';
%     nn = length(symbolsignal_2);
%     n = 1:length(symbolsignal_2);
%     for n_tau = 0:6
%         for n_fd = 1:length(fd)
%             ref = [zeros(1,n_tau),symbolsignal_2];
%             ref = ref(1:nn);
%             COR(n_tau+1,n_fd) = sum(conj(ref).*sportsignal2_2.*exp(-1j*2*pi*n*Ts*fd(n_fd)));
%         end
%     end
%     [maxValue, maxIndex] = max(COR(:));
%     [row, col] = ind2sub(size(COR), maxIndex);
%     doppler_time(floor(index/delta_n)+1,:) = COR(row,:)/maxValue;
% end
% % mesh(lamda*fd,3*10^8*Ts*(0:8),abs(COR),'CDataMapping','scaled');
% mesh(Ts*(1:delta_n:n_total),lamda*fd,abs(doppler_time.'),'CDataMapping','scaled');
% colorbar
% set(gca,'YDir','normal');
% 
% % pcolor(Ts*((1:delta_n:n_total)+delta_n),lamda*fd,abs(doppler_time.')),colorbar
% % yticks(lamda*fd);
fs = 60e6;
load('sursignal_002.mat')
load('symbolsignal_017.mat')




n = 1:length(sursignal);
rs1=fftshift(20*log10(abs(fft(sursignal))));
f1 = linspace(0,fs,length(rs1))-fs/2;
figure(1);
subplot(2,1,1)
plot(f1,rs1)
subplot(2,1,2)
plot(n,sursignal);

rs2=fftshift(20*log10(abs(fft(symbolsignal))));
f2 = linspace(0,fs,length(rs2))-fs/2;
n = 1:length(symbolsignal);
figure(2);
subplot(2,1,1)
plot(f2,rs2)
subplot(2,1,2)
plot(n,sursignal);


% beacon_signal = sursignal(14859200:14861300);
% beacon_signal = sursignal(8557620:8559860);
% beacon_signal = sursignal(7904740:7921650);
beacon_signal = sursignal(15075000:15111000);
rs2=fftshift(20*log10(abs(fft(beacon_signal))));
f2 = linspace(0,fs,length(rs2))-fs/2;
n = 1:length(beacon_signal);
figure(3);
subplot(2,1,1)
plot(f2,rs2)
subplot(2,1,2)
plot(n,beacon_signal);


sursignal_adjusted = sursignal;
for i = 1:10
    disp(i);
    corr = xcorr(sursignal_adjusted,beacon_signal);
    [max_corr, max_corr_index] = maxk(corr,1);
    original_position = max_corr_index - length(sursignal) + 1;
    sursignal_adjusted(original_position:original_position+length(beacon_signal)-1) = 0;
end
rs2=fftshift(20*log10(abs(fft(sursignal_adjusted))));
f2 = linspace(0,fs,length(rs2))-fs/2;
n = 1:length(sursignal_adjusted);
figure(3);
subplot(2,1,1)
plot(f2,rs2)
subplot(2,1,2)
plot(n,sursignal_adjusted);

rs2=fftshift(20*log10(abs(fft(corr))));
f2 = linspace(0,fs,length(rs2))-fs/2;
n = 1:length(corr);
figure(4);
subplot(2,1,1)
plot(f2,rs2)
subplot(2,1,2)
plot(n,corr);
%% IAA
M = length(sursignal);
K = 160;
A = zeros(M,K);
f = linspace(-80, 80, K);
s = zeros(1,K);
for k = 1:K
    A(:,k) = doppler_steering_vector(M,f(k),fs);
end
R = eye(M);
temp = zeros(K,K);
diff = 10;
while(diff>1)
    P = zeros(K,K);
    for k = 1:K
        s(k) = (((A(:,k)'/R)*y))/(((A(:,k)'/R)*A(:,k)));
        P(k,k) = abs(s(k))^2;
    end
    R = A*P*A';
    diff = trace(abs(P-temp));
    temp = P;
end
%% IAA test
signal = conj(ref).*sportsignal_2;
Ml = length(signal);
fs = 30e6;
K = 50;
A = zeros(Ml,K);
f = linspace(-20, 20, K);
s = zeros(1,K);
out = zeros(1,K);
for k = 1:K
    A(:,k) = doppler_steering_vector(Ml,f(k),fs);
end
R = eye(Ml);
temp = zeros(K,K);
diff = 10;
while(diff>1)
    P = zeros(K,K);
    for k = 1:K
        s(k) = (((A(:,k)'/R)*signal))/(((A(:,k)'/R)*A(:,k)));
        P(k,k) = abs(s(k))^2;
    end
    R = A*P*A';
    diff = trace(abs(P-temp));
    temp = P;
end
for o = 1:K
    out(o) = sqrt(P(o,o));
end


    P = 10;
    V = zeros(n_total,P);
    for p = 1:P
        M = [zeros(1,p-1),ref_adjusted];
        V(:,p) = M(1:n_total);
    end
    W = (V'*V)\(V'*sur_adjusted);
    U = sur_adjusted-V*W;
