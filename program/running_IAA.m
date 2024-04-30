
figure;
signal_number = 20;
K = 50;
fd = linspace(-40, 40, K);
fs = 60e6;
Ts = 1/fs;
doppler_time = zeros(signal_number*2000,length(fd));
for order = 1:signal_number
    disp(order);
    filename_sursignal = sprintf('sursignal_%03d.mat', order);
    filename_symbolsignal = sprintf('symbolsignal_%03d.mat', order);
    load(filename_sursignal);
    load(filename_symbolsignal);
    n_total = length(sursignal);
%     fd2 = randi([-40, 40], 1, n_total);
%     w = 0:n_total-1;
%     array = exp(1j * 2 * pi * fd2 .* w * Ts);
    sursignal_adjusted = sursignal;
    symbol_adjusted = symbolsignal;
%     symbol_adjusted = symbolsignal.*array.';

%     for i = 1:10
%         sur_corr = xcorr(sursignal_adjusted,beacon_signal);
%         symbol_corr = xcorr(symbol_adjusted,beacon_signal);
%         [max_corr_sur, max_corr_index_sur] = maxk(sur_corr,1);
%         [max_corr_symbol, max_corr_index_symbol] = maxk(symbol_corr,1);
%         original_position_sur = max_corr_index_sur - length(sursignal) + 1;
%         original_position_symbol = max_corr_index_symbol - length(symbolsignal) + 1;
%         sursignal_adjusted(original_position_sur:original_position_sur+length(beacon_signal)-1) = 0;
%         symbol_adjusted(original_position_symbol:original_position_symbol+length(beacon_signal)-1) = 0;
%     end

    P = 10;
    V = zeros(n_total,P);
    for p = 1:P
        M = [zeros(1,p-1),sursignal_adjusted.'];
        V(:,p) = M(1:n_total);
    end
    W = (V'*V)\(V'*symbol_adjusted);
    U = symbol_adjusted-V*W;
    lamda = 3/24;
    CIT = 0.5;
    COR = zeros(12,K);
    symbolsignal_1 = sursignal_adjusted;
    sportsignal_1 = U;
    symbolsignal_1 = symbolsignal_1.';
    sportsignal_1 = sportsignal_1.';

    for segement  = 1:2000
        disp(segement);
        symbolsignal_2 = symbolsignal_1(floor((segement-1)*n_total/2000)+1:floor(segement*n_total/2000));
        sportsignal_2 = sportsignal_1(floor((segement-1)*n_total/2000)+1:floor(segement*n_total/2000));
        nn = length(symbolsignal_2);
        n = 1:length(symbolsignal_2);
        % for n_tau = 0:6
        %     for n_fd = 1:length(fd)
        %         ref = [zeros(1,n_tau),symbolsignal_2];
        %         ref = ref(1:nn);
        %         COR(n_tau+1,n_fd) = sum(conj(ref).*sportsignal_2.*exp(-1j*2*pi*n*Ts*fd(n_fd)));
        %     end
        % end
        
        for n_tau = 0:12
            ref = [zeros(1,n_tau),symbolsignal_2];
            ref = ref(1:nn);
    %         FFT(n_tau+1,:) = abs(fftshift(fft(conj(ref).*sportsignal_2)));
            COR(n_tau+1,:) = IAA(conj(ref).*sportsignal_2);
        end
%         COR=FFT(:,(floor(-40+nn/2+1):floor(40+nn/2+1)));
        
    %     [maxValue, maxIndex] = max(COR(1,:));
        [maxValue, maxIndex] = max(COR(:));
        [row, col] = ind2sub(size(COR), maxIndex);
        temp = COR(row,:);
    %     temp(temp<0.8*maxValue)=0;
        doppler_time((order-1)*5+segement,:) = temp;
        colorbar
        set(gca,'YDir','normal');
    end
end
surf(CIT*(1:signal_number*2000),lamda*fd,abs(doppler_time.'),'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');