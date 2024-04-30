figure;
signal_number = 20;
fd = -80:2:80;
fs = 60e6;
Ts = 1/fs;
delta_t = 0.1;
segement_number = 0.5/delta_t;
doppler_time = zeros((signal_number-1)*segement_number,length(fd));
for order = 1:signal_number-1
    disp(order);
    filename_sursignal = sprintf('sursignal_%03d.mat', order);
    filename_symbolsignal = sprintf('symbolsignal_%03d.mat', order);
    filename_sursignal_2 = sprintf('sursignal_%03d.mat', order+1);
    filename_symbolsignal_2 = sprintf('symbolsignal_%03d.mat', order+1);
    load(filename_sursignal);
    refsignal1 = sursignal;
    load(filename_sursignal_2);
    refsignal2 = sursignal;
    load(filename_symbolsignal);
    sursignal1 = symbolsignal;
    load(filename_symbolsignal_2);
    sursignal2 = symbolsignal;
    n_total = length(refsignal1);
    ref_adjusted1 = refsignal1.';
    ref_adjusted2 = refsignal2.';
    sur_adjusted1 = sursignal1;
    sur_adjusted2 = sursignal2;

    for i = 1:10
        sur_corr1 = xcorr(sur_adjusted1,beacon_signal);
        sur_corr2 = xcorr(sur_adjusted2,beacon_signal);
        ref_corr1 = xcorr(ref_adjusted1,beacon_signal);
        ref_corr2 = xcorr(ref_adjusted2,beacon_signal);
        [max_corr_sur1, max_corr_index_sur1] = maxk(sur_corr1,1);
        [max_corr_sur2, max_corr_index_sur2] = maxk(sur_corr2,1);
        [max_corr_ref1, max_corr_index_ref1] = maxk(ref_corr1,1);
        [max_corr_ref2, max_corr_index_ref2] = maxk(ref_corr2,1);

        original_position_sur1 = max_corr_index_sur1 - n_total + 1;
        original_position_sur2 = max_corr_index_sur2 - n_total + 1;
        original_position_ref1 = max_corr_index_ref1 - n_total + 1;
        original_position_ref2 = max_corr_index_ref2 - n_total + 1;

        sur_adjusted1(original_position_sur1:original_position_sur1+length(beacon_signal)-1) = 0;
        sur_adjusted2(original_position_sur2:original_position_sur2+length(beacon_signal)-1) = 0;
        ref_adjusted1(original_position_ref1:original_position_ref1+length(beacon_signal)-1) = 0;
        ref_adjusted2(original_position_ref2:original_position_ref2+length(beacon_signal)-1) = 0;
    end

    P = 8;
    V1 = zeros(n_total,P);
    V2 = zeros(n_total,P);
    for p = 1:P
        M1 = [zeros(1,p-1),ref_adjusted1];
        M2 = [zeros(1,p-1),ref_adjusted2];
        V1(:,p) = M1(1:n_total);
        V2(:,p) = M2(1:n_total);
    end
    W1 = (V1'*V1)\(V1'*sur_adjusted1);
    W2 = (V2'*V2)\(V2'*sur_adjusted2);
    U1 = sur_adjusted1-V1*W1;
    U2 = sur_adjusted2-V2*W2;
    sur_adjusted1 = U1.';
    sur_adjusted2 = U2.';

    lamda = 3/24;
    CIT = 0.5;
    FFT = zeros(12,n_total);

    n = 1:length(sur_adjusted2);
    for segement = 0:segement_number-1
        refsignal_s = [ref_adjusted1(segement*delta_t*fs+1:n_total),ref_adjusted2(1:segement*delta_t*fs)];
        sursignal_s = [sur_adjusted1(segement*delta_t*fs+1:n_total),sur_adjusted2(1:segement*delta_t*fs)];
        for n_tau = 0:11
            ref = [zeros(1,n_tau),refsignal_s];
            ref = ref(1:n_total);
            FFT(n_tau+1,:) = abs(fftshift(fft(conj(ref).*sursignal_s)));
    %         FFT(n_tau+1,:) = IAA(conj(ref).*sportsignal_2);
        end
        COR=FFT(:,(floor(-40+n_total/2+1):floor(40+n_total/2+1)));
        
    %     [maxValue, maxIndex] = max(COR(1,:));
        [maxValue, maxIndex] = max(COR(:));
        [row, col] = ind2sub(size(COR), maxIndex);
        temp = COR(1,:);
    %     temp(temp<0.8*maxValue)=0;
        doppler_time((order-1)*segement_number+segement+1,:) = temp;
    end
end
surf(delta_t*(1:(signal_number-1)*segement_number),lamda*fd,abs(doppler_time.'),'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');