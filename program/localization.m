% load("signal2.mat");
% X0=rxmimo2x2;% the received signal
% load("signal3.mat")
% susrp=rxmimo2x2;
% load("signal1.mat")
% X1=rxmimo2x2;
% X=X0-susrp;
% % X = sqrt(sum(abs(X).^2));
% X_tr=X1-susrp;
% sum(abs(X_tr).^2)
% sum(abs(X).^2)
% 10*log10(sum(abs(X_tr).^2)/sum(abs(X).^2))
% fs = 56e6;
% rs=fftshift(abs(fft(X)));
% f2 = linspace(0,fs,length(rs))-fs/2;
% plot(f2,rs)
% The number of the grid 
r = 800; % row
c = 800; % column
L = 6; % the length of the room
W = 4; % the width of the room
long = L/r; % the length of each grid
width = W/c; % the width of each grid
% The transmitted power in decibel
Pt = 0.687690664741878;
gama = 2;
PathLoss0 = 16.081463599820760;
d0 = 1.85;
X = [sqrt(0.0075),sqrt(0.0132),sqrt(0.0089),sqrt(0.00751),sqrt(0.00708),sqrt(0.00638)];
coord_x = 0.01*[21,21+122,21+122+78,21+122+78+88,21+122+78+88+39,21+122+78+88+39+63];
% The number of samples
sum_error = zeros(c,r);
fitting_error2 = zeros(c,r);
for li = 1:r
    for ci = 1:c
        d = sqrt((ci*width-coord_x).^2+((li*long+0.8))^2);
        Pr = 10*log10(Pt)-PathLoss0-10*gama*log10(d/d0);
        Yk = sqrt(10.^(Pr/10));
        ak = linspace(0.05,1.1,1000);
        fitting_error1 = zeros(1,length(ak));
        for j = 1:length(ak)
            fitting_error1(j) = 1/6*sum(abs(X-ak(j).*Yk).^2);
        end
        [m_a,k_a] = min(fitting_error1);
        a = ak(k_a);
        fitting_error2(li,ci) = 1/6*sum(abs(X-a*Yk).^2);
    end
end
step_x = width*(1:c);
step_y = long*(1:r);
[minValue, minIndex] = min(fitting_error2(:));
[row, col] = ind2sub(size(fitting_error2), minIndex);
figure;
mesh(step_x,step_y,fitting_error2,'CDataMapping','scaled');
colorbar
set(gca,'YDir','normal');
hold on
scatter3(1.7,d0,0.001,120,'x', 'y', 'LineWidth', 3)
hold on
scatter3(row*long,col*width,0.001,120,'o', 'y', 'LineWidth', 3)
hold off
