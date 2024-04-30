function out = IAA(signal)
Ml = length(signal);
fs = 60e6;
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
        s(k) = (((A(:,k)'/R)*signal.'))/(((A(:,k)'/R)*A(:,k)));
        P(k,k) = abs(s(k))^2;
    end
    R = A*P*A';
    diff = trace(abs(P-temp));
    temp = P;
end
for o = 1:K
    out(o) = sqrt(P(o,o));
end
end