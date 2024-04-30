function a = doppler_steering_vector(M,f,fs)
    m = 0:M-1;
    a = exp(1j*2*pi*f/fs*m);
end