.clc;
clear;
close;

%% Discover Radio
connectedRadios = findsdru;
if strncmp(connectedRadios(1).Status, 'Success', 7)
    radioFound = true;
    platform = connectedRadios(1).Platform;
    switch connectedRadios(1).Platform
        case {'B200','B210'}
            address = connectedRadios(1).SerialNum;
        case {'N200/N210/USRP2','X300','X310'}
            address = connectedRadios(1).IPAddress;
    end
else
    radioFound = false;
    address = '192.168.1.4';
    platform = 'N200/N210/USRP2';
end
a
if(radioFound == true)
    disp('discover radio!');
else
    error('no discover radio!');
end

SYS.Platform = platform;
SYS.Address  = address;

%% USRP transmitter parameters
receive_time = 0.5;
SYS.MasterClockRate        = 60e6;% Hz
SYS.USRPCenterFrequency    = 2.44e9;
SYS.USRPGain               = 20;
SYS.USRPFrontEndSampleRate = 60e6;
SYS.USRPDecimationFactor   = SYS.MasterClockRate/SYS.USRPFrontEndSampleRate;
SYS.USRPFrameLength        = 45000;
% Experiment Parameters
SYS.numRxFrame = ceil(receive_time*SYS.MasterClockRate/SYS.USRPFrameLength);

%% USRP initation
    radio = comm.SDRuReceiver(...
            'Platform',             SYS.Platform, ...
            'SerialNum',            SYS.Address, ...
            'MasterClockRate',      SYS.MasterClockRate, ...
            'CenterFrequency',      SYS.USRPCenterFrequency, ...
            'Gain',                 SYS.USRPGain, ...
            'DecimationFactor',     SYS.USRPDecimationFactor, ...
            'SamplesPerFrame',      SYS.USRPFrameLength, ...
            'OutputDataType',       'double', ...
            'IPAddress',            '192.168.1.4');

radio.ChannelMapping = 1;     % Use both TX channels
radio.OverrunOutputPort = true;

%% Initialize variables
len = uint32(0);
rcvdSignal = complex(zeros(SYS.USRPFrameLength,2));

disp(SYS);
disp('数据采集进行中...');

hlog = dsp.SignalSink;
for idx = 1:SYS.numRxFrame
    % Keep accessing the SDRu System object output until it is valid
    while len <= 0
        [rcvdSignal, len] = step(radio);
    end
    % When the SDRu System object output is valid, decode the received message
    if (len > 0)
        hlog(rcvdSignal);
    end
    % update
    len = uint32(0);
end
disp('数据采集结束！');
fs = 60e6;
sursignal4 = hlog.Buffer;
n = 1:length(sursignal4);
rs=fftshift(20*log10(abs(fft(sursignal4))));
f2 = linspace(0,fs,length(rs))-fs/2;
subplot(2,1,1)
plot(f2,rs)
subplot(2,1,2)
plot(n,sursignal4);
save('sursignal4','sursignal4')       
    % save message and decode
