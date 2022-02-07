clc

close all

clear all

#-------------------Read in Data------------------------------

Data = load('student.mat');



output = Data.cn;
input  = Data.r;
t      = Data.t;

#-------------------Plot Data------------------------------

original_data = figure;
plot(t,input);
grid on;
xlabel('Time ')
ylabel ('Magnitude')


hold on

plot(t,output);

#-------------------Fitted model------------------------------

zeta        = 0.15;
natural_f   = 5;
theta       = acos(zeta);
damp_f      = natural_f*sqrt(1-zeta^2);

model       = (1 - (1/sqrt(1-zeta^2))*exp(-zeta*natural_f*t).*sin(damp_f*t + theta));

hold on

plot (t, model, 'g');


heading = strcat('Fitted Model Data ');
title(heading);

legend('Input', 'Output', 'Fitted');

hold off

#-------------------Fitted model Transfer Function------------------------------
K           = 1;
numerator   = K* natural_f;
denominator = [1 , 2*zeta*natural_f ,natural_f^2 ];

sys         = tf(numerator , denominator)

#-------------------Frequency Spectra------------------------------------

L                       =   length(output);

frequency_spectra       =   figure;

Measured                =   fft(output);

P2 = abs(Measured/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = 0.01*(0:(L/2))/L;


plot(f , P1, 'r');


hold on


Model  = fft(model);
P2     = abs(Model/L);
P1     = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = (0.01)*(0:(L/2))/L;

plot(f , P1, 'g');

grid on
heading = strcat('Frequency Power Spectra');
title(heading);
legend('Measured', 'Model');
xlabel('Frequency ')
ylabel ('Magnitude')

hold off

#------------------Naive Approach----------------------------



#------Expected Noise (Unknown vs Model) - Time Domain Analysis - Plot----------

noise_signal = figure;

#difference of signals
noise = (output - model);

plot(t,noise);
grid on;
xlabel('Time ')
ylabel ('Magnitude')

heading = strcat('Noise of Data in Time Domain ');

title(heading);

#-------------------Expected Noise Power - Time Domain Analysis - Calculation-----------

L=length(noise);


rms = 1000*std (noise);



power_of_noise = 20*log10(rms/1000);

#power_of_noise  = 20*log10((noise)^2)/L;
#power_of_noise  = std(noise);

sprintf('Power of noise signal Time Domain analysis %f',power_of_noise)


#-------------------Expected Noise Power - Frequency Domain Analysis---------
noise_signal = figure;

X   = fft(noise);
Px  = mag2db(sum(X.*conj(X))/(L^2));

stem(Px);


sprintf('Power of noise signal Frequency Domain Analysis %f',Px)

heading = strcat('Power of Noise in the Unknown System ');
title(heading);
xlabel('Frequency (rad/sec) ')
ylabel ('Magnitude')
grid on


#-------------------SNR-Calculation--------
#
#SNR = 20*log10(rms(model)/rms(noise));
#
#sprintf('Signal to Noise Ratio %f',SNR)




#----------------Expected Noise Power - Using Model Only-------------------------------

fs = 1000;                 % sampling frequency, Hz
T = 6;                     % signal duraton, s
N = round(T*fs);           % number of signal samples



wlen = round(N/100);
win = hamming(wlen, 'periodic');
noverlap = 0;
nfft = 2*wlen;








#remove DC
noise = model - mean(output);



Noise = figure;



plot(t, noise, 'r')
grid on
hold on



A_weighted_noise = filterC(model, fs);

plot(t, A_weighted_noise, 'b')


hold on



B_weighted_noise = filterA(model, fs);

plot(t, B_weighted_noise, 'g')

xlabel('Time(s)')
ylabel('Amplitude')

title('Noise signal in the time domain')

legend('Original signal', 'A - Weighted signal','B - Weighted signal' )


L   = length(output);

power_of_noise  = (norm(noise)^2)/L;
sprintf('Power of noise signal unweighted Time Domain analysis %f',power_of_noise)


power_of_noise =(norm(A_weighted_noise)^2)/L;
sprintf('Power of noise signal A-weighted Time Domain analysis %f',power_of_noise)



power_of_noise =(norm(B_weighted_noise)^2)/L;
sprintf('Power of noise signal B-weighted Time Domain analysis %f',power_of_noise)


#----------------Power Spectral Density-------------------------------



#--------Average Power Spectral Density before  weighting-------


[PSD, f] = pwelch(noise, win, noverlap, nfft, fs, 'onesided');
PSDdB = 20*log10(PSD);


Power_Spectra_Noise = figure;

plot(f, PSDdB, 'k')


hold on
grid on



power_of_noise = mean(PSDdB);
sprintf('Average Power Spectral Density before the weighting %f',power_of_noise)


#-------Average Power Spectral Density after  weighting----------

[PSD, f] = pwelch(A_weighted_noise, win, noverlap, nfft, fs, 'onesided');
PSDdB = 20*log10(PSD);
plot(f, PSDdB, 'r')
hold on

power_of_noise = mean(PSDdB);
sprintf('Average Power Spectral Density after A-weighting %f',power_of_noise)




[PSD, f] = pwelch(B_weighted_noise, win, noverlap, nfft, fs, 'onesided');
PSDdB = 20*log10(PSD);

plot(f, PSDdB, 'g')



power_of_noise = mean(PSDdB);
sprintf('Average Power Spectral Density after B-weighting %f',power_of_noise)


xlabel('Frequency')
ylabel('Magnitude')
title('Power Spectral Density of the signal')
legend('Unweighted', 'A-weighted', 'B-weighted');


hold off

#----------------------SNR-----------------------


[PSD, f] = pwelch(model, win, noverlap, nfft, fs, 'onesided');
PSDdB = 20*log10(PSD);


vrms = 1000*std (noise);



true_noise = 20*log10(vrms/1000);
sprintf('True Noise %f',true_noise)

SNR = (PSDdB/true_noise);

plot (f,SNR);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
title('Signal to Noise Ration')
#legend('Unweighted', 'A-weighted', 'B-weighted');

#sprintf('Signal to Noise Ratio %f',SNR)
