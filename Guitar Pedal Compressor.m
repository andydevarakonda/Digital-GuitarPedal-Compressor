
﻿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital Guitar Pedal Compressor Filter
% Name: Anirudh Devarakonda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global freq_filter,
%global norm_filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import = importdata('');%paste txt file of any sound data 
% fs = 48000;                             % sampling freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play sound of raw audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%soundsc(y,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts = 1/fs;                              % sampling period      
Length_y = length(import(:,1));              % length of signal
time = (0:Length_y-1)*ts;               % time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph with DFT/FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = fft(import);                             % Discrete Fourier transform
F1 = abs(Y/Length_y);                   % frequency
F2 = F1(1:Length_y/2+1);                % half of frequency
F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
f = fs*(0:(Length_y/2))/Length_y;       % freq vector [Hz]
f_kHz = f/1000;                         % freq vector [kHz]

figure(1) ;                             
subplot(1,2,1)
plot(time,import)
axis tight
title('y(t) vs t (original)')           % label for the graph title
xlabel('t (s)')                         % label for x axis
ylabel('y(t)')                          % label for y axis
subplot(1,2,2)
plot(f_kHz,F2) 
axis([0 5  0 0.06])
title('Y(F) vs F ')                           % label for the graph title                 
xlabel('F (kHz)')                              % label for x axis                     
ylabel('Y(F)')                                 % label for y axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D major chord frequencies for notes D3, A3, D4, F#4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D3 = 146.83;                            % freq of note D3 [Hz]
D3_int = round(D3/max(f)*length(f))     % associated integer to above freq
A3 = 220.00;                            % freq of note A3 [Hz]
A3_int = round(A3/max(f)*length(f))     % associated integer to above freq
D4 = 293.66;                            % freq of note D4 [Hz]
D4_int = round(D4/max(f)*length(f))     % associated integer to above freq
F_sharp_4 = 369.99;                     % freq of note F#4 [Hz]
F_sharp_4_int = round(F_sharp_4/max(f)*length(f))   % associated integer to above freq

note_freq = [D3 A3 D4 F_sharp_4];       % vector of all note freqs
note_freq_int = [D3_int A3_int D4_int F_sharp_4_int]; % vector of all int note freqs


    delta = 3 %cutoff frequency paramter
    order = 6 %order of the filter
    ripple = 1 %ripple frequency
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this for loop is the filter bank and filters the frequencies based on the
%note_freq frequencies
   
for z = 1:length(note_freq)
    c_freq = round(note_freq(z));
    [x1, y1] = cheby1(order, ripple, (c_freq-delta)/(fs/2), 'High');%high pass
    [x2, y2] = cheby1(order, ripple, (c_freq+delta)/(fs/2), 'Low');%low pass
    if(z == 1)                                      
        D3_Z = filter(x2, y2, filter(x1, y1, import)); %creates a bandpass  
    end
    if(z == 2)                                      
        A3_Z= filter(x2, y2, filter(x1, y1, import));
         
    end
    if(z == 3)                                      
        D4_Z = filter(x2, y2, filter(x1, y1, import));
        
    end
    if(z == 4)                                     
        FS4_Z = filter(x2, y2, filter(x1, y1, import));
        
    end
   
   
end
   
filter_bank = [D3_Z, A3_Z,D4_Z,FS4_Z] %creates an array of filtered signals
total_filter = D3_Z + A3_Z+ D4_Z+ FS4_Z %adds the parallel bandpass filter output
                                        %in the filter bank to
                                        %provide the reconstructed
                                        %signal that which has the required
                                        %isolated frequencies
                                        
%Taking the fourier transform of the filtered signals

D3_ZF = fft(D3_Z);                  %D3 Filtered signal FT
D3F1 = abs(D3_ZF/Length_y);
D3F2 = D3F1(1: Length_y/2+1);
D3F2(2:end-1) = 2*D3F2(2:end-1);

A3_ZF = fft(A3_Z);                   %A3 Filtered signal FT
A3F1 = abs(A3_ZF/Length_y);
A3F2 = A3F1(1: Length_y/2+1);
A3F2(2:end-1) = 2*A3F2(2:end-1);

D4_ZF = fft(D4_Z);                   %D4 Filtered signal FT
D4F1 = abs(D4_ZF/Length_y);
D4F2 = D4F1(1: Length_y/2+1);
D4F2(2:end-1) = 2*D4F2(2:end-1);

FS4_ZF = fft(FS4_Z);                %FS4 Filtered signal FT
FS4F1 = abs(FS4_ZF/Length_y);
FS4F2 = FS4F1(1: Length_y/2+1);
FS4F2(2:end-1) = 2*FS4F2(2:end-1);


%Scale and normalize the fourier transformed signals and combine them

max_amplitude = [max(D3F2), max(A3F2), max(D4F2), max(FS4F2)] % calculate and record the peak amplitude of the fft of each filtered signal

norm_factor = max(max_amplitude)./max_amplitude %scaling factor for the normalization of the filtered signals

D3F2 = D3F2.*norm_factor(1) %normalizing the signals one by one
A3F2 = A3F2.*norm_factor(2)
D4F2 = D4F2.*norm_factor(3)
FS4F2 = FS4F2.*norm_factor(4)

total_filter_fft = D3F2 + A3F2 + D4F2 + FS4F2
total_filter_norm = D3_Z*norm_factor(1)+A3_Z*norm_factor(2)+D4_Z*norm_factor(3)+FS4_Z*norm_factor(4) %combined normalized fft of for the 4 isolated frequencies 

subplot(4,1,1);    %plot the frequency response of each of the filtered signals according to the 4 frequencies
[h,w] = freqz(D3_Z);%Frequency response of D3 filtered signal
plot(w/pi,20*log10(abs(h)))

ylabel('Magnitude (dB)')
subplot(4,1,2);
[h1,w1] = freqz(A3_Z); %Frequency response of A3 filtered signal
plot(w1/pi,20*log10(abs(h1)))

ylabel('Magnitude (dB)')
subplot(4,1,3);
[h2,w2] = freqz(D4_Z);%Frequency response of D4 filtered signal
plot(w2/pi,20*log10(abs(h2)))

ylabel('Magnitude (dB)')
subplot(4,1,4);
[h3,w3] = freqz(FS4_Z);%Frequency response of F#4 filtered signal
plot(w3/pi,20*log10(abs(h3)))
ax = gca;
ax.YLim = [-100 20];
ax.XTick = 0:.5:2;
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
title('Frequency Response vs F ')

%normalized fft amplitude vs F frequency signal
figure(2);
hold on
plot(f_kHz,D3F2) %Plot the FFT filtered signals and combine them into one plot
plot(f_kHz,A3F2)
plot(f_kHz,D4F2)
plot(f_kHz,FS4F2)
axis([0 0.5 0 0.05])
title('Normazlized and Filtered FFT vs F')                           
xlabel('F (kHz)')                       
ylabel('Y(F)')
hold off

%filtered and normalized time signal
figure(3);
subplot(2,1,1);
plot(time, total_filter);
axis([0 max(time) -.65 .65]);
title('y(t) vs t (Filtered)')                           
xlabel('t (s)')                       
ylabel('y(t)')
subplot(2,1,2);
plot(time, total_filter_norm);
axis([0 max(time) -2.2 2.2]);
title('y(t) vs t (Normalized)')                           
xlabel('t (s)')                       
ylabel('y(t)')

soundsc(total_filter_norm,fs)
