%% BMEG 321 - Lab 3
% Author: Lauryn Peters
% Date: February 02, 2024

clear all; clc

%% Loading data

bioRadioData_raw = readtable('EMG_raw.csv');
bioRadioData_filt = readtable('EMG_filtered.csv');

%% Processing Signals in the Time Domain

% removing DC offset
dcOffset_filt = bioRadioData_filt-mean(bioRadioData_filt.EMGBICEP);

fs = 2000; % sample rate
[b,a] = butter(4,[10 450]/(fs/2)); % define Butterworth filter parameters

% full wave rectification
data_filt = filtfilt(b,a, dcOffset_filt.EMGBICEP);
signal_rect = abs(data_filt);

% low pass filter - envelope
fpass = 2.8; 
[d,c] = butter(4,fpass/(fs/2)); % define Butterworth filter parameters
signal_env = filtfilt(d,c, signal_rect);

%% Plotting

% convert raw data to mV because raw data output in V while filtered data 
% output in mV

bioRadioData_raw_mV = bioRadioData_raw{:,2}.*1000;
time = (0:inv(fs):54.7210)'; % time vector 

figure(1)
plot(time, bioRadioData_raw_mV) % raw signal
title("EMG Signals")
xlabel("Time (s)")
ylabel("Voltage (mV)")

hold on
plot(time, bioRadioData_filt{:,2}) % filtered signal
hold on 
plot(time, signal_rect) % processed signal
hold on
plot(time, signal_env) % enveloped signal
hold off

legend('Raw EMG','Filtered EMG', 'Processed EMG','EMG envelope')

%% FFTs of Signals to Isolate Noise

% specify params 
T = 1/fs; % period of signal 
L = 109443; % length of signal vector  

% raw signal
ft_raw = fft(bioRadioData_raw_mV); % fourier transform of raw data 

% convert to single sided spectrum
P2_raw = abs(ft_raw/L);
P1_raw = P2_raw(1:L/2+1);
P1_raw(2:end-1) = 2*P1_raw(2:end-1);

% filtered signal
ft_filt = fft(bioRadioData_filt{:,2}); % fourier transform of filtered data 

% convert to single sided spectrum
P2_filt = abs(ft_filt/L);
P1_filt = P2_filt(1:L/2+1);
P1_filt(2:end-1) = 2*P1_filt(2:end-1);

% processed signal
ft_proc = fft(signal_rect); % fourier transform of processed data 

% convert to single sided spectrum
P2_proc = abs(ft_proc/L);
P1_proc = P2_proc(1:L/2+1);
P1_proc(2:end-1) = 2*P1_proc(2:end-1);

% noisy signal
ft_noise = fft(bioRadioData_raw_mV - signal_rect); % fourier transform of noise 

% convert to single sided spectrum
P2_noise = abs(ft_noise/L);
P1_noise = P2_noise(1:L/2+1);
P1_noise(2:end-1) = 2*P1_noise(2:end-1);

% plotting tiles
f = fs/L*(0:(L/2));

figure(2)
tiledlayout(2,2);
sgtitle("Single-Sided Amplitude Spectrums");

% Tile 1
nexttile
plot(f,P1_raw) 
title("Raw Signal")
xlabel("f (Hz)")
ylabel("Amplitude")
ylim([0 0.03])

% Tile 2
nexttile
plot(f,P1_filt) 
title("Filtered Signal")
xlabel("f (Hz)")
ylabel("Amplitude")
ylim([0 0.03])

% Tile 3
nexttile
plot(f,P1_proc) 
title("Processed Signal")
xlabel("f (Hz)")
ylabel("Amplitude")
ylim([0 0.03])

% Tile 4
nexttile
plot(f,P1_noise) 
title("Noise")
xlabel("f (Hz)")
ylabel("Amplitude")
ylim([0 0.03])

%% Contraction EMG Values

idx1 = find(time > 4.1815 & time < 8.309); % indices of first contraction
contrac1 = cumtrapz(signal_env(idx1)); % integration vector 
emg1 = contrac1(end); % final estimated EMG value 

idx2 = find(time > 11.7395 & time < 16.706); % indices of second contraction
contrac2 = cumtrapz(signal_env(idx2)); % integration vector 
emg2 = contrac2(end); % final estimated EMG value 

idx3 = find(time > 21.572 & time < 27.9925); % indices of third contraction
contrac3 = cumtrapz(signal_env(idx3)); % integration vector 
emg3 = contrac3(end); % final estimated EMG value 

idx4 = find(time > 31.7205 & time < 37.8485); % indices of fourth contraction
contrac4 = cumtrapz(signal_env(idx4)); % integration vector 
emg4 = contrac4(end); % final estimated EMG value 

idx5 = find(time > 41.2575 & time < 47.3865); % indices of fifth contraction
contrac5 = cumtrapz(signal_env(idx5)); % integration vector 
emg5 = contrac5(end); % final estimated EMG value 
