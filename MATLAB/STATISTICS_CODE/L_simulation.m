% Runs the whole L-band receiver pipeline. If desired,
%   - Plots the each step
%   - Saves each plot
%   - Allows to choose between two different filters
%   - Allows to choose whether to apply the 2nd filter, which is redundant


clear all; close all;

display = true;     %display = true - all figures are displayed
                    %display = false - no figure is generated

save = true;    %save = true - images are saved
                %save = false - images are not saved

filter = 1; %filter = 1 for FFT&IFFT            
            %filter = 2 for designfilt
            
LP2applied = true; %LP2applied = true - the RSSI output filter is applied
                    %LP2applied = false - the RSSI output filter is not applied
%%
fs_adc = 3200;
N_adc = 2000;
n_adc = 1:N_adc;
fact_decimate = 25000;
fs = fs_adc * fact_decimate;
N = N_adc * fact_decimate;
n = 1:N;
Nfft = N;

fpass1 = 2.5e6; %2.5 MHz for the L-band
                %10 MHz for the K-band
fpass2 = 6e6;
fpass3 = 1600;
if(filter==2)
    d1 = designfilt('lowpassfir', 'FilterOrder', 200, 'CutoffFrequency', fpass1, 'SampleRate', fs);
    d2 = designfilt('lowpassfir', 'FilterOrder', 150, 'CutoffFrequency', fpass2, 'SampleRate', fs);
    d3 = designfilt('lowpassfir', 'FilterOrder', 20000, 'CutoffFrequency', fpass3, 'SampleRate', fs);
end
%fvtool(d2)
%%

tic

INR = -5;
DC = 0.4;
A0 = 1;
noise = A0*normrnd(0,1,[1,N]);
%A1 = sqrt((2*10.^(INR/10))/DC);
A1 = sqrt(2*10^(INR/10));
fo = 2e6;
theta = 0;
Mmin = 10;
Mmax = 50;
%M = N/(Mmin + (Mmax-Mmin)*rand());
M = N/10;
offset = rand()*M;
i1 = A1*cos(2*pi*(fo/fs)*n + theta).*pulsetrain(n,M,DC,offset);

%% Figure 1: signal + noise
x = i1 + noise;

if(display)
    Nplot = 5e6;
    
    fig = figure('Name', 'signal + noise', 'Position', [10 250 1500 500]);
    
    subplot(1,2,1)
    plot(n(1:Nplot)/fs*1e3, x(1:Nplot));
    title(['Time domain (', num2str(Nplot), ' samples out of ', num2str(N),')'])
    xlabel('Time [ms]')
    ylabel('Amplitude')
    xlim([0 Nplot/fs*1e3])
    grid on
    grid minor
    
    subplot(1,2,2)
    plot((1:Nfft)/Nfft*fs/1e6, log10(abs(fft(x,Nfft))))
    title(['Frequency domain (', num2str(Nfft), ' samples)'])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude (logarithmic scale)')
    grid on
    grid minor
    
    sgtitle("signal + noise")
    if(save)
        saveas(fig, 'signal+noise', 'png');
    end
end
%% Figure 2: lowpass(signal + noise)

if(filter==1)
    y = lowpass_FFT(x, fpass1, fs, Nfft);
elseif(filter==2)
    y = filtfilt(d1, x);
end

if(display)
    Nplot = 5e6;

    fig = figure('Name', 'lowpass(signal + noise)', 'Position', [10 250 1500 500]);
    
    subplot(1,2,1)
    plot(n(1:Nplot)/fs*1e3, y(1:Nplot));
    title(['Time domain (', num2str(Nplot), ' samples out of ', num2str(N),')'])
    xlabel('Time [ms]')
    ylabel('Amplitude')
    xlim([0 Nplot/fs*1e3])
    grid on
    grid minor
    
    subplot(1,2,2)
    plot((1:Nfft)/Nfft*fs/1e6, log10(abs(fft(y,Nfft))))
    title(['Frequency domain (', num2str(Nfft), ' samples)'])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude (logarithmic scale)')
    grid on
    grid minor
    
    sgtitle("lowpass(signal + noise)")
    if(save)
        saveas(fig, 'lowpass(signal+noise)', 'png');
    end
end
%% Figure 3: (lowpass(signal + noise))^2

z = y.^2;

if(display)
    Nplot = 5e6;

    fig = figure('Name', '(lowpass(signal + noise))^2', 'Position', [10 250 1500 500]);
    
    subplot(1,2,1)
    title(['Time domain (first ', Nplot, ' samples']);
    plot(n(1:Nplot)/fs*1e3, z(1:Nplot));
    title(['Time domain (', num2str(Nplot), ' samples out of ', num2str(N),')'])
    xlabel('Time [ms]')
    ylabel('Amplitude')
    xlim([0 Nplot/fs*1e3])
    grid on
    grid minor
    
    subplot(1,2,2)
    plot((1:Nfft)/Nfft*fs/1e6, log10(abs(fft(z,Nfft))))
    title(['Frequency domain (', num2str(Nfft), ' samples)'])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude (logarithmic scale)')
    grid on
    grid minor
    
    sgtitle("(lowpass(signal + noise))^2")
    if(save)
        saveas(fig, '(lowpass(signal+noise))^2', 'png');
    end
end
%% Figure 4: 10*log10((lowpass(signal + noise))^2)

s = 10*log10(z);

if(display)
    Nplot = 5e6;

    fig = figure('Name', '10*log10((lowpass(signal + noise))^2)', 'Position', [10 250 1500 500]);
    
    subplot(1,2,1)
    plot(n(1:Nplot)/fs*1e3, s(1:Nplot));
    title(['Time domain (', num2str(Nplot), ' samples out of ', num2str(N),')'])
    xlabel('Time [ms]')
    ylabel('Amplitude')
    xlim([0 Nplot/fs*1e3])
    grid on
    grid minor
    
    subplot(1,2,2)
    plot((1:Nfft)/Nfft*fs/1e6, log10(abs(fft(s,Nfft))))
    title(['Frequency domain (', num2str(Nfft), ' samples)'])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude (logarithmic scale)')
    grid on
    grid minor
    
    sgtitle("10*log10((lowpass(signal + noise))^2)")
    if(save)
        saveas(fig, '10log10((lowpass(signal + noise))^2)', 'png');
    end
end
%% Figure 5: lowpass(10*log10((lowpass(signal + noise))^2))

if(LP2applied)
    if(filter==1)
        r = lowpass_FFT(s, fpass2, fs, Nfft);
    elseif(filter==2)
        r = filtfilt(d2, s);
    end
else
    r = s;
end

if(display && LP2applied)
    Nplot = 5e6;

    fig = figure('Name', 'lowpass(10*log10((lowpass(signal + noise))^2))', 'Position', [10 250 1500 500]);
    
    subplot(1,2,1)
    plot(n(1:Nplot)/fs*1e3, r(1:Nplot));
    title(['Time domain (', num2str(Nplot), ' samples out of ', num2str(N),')'])
    xlabel('Time [ms]')
    ylabel('Amplitude')
    xlim([0 Nplot/fs*1e3])
    grid on
    grid minor
    
    subplot(1,2,2)
    plot((1:Nfft)/Nfft*fs/1e6, log10(abs(fft(r,Nfft))))
    title(['Frequency domain (', num2str(Nfft), ' samples)'])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude (logarithmic scale)')
    grid on
    grid minor
    
    sgtitle("lowpass(10*log10((lowpass(signal + noise))^2))")
    if(save)
        saveas(fig, 'lowpass(10log10((lowpass(signal + noise))^2))', 'png');
    end
end

%% Figure 6: lowpass(lowpass(10*log10((lowpass(signal + noise))^2)))

if(filter==1)
    t = lowpass_FFT(r,fpass3,fs,Nfft);
elseif(filter==2)
    t = filtfilt(d3, r);
end

if(display)
    Nplot = 5e6;

    fig = figure('Name', 'lowpass(lowpass(10*log10((lowpass(signal + noise))^2)))', 'Position', [10 250 1500 500]);
    
    subplot(1,2,1)
    plot(n(1:Nplot)/fs*1000, t(1:Nplot));
    title(['Time domain (', num2str(Nplot), ' samples out of ', num2str(N),')'])
    xlabel('Time [ms]')
    ylabel('Amplitude')
    xlim([0 Nplot/fs*1000])
    grid on
    grid minor
    
    subplot(1,2,2)
    plot((1:Nfft)/Nfft*fs/1e6, log10(abs(fft(t,Nfft))))
    title(['Frequency domain (', num2str(Nfft), ' samples)'])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude (logarithmic scale)')
    grid on
    grid minor
    
    sgtitle("lowpass(lowpass(10*log10((lowpass(signal + noise))^2)))")
    if(save)
        saveas(fig, 'lowpass(lowpass(10log10((lowpass(signal + noise))^2)))', 'png');
    end
end

%% Figure 7: sampling(lowpass(lowpass(10*log10((lowpass(signal + noise))^2))))

ind = 1:fact_decimate:N;
u = t(ind);

if(display)
    Nplot = 200;
    %Nplot = N_adc;
    fig = figure('Name', 'sampling(lowpass(lowpass(10*log10((lowpass(signal + noise))^2))))', 'Position', [10 250 1500 500]);
    
    subplot(1,2,1)
    plot(n_adc(1:Nplot)/fs_adc*1000, u(1:Nplot));
    title(['Time domain (', num2str(Nplot), ' samples out of ', num2str(N_adc),')'])
    xlabel('Time [ms]')
    ylabel('Amplitude')
    xlim([0 Nplot/fs_adc*1000])
    grid on
    grid minor
    
    subplot(1,2,2)
    plot((1:Nfft)/Nfft*fs_adc, log10(abs(fft(u,Nfft))))
    title(['Frequency domain (', num2str(Nfft), ' samples)'])
    xlim([0 fs_adc])
    xlabel('Frequency [MHz]')
    ylabel('Amplitude (logarithmic scale)')
    grid on
    grid minor
    
    sgtitle("sampling(lowpass(lowpass(10*log10((lowpass(signal + noise))^2))))")
    if(save)
        saveas(fig, 'sampling(lowpass(lowpass(10log10((lowpass(signal + noise))^2))))', 'png');
    end
end

%% Figure 8: Histogram
toc

fig = figure('Name', 'Histogram', 'Position', [10 250 1500 500]);
histogram(real(u),250)
ku = kurtosis(real(u));
sk = skewness(real(u));
dim = [0.02 0 0.3 0.3];
str = {['Kurtosis: ', num2str(ku)],['Skewness: ', num2str(sk)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title("Histogram of sampling(lowpass(10*log10(lowpass((lowpass(signal + noise))^2))))")
grid on
grid minor
if(save)
    saveas(fig, 'hist_sampling(lowpass(lowpass(10log10((lowpass(signal + noise))^2))))', 'png');
end