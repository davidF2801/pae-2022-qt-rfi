clear all, close all

% Set the center frequency and bandwidth of the signal
%For L band
center_frequency = 2.5e6; % 20 MHz
bandwidth = 5e6; % 20 MHz

%For K band
% center_frequency = 20e6; 
% bandwidth = 20e6;
% Generate the base signal (white noise)


total_bandwidth = 1e9;
duration = 0.7;
%For L band
fpass11 = 0e6;
fpass12 = 5e6;
%For K band
% fpass11 = 10e6;
% fpass12 = 30e6;
fs_adc = 3200;
N_adc = 2000;
fact_decimate = fix(2*(bandwidth+center_frequency)/fs_adc);
fs = fs_adc * fact_decimate+1;
t = 0:1/fs:duration;
N = N_adc * fact_decimate;
n_vec = 1:N;
Nfft = fix(N);
ind = 1:fact_decimate:N;
fpass1 = 2.5e6;
fpass2 = 6e6;
fpass3 = 1600;
NBINS = total_bandwidth/bandwidth;
nb = 1:NBINS;
% Iterate 50 times
power = zeros(1,NBINS);
INRNUM = 6;
INR = [-30 -20 -10 -5 0 5];
for j = 1:INRNUM
    INRref = INR(j);
    rfi_location = zeros(1,NBINS);
    rfi_amplitude = rand()*sqrt(2*10^(INRref/10));
    
    for i = 1:NBINS
       disp(i)
       % Generate RFI signal
       base_signal = randn(size(t));
       rfi_frequency = rand() * bandwidth + (center_frequency - bandwidth/2); % Random frequency within bandwidth
       rfi_signal = rfi_amplitude*cos(2*pi*rfi_frequency*real(t));
       t1 = t(ind);
       % Add RFI to base signal with probability 0.1
       if rand() >= 0.7
           signal = base_signal + rfi_signal;
           rfi_location(i) = 1;
           signal1 = signal(ind);
       else
           signal = base_signal;
       end
       
       
     
       %y = bandpass(real(signal), [fpass11 fpass12], fs); %For K band
       y = lowpass_MXF(real(signal),fpass12,fs,Nfft); %For L band
       z = y.^2;
       s = 10*log10(z);
       r = lowpass_MXF(s, fpass2, fs, Nfft);
       t = lowpass_MXF(r,fpass3,fs,Nfft);
       u = t(ind);
       power(i) = mean(u);
       
    
       % Do something with the signal here
       % ...
    end
    
    th = zeros(1,NBINS);
    pd = zeros(1,NBINS);
    pd(1) = 0;
    sigma = std(power);
    mn = mean (power);
    beta = 1;
    ndet = 0;
    nfa = 0;
    for i = 1:NBINS
        th(i) = beta*sigma + mn;
        if power(i) < th(i)
            pd(i) = 0; % no pulse detected
        else
            pd(i) = 1; % pulse detected
            if rfi_location(i) == 1
                ndet = ndet +1;
            else
                nfa = nfa +1;
            end
            
        end
    end
    pdet = ndet/sum(rfi_location == 1);
    pfa = nfa/sum(rfi_location == 0);
    
    disp(pdet)
    disp(pfa)
    fig1 = figure(1);
    subplot(211)
    bar(nb,power,'b')
    hold on
    plot(nb, th, 'r', 'LineWidth', 2)
    grid
    hold off
    xlabel('Bin Index'), ylabel('Power (dBm)')
    title(['Signal power (blue), Threshold = \beta\sigma (red), \beta^2=' num2str(beta)])
    title_string = sprintf('Signal power (blue), Threshold = \\beta\\sigma (red), pdet=%.2f,pfa=%.2f, INRref=%.2f', pdet,pfa, INRref);
    title(title_string);
    subplot(212),bar(nb,pd),grid
    xlabel('Sample Index'), ylabel('Blank on/off')
    title('Pulse Detection')
    path = 'C:\Users\dforn\Documents\PAE MATLAB\Matlab scripts\Matlab scripts\FB\';
    filename  = strcat(path,'FBpng',num2str(j));
    saveas(fig1,filename,'png');
    filename  = strcat(path,'FBfig',num2str(j));
    saveas(fig1,filename,'fig');

end