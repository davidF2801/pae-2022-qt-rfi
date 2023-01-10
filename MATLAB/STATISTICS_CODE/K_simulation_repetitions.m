% Runs the whole K-band receiver pipeline over R repetitions.
% Instead of plotting only with one particular INR or DC, R repetitions can be made on various INR-DC pairs defined in the corresponding vectors (designed like this
% to have all desired results at once)
% Results are also saved


fs_adc = 3200;
N_adc = 2000;
fact_decimate = 50000;
fs = fs_adc * fact_decimate;
N = N_adc * fact_decimate;
n = 1:N;
Nfft = N;

fpass11 = 20e6;
fpass12 = 40e6;                
fpass2 = 6e6;
fpass3 = 1600;

INR = [-21 -21];
DC = [0.65 0.79];

f0min = fpass11;
f0max = fpass12;

Mmin = 10; 
Mmax = 50;

R = 1000;
ku = ones(1,R);
sk = ones(1,R);

ind = 1:fact_decimate:N;

for i = 1:length(INR)
    for r=1:R
        disp(strcat("DC = ", num2str(DC(i)), "   INR = ", num2str(INR(i)), "   R = ", num2str(r)))
        A0 = 1;
        noise = A0*normrnd(0,1,[1,N]);
        f0 = f0min + (f0max-f0min)*rand();
        theta = 2*pi*rand();
        M = N/(Mmin + (Mmax-Mmin)*rand());
        offset = rand()*M;
        t = lowpass_FFT(10*log10(bandpass((sqrt(2*10^(INR(i)/10))*cos(2*pi*(f0/fs)*n + theta).*pulsetrain(n,M,DC(i),offset)) + noise, [fpass11 fpass12], fs).^2),fpass3,fs,Nfft);
        u = t(ind);
        ku(r) = kurtosis(real(u));
        sk(r) = skewness(real(u));
    end
    ku_mean = mean(ku);
    sk_mean = mean(sk);
    str = strcat("K_rep_save_DC",num2str(DC(i)*100),"_INR",num2str(abs(INR(i))));
    save(str, 'ku_mean', 'sk_mean', 'ku', 'sk', 'DC', 'INR', 'i', 'R','N_adc','fs_adc','fact_decimate');
end
