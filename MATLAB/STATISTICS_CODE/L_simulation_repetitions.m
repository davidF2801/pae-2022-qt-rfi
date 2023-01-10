% Runs the whole L-band receiver pipeline over R repetitions.
% Instead of plotting only with one particular INR or DC, R repetitions can be made on various INR-DC pairs defined in the corresponding vectors (designed like this
% to have all desired results at once)
% Results are also saved


fs_adc = 3200;
N_adc = 2000;
n_adc = 1:N_adc;
fact_decimate = 25000;
fs = fs_adc * fact_decimate;
N = N_adc * fact_decimate;
n = 1:N;
Nfft = N;
 
fpass1 = 2.5e6;
fpass2 = 6e6;
fpass3 = 1600;

INR = [-23 -23];
DC = [0.35 0.65];

R = 1000;
ku = ones(1,R);
sk = ones(1,R);

ind = 1:fact_decimate:N;

for i = 1:length(INR)
    for r=1:R
        disp(strcat("DC = ", num2str(DC(i)), "   INR = ", num2str(INR(i)), "   R = ", num2str(r)))
        A0 = 1;
        noise = A0*normrnd(0,1,[1,N]);
        fo = rand()*fpass1;
        theta = 2*pi*rand();
        Mmin = 10;
        Mmax = 50;
        M = N/(Mmin + (Mmax-Mmin)*rand());
        offset = rand()*M;
        t = lowpass_FFT(10*log10(lowpass_FFT((sqrt(2*10^(INR(i)/10))*cos(2*pi*(fo/fs)*n + theta).*pulsetrain(n,M,DC(i),offset)) + noise, fpass1, fs, Nfft).^2),fpass3,fs,Nfft);
        u = t(ind);
        ku(r) = kurtosis(real(u));
        sk(r) = skewness(real(u));
    end
    ku_mean = mean(ku);
    sk_mean = mean(sk);
    str = strcat("rep_save_DC",num2str(DC(i)*100),"_INR",num2str(abs(INR(i))));
    save(str, 'ku_mean', 'sk_mean', 'ku', 'sk', 'DC', 'INR', 'i', 'R','N_adc','fs_adc','fact_decimate');
end
