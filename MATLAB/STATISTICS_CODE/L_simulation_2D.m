% For the L-band receiver
% Allows to compute R repetitions from the INR-DC defined pairs, thus defining a 2D measurement plane
% In the end, results are also saved


fs_adc = 3200;
N_adc = 2000;
fact_decimate = 25000;
fs = fs_adc * fact_decimate;
N = N_adc * fact_decimate;
n = 1:N;
Nfft = N;

fpass1 = 2.5e6;
fpass2 = 6e6;
fpass3 = 1600;

INRmax = 0;
INRmin = -30;
INRnum = 16;
INR = linspace(INRmin, INRmax, INRnum);

DCmax = 1;
DCmin = 0;
DCnum = 21;
DC = linspace(DCmin,DCmax,DCnum);

R = 5;
ku = ones(INRnum,DCnum,R);
sk = ones(INRnum,DCnum,R);

ind = 1:fact_decimate:N;

for r=1:R
    disp(strcat("-->R = ", num2str(r)))
    A0 = 1;
    noise = A0*normrnd(0,1,[1,N]);
    f0 = fpass1*rand();
    theta = 2*pi*rand();
    Mmin = 10;
    Mmax = 50;
    M = N/(Mmin + (Mmax-Mmin)*rand());
    offset = rand()*M;
    for i=1:INRnum
        disp(strcat("*** R = ",num2str(r)," INR = ", num2str(INR(i))));
        for j=1:DCnum
            disp(strcat("R = ",num2str(r)," INR = ", num2str(INR(i))," DC = ", num2str(DC(j))));
            t = lowpass_FFT(10*log10(lowpass_FFT((sqrt(2*10^(INR(i)/10))*cos(2*pi*(f0/fs)*n + theta).*pulsetrain(n,M,DC(j),offset)) + noise, fpass1, fs, Nfft).^2),fpass3,fs,Nfft);
            u = t(ind);
            ku(i,j,r) = kurtosis(real(u));
            sk(i,j,r) = skewness(real(u));
        end
    end
end
ku_mean = mean(ku,3);
sk_mean = mean(sk,3);
save('L_save.mat', 'ku_mean', 'sk_mean', 'ku', 'sk', 'DC', 'INR', 'N_adc','fs_adc','fact_decimate');