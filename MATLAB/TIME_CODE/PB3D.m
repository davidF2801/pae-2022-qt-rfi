clear all; close all;
fs_adc = 3200;
N_adc = 2000;
n_adc = 1:N_adc;
fact_decimate = 5000;
fs = fs_adc * fact_decimate;
N = N_adc * fact_decimate;
n = 1:N;
Nfft = N;
beta = 0.75;
fpass1 = 5e6;
fpass2 = 1e7;
fpass3 = 8000;

INRmax = 0;
INRmin = -30;
INRnum = 16;
%INRnum = 4;
INR = linspace(INRmin, INRmax, INRnum);

DCmax = 1;
DCmin = 0;
DCnum = 21;
%DCnum = 6;
DC = linspace(DCmin,DCmax,DCnum);

R = 5;
%R = 1;
%R = 3;
pdet = ones(INRnum,DCnum,R);
pfa = ones(INRnum,DCnum,R);

ind = 1:fact_decimate:N;

for r=1:R
    disp(strcat("-->R = ", num2str(r)))
    A0 = 1;
    noise = A0*normrnd(0,1,[1,N]);
    fo = rand()*fpass1;
    theta = 2*pi*rand();
    M = N/(rand()*100);
    for i=1:INRnum
        disp(strcat("*** R = ",num2str(r)," INR = ", num2str(INR(i))));
        for j=1:DCnum
            disp(strcat("R = ",num2str(r)," INR = ", num2str(INR(i))," DC = ", num2str(DC(j))));
            if(DC(j)==0)
                A1 = 0;
            else
                A1 = sqrt((2*10.^(INR(i)/10))/DC(j));
            end
            pulses = pulsetrain(n,M,DC(j),rand()*M);
            i1 = A1*cos(2*pi*(fo/fs)*n + theta).*pulses;
            pulses_dec = pulses(ind);
            
            t = lowpass_MXF(10*log10(lowpass_MXF(i1 + noise, fpass1, fs, Nfft).^2),fpass3,fs,Nfft);
            u = t(ind);
           
            
            th = zeros(1,N_adc);
            pd = zeros(1,N_adc);
            sigma = std(u);
            offset = mean(u);
            
            for k = 1:N_adc
                th(k) = beta*sigma + offset;
                if u(k) < th(k)
                    pd(r) = 0; % no pulse detected
                    detected = 0;
                else
                    pd(k) = 1; % pulse detected
                    detected = 1;
                end
            end
            ndet = 0;
            nfa = 0;
            pulses_dec = pulses(ind);
            for b=1:N_adc
                if(pulses_dec(b)==1 && pd(b)==1)
                    ndet=ndet+1;
                
                elseif(pulses_dec(b)==0 && pd(b)==1)
                    nfa = nfa+1;
           
                end
            end
            if(DC(j)==0)
                pdet(i,j,r) = 0;
                pfa(i,j,r) = nfa/(N_adc*(1-DC(j)));
            elseif(DC(j) == 1)
                
                pdet(i,j,r) = ndet/(DC(j)*N_adc);
                pfa(i,j,r) = nfa/(N_adc);
            else
                pdet(i,j,r) = ndet/(DC(j)*N_adc);
                pfa(i,j,r) = nfa/(N_adc*(1-DC(j)));
            end
            
            
        end
    end
end
pdet_mean = mean(pdet,3);
pfa_mean = mean(pfa,3);
figure(1);
surf(DC, INR, pdet_mean);
xlabel('DC');
ylabel('INR');
zlabel('Detection probability');
title('Detection probability for \beta = 0.75 vs DC and INR');
colorbar

figure(2);
surf(DC, INR, pfa_mean);
xlabel('DC');
ylabel('INR');
zlabel('False alarm probability');
title('False alarm probability for \beta = 0.75 vs DC and INR');
colorbar
save('prova_save.mat', 'th', 'DC', 'INR', 'N_adc','fs_adc','fact_decimate');

