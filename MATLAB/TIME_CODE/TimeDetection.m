clear all
close all
newcolors = [0.000 0.447 0.741
             0.850 0.325 0.098
             0.929 0.694 0.125
             0.494 0.184 0.556
             0.466 0.674 0.188
             0.301 0.745 0.933
             0.635 0.078 0.184
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54
             0.00 0.00 1.00
             0.00 0.50 0.00
             1.00 0.00 0.00
             0.00 0.75 0.75
             0.75 0.00 0.75
             0.75 0.75 0.00
             0.25 0.25 0.25];
         
N = 12; % input data x resolution
D = 4; % decimation factor
L = 12; % time constant mumean resolution
S = 2000; % initial number of x samples
s_vec = (1:S); % sample index % blank inhibit (0=always update, 1=selective update)
A0 = 1;
BW = 20e6;
f0 = 869e6;
fin = f0-BW/2;
ffin = f0+BW/2;
INRrefv = -20:4:0;
fs = 5e6;
U = 20;
I = S*0.2;
N_avg = 150;
b = linspace(0.2,2,U);
C = 1.3339;
det_avg = zeros(U,N_avg);
pfa_avg = zeros(U,N_avg);
betavspd = ones(length(INRrefv),U);
betavspfa = ones(length(INRrefv),U);
%betavspd_avg = ones(length(INRrefv),U,M);
pfa = 0;
pd = 0;


    fs_adc = 3200;
    N_adc = 2000;
    fact_decimate = 25000;
    fs = fs_adc * fact_decimate;
    N = N_adc * fact_decimate;
    n_vec = 1:N;
    Nfft = N;
    ind = 1:fact_decimate:N;
    fpass1 = 2.5e6;
    fpass2 = 6e6;
    fpass3 = 1600;
for l = 1:length(INRrefv)
    INRref = INRrefv(l);
    
 
    for n=1:N_avg
        
    A0 = 1;
    noise = A0*normrnd(0,1,[1,N]);
    fo = rand()*fpass1;
    theta = 2*pi*rand();
    Mmin = 10;
    Mmax = 50;
    M = N/(Mmin + (Mmax-Mmin)*rand());
    offset = rand()*M;
    DC = 0.5;
   pulses = pulsetrain(n_vec,M,DC,offset);
   i1 = sqrt(2*10^(INRrefv(l)/10))*cos(2*pi*(fo/fs)*n_vec + theta).*pulses;
   r1 = noise;
   t = noise + i1;
   u = lowpass_MXF(10*log10(abs(lowpass_MXF(t, fpass1, fs, Nfft)).^2),fpass3,fs,Nfft);


    pq = u(ind); 
    offset = mean(pq);
    sigma = std(pq);
        for k=1:U
           beta = b(k);

%% 
            th = zeros(1,S);
            pd = zeros(1,S);
            pd(1) = 0;
            for i = 2:S
                th(i) = beta*sigma + offset;
                if pq(i) < th(i)
                    pd(i) = 0; % no pulse detected
                    detected = 0;
                else
                    pd(i) = 1; % pulse detected
                    detected = 1;
                end
            end
            ndet = 0;
            nfa = 0;
            pulses_dec = pulses(ind);
            for j=1:S
                if(pulses_dec(j)==1 && pd(j)==1)
                    ndet=ndet+1;
                
                elseif(pulses_dec(j)==0 && pd(j)==1)
                    nfa = nfa+1;
           
                end
            end
            
            det_avg(k,n) = ndet/(DC*N_adc);

            pfa_avg(k,n) = nfa/(N_adc*(1-DC));

        end
        
    end
    for u = 1:U
        det(u)= mean(det_avg(u));
        pfa(u) = mean(pfa_avg(u));
    end
  betavspd(l,:) = det;
  betavspfa(l,:) = pfa;
  disp(l)
end
    


figure('Name', 'Pulse Blanking performance vs aggressiveness ')
colororder(newcolors)
plot(b, betavspd, 'LineWidth', 2);
grid on
grid minor
title('PB performance vs aggressiveness for different INR values');
%legend('0 dB', '1 dB', '2 dB', '3 dB', '4 dB', '5 dB','Location','NorthEastOutside');
%legend('0 dB', '2 dB', '4 dB', '6 dB', '8 dB', '10 dB', '12 dB','14 dB','16 dB','18 dB','Location','NorthEastOutside');
legend('-20 dB', '-16 dB', '-12 dB', '-8 dB', '-4 dB', '0 dB','Location','NorthEastOutside');
xlabel('Beta (Aggressiveness)','FontWeight','bold');
xlim([0.2 2]);
ylabel('Detection probability','FontWeight','bold');
ylim([0 1]);

figure('Name', 'Pulse Blanking PFA vs aggressiveness ')
colororder(newcolors)
plot(b, betavspfa, 'LineWidth', 2);
grid on
grid minor
title('PFA vs aggressiveness for different INR values');
legend('-20 dB', '-16 dB', '-12 dB', '-8 dB', '-4 dB', '0 dB','Location','NorthEastOutside');
xlabel('Beta (Aggressiveness)','FontWeight','bold');
xlim([0.2 2]);
ylabel('False alarm probability','FontWeight','bold'), 
ylim([0 0.05]);

figure(4)
subplot(211)
plot(s_vec,pq,'b')
hold on
plot(s_vec,th,'r')
grid
hold off
xlabel('Sample Index'), ylabel('Power (dBm)')
title(['Signal power (blue), Threshold = \beta\sigma (red), \beta=' num2str(beta)])
subplot(212),plot(s_vec,pd),grid
xlabel('Sample Index'), ylabel('Blank on/off')
title('Pulse Detection')