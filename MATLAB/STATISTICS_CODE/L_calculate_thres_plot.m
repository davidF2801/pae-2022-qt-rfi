% For the L-band receiver
% Generates thresholds for kurtosis and skewness given a certain probability of false alarm, Pfa. It also plots such thresholds together with the DC = 0% PDF
% ku_mean, ku, sk_mean and sk are generated from "L_simulation_repetitions", being ku and sk 1000-posiitons-long vectors holding ku and sk from each repetition


close all;
N = 1e7;
ku_noise_mean = 3.0012;
ku_noise_std = 0.10812;
ku_noise = normrnd(ku_noise_mean, ku_noise_std,[1,N]);

ku_ref_mean = ku_mean;
ku_ref_std = std(ku);
ku_ref = normrnd(ku_mean,ku_ref_std,[1,N]);

sk_noise_mean = -0.025487;
sk_noise_std = 0.05581;
sk_noise = normrnd(sk_noise_mean,sk_noise_std,[1,N]);

sk_ref_mean = sk_mean;
sk_ref_std = std(sk);
sk_ref = normrnd(sk_mean,sk_ref_std,[1,N]);

Pfa = 5e-4;
Pfa_kusk = (1-sqrt(1-4*Pfa))/2;
Pfa_kusk_inf = Pfa_kusk/2;
Pfa_kusk_sup = 1-Pfa_kusk/2;

ku_low_thres = icdf("Normal", Pfa_kusk_inf, ku_noise_mean, ku_noise_std);
ku_upp_thres = icdf("Normal", Pfa_kusk_sup, ku_noise_mean, ku_noise_std);
sk_low_thres = icdf("Normal", Pfa_kusk_inf, sk_noise_mean, sk_noise_std);
sk_upp_thres = icdf("Normal", Pfa_kusk_sup, sk_noise_mean, sk_noise_std);



ku_min = (mean(ku_noise)-4*std(ku_noise));
ku_max = (mean(ku_noise)+4*std(ku_noise));
edges_ku = ku_min:0.007:ku_max;
fig1=figure(1);
histogram(ku_noise,edges_ku,'Normalization','probability');
grid on;
grid minor;
hold on
xline(ku_upp_thres,'LineWidth',1,'Color',[1 0.1 0.1]);
xline(ku_low_thres,'Linewidth',1,'Color',[0.1 1 0.1]);
title('Kurtosis PDF with the corresponding decision thresholds')
dim = [0.22 0.62 0.3 0.3];
str = {strcat('Upper threshold = ',num2str(ku_upp_thres)), strcat('Lower threshold = ',num2str(ku_low_thres)), strcat('mean = ', num2str(ku_noise_mean)), strcat('stdv = ', num2str(ku_noise_std))};
annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
saveas(fig1,strcat('L_ku_hist_thres.fig'));
saveas(fig1,strcat('L_ku_hist_thres.png'));

sk_min = (mean(sk_noise)-4*std(sk_noise));
sk_max = (mean(sk_noise)+4*std(sk_noise));
edges_sk = sk_min:0.003:sk_max;
fig2=figure(2);
histogram(sk_noise,edges_sk,'Normalization','probability');
grid on;
grid minor;
hold on
xline(sk_upp_thres,'LineWidth',1,'Color',[1 0.1 0.1]);
xline(sk_low_thres,'Linewidth',1,'Color',[0.1 1 0.1]);
title('Skewness PDF with the corresponding decision thresholds')
dim = [0.22 0.62 0.3 0.3];
str = {strcat('Upper threshold = ',num2str(sk_upp_thres)), strcat('Lower threshold = ',num2str(sk_low_thres)), strcat('mean = ', num2str(sk_noise_mean)), strcat('stdv = ', num2str(sk_noise_std))};
annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
saveas(fig2,strcat('L_sk_hist_thres.fig'));
saveas(fig2,strcat('L_sk_hist_thres.png'));
