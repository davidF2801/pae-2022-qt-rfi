% For K-band receiver
% Plots the PDFs of the DC=0% and the one coming from ku_mean, sk_mean, ku and sk, generated from "K_simulation_repetitions".
% The thresholds are also plotted
% It calculates the probability of detection
% Saves the results


close all

% For Pfa = 5e-4
ku_upp_thres = 3.381210474816579;
ku_low_thres = 2.619789525183424;
sk_upp_thres = 0.181464536987469;
sk_low_thres = -0.200359736987467;

ku_noise = normrnd(3.0005, 0.10938,[1,1000000]);
ku_ref = normrnd(ku_mean,std(ku),[1,1000000]);

sk_noise = normrnd(-0.0094476,0.05485,[1,1000000]);
sk_ref = normrnd(sk_mean,std(sk),[1,1000000]);

Pd_ku = cdf("Normal", ku_low_thres, ku_mean, std(ku)) + (1 - cdf("Normal", ku_upp_thres, ku_mean, std(ku))) 
Pd_sk = cdf("Normal", sk_low_thres, sk_mean, std(sk)) + (1 - cdf("Normal", sk_upp_thres, sk_mean, std(sk)))
Pd = Pd_ku + Pd_sk - Pd_ku*Pd_sk

ku_min = min([(mean(ku_noise)-4*std(ku_noise)) (mean(ku_ref)-4*std(ku_ref))]);
ku_max = max([(mean(ku_noise)+4*std(ku_noise)) (mean(ku_ref)+4*std(ku_ref))]);
edges_ku = ku_min:0.008:ku_max;
fig1=figure(1);
histogram(ku_noise,edges_ku,'Normalization','probability');
hold on
histogram(ku_ref,edges_ku,'Normalization', 'probability');
grid on
grid minor
hold on
xline(ku_upp_thres,'LineWidth',1,'Color',[1 0.1 0.1]);
xline(ku_low_thres,'Linewidth',1,'Color',[0.1 1 0.1]);
title(strcat('Kurtosis distr. of noise compared to DC = ', num2str(DC(i)*100),'%, INR = ',num2str(INR(i)),' dB compared'))
legend('Noise kurtosis distribution', strcat('Noise + interference (DC=',num2str(DC(i)*100),'%, INR=',num2str(INR(i)),'dB) kurtosis distribution'), 'Location','southoutside');
saveas(fig1,strcat('K_ku_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'_compare.fig'));
saveas(fig1,strcat('K_ku_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'_compare.png'));


sk_min = min([(mean(sk_noise)-4*std(sk_noise)) (mean(sk_ref)-4*std(sk_ref))]);
sk_max = max([(mean(sk_noise)+4*std(sk_noise)) (mean(sk_ref)+4*std(sk_ref))]);
edges_sk = sk_min:0.005:sk_max;
fig2=figure(2);
histogram(sk_noise,edges_sk,'Normalization','probability');
hold on
histogram(sk_ref,edges_sk,'Normalization', 'probability');
grid on
grid minor
hold on
xline(sk_upp_thres,'LineWidth',1,'Color',[1 0.1 0.1]);
xline(sk_low_thres,'Linewidth',1,'Color',[0.1 1 0.1]);
title(strcat('Skewness distr. of noise compared to DC = ', num2str(DC(i)*100),'%, INR = ',num2str(INR(i)),' dB compared'))
legend('Noise skewness distribution', strcat('Noise + interference (DC=',num2str(DC(i)*100),'%, INR=',num2str(INR(i)),'dB) skewness distribution'), 'Location','southoutside');
saveas(fig2,strcat('K_sk_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'_compare.fig'));
saveas(fig2,strcat('K_sk_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'_compare.png'));

save(strcat('K_prob_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.mat'), 'ku_upp_thres', 'ku_low_thres', 'sk_upp_thres', 'sk_low_thres', 'Pfa', 'Pfa_kusk', 'Pd', 'Pd_ku', 'Pd_sk');