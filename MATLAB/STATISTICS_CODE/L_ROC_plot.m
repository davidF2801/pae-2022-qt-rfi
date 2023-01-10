% For L-band receiver
% Plots the ROC from the DC=0% PDF together with the one coming from ku_mean, sk_mean, ku and sk, generated from "L_simulation_repetitions".
% Saves the plots


close all

line_ku = true;
line_sk = true;

ku_noise_mean = 3.0012;
ku_noise_std = 0.10812;

ku_ref_mean = ku_mean;
ku_ref_std = std(ku);

sk_noise_mean = -0.025487;
sk_noise_std = 0.05581;

sk_ref_mean = sk_mean;
sk_ref_std = std(sk);

ku_thres_min = 2.6;
ku_thres_max = 3.8;
ku_thres = ku_thres_min:0.001:ku_thres_max;

if(ku_mean<3.0012) %lower kurtosis threshold
    ku_pd = cdf("Normal", ku_thres, ku_ref_mean, ku_ref_std);
    ku_pfa = cdf("Normal", ku_thres, ku_noise_mean, ku_noise_std);
else %upper kurtosis threshold
    ku_pd = 1-cdf("Normal", ku_thres, ku_ref_mean, ku_ref_std);
    ku_pfa = 1-cdf("Normal", ku_thres, ku_noise_mean, ku_noise_std);
end

fig1 = figure(1);
plot(ku_pfa,ku_pd,'LineWidth',1.2,'Color','blue')
if(line_ku)
    hold on;
    plot(ku_pfa,ku_pfa,'LineWidth',1.2,'Color',[0.6 0.6 1],'LineStyle',':')
end
% xlim([0 max(ku_pfa)])
% ylim([min(ku_pd) 1])
xlim([0 1])
ylim([0 1])
grid on;
grid minor;
title(strcat('ROC curve on kurtosis of a pulsed sine RFI with DC=',num2str(DC(i)*100),'% and INR=',num2str(INR(i)),'dB'));
xlabel('Probability of false alarm (P_f_a)');
ylabel('Probability of detection (P_d)')
saveas(fig1,strcat('L_ku_ROC_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.fig'));
saveas(fig1,strcat('L_ku_ROC_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.png'));


sk_thres_min = -0.4;
sk_thres_max = 0.2;
sk_thres = sk_thres_min:0.001:sk_thres_max;

if(sk_mean<-0.025487) %lower skewness threshold
    sk_pd = cdf("Normal", sk_thres, sk_ref_mean, sk_ref_std);
    sk_pfa = cdf("Normal", sk_thres, sk_noise_mean, sk_noise_std);
else %upper skewness threshold
    sk_pd = 1-cdf("Normal", sk_thres, sk_ref_mean, sk_ref_std);
    sk_pfa = 1-cdf("Normal", sk_thres, sk_noise_mean, sk_noise_std);
end

fig2 = figure(2);
plot(sk_pfa,sk_pd,'LineWidth',1.2,'Color','red')
if(line_sk)
    hold on;
    plot(sk_pfa,sk_pfa,'LineWidth',1.2,'Color',[1 0.6 0.6],'LineStyle',':')
end
% xlim([0 max(sk_pfa)])
% ylim([min(sk_pd) 1])
xlim([0 1])
ylim([0 1])
grid on;
grid minor;
title(strcat('ROC curve on skewness of a pulsed sine RFI with DC=',num2str(DC(i)*100),'% and INR=',num2str(INR(i)),'dB'));
xlabel('Probability of false alarm (P_f_a)');
ylabel('Probability of detection (P_d)');
saveas(fig2,strcat('L_sk_ROC_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.fig'));
saveas(fig2,strcat('L_sk_ROC_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.png'));
