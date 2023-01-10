% Can be used for both K and L-band
% The results from "L_simulation_repetitions" and "K_simulation_repetitions" of are used to plot comparative graphs of both the histogram with 1000 samples and the PDF calculated from such samples


close all

ku_new = normrnd(mean(ku),std(ku),[1,1000000]);
sk_new = normrnd(mean(sk),std(sk),[1,1000000]);

edges_ku = (mean(ku)-4*std(ku)):0.02:(mean(ku)+4*std(ku));
fig1=figure(1);
histogram(ku,edges_ku,'Normalization','probability');
hold on
histogram(ku_new,edges_ku,'Normalization', 'probability');
grid on
grid minor
title(strcat('Kurtosis distribution for DC = ',num2str(DC(i)*100),'%, INR = ',num2str(INR(i)),' dB'))
legend('Real raw samples, 1000 values', 'Generated from std and mean, 1e6 values', 'Location','southoutside');
dim = [0.2 0.5 0.3 0.3];
str = {strcat('mean = ', num2str(mean(ku))),strcat('stdv = ', num2str(std(ku)))};
annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
saveas(fig1,strcat('K_ku_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.fig'));
saveas(fig1,strcat('K_ku_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.png'));

edges_sk = (mean(sk)-4*std(sk)):0.008:(mean(sk)+4*std(sk));
fig2=figure(2);
histogram(sk,edges_sk,'Normalization','probability');
hold on
histogram(sk_new,edges_sk,'Normalization','probability');
grid on
grid minor
title(strcat('Skewness distribution for DC = ',num2str(DC(i)*100),'%, INR = ',num2str(INR(i)),' dB'))
legend('Real raw samples, 1000 values', 'Generated from std and mean, 1e6 values', 'Location','southoutside');
dim = [0.2 0.5 0.3 0.3];
str = {strcat('mean = ', num2str(mean(sk))),strcat('stdv = ', num2str(std(sk)))};
annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
saveas(fig2,strcat('K_sk_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.fig'));
saveas(fig2,strcat('K_sk_hist_DC',num2str(DC(i)*100),'_INR',num2str(abs(INR(i))),'.png'));
