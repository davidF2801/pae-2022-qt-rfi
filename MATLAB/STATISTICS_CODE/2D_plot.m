% Obtain 2D plots from the results of K_simulation_2D or L_simulation_2D 


fig1 = figure('Name', 'Kurtosis value of PULSED SIN + AWGN over INR and DC', 'Position', [10 100 625 625]);
imagesc([0 1],[-30 -0],ku_mean)
%imagesc([0 1],[-30 -0], std(ku_comb12345,0,3))
c = colorbar;
c.Label.String = 'Kurtosis value';
title('Kurtosis value of PULSED SIN + AWGN over INR and DC')
xlabel('DC')
ylabel('INR [dB]')
xticks(0:0.05:1)
yticks(-30:2:0)
%colormap(fig1,"hot")
saveas(fig1,'K_2D_sweep_ku_interp.fig');
saveas(fig1,'K_2D_sweep_ku_interp.png');

fig2 = figure('Name', 'Skewness value of PULSED SIN + AWGN over INR and DC', 'Position', [635 100 625 625]);
imagesc([0 1],[-30 -0],sk_mean)
%imagesc([0 1],[-30 -0],std(sk_comb12345,0,3))
c = colorbar;
c.Label.String = 'Skewness value';
title('Skewness value of PULSED SIN + AWGN over INR and DC')
xlabel('DC')
ylabel('INR [dB]')
xticks(0:0.05:1)
yticks(-30:2:0)
%colormap(fig2,"cool")
saveas(fig2,'K_2D_sweep_sk_interp.fig');
saveas(fig2,'K_2D_sweep_sk_interp.png');