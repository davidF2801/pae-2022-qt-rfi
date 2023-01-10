% Obtain 3D plots from the results of K_simulation_2D or L_simulation_2D 


DC = 0:0.05:1;
INR = -30:2:0;

fig1 = figure('Name', 'Kurtosis value of PULSED SIN + AWGN over INR and DC', 'Position', [10 100 625 625]);
surf(DC, INR, log10(ku_mean),'EdgeColor','none')
% colormap(fig1,"winter")
hold on
surf(DC, INR, log10(ku_upp_thres)*ones(length(INR), length(DC)), 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceColor','red');
hold on
surf(DC, INR, log10(ku_low_thres)*ones(length(INR), length(DC)), 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceColor','magenta');
title('Kurtosis value of PULSED SIN + AWGN over INR and DC')
xlabel('DC')
ylabel('INR [dB]')
zlabel('Kurtosis value (logarithmic scale)')
xticks(DC)
yticks(INR)
view(35,10)
saveas(fig1,'K_3D_sweep_ku_interp.fig');
saveas(fig1,'K_3D_sweep_ku_interp.png');

fig2 = figure('Name', 'Skewness value of PULSED SIN + AWGN over INR and DC', 'Position', [635 100 625 625]);
surf(DC, INR, sk_mean,'EdgeColor','none')
% colormap(fig2,"spring")
hold on
surf(DC, INR, sk_upp_thres*ones(length(INR), length(DC)), 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceColor','red');
hold on
surf(DC, INR, sk_low_thres*ones(length(INR), length(DC)), 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceColor','magenta');
title('Skewness value of PULSED SIN + AWGN over INR and DC')
xlabel('DC')
ylabel('INR [dB]')
zlabel('Skewness value')
xticks(DC)
yticks(INR)
view(35,10)
saveas(fig2,'K_3D_sweep_sk_interp.fig');
saveas(fig2,'K_3D_sweep_sk_interp.png');