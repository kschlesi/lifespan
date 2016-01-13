%%% sydney

sydney = [10, 250 ;15, 150 ;20, 100;23, 75;29, 50;34, 37 ;39, 25 ;42, 20 ;47, 15 ;52, 10];
extrapx = [49; 50];

% plot original table data
figure; scatter(sydney(:,1),sydney(:,2));
hold on;

% n = 1 polynomial fit
p = polyfit(sydney(:,1),sydney(:,2),1);
plot(sydney(:,1),p(1)*sydney(:,1)+p(2));
plot(extrapx,p(1)*extrapx+p(2),'bd');

% n = 2 polynomial fit
p2 = polyfit(sydney(:,1),sydney(:,2),2);
plot(sydney(:,1),p2(1)*sydney(:,1).^2+p2(2)*sydney(:,1)+p2(3));
plot(extrapx,p2(1)*extrapx.^2+p2(2)*extrapx+p2(3),'md')

% exponential fit
efit = fit(sydney(:,1),sydney(:,2),'exp1');
plot(efit);
plot(extrapx,efit.a*exp(efit.b.*extrapx),'kd');

title('data extrapolation, UA biotech class');
xlabel('band (mm)'); ylabel('kD');

legend('table data', 'best fit line (n=1)', 'measured data (n=1)',...
        'best fit curve (n=2)', 'measured data (n=2)',...
        'best fit curve (exp)', 'measured data (exp)');


