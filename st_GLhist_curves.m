% create plots of spatiotemporal dissimilarity curves
% for GL history over gamma and omega

% first load all for age data.
% set parameters:
subj = 10;
%gammas = [0.85,1,1.15,1.3,1.45];
gamma = goRange(1);
%gamma = 0;
omegas = [0,0.005,0.05,0.1,0.2,0.3,0.4,0.5,1];
%omega = goRange(2);
%omega = 0;
p = 10;

% grab data and create plot.
matcell = A(subj).adj;
%for gamma = gammas
for omega = omegas
    [C,Q,C_hist,mnp] = genlouvainREPs(matcell,p,gamma,omega);
%%%%% or if omega=0
%     if ~omega
%         C = zeros(p,n,t);
%         C_hist = zeros(p,n,t,7);
%         maxp = 0;
%         for T=1:t
%             [C(:,:,T),~,Chi] = genlouvainREPs(matcell{T},p,gg,omega);
%             C_hist(:,:,T,1:size(Chi,3)) = Chi;
%             maxp = max(maxp,size(Chi,3));
%         end
%         C_hist = C_hist(:,:,:,1:maxp);
%     else
%         [C,~,C_hist,mnp] = genlouvainREPs(matcell,p,gg,omega);
%     end    
    figure(1);
    hold on; 
    [sC, tC] = st_GLhist(C_hist,mnp,'ExtFigureCmd');
    hold off;
    figure(2);
    hold on; 
    [sCm, tCm] = st_GLhist(C_hist,mnp,'ExtFigureCmd','MeanTraj');
    hold off;

%     figure; 
%     for T=1:t
%         subplot(1,t,T); bcolor(C(:,:,T)'); colorbar;
%      %   title(['\gamma=' num2str(gamma)])
%         title(['\omega=' num2str(omega)])
%     end
%     %suptitle(['communities, subject ' num2str(subj) ', \omega=' num2str(omega)]);
%     suptitle(['communities, subject ' num2str(subj) ', \gamma=' num2str(gamma)]);

    figure; 
    for mm=1:size(C_hist,4)-1
        subplot(1,size(C_hist,4)-1,mm); bcolor(squeeze(C_hist(1,:,:,mm))); colorbar;
        title(['step' num2str(mm)]);
    end
    suptitle(['GL history, subject ' num2str(subj) ...
        ', \gamma=' num2str(gamma) ', \omega=' num2str(omega)]);
end

figure(1);
hold on;
title(['subject ' num2str(subj) ', 316tw, \omega=' num2str(omega)]);
hold off;

figure(2);
hold on;
%title(['subject ' num2str(subj) ', 316tw, \omega=' num2str(omega)]);
%legend(num2str(gammas'),'Location','NorthWest');
title(['subject ' num2str(subj) ', 316tw, \gamma=' num2str(gamma)]);
legend(num2str(omegas'),'Location','NorthWest');
hold off;
