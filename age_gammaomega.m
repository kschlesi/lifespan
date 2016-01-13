% age gamma/omega chooser

inpath = '/Users/kimberly/Documents/lifespan/CD/';

%goRange = [1,0;1.1,0;1.15,0;1.2,0;1.25,0];  % for novis
%goRange = [1.15,0.001;1.15,0.005;1.15,0.01;1.15,0.05;1.15,0.1;1.15,0.15;1.15,0.2];   % for novis
 goRange = [0.8,0;...
          1.0,0;...
          1.05,0;...    
          1.1,0;...
           1.15,0;...
          1.2,0;...
            1.25,0;...
%            1.3,0;...
%           1.4,0;...
%            1.0, 0.005;
%            1.10,0.005;
%            1.15,0.005;
%           1.20,0.005;
%            1.10,0.01;...
%            1.1,0.01;...
%            1.15,0.01;...
%           1.2,0.01;...
%            1.0,0.05;...
%            1.1,0.05;...
%            1.15,0.05;...
%            1.2,0.05;
%            1.2,0.1;
%            1.2,0.25;
%            1.2,0.5;
%            1.2,0.75;
            ];
overgos = 1;  % 1 if plotting over gammas, 2 if plotting over omegas
%kplot = 1:105;
%kplot = 1:4;
kplot = 1:numel(ib);
fplot = {'Zavg','Zvar','Qavg','Qvar','NCavg','Fxavg','Fxvar'};
% Zavg Zvar Qavg Qvar NCavg NCvar Fxavg Fxvar

nsubjs = 108;
missing_subjs = [46;64;81;82];
%nsubjs = 4;
ts = 316;
ts_run = 316;
nruns = 3;
inc_leg = 0;
tag = []; assert(n==194); nodeix = 1:n;
%tag = '_nv'; nodeix = ~ismemvar(regClass,'Visual');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = size(goRange,1);
[t,n,ib,Cplot,Cplotall,Zvg,Qvg,Ncvg] = load_CD_results(nsubjs,...
                        missing_subjs, nruns, ts, ts_run, goRange, tag, inpath);
p = size(Cplot,2)/numel(ib);
pplot = ib(kplot);
flexS = zeros(g,numel(ib),n); % compute flexibilities
for goix=1:g
    for k=1:numel(ib)
        partns = squeeze(mode(squeeze(Cplotall(goix,p*(k-1)+1:p*k,:,:)),1));
        if goRange(goix,2)==0
            
        end
        flexS(goix,k,:) = shiftdim(flexibility(partns',1),-1);
    end
end                    
                    
if ~any(goRange(:,2))
    sliceix = 1:t;
else
    sliceix = 1;
end

    
figure;
bcolor(squeeze(mode(squeeze(Cplot(1,p*(k-1)+1:p*k,:,:)),1)));
title(['\gamma = ' num2str(goRange(1,1)) ', \omega = ' ...
        num2str(goRange(1,2))]);
    colorbar('Location','EastOutside');
for goix=2:g;
    figure; 
    bcolor(squeeze(mode(squeeze(Cplotall(goix,p*(k-1)+1:p*k,:,:)),1)));
    title(['\gamma = ' num2str(goRange(goix,1)) ', \omega = ' ...
        num2str(goRange(goix,2))]);
    colorbar('Location','EastOutside');
end;

%%  plotting

    for i=overgos
    %%% average and variance in z-scores, modularity, NComms per subject
    % zscore
        % average
        if any(ismember(fplot,'Zavg'))
        figure;
        plot(goRange(:,i),mean(Zvg(:,kplot,1,sliceix),4),':');
        hold on;
        plot(goRange(:,i),mean(mean(Zvg(:,kplot,1,sliceix),4),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('rand z-score')
        title(['average rand z-score over ' num2str(p) ' genlouvain runs, orig.'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        if any(ismember(fplot,'Zavg'))
        figure;
        plot(goRange(:,i),mean(Zvg(:,kplot,2,sliceix),4),':');
        hold on;
        plot(goRange(:,i),mean(mean(Zvg(:,kplot,2,sliceix),4),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('rand z-score')
        title(['average rand z-score over ' num2str(p) ' genlouvain runs, CC'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        % variance
        if any(ismember(fplot,'Zvar'))
        figure;
        plot(goRange(:,i),mean(Zvg(:,kplot,3,sliceix),4),':');
        hold on;
        plot(goRange(:,i),mean(mean(Zvg(:,kplot,3,sliceix),4),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('rand z-score')
        title(['variance in rand z-score over ' num2str(p) ' genlouvain runs, orig.'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        figure;
        zorig = mean(mean(Zvg(3:end,kplot,1,sliceix),4),2);
        csize = (n./mean(Ncvg(3:end,kplot,1,end),2));
        plot(goRange(3:end,i),zorig.*csize)%,goRange(3:end,i),csize)
        
        if any(ismember(fplot,'Zvar'))
        figure;
        plot(goRange(:,i),mean(Zvg(:,kplot,4,sliceix),4),':');
        hold on;
        plot(goRange(:,i),mean(mean(Zvg(:,kplot,4,sliceix),4),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('rand z-score')
        title(['variance in rand z-score over ' num2str(p) ' genlouvain runs, CC'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
    
    % modularity
        % average
        if any(ismember(fplot,'Qavg'))
        figure;
        plot(goRange(:,i),mean(Qvg(:,kplot,2,sliceix),4),':');
        hold on;
        plot(goRange(:,i),mean(mean(Qvg(:,kplot,2,sliceix),4),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('modularity')
        title(['average modularity over ' num2str(p) ' genlouvain runs'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        % variance
        if any(ismember(fplot,'Qvar'))
        figure;
        plot(goRange(:,i),mean(Qvg(:,kplot,4,sliceix),4),':');
        hold on;
        plot(goRange(:,i),mean(mean(Qvg(:,kplot,4,sliceix),4),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('modularity')
        title(['variance in modularity over ' num2str(p) ' genlouvain runs'])
        if inc_leg
        legend(num2str(pplot))
        end
        end

    % number of communities
        % average
        if any(ismember(fplot,'NCavg'))
        figure;
        plot(goRange(:,i),Ncvg(:,kplot,1,end),':');
        hold on;
        plot(goRange(:,i),mean(Ncvg(:,kplot,1,end),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('number of communities')
        title(['average number communities over ' num2str(p) ' genlouvain runs, orig'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        if any(ismember(fplot,'NCavg'))
        figure;
        plot(goRange(:,i),Ncvg(:,kplot,2,end),':');
        hold on;
        plot(goRange(:,i),mean(Ncvg(:,kplot,2,end),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('number of communities')
        title(['average number communities over ' num2str(p) ' genlouvain runs, CC'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        % variance
        if any(ismember(fplot,'NCvar'))
        figure;
        plot(goRange(:,i),Ncvg(:,kplot,3,end),':');
        hold on;
        plot(goRange(:,i),mean(Ncvg(:,kplot,3,end),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('number of communities')
        title(['variance in number communities over ' num2str(p) ' genlouvain runs, orig'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        if any(ismember(fplot,'NCvar'))
        figure;
        plot(goRange(:,i),Ncvg(:,kplot,4,end),':');
        hold on;
        plot(goRange(:,i),mean(Ncvg(:,kplot,4,end),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('number of communities')
        title(['variance in number communities over ' num2str(p) ' genlouvain runs, CC'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
    %%% average and variance in flexibility per subject
        % average
        if any(ismember(fplot,'Fxavg'))
        figure;
        plot(goRange(:,i),mean(flexS(:,kplot,~~nodeix),3),':');
        hold on;
        plot(goRange(:,i),mean(mean(flexS(:,kplot,~~nodeix),3),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('flexibility')
        title(['average flexibility over ' num2str(n) ' nodes'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
        % variance
        if any(ismember(fplot,'Fxvar'))
        figure;
        plot(goRange(:,i),var(flexS(:,kplot,~~nodeix),[],3),':');
        hold on;
        plot(goRange(:,i),mean(var(flexS(:,kplot,~~nodeix),[],3),2),'k');
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('flexibility')
        title(['variance in flexibility over ' num2str(n) ' nodes'])
        if inc_leg
        legend(num2str(pplot))
        end
        end
        
    end  % loop over overgos
    
    
% goRange = [1.15,0;
%            1.15,0.0001;
%            1.15,0.0005;
%            1.15,0.001;
%            1.15,0.005;
%            1.15,0.01;
% %            1.15,0.0125;
% %            1.15,0.025;
% %            1.15,0.0375;
% %            1.15,0.05;
% %            1.15,0.0625;
% %            1.15,0.1;
% %            1.15,0.15;
% %            1.15,0.2;
% %            1.15,0.25;
% %            1.15,0.3;
% %            1.15,0.5
%            ];
% goRange = [0.8,0;
%            0.9,0;
%            1.0,0;
%            1.05,0;
%            1.1,0;
%            1.15,0;
%            1.2,0;
%            1.25,0;
%            1.3,0;
%            1.4,0;
%            1.65,0
%           ];

    