% age gamma/omega chooser

inpath = '/Users/kimberly/Documents/lifespan/CD/';

goRange = [1,0.05];
overgos = 2;  % 1 if plotting over gammas, 2 if plotting over omegas
kplot = [1 2 3 4 5];
fplot = {'Zavg','Zvar'};
% Zavg Zvar Qavg Qvar NCavg NCvar Fxavg Fxvar

nsubjs = 108;
missing_subjs = [46;81;82];
ts = 316;
ts_run = 316;
nruns = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = size(goRange,1);
[t,n,ib,Cplot,Cplotall,Zvg,Qvg,Ncvg] = load_CD_results(nsubjs,...
                        missing_subjs, nruns, ts, ts_run, goRange, inpath);
pplot = ib(kplot);
flexS = zeros(g,numel(ib),n); % compute flexibilities
for goix=1:g
    for k=1:numel(ib)
        partns = squeeze(mode(squeeze(Cplotall(goix,p*(k-1)+1:p*k,:,:)),1));
        flexS(goix,k,:) = shiftdim(flexibility(partns',1),-1);
    end
end                    
                    
if ~omega
    sliceix = 1:t;
else
    sliceix = 1;
end
    
    for i=overgos
    %%% average and variance in z-scores, modularity, NComms per subject
    % zscore
        % average
        if any(ismember(fplot,'Zavg'))
        figure;
        plot(goRange(:,i),mean(Zvg(:,kplot,2,sliceix),4));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('rand z-score')
        title(['average rand z-score over ' num2str(p) ' genlouvain runs'])
        legend(num2str(pplot))
        end
        
        % variance
        if any(ismember(fplot,'Zvar'))
        figure;
        plot(goRange(:,i),mean(Zvg(:,kplot,4,sliceix),4));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('rand z-score')
        title(['variance in rand z-score over ' num2str(p) ' genlouvain runs'])
        legend(num2str(pplot))
        end
    
    % modularity
        % average
        if any(ismember(fplot,'Qavg'))
        figure;
        plot(goRange(:,i),mean(Qvg(:,kplot,2,sliceix),4));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('modularity')
        title(['average modularity over ' num2str(p) ' genlouvain runs'])
        legend(num2str(pplot))
        end
        
        % variance
        if any(ismember(fplot,'Qvar'))
        figure;
        plot(goRange(:,i),mean(Qvg(:,kplot,4,sliceix),4));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('modularity')
        title(['variance in modularity over ' num2str(p) ' genlouvain runs'])
        legend(num2str(pplot))
        end

    % number of communities
        % average
        if any(ismember(fplot,'NCavg'))
        figure;
        plot(goRange(:,i),mean(Ncvg(:,kplot,2,end),4));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('number of communities')
        title(['average number communities over ' num2str(p) ' genlouvain runs'])
        legend(num2str(pplot))
        end
        
        % variance
        if any(ismember(fplot,'NCvar'))
        figure;
        plot(goRange(:,i),mean(Ncvg(:,kplot,4,end),4));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('number of communities')
        title(['variance in number communities over ' num2str(p) ' genlouvain runs'])
        legend(num2str(pplot))
        end
        
    %%% average and variance in flexibility per subject
        % average
        if any(ismember(fplot,'Fxavg'))
        figure;
        plot(goRange(:,i),mean(flexS(:,kplot,:),3));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('flexibility')
        title(['average flexibility over ' num2str(n) ' nodes'])
        legend(num2str(pplot))
        end
        
        % variance
        if any(ismember(fplot,'Fxvar'))
        figure;
        plot(goRange(:,i),var(flexS(:,kplot,:),[],3));
        if i==1, xlabel('\gamma');
        else xlabel('\omega'); end
        ylabel('flexibility')
        title(['variance in flexibility over ' num2str(n) ' nodes'])
        legend(num2str(pplot))
        end
        
    end  % loop over overgos
    