% wtc test

nregs = 60; % number of cortical brain regions per hemisphere (choose atlas)
tothresh = 0;
subjects = 22; % number of subjects
missing_subjs = [2;5;8;14;21];
nruns = 9;

% ib = vector containing indices of subjects with full data for all tasks
ib = removeval((1:subjects)',missing_subjs);
ib = ib(1:4);
disp(ib);

% set number of nodes (brain regions) & timesteps per run
    missing = zeros(subjects,nruns,2);
    switch nregs
        case 33, n = 83;  
        case 60, n = 129; missing(7,:,1) = [110,110,110,110,0,0,0,110,0];
                          missing(12,2,1) = 68;
                          missing(19,:,1) = [46,46,46,46,46,46,46,46,46];
                          %missing(21,1,:) = [ALL];
        case 125, n = 234; missing(1,:,1) = [208,0,0,92,92,92,0,92,92];
    end
    ts = 145; 
    TR = 2; % in seconds
    bandpass = [0.06,0.125]; % frequency band in Hz
%     tvec = (1:ts).*TR;

for k=1:1%numel(ib)
    person = ib(k);
    disp(['SUBJECT ' num2str(person)]);
     
    for runnum=1:1%nruns
        disp(['run ' num2str(runnum)]);
        %inrun = Linfo(:,rNum>0)==runnum; % rows of reward data to use in run
        
        % load data (timeseries)
        Lseries = zeros(n,ts);
        for reg=removeval(1:n,squeeze(missing(person,runnum,:))) % list of regions to include in this run
            Lseries(reg,:) = dlmread(['/Users/kimberly/Google Drive/choking/' num2str(person+100) ...
                '_ts/run' num2str(runnum) '/laus' num2str(nregs) '_ts' num2str(reg) '.txt']);
            if size(Lseries,2)~=ts
                disp([person,runnum,reg,size(Lseries,2)]);
                error('Timeseries do not have consistent lengths!');
            end
        end

        adj = zeros(n,n);
        tic;
        for i=1:n
            pctdone=i/n
            adj(i,i) = 0;
            for j=i+1:n
                [Rsq,period,scale] = wtc(Lseries(i,:)',Lseries(j,:)','mcc',0);
%                 [w1,period1,scale1] = wt(Lseries(i,:)');
%                 [w2,period2,scale2] = wt(Lseries(j,:)');
%                 [Xwt,periodx,scalex] = xwt(Lseries(i,:)']);
                freq = 1./(period*TR);
                    if j==2 && i==1
%                     figure
%                     plot(w1)
%                     title('wavelet spectrum 1')
%                     figure
%                     plot(w2)
%                     figure
%                     plot(Xwt)
                        
                    figure
                    h = pcolor(Rsq);
                    set(h,'EdgeColor','none')
                    xlabel('time (sec)')
                    set(gca,'XTickLabel',(1:size(Rsq,2)).*TR)
                    ylabel('frequency (1/sec)')
                    set(gca,'YTickLabel',freq)
                    end
                adj(i,j) = mean(mean(Rsq((freq<bandpass(2)&freq>bandpass(1)),:)));
            end
        end
        adj = adj + adj';
        toc;
        
    end
    
    figure
    h = pcolor(Rsq);
    set(h,'EdgeColor','none')
    
end

%                 [Rsq,period,scale] = wtc([tvec',Lseries(i,:)'],...
%                     [tvec',Lseries(j,:)'],'mcc',0);
%                 freq = 1./(period);
%                 adj(i,j) = mean(mean(Rsq((freq<bandpass(2)&freq>bandpass(1)),:)));

%                 [w1,period1,scale1] = wt([tvec',Lseries(i,:)'],'S0',4*tslength,'J1',0,'Dj',1);
%                 [w2,period2,scale2] = wt([tvec',Lseries(j,:)'],'S0',4*tslength,'J1',0,'Dj',1);
%                 [Xwt,periodx,scalex] = xwt([tvec',Lseries(i,:)'],...
%                 [tvec',Lseries(j,:)'],'S0',4*tslength,'J1',0,'Dj',1);