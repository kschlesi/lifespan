function output = age_nodals_chunks(person,inpath,outpath,nruns,ts,...
                                    TR,bandpass,edgedef,pval)

%%% CREATE AND SAVE RUNWISE NODE-NODE CORRELATION MATRICES FOR CHOKING DATA

% %%%%%%%%%%%%%%%%%% MODIFY THIS PART %%%%%%%%%%%%%%%%%%%%%
% inpath = '/Users/kimberly/Documents/lifespan/'; % should end with /
% outpath = '/Users/kimberly/Documents/lifespan/NNs/'; % should end with /
% nruns = 3; % number of functional blocks / runs
% ts = 316; % number of TRs per time window
% TR = 2;  % repetition time (TR) in seconds
% bandpass = [0.06,0.125]; % in Hz
% edgedef = 'cohr'; % 'corr' if using correlation, 'cohr' if wavelet coherence
% pval = 1;  % whether to FDR-threshold (no threshold if pval = 1 or if edgedef = 'cohr')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load timeseries data struct
load([inpath 'age_timeseries.mat']);
timeseries = age_timeseries;
clear age_timeseries;
%n = size(timeseries(1).run1,1); % number of nodes (brain regions)
%nsubjs = length(timeseries);
%ib = (1:nsubjs)';  % vector containing indices of subjects to use
%disp(ib);

%%%%%%% for each subject k %%%%%%%

    disp(['SUBJECT ' num2str(person)]);
    
    for runnum=1:nruns
        disp(['run ' num2str(runnum)]);
        
        % load data (timeseries)
        eval(['Lseries = timeseries(person).run' num2str(runnum) ';']);

        % set Nts = number of timesteps in this functional run (task)
        % and t = number of chunks (time windows) in this run
        Nts = size(Lseries,2);
        t = floor(Nts/ts);
        
        % loops over each chunk within the functional run
        for i=1:t
            disp(['run ' num2str(runnum) ', ' num2str(i) ' of ' num2str(t)]);
            chunk = Lseries(:,((i-1)*ts)+1:i*ts);

            % send chunked matrices into waveletter
            [Nseries,~] = waveletter(chunk,edgedef,bandpass,TR,pval);

            % save nodenode matrix
            dlmwrite([outpath 'ageNN' num2str(person) '_run' num2str(runnum) '-' num2str(i) '_tw' num2str(ts) '.txt'],Nseries);
            clear Nseries;            

        end % end loop over chunks
    end % end loop over runs
    
    output = [num2str(person) ' matrix created'];
    
end  % function
