%%% CREATE AND SAVE NODE_NODE CORRELATION MATRICES FOR OFFICER DATA
%%% (ONE PER TASK, PER SUBJECT)

load('age_timeseries.mat'); % struct array containing all data
timeseries = age_timeseries;
clear age_timeseries;
n = size(timeseries(1).run1,1); % number of nodes (brain regions)
nsubjs = length(timeseries);
nruns = 3;
tw = size(timeseries(1).run1,2); % number of TRs per time window
ib = (1:nsubjs)';
bandpass = [0.06,0.125]; % in Hz
file_path = '/Users/kimberly/Documents/lifespan/NNs';% string giving 
                                    % location to save node-node matrices

% OVERRIDE ib if desired to restrict to only a few subjects
ib = [2;68];

for k=1:numel(ib)
    person = ib(k);

    for r = 1:nruns
    Lrun1 = timeseries(k).run1(:,:);
    Lrun1 = timeseries(k).run1(:,:);
    Lrun1 = timeseries(k).run1(:,:);

disp([num2str(person) ' REST']); tic;
    [Nrest,Prest] = waveletter(Lrest,'cohr',bandpass,2,1); % 2 = TR in sec; 1 = p-val threshold
    dlmwrite([file_path 'Nrest_' num2str(person) '.dat'],Nrest);
    dlmwrite([file_path 'Prest_' num2str(person) '.dat'],Prest);
    toc;
disp([num2str(person) ' ATTN']); tic;
    [Nattn,Pattn] = waveletter([Latt1,Latt2],'cohr',bandpass,2,1); % 2 = TR in sec
    dlmwrite([file_path 'Nattn_' num2str(person) '.dat'],Nattn);
    dlmwrite([file_path 'Pattn_' num2str(person) '.dat'],Pattn);
    toc;
disp([num2str(person) ' WORD']); tic;
    [Nword,Pword] = waveletter(Lword,'cohr',bandpass,2.5,1); % 2.5 = TR in sec
    dlmwrite([file_path 'Nword_' num2str(person) '.dat'],Nword);
    dlmwrite([file_path 'Pword_' num2str(person) '.dat'],Pword);
    toc;
disp([num2str(person) ' FACE']); tic;
    [Nface,Pface] = waveletter(Lface,'cohr',bandpass,2.5,1); % 2.5 = TR in sec
    dlmwrite([file_path 'Nface_' num2str(person) '.dat'],Nface);
    dlmwrite([file_path 'Pface_' num2str(person) '.dat'],Pface);
    toc;
    
end  % end loop over subjects
