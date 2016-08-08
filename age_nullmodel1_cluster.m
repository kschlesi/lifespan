% cluster_age_nullmodel1

% define cdINPATH AND nnINPATH.....
myPwd = pwd;
totalsubjs = 108;
missing_subjs = [46;64;81;82]; % final list
nruns = 3;
goRange = [1.2, 0.05]; ts = 52; stats_tag = '52';
ts_run = 316; tag = [];
NNinpath = [myPwd 'NNs/'];
CDinpath = [myPwd 'CD/'];

% load all connectivity matrices
p = 100;
[t,n,ib] = load_CD_results(totalsubjs,missing_subjs,nruns,ts,ts_run,goRange,tag,CDinpath);
A = load_adjs(t,ib,subj_IDs,ts,ts_run,NNinpath);
ix = find(triu(ones(n),1));     % find triu indices

tic;
matlabpool OPEN 13;
numtasks = numel(ib)*p;
outputarray = cell(numtasks,1);
parfor tasknum = 1:numtasks
    k = floor((tasknum-1)/p)+1;
    px = mod(tasknum-1,p)+1;
    output = nullmodel1_wrapper(myPwd,px,k,totalsubjs,missing_subjs,nruns,ts,ts_run,goRange,tag,p);
end
matlabpool CLOSE;