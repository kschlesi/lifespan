function [rh_set,lh_set,rh_map,lh_map,dfilledR,dfilledL] = ME2dset(MEdata,tag)
% this function converts a vector of mary-ellen hybrid atlas data to a
% .dset files usable in AFNI Suma to plot the data voxel-wise on the brain.
%
% INPUT:  MEdata, a 194-element vector of numerical color axis values
%                 to be plotted
%         tag,    a string to be used in the .dset filename
% OUTPUT: lh_set, a 155029x1 vector containing the voxelwise-mapped values 
%                 for the left hemisphere
%         rh_set, a 155100x1 vector containing the voxelwise-mapped values 
%                 for the right hemisphere
%
% NOTE: this function will create new files with the following names in the
%       runtime directory: 'MEhybrid.lh.<tag>.dset'
%                          'MEhybrid.rh.<tag>.dset'
%       this function requires the following files to be present on the
%       matlab path: 'hybrid_data.txt'
%                    'Lausanne2008.lh.scale1251D.dset'
%                    'Lausanne2008.lh.scale60.1D.dset'
%                    'Lausanne2008.lh.scale36.1D.dset'
%                    'Lausanne2008.rh.scale1251D.dset'
%                    'Lausanne2008.rh.scale60.1D.dset'
%                    'Lausanne2008.rh.scale36.1D.dset'

% check/fix dimensions of MEdata
assert(numel(MEdata)==194); 
MEdata = MEdata(:);

% read in voxel mappings to Lausanne2008 regions at relevant scales
laus_125_lh = importdata('Lausanne2008.lh.scale125.1D.dset');
laus_125_rh = importdata('Lausanne2008.rh.scale125.1D.dset');
laus_60_lh = importdata('Lausanne2008.lh.scale60.1D.dset');
laus_60_rh = importdata('Lausanne2008.rh.scale60.1D.dset');
laus_33_lh = importdata('Lausanne2008.lh.scale36.1D.dset');
laus_33_rh = importdata('Lausanne2008.rh.scale36.1D.dset');

% import map info and initialize lh_set, rh_set 
map = dlmread('hybrid_data.txt',' ');
rh_map = zeros(155100,3);
lh_map = zeros(155029,3);
rh_set = zeros(155100,1);
lh_set = zeros(155029,1);
%nfilled = zeros(size(MEdata));
%dupfilled = zeros(size(MEdata));

% map right-hemisphere ME regions
for i=1:90  
    res = map(map(:,3)==i,2);
    regID = map(map(:,3)==i,1);
  % apply subcortical region cutoff
    %if regID <= cutoff(res,1)%; regID = cutoff(res,1); end;
      % map into appropriately labeled entries of rh_map
        %eval(['nfilled(i) = sum(laus_' num2str(res) '_rh==regID);']);
        eval(['dfilledR(i) = sum(sum(rh_map(laus_' num2str(res) '_rh==regID,:),2)~=0);']);
        eval(['rh_map(laus_' num2str(res) '_rh==regID,resix(res)) = MEdata(i);']);
    %end
end

% map left-hemisphere ME regions
for i=97:187 
    res = map(map(:,3)==i,2);
    regID = map(map(:,3)==i,1) - offset(res);
  % apply subcortical region cutoff
    %if regID <= cutoff(res,2)%; regID = cutoff(res,2); end;
      % apply hemisphere offset
        %regID = regID - offset(res); 
      % map into appropriately labeled entries of lh_set
        %eval(['nfilled(i) = sum(laus_' num2str(res) '_lh==regID);']);
        eval(['dfilledL(i) = sum(sum(lh_map(laus_' num2str(res) '_lh==regID,:),2)~=0);']);
        eval(['lh_map(laus_' num2str(res) '_lh==regID,resix(res)) = MEdata(i);']);
    %end
end

% resolve conflicts
twoixR = find(sum((rh_map>0),2)==2);
twoixL = find(sum((lh_map>0),2)==2);
%dfilled = numel(twoixR) + numel(twoixL);
oneixR = find(sum((rh_map>0),2)==1);
oneixL = find(sum((lh_map>0),2)==1);
%nfilled = dfilled + numel(oneixR) + numel(oneixL);
rh_set(twoixR) = rh_map(twoixR,2);
rh_set(oneixR) = sum(rh_map(oneixR,:),2);
lh_set(twoixL) = lh_map(twoixL,2);
lh_set(oneixL) = sum(lh_map(oneixL,:),2);

% save as .dset file
dlmwrite(['MEhybrid.rh.' tag '.dset'], rh_set);
dlmwrite(['MEhybrid.lh.' tag '.dset'], lh_set);

end

function offset = offset(res)
% set ID offset for left hemisphere Lausanne regions
    switch res
        case 33,  offset = 41;
        case 60,  offset = 64;
        case 125, offset = 115;
    end
end

function cutoff = cutoff(res,hem)
% set ID cutoff for subcortical Lausanne regions 
% hem=1 -> right, hem=2 -> left
    switch res
        case 33,  cutoff = [35,76];
        case 60,  cutoff = [58,122];
        case 125, cutoff = [109,227];
    end
    cutoff = cutoff(hem);
end

function ix = resix(res)
% give index of resolution atlas
    switch res
        case 33,  ix = 1;
        case 36,  ix = 1;
        case 60,  ix = 2;
        case 125, ix = 3;
    end
end

%[map(d_test~=0,[1,2]),n_test(d_test~=0),d_test(d_test~=0)]