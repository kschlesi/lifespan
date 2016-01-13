% resolving conflicts...
hybrid_labels_long = importdata('hybrid_labels_long.txt');
hybrid_labels = importdata('hybrid_labels_long.txt');
% re-format hybrid labels to match Lausanne labels
for i=1:size(hybrid_labels, 1);
    if ~isletter(hybrid_labels{i}(end))
        hybrid_labels{i} = hybrid_labels{i}(1:end-3);
    else
        hybrid_labels{i} = [hybrid_labels{i} '_1'];
    end
end

% load Lausanne atlas labels
laus_125_lh = importdata('Lausanne2008.lh.scale125.1D.dset');
laus_125_rh = importdata('Lausanne2008.rh.scale125.1D.dset');
laus_60_lh = importdata('Lausanne2008.lh.scale60.1D.dset');
laus_60_rh = importdata('Lausanne2008.rh.scale60.1D.dset');
laus_33_lh = importdata('Lausanne2008.lh.scale36.1D.dset');
laus_33_rh = importdata('Lausanne2008.rh.scale36.1D.dset');
[~,laus_labels_60] = read_regions('ParcellationLausanne2008.csv',60);
[~,laus_labels_125] = read_regions('ParcellationLausanne2008.csv',125);
[~,laus_labels_33] = read_regions('ParcellationLausanne2008.csv',33);

% load map info
map = dlmread('hybrid_data.txt',' ');
[rh_set,lh_set,rh_map,lh_map] = ME2dset(1:194,'test');

% identify voxels with conflicts
twoixR = find(sum((rh_map>0),2)==2);
twoixL = find(sum((lh_map>0),2)==2);
oneixR = find(sum((rh_map>0),2)==1);
oneixL = find(sum((lh_map>0),2)==1);

% identify Lausanne regions & labels of those voxels
r125R = laus_125_rh(twoixR);
r60R = laus_60_rh(twoixR);
r125L = laus_125_lh(twoixL) + offset(125);
r60L = laus_60_lh(twoixL) + offset(60);
b125R = laus_labels_125(r125R);
b60R = laus_labels_60(r60R);
b125L = laus_labels_125(r125L);
b60L = laus_labels_60(r60L);

% create final display case
regnsR = [r60R, r125R, zeros(size(r60R)), zeros(size(r60R))];
lablsR = [b60R, b125R, cell(size(b60R)), cell(size(b60R))];
regnsL = [r60L, r125L, zeros(size(r60L)), zeros(size(r60L))];
lablsL = [b60L, b125L, cell(size(b60L)), cell(size(b60L))];

% find conflicting hybrid regions & labels
for i=1:numel(r60R)
    h60 = map(~~((map(:,1)==regnsR(i,1)).*(map(:,2)==60)),3);
    h125 = map(~~((map(:,1)==regnsR(i,2)).*(map(:,2)==125)),3);
    regnsR(i,3) = h60; lablsR(i,3) = hybrid_labels(h60);
    regnsR(i,4) = h125; lablsR(i,4) = hybrid_labels(h125);
    h60 = map(~~((map(:,1)==regnsL(i,1)).*(map(:,2)==60)),3);
    h125 = map(~~((map(:,1)==regnsL(i,2)).*(map(:,2)==125)),3);
    regnsL(i,3) = h60; lablsL(i,3) = hybrid_labels(h60);
    regnsL(i,4) = h125; lablsL(i,4) = hybrid_labels(h125);
end



voxR = zeros(numel(laus_60_rh),2);
voxL = zeros(numel(laus_60_lh),2);
for v=1:numel(laus_60_rh)
    % is the voxel assigned at all
    voxR(v,1) = ~~laus_60_rh(v);
    % is the voxel's assignment used in map (& how many times)
    voxR(v,2) = ismember(laus_33_rh(v),map(map(:,2)==33,1)) + ...
                ismember(laus_60_rh(v),map(map(:,2)==60,1)) + ...
                ismember(laus_125_rh(v),map(map(:,2)==125,1));
  if v<=numel(laus_60_lh)
    voxL(v,1) = ~~laus_60_lh(v);
    if voxL(v,1)
    voxL(v,2) = ismember(laus_33_lh(v)+offset(33),map(map(:,2)==33,1)) + ...
                ismember(laus_60_lh(v)+offset(60),map(map(:,2)==60,1)) + ...
                ismember(laus_125_lh(v)+offset(125),map(map(:,2)==125,1));
    end
  end
end

empty=sum(voxR(:,2)==0)+sum(voxL(:,2)==0);
empty2=sum(lh_set==0)+sum(rh_set==0);