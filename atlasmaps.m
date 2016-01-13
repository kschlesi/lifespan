% load mary-ellen atlas labels
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
[~,laus_labels_60] = read_regions('ParcellationLausanne2008.csv',60);
[~,laus_labels_125] = read_regions('ParcellationLausanne2008.csv',125);
[~,laus_labels_33] = read_regions('ParcellationLausanne2008.csv',33);

% load map info
mapinfo = dlmread('hybrid_data.txt',' ');

% check whether mapping is accurate
mapcheck = zeros(length(hybrid_labels),1);
for i=1:length(hybrid_labels)
    eval([ 'reslabels = laus_labels_' num2str(mapinfo(mapinfo(:,3)==i,2)) ';' ]);
    mapcheck(i) = strcmp(hybrid_labels{i},reslabels(mapinfo(mapinfo(:,3)==i,1)));
    if ~mapcheck(i)
        disp([hybrid_labels{i},reslabels(mapinfo(mapinfo(:,3)==i,1))]);
    end
end