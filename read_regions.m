function [theRegions,theLabels,nRightHem] = read_regions(filename,atlasID)

fid = fopen(filename);

colHeads = textscan(fid,'%s',20,'HeaderLines',1,'Delimiter',',');
colHeads = colHeads{:};
idix = zeros(size(colHeads));
labelix = zeros(size(colHeads));
for i=1:length(colHeads)
    idix(i) = strcmp(colHeads{i},['scale' num2str(atlasID) ' ID']);
    labelix(i) = strcmp(colHeads{i},['scale' num2str(atlasID) ' labels']);
end

theSheet = textscan(fid,repmat('%n %s ',1,10),'Delimiter',',');
fclose(fid);

lHemLabels = theSheet{find(labelix>0,1,'last')};
lHemLabels = lHemLabels(~strcmp(lHemLabels,''));
rHemLabels = theSheet{find(labelix>0,1,'first')};
rHemLabels = rHemLabels(~strcmp(rHemLabels,''));

lHemIDs = theSheet{find(idix>0,1,'last')};
lHemIDs = lHemIDs(~isnan(lHemIDs));
rHemIDs = theSheet{find(idix>0,1,'first')};
rHemIDs = rHemIDs(~isnan(rHemIDs));
nRightHem = size(rHemIDs,1);

theRegions = [rHemIDs;lHemIDs];
theLabels = cell(size(theRegions));
 for i=1:nRightHem
     theLabels{i} = ['R_' rHemLabels{i}];
 end
 for i=nRightHem+1:size(theRegions,1)
     theLabels{i} = ['L_' lHemLabels{i-nRightHem}];
 end
 for i=1:size(theLabels,1)
    if isletter(theLabels{i}(end))
        theLabels{i} = [theLabels{i} '_1'];
    end
 end
 theLabels = deblank(theLabels);
 
% theRegions = [rHemIDs,lHemIDs];
% theLabels = [rHemLabels,lHemLabels];

end

function new_str = text_replace(in_str,find_exp,sub_exp)



end