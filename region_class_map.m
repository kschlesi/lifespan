function [regClass,regCN] = region_class_map(regName,filepath)
% assigns region classes (both named and numerical) to each region in the
% nx1 cell of strings 'regNames' that is mappable by name to a
% Lausanne-identified region from the 83 atlas. region classes are taken
% from http://arxiv.org/pdf/1406.5197v1.pdf

[~,classMap,~] = xlsread([filepath 'ParcellationLausanne2008.xls'],'CLASS','A1:B43');
classNum = xlsread([filepath 'ParcellationLausanne2008.xls'],'CLASS','C1:C43');
regClass = cell(size(regName));
regCN = zeros(size(regName));
for i=1:size(classMap,1)
    ix = find(ismemvar(strfind(regName,classMap{i,1}),3));
    if numel(ix)
      for j=1:length(ix)
        regClass{ix(j)} = classMap{i,2};
        regCN(ix(j)) = classNum(i);
      end
    end
end

end
