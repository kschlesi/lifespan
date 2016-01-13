% age data processing: read in raw denoised .mat files and create
% timeseries struct (similar to officer timeseries struct)

inpath = '/Users/kimberly/Documents/age_data/';
outpath = '/Users/kimberly/Documents/lifespan/';

% read in .mat files from inpath
fid = fopen([inpath 'age_files.txt']);
files = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
files = files{:};
nsubjs = length(files);

% create struct with 'run' fields and subject ID field
field0 = 'ID'; value = cell(nsubjs,1);
field1 = 'run1'; field2 = 'run2'; field3 = 'run3'; 
age_timeseries = struct(field0,value,field1,value,field2,value,...
                    field3,value);

% useful numbers & preallocations       
load([inpath files{1}]);
nruns = length(tsdata);
[n,ts] = size(tsdata{1});
for k=1:nsubjs  % loop over subject files
    load([inpath files{k}]);
    for r=1:nruns
        % assign timeseries to struct
        eval(['age_timeseries(k).run' num2str(r) ' = tsdata{r};']);
    end
    age_timeseries(k).ID = files{k}(end-8:end-4);
end

save([outpath 'age_timeseries1.mat'],'age_timeseries');