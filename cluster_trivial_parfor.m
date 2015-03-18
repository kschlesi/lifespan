% A simple test script that uses the default parallel configuration.
%
%  Using matlabpool and parfor
%
tic;
matlabpool OPEN 13;
numtasks = 25;
outputarray = cell(numtasks,1);
parfor tasknum = 1:numtasks
    [status, outputarray{tasknum}] = system('hostname');
    system('sleep 10');
end
matlabpool CLOSE;
for tasknum = 1:numtasks
    disp(outputarray{tasknum});
end
toc;
