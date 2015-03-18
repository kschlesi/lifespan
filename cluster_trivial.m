% A simple test script that uses the default parallel configuration.
tic;
job = createJob('Name','rand_try');
% The following adds the current directory to the Matlab path for
% cluster workers:
my_pwd = pwd;
set(job,'PathDependencies',{my_pwd});
%
numtasks = 10;
for tasknum = 1:numtasks
  createTask(job,@rand_try,1,{4,4});
end
submit(job);
waitForState(job,'finished');
results = getAllOutputArguments(job);
for tasknum = 1:numtasks
  disp(results{tasknum});
end
destroy(job);
toc;
