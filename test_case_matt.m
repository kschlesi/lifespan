function pmat = test_case_matt()

nlayers = 2;
ndims = 2;
nstream = 20;
pmat = zeros(nlayers,ndims,2,nstream);

% add all specific data
pmat(:,:,1,1) = [3,1;2,3]; pmat(:,:,2,1) = [1,11;5,7];
pmat(:,:,1,2) = [2,1;1,4]; pmat(:,:,2,2) = [1,10;4,7];
pmat(:,:,1,3) = [2,2;2,4]; pmat(:,:,2,3) = [1,8;4,6];
pmat(:,:,1,4) = [2,3;2,4]; pmat(:,:,2,4) = [1,8;4,6];
pmat(:,:,1,5) = [3,2;3,4]; pmat(:,:,2,5) = [1,8;4,6];
pmat(:,:,1,6) = [3,3;3,4]; pmat(:,:,2,6) = [2,9;5,6];
pmat(:,:,1,7) = [3,1;2,3]; pmat(:,:,2,7) = [2,10;5,6];
pmat(:,:,1,8) = [3,2;3,3]; pmat(:,:,2,8) = [3,8;5,6];
pmat(:,:,1,9) = [4,4;3,4]; pmat(:,:,2,9) = [10,9;7,4];
pmat(:,:,1,10) = [4,4;4,4]; pmat(:,:,2,10) = [10,8;7,4];
pmat(:,:,1,11) = [5,5;4,4]; pmat(:,:,2,11) = [10,7;7,3];
pmat(:,:,1,12) = [6,5;5,4]; pmat(:,:,2,12) = [11,7;7,3];
pmat(:,:,1,13) = [5,4;4,4]; pmat(:,:,2,13) = [10,1;4,1];
pmat(:,:,1,14) = [6,4;4,3]; pmat(:,:,2,14) = [11,1;4,1];
pmat(:,:,1,15) = [6,5;4,3]; pmat(:,:,2,15) = [11,2;5,1];
pmat(:,:,1,16) = [10,2;5,1]; pmat(:,:,2,16) = [3,9;5,6];
pmat(:,:,1,17) = [9,2;4,1]; pmat(:,:,2,17) = [3,10;5,6];
pmat(:,:,1,18) = [10,3;5,1]; pmat(:,:,2,18) = [3,10;6,6];
pmat(:,:,1,19) = [11,3;5,1]; pmat(:,:,2,19) = [3,11;6,6];
pmat(:,:,1,20) = [11,8;7,4]; pmat(:,:,2,20) = [3,11;6,6];


end