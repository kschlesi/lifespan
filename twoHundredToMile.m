function [mMin,mSec,tSec] = twoHundredToMile(time)

tSec = time*8;
mMin = floor(tSec/60);
mSec = mod(tSec,60);

end