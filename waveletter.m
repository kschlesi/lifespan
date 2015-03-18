function [Nseries,Pseries,Wseries] = waveletter(chunkmatrix,adjscheme,bandpass,TR,pcutoff)
% extracts wavelet coefficients and creates a node-node adjacency matrix
% INPUT:  chunkmatrix, a nxts matrix: n raw timeseries of ts timesteps each
%         adjscheme, a string: if 'corr', use Pearson correlation
%                              if 'cohr', use average wavelet coherence
%         bandpass, a 2-element vector: frequency range of interest, in Hz
%         TR, a scalar: length of sampling period (repetition time) in sec
%         pcutoff, a scalar: threshold p-value, no threshold if = 1
% OUTPUT: Nseries, a symmetric nxn matrix of node-node correlations
%         Pseries, a symmetric nxn matrix of p-values for those correlations
%         Wseries, an nx(unknown) matrix of 2nd-level wavelet coefficients
%                   for each region
%
% uses Grinsted et al. "wtc" package (http://noc.ac.uk/using-science/crosswavelet-wavelet-coherence)
% uses Nichols "FDR" software (http://www-personal.umich.edu/~nichols/FDR)

[n,ts]=size(chunkmatrix);

if nargin<5
    pcutoff = 1; % default: no thresholding
end

if bandpass(1)>=bandpass(2) || numel(bandpass)~=2
    error('wrongly formatted frequency band input!');
end

switch adjscheme
    case 'corr', % uses Pearson correlation
        
        % code is not modified to take argument for wavelet scale
        warning('default: scale two wavelet coefficients');
        
        % extract wavelet coefficients (level 2, daubechies 4 discrete):
        disp('wavelet extraction...');
        tic;
        for i=1:n % (needs MATLAB wavelet toolbox)
            w = ndwt(chunkmatrix(i,:)-mean(chunkmatrix),3,'db4');
            Wseries(i,:) = w.dec{3}; 
            clear w;
        end
        toc;
        
        % compute node-node Pearson correlations & p-values
        disp('node-node correlation matrices...');
        tic;
        [Nseries,Pseries] = corrcoef(Wseries');
        % force symmetric matrix
        Nseries = triu(Nseries)+triu(Nseries,1)'; 
        Pseries = triu(Pseries)+triu(Pseries,1)';
        toc;
        
        % threshold with FDR (needs FDR package)
        if pcutoff<1 && pcutoff >=0
            [pID,~] = FDR(reshape(Pseries,n*n,1),pcutoff);
            Nseries = Nseries.*(Pseries<pID);
        else
            disp('No FDR threshold applied');
        end
        
    case 'cohr', % uses average wavelet coherence
        
        % extracting wavelet coherence (needs Grinsted wtc package)
        disp('wavelet extraction...');
        Nseries = zeros(n,n);
        tic;
        for i=1:n
            disp(i/n);
            Nseries(i,i) = 0;
            for j=i+1:n
                [Rsq,period] = wtc(chunkmatrix(i,:)',chunkmatrix(j,:)','mcc',0);
                freq = 1./(period*TR);
                % average over time and over frequencies of interest
                Nseries(i,j) = mean(mean(Rsq((freq<bandpass(2)&freq>bandpass(1)),:)));
            end
        end
        Nseries = Nseries + Nseries'; % force symmetric matrix
        toc;
        Wseries = [];
        
        % threshold (would need to compute significance with monte carlo)
        Pseries = [];
        disp('No FDR threshold applied');
        
end

end