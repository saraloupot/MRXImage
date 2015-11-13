%% Calculate the B field
if(~exist('B','var'))
    B_ = forwardproblem(sourcestrength,sourcelocation,rdet);
    B = B_(:,3);
end
%% define snr
if(~exist('MRXsnr','var'))
    MRXsnr = [];
end
%% Add noise
if(~isempty(MRXsnr));
    B= addAWGN(B,MRXsnr);
end
%
%% exclude data
if(length(B)>size(rdet,1))
    B = B(use);
end
Bfake = [B; 0; 0; 0; 0];
