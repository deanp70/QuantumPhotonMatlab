function fidelity = showFidelity(targetMat,compareMat)
%% Show Fidelity
% 
% SHOWFIDELITY shows how close the compareMat is to the targetMat, by
% showing the trace of their scalar multiplication.
%
%   SHOWFIDELITY(TARGETMAT,COMPAREMAT) - compare the compare matrix to the
%   target matrix.
    if (size(targetMat) ~= size(compareMat))
        error('Illegal input, matrix sizes must be same');
    end
    fidelity = trace(targetMat*compareMat);
end