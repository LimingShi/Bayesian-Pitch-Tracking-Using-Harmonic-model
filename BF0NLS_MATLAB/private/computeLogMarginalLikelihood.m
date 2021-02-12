function logMarginalLikelihood = ...
        computeLogMarginalLikelihood(logPitchLikelihood, fullPitchGrid)
    % compute the unnormalised log marginal likelihood
    deltaFreq = 1;
    maxModelOrder = size(logPitchLikelihood,1);
    logMarginalLikelihood = ones(1,maxModelOrder);
    for iOrder = 1:maxModelOrder
        % find the non-NaN indices
        indices = find(~isnan(logPitchLikelihood(iOrder,:)));
        % find the maximum value
        maxLogMarginalLikelihood = max(logPitchLikelihood(iOrder,indices));
        logMarginalLikelihood(iOrder) = maxLogMarginalLikelihood+...
            log(deltaFreq*sum(exp(...
            logPitchLikelihood(iOrder,indices)-maxLogMarginalLikelihood)));
        
%         log(trapz(exp(...
%             logPitchLikelihood(iOrder,indices)-maxLogMarginalLikelihood)))+log(deltaFreq)+maxLogMarginalLikelihood
%         log(sum(  (exp(...
%             logPitchLikelihood(iOrder,indices)-maxLogMarginalLikelihood))))+log(deltaFreq)+maxLogMarginalLikelihood
    end
end