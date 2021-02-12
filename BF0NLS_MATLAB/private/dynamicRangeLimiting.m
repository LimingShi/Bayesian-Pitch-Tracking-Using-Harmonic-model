%% dynamicRangeLimiting
% Limit the dynamic range of nonnegative data
% 
%% Syntax:
% limitedNonnegativeData = ...
%     dynamicRangeLimiting(nonnegativeData, maxRangeDb)
%
%% Description:
% Limit the dynamic range of nonnegative data. This is usefull for plotting
% purposes, e.g., the the spectrogram. The data in the returned vector
% satifies that
% $$maxRangeDb >= 10*log10(max(limitedNonnegativeData)/...
%     min(limitedNonnegativeData))$$
% 
% * nonnegativeData: A vector of nonnegative data
% * maxRangeDb: The maximum dynamic range in dB
% * limitedNonnegativeData: The input vector with a limited dynamic range
%
%% Examples:
% maxRangeDb = 60 dB;
% limitedMagnitudeSpectrum = ...
%     dft2SymmetricDft(abs(fft(data)).^2, maxRangeDb);
%
function limitedNonnegativeData = ...
        dynamicRangeLimiting(nonnegativeData, maxRangeDb)
    if min(nonnegativeData) < 0
        error('dynamicRangeLimiting:argChk', ...
            'All data must be nonnegative');
    end
    logPowerSpectrum = 10*log10(nonnegativeData);
    limitedLogPowerSpectrum = max(logPowerSpectrum, ...
            max(max(logPowerSpectrum)-abs(maxRangeDb)));
    limitedNonnegativeData = 10.^(limitedLogPowerSpectrum/10);
end