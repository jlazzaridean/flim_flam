function [result] = evalGaussFunc_noIntegral(tRange,tauBar,tauSig)
%log normal distribution
%1/sqrt(tau.^2*tauSig.^2*pi*2) * exp(-(log(tau)-tauBar).^2/(2*tauSig.^2))

% Make 1000 equally spaced lifetimes of each component
ts1 = linspace(1/1000,10,1000)'; %space of taus to evaluate over
%set up a log-normal distribution with these taus
%anonF = @(sig, bar, tau) (1/sqrt(tau.^2*sig.^2*pi*2) * exp(-(log(tau)-bar).^2/(2*sig.^2)));
P1 = lognpdf(ts1,tauBar,tauSig);

result = zeros(size(tRange,1),1);
%for each point in time, determine the contribution from all taus
for i=1:size(tRange,1)
    for j=1:size(ts1,1)
        result(i,1) = result(i,1) + exp(-tRange(i,1)/ts1(j,1))*P1(j,1);
    end
end

end

