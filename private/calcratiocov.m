function V=calcratiocov(composition,errormodel,di,isonorm,prop)
% Calculate the covariance matrix of the ratios based on the given error model and composition
% di is the isotope with which to denominator
% isonorm are the isotopes to use in the normalisation
% prop is the proportion of spike in the spike-sample mix

if (nargin<5) || isempty(isonorm)
	prop = 0.0  % prop is only used when the error model is fixed-sample
end

if (nargin<4) || isempty(isonorm)
	isonorm=1:size(composition,2);   % by default use all isotopes for normalisation (same amount of element run each time)
end

% first normalise composition so it is really a composition (unit sum)
composition=composition./sum(composition);

meanbeams=composition.*errormodel.intensity./sum(composition(isonorm));

if isequal(errormodel.type, 'fixed-sample')
   meanbeams = meanbeams/(1.0 - prop);
end

covbeams=calcbeamcov(meanbeams,errormodel);
V=covbeamtoratio(meanbeams,covbeams,di);

function beamcov=calcbeamcov(meanbeams,errormodel)
% the beam covariance matrix
beamvar=errormodel.a + meanbeams.*errormodel.b + (meanbeams.^2).*errormodel.c;
beamcov=diag(beamvar);

function V=covbeamtoratio(meanbeams,covbeams,di)
% converts a covariance matrix for beams to one for ratios
% di is the isotope to denominator with
% assumes last row and column of M correspond to denominator
isonums=1:length(meanbeams);
ni=isonums(isonums~=di);

n=meanbeams(ni);
d=meanbeams(di);

M=covbeams([ni di],[ni di]);  % move denominator to end
%A=[diag(repmat(1/d,1,length(n))) -n'./(d^2)];
%A=[(1/d).*eye(length(n)) -n'./(d^2)];
A=[diag((1/d).*ones(1,length(n))) -n'./(d^2)];
V=(A*M)*(A');
