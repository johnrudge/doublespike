function [measuredv standardv spikev]=monterun(element,prop,spike,alpha,beta,n)
%MONTERUN Generate a fake mass spectrometer run by Monte-Carlo simulation
%    [measuredv standardv spikev]=MONTERUN(element,prop,spike,alpha,beta,n)
%             element -- element used in double spike, e.g. 'Fe'
%             prop -- proportion of spike in double spike-sample mix.
%                A vector of values can be specified if desired, to reflect changes over a run.
%             spike -- the isotopic composition of the double spike e.g. [0 0.5 0 0.5]
%                corresponds to a 50-50 mixture of the 2nd and 4th isotopes
%                (56Fe and 58Fe) in the case of Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (instrumental and natural, or alpha and beta). Values of beta and
%                alpha can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%                A vector of values can be specified if desired, to reflect changes over a run.
%             n -- number of Monte-Carlo samples to take. Default is 1000.
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% This function produces a fake mass spectrometer run using the error model specified in 
% ISODATA.(element).errormodel. The output is given as ion beam intensities for measured,
% standard, and double spike.
%
% Example
%    measured=monterun('Fe',0.5,[0 0 0.5 0.5]);
%
% See also dsinversion
global ISODATA

% Set some defaults
if isempty(ISODATA)
	dsstartup;
end
if (nargin<6) || (isempty(n))
	n=1000;
end
if (nargin<5) || (isempty(beta))
	beta=0;
end
if (nargin<4) || (isempty(alpha))
	alpha=0;
end

rawdata=ISODATA.(element);
standard=rawdata.standard;
mass=rawdata.mass;
emodel=rawdata.errormodel;
spike=spike/sum(spike);
nisos=ISODATA.(element).nisos;

if (length(alpha)==1)
	alpha=repmat(alpha,[n 1]);
end
if (length(beta)==1)
	beta=repmat(beta,[n 1]);
end
if (length(prop)==1)
	prop=repmat(prop,[n 1]);
end

% This code needs vectorising...
% calculate expected sample composition
sample=zeros(n,nisos);
mixture=zeros(n,nisos);
measured=zeros(n,nisos);
for i=1:n
	sample(i,:)=standard.*exp(-log(mass).*alpha(i));
	sample(i,:)=sample(i,:)./sum(sample(i,:));
	
	mixture(i,:)=prop(i).*spike + (1-prop(i)).*sample(i,:);
	
	measured(i,:)=mixture(i,:).*exp(log(mass).*beta(i));
	measured(i,:)=measured(i,:)/sum(measured(i,:));

	measuredi(i,:)=measured(i,:).*emodel.measured.intensity; % scale up to V
	measuredivar(i,:)=emodel.measured.a + emodel.measured.b.*measuredi(i,:) + emodel.measured.c.*(measuredi(i,:).^2);
	measuredicov(:,:,i)=diag(measuredivar(i,:));

	standardi(i,:)=standard.*emodel.standard.intensity; % scale up to V
	standardivar(i,:)=emodel.standard.a + emodel.standard.b.*standardi(i,:) + emodel.standard.c.*(standardi(i,:).^2);
	standardicov(:,:,i)=diag(standardivar(i,:));

	spikei(i,:)=spike.*emodel.spike.intensity; % scale up to V
	spikeivar(i,:)=emodel.spike.a + emodel.spike.b.*spikei(i,:) + emodel.spike.c.*(spikei(i,:).^2);
	spikeicov(:,:,i)=diag(spikeivar(i,:));
end

measuredv=mvnrnd(measuredi,measuredicov,n);
standardv=mvnrnd(standardi,standardicov,n);
spikev=mvnrnd(spikei,spikeicov,n);
