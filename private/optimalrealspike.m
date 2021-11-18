function [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalrealspike(element,beta,alpha,errorratio,isospike,isoinv)
%OPTIMALREALSPIKE    Finds the best real spike
%    OPTIMALREALSPIKE(element,beta,alpha,errorratio,isospike,isoinv)
%             element -- the particular element
%             beta -- instrumental fractionation
%             alpha -- natural fractionation
%             errorratio -- the ratio whose error we are targeting
%             isospike -- the isotopes to spike
%             isoinv -- the isotopes used in the inversion
global ISODATA
rawdata=ISODATA.(element);

% Have some default arguments
if (nargin<6) || isempty(isoinv)
	isoinv=[];
end
if (nargin<5) || isempty(isospike)
	isospike=[];
end
if (nargin<4) || isempty(errorratio)
	errorratio=[];
end
if (nargin<3) || isempty(alpha)
	alpha=0;
end
if (nargin<2) || isempty(beta)
	beta=0;
end

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isospike=rawdata.isoindex(isospike);
isoinv=rawdata.isoindex(isoinv);

if (isempty(isoinv))
	isoinv=combnk(1:rawdata.nisos,4);
end

if isempty(isospike)
	isospikev=combnk(1:rawdata.nspikes,2);
else
	isospikev=isospike;
end

isoinvvals=[];            % inverse values to check
isospikevals=[];          % spike values to check
for (i=1:size(isoinv,1))
	isospikevals=[isospikevals; isospikev];
	isoinvvals=[isoinvvals; repmat(isoinv(i,:),size(isospikev,1),1)];
end

for (i=1:size(isoinvvals,1))
	[optspike(i,:),optprop(i,:),opterr(i,:),optspikeprop(i,:),optppmperamu(i,:)]=singlerealoptimalspike(element,beta,alpha,errorratio,isospikevals(i,:),isoinvvals(i,:));
end
optisoinv=isoinvvals;

% restrict to the top six to save masses of output...
indx=[];
for (i=1:size(isoinv,1))
	nspikechoice=size(isospikev,1);
	oe=opterr(nspikechoice*(i-1)+1:nspikechoice*i,:);
	[oe,ix]=sort(oe);
	if (nspikechoice>6)
		ix=ix(1:6,:);
	end
	indx=[indx; (nspikechoice*(i-1)+ix)];
end
optspike=optspike(indx,:);
optppmperamu=optppmperamu(indx,:);
optprop=optprop(indx,:);
opterr=opterr(indx,:);
optspikeprop=optspikeprop(indx,:);
optisoinv=optisoinv(indx,:);

% Sort in ascending order of error
[opterr,ix]=sort(opterr);
optppmperamu=optppmperamu(ix,:);
optspike=optspike(ix,:);
optprop=optprop(ix,:);
optisoinv=optisoinv(ix,:);
optisoinv=rawdata.isonum(optisoinv);
optspikeprop=optspikeprop(ix,:);

function [optspike,optprop,opterr,optspikeprop,optppmperamu]=singlerealoptimalspike(element,beta,alpha,errorratio,isospike,isoinv)
% Calculate the composition of the optimal double spike given the isotopes used in the inversion
% and of those the isotopes we are spiking
global ISODATA
rawdata=ISODATA.(element);

spikevector1=rawdata.rawspike(isospike(1),:);
spikevector2=rawdata.rawspike(isospike(2),:);

if (verLessThan('optim','4.0'))
	options=optimset('Display','notify','TolX',1e-8,'TolFun',1e-10,'TolCon',1e-6,'LargeScale','off','MaxFunEvals',10000);
else
	options=optimset('Display','notify','TolX',1e-8,'TolFun',1e-10,'TolCon',1e-6,'Algorithm','active-set','MaxFunEvals',10000);
end

tol=1e-5;    % how close can we be to the edges of the domain?
lb=[tol;
	tol];

ub=[1-tol;
	1-tol];

y0=[0.5 0.5]';

% Helpful to rescale the error, to make everything roughly order 1 for the optimiser
initialerror=errorestimate(element,y0(1),y0(2).*spikevector1+(1-y0(2)).*spikevector2,isoinv,errorratio,beta,alpha);

[y,opterr,exitflag] = fmincon(@(y) errorestimate(element,y(1),y(2).*spikevector1 +(1-y(2)).*spikevector2,isoinv,errorratio,beta,alpha)./initialerror,y0,[],[],[],[],lb,ub,[],options);
[opterr,optppmperamu]=errorestimate(element,y(1),y(2).*spikevector1 +(1-y(2)).*spikevector2,isoinv,errorratio,beta,alpha);

optprop=y(1);
optspike=y(2)*spikevector1 +(1-y(2))*spikevector2;

optspikeprop=zeros(rawdata.nspikes,1);
optspikeprop(isospike(1))=y(2);
optspikeprop(isospike(2))=1-y(2);
%optspikeprop=[y(2) 1-y(2)];