function [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalpurespike(element,beta,alpha,errorratio,isospike,isoinv)
%OPTIMALPURESPIKE    Finds the best pure spike
%    OPTIMALPURESPIKE(rawdata,beta,alpha,errorratio,isospike,isoinv)
%             rawdata -- data about a particular element
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

isoinvvals=[];
isospikevals=[];
for i=1:size(isoinv,1)
	if isempty(isospike)
		isospikev=combnk(isoinv(i,:),2);
	else
		if length(intersect(isospike,isoinv(i,:)))==2
			isospikev=isospike;
		else
			isospikev=[];
		end
	end
	isospikevals=[isospikevals; isospikev];
	isoinvvals=[isoinvvals; repmat(isoinv(i,:),size(isospikev,1),1)];
end

for i=1:size(isoinvvals,1)
	[optspike(i,:),optprop(i,:),opterr(i,:),optppmperamu(i,:)]=singlepureoptimalspike(element,beta,alpha,errorratio,isospikevals(i,:),isoinvvals(i,:));
end
optisoinv=isoinvvals;

% Sort in ascending order of error
[opterr,ix]=sort(opterr);
optppmperamu=optppmperamu(ix,:);
optspike=optspike(ix,:);
optprop=optprop(ix,:);
optisoinv=optisoinv(ix,:);
optisoinv=rawdata.isonum(optisoinv);
%optspikeprop=optspikeprop(ix,:);
optspikeprop=optspike;

function [optspike,optprop,opterr,optppmperamu]=singlepureoptimalspike(element,beta,alpha,errorratio,isospike,isoinv)
% Calculate the composition of the optimal double spike given the isotopes used in the inversion
% and of those the isotopes we are spiking
global ISODATA
rawdata=ISODATA.(element);

spikevector1=zeros(1,rawdata.nisos);
spikevector1(isospike(1))=1;
spikevector2=zeros(1,rawdata.nisos);
spikevector2(isospike(2))=1;

if (verLessThan('optim','4.0'))
	options=optimset('Display','notify','TolX',1e-8,'TolFun',1e-10,'TolCon',1e-6,'LargeScale','off','MaxFunEvals',10000);
else
	options=optimset('Display','notify','TolX',1e-8,'TolFun',1e-10,'TolCon',1e-6,'Algorithm','active-set','MaxFunEvals',10000);
end

tol=2e-5;    % how close can we be to the edges of the domain?
lb=[tol;
	tol];

ub=[1-tol;
	1-tol];

y0=[0.5 0.5]';

% Helpful to rescale the error, to make everything roughly order 1 for the optimiser
initialerror=errorestimate(element,y0(1),y0(2).*spikevector1+(1-y0(2)).*spikevector2,isoinv,errorratio,beta,alpha);

[y,opterr] = fmincon(@(y) errorestimate(element,y(1),y(2).*spikevector1 +(1-y(2)).*spikevector2,isoinv,errorratio,beta,alpha)./initialerror,y0,[],[],[],[],lb,ub,[],options);
[opterr,optppmperamu]=errorestimate(element,y(1),y(2).*spikevector1 +(1-y(2)).*spikevector2,isoinv,errorratio,beta,alpha);
optprop=y(1);
optspike=y(2)*spikevector1 +(1-y(2))*spikevector2;