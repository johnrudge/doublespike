function [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=errorcurveoptimalspike(element,type,isospike,isoinv,errorratio,alpha,beta,plottype,varargin)
%ERRORCURVEOPTIMALSPIKE    Find the optimal double spike compositions and plot the corresponding error curves
% [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]
%  =ERRORCURVEOPTIMALSPIKE(element,type,isoinv,isospike,errorratio,alpha,beta,...)
%             element -- element used in double spike, e.g. 'Fe'
%                This is the only mandatory argument.
%             type -- type of spike, 'pure' or 'real'. Real spikes, such as those from
%                Oak Ridge National Labs, contain impurities (see 'maininput.csv'
%                or ISODATA.(element).rawspike) for their assumed compositions.
%                By default pure spikes are used.
%             isospike -- the isotopes used in the double spike e.g. [54 57].
%                By default all choices of 2 isotopes are tried.
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default the first four isotopes are used.
%             errorratio -- by default, the optimal spike is chosen as that which
%                minimises the error on the natural fractionation factor (known as
%                alpha). Instead, the optimiser can be told to minimise the
%                error on a particular ratio by setting errorratio. e.g.
%                setting errorratio=[58 56] will minimise the error on 58Fe/56Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (natural and instrumental). Values of alpha and
%                beta can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%             plottype -- by default, the error is plotted. By setting this to 'ppmperamu'
%                an estimate of the ppm per amu is plotted instead.
%             ... -- additional arguments are passed to the plot command
%
% Outputs are the same as those of optimalspike.
%
% Example
%   errorcurveoptimalspike('Fe')
%
% See also optimalspike, errorcurve

global ISODATA

% Have some default arguments
if isempty(ISODATA)
	dsstartup;
end
if (nargin<8) || isempty(plottype)
	plottype='default';
end
if (nargin<7) || isempty(beta)
	beta=0;
end
if (nargin<6) || isempty(alpha)
	alpha=0;
end
if (nargin<5) || isempty(errorratio)
	errorratio=[];
end
if (nargin<4) || isempty(isoinv)
	isoinv=[1 2 3 4];
end
if (nargin<3) || isempty(isospike)
	isospike=[];
end
if (nargin<2) || isempty(type)
	type='pure';
end
rawdata=ISODATA.(element);


% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isospike=rawdata.isoindex(isospike);
isoinv=rawdata.isoindex(isoinv);

% Find the optimal spikes
[optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalspike(element,type,isospike,isoinv,errorratio,alpha,beta);

cols=repmat('brgcmky',1,size(optspike,1));

for j=1:size(optspike,1)
	if strcmp(type,'pure')
		spiked=find(optspike(j,:)>0);
		leglabel=[rawdata.isolabel{spiked(1)} '-' rawdata.isolabel{spiked(2)}];
	else
		spiked=find(optspikeprop(j,:)>0);
		leglabel=[rawdata.rawspikelabel{spiked(1)} '-' rawdata.rawspikelabel{spiked(2)}];
		%leglabel=['rawspikes ' num2str(spiked(1)) '-' num2str(spiked(2))];
%  		if (size(optspikeprop,2)==2)
%  			leglabel=['rawspikes ' num2str(isospike(1)) '-' num2str(isospike(2))];
%  		else
%  			leglabel=['best real spike'];
%  		end
	end
	%errorcurve(element,alpha,beta,optspike(j,:),errorratio,isoinv,cols(j),'DisplayName',leglabel);
	errorcurve(element,optspike(j,:),isoinv,errorratio,alpha,beta,plottype,cols(j),varargin{:},'DisplayName',leglabel);
	hold on;
end
if isequal(plottype,'ppmperamu')
	ylim([0 5*min(optppmperamu)]);
else
	ylim([0 5*min(opterr)]);
end
hold off;

legend('show');