function errorcurve2(element,type,prop,isospike,isoinv,errorratio,alpha,beta,plottype,varargin)
%ERRORCURVE2    A plot of error as a function of double spike proportions for a given double spike-sample proportion
%  ERRORCURVE2(element,type,prop,isospike,isoinv,errorratio,alpha,beta,...)
%             element -- element used in double spike, e.g. 'Fe'
%             type -- type of spike, 'pure' or 'real'. Real spikes, such as those from
%                Oak Ridge National Labs, contain impurities (see 'data/maininput.csv'
%                or ISODATA.(element).rawspike) for their assumed compositions.
%                By default pure spikes are used.
%             prop -- proportion of double spike in double spike-sample mixture e.g. 0.5.
%             isospike -- the isotopes used in the double spike e.g. [54 57].
%                By default the first two isotopes are chosen.
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default the first four isotopes are chosen.
%             errorratio -- by default, the error on the natural fractionation
%                factor (known as alpha or alpha) is given. Instead, the
%                error on a particular ratio can be given by setting errorratio. e.g.
%                setting errorratio=[58 56] will give the error on 58Fe/56Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (natural and instrumental). Values of beta and
%                alpha can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%             plottype -- by default, the error is plotted. By setting this to 'ppmperamu'
%                an estimate of the ppm per amu is plotted instead.
%             ... -- additional arguments are passed to the plot command.
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%    errorcurve2('Fe','real',0.5,[54 57])
%
% See also errorestimate, errorcurve
global ISODATA

if isempty(ISODATA)
	dsstartup;
end
if (nargin<9) || isempty(plottype)
	plottype='default';
end
if (nargin<8) || isempty(beta)
	beta=0;
end
if (nargin<7) || isempty(alpha)
	alpha=0;
end
if (nargin<6) || isempty(errorratio)
	errorratio=[];
end
if (nargin<5) || isempty(isoinv)
	isoinv=[1 2 3 4];
end
if (nargin<4) || isempty(isospike)
	isospike=[1 2];
end
if (nargin<3) || isempty(prop)
	prop=0.5;
end
if (nargin<2) || isempty(type)
	type='pure';
end
rawdata=ISODATA.(element);
rawspike=rawdata.rawspike;

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isoinv=rawdata.isoindex(isoinv);
isospike=rawdata.isoindex(isospike);

qvals=linspace(0.001,0.999,1000);
errvals=zeros(size(qvals));
ppmperamuvals=zeros(size(qvals));
for i=1:length(qvals)
	spike=(qvals(i).*rawspike(isospike(1),:))+((1-qvals(i)).*rawspike(isospike(2),:));
	[errvals(i) ppmperamuvals(i)]=errorestimate(element,prop,spike,isoinv,errorratio,alpha,beta);
end

if isequal(plottype,'ppmperamu')
	plotvals=ppmperamuvals;
else
	plotvals=errvals;
end

plot(qvals,plotvals,varargin{:});
mine=min(plotvals);
xlim([0 1]);
ylim([0 5*mine]);

if strcmp(type,'pure')
	xlabel(['proportion of ' rawdata.isolabel{isospike(1)} ' in ' rawdata.isolabel{isospike(1)} '-' rawdata.isolabel{isospike(2)} ' double spike']);
else
	xlabel('proportion of first rawspike in double spike');
end

if isempty(errorratio)
	ylabel('Error in \alpha (1SD)');
else
	ylabel(['Error in ' rawdata.isolabel{errorratio(1)} '/' rawdata.isolabel{errorratio(2)} ' (1SD)']);
end
title([rawdata.isolabel{isoinv(1)} ', ' rawdata.isolabel{isoinv(2)} ', ' rawdata.isolabel{isoinv(3)} ', ' rawdata.isolabel{isoinv(4)} ' inversion' ])