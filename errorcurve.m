function errorcurve(element,spike,isoinv,errorratio,alpha,beta,plottype,varargin)
%ERRORCURVE    A plot of error as a function of double spike-sample proportions for a given double spike composition
%  ERRORCURVE(element,spike,isoinv,errorratio,alpha,beta,...)
%             element -- element used in double spike, e.g. 'Fe'
%             spike -- the composition of the double spike as a composition vector e.g. [0 0 0.5 0.5]
%                represents a 50-50 mixture of the third and fourth isotopes (57-58 for Fe).
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default the first four isotopes are chosen.
%             errorratio -- by default, the error on the natural fractionation
%                factor (known as alpha) is given. Instead, the error on a
%                 particular ratio can be given by setting errorratio. e.g.
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
%    errorcurve('Fe',[0 0 0.5 0.5])
%
% See also errorestimate, errorcurve2
global ISODATA

if isempty(ISODATA)
	dsstartup;
end
if (nargin<7) || isempty(plottype)
	plottype='default';
end
if (nargin<6) || isempty(beta)
	beta=0;
end
if (nargin<5) || isempty(alpha)
	alpha=0;
end
if (nargin<4) || isempty(errorratio)
	errorratio=[];
end
if (nargin<3) || isempty(isoinv)
	isoinv=[1 2 3 4];
end
rawdata=ISODATA.(element);
spike=spike/sum(spike);

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isoinv=rawdata.isoindex(isoinv);

pvals=linspace(0.001,0.999,1000);

errvals=zeros(size(pvals));
ppmperamuvals=zeros(size(pvals));
for i=1:length(pvals)
	[errvals(i) ppmperamuvals(i)]=errorestimate(element,pvals(i),spike,isoinv,errorratio,alpha,beta);
end

if isequal(plottype,'ppmperamu')
	plotvals=ppmperamuvals;
else
	plotvals=errvals;
end

plot(pvals,plotvals,varargin{:});
mine=min(plotvals);
xlim([0 1]);
ylim([0 5*mine]);
xlabel('proportion of double spike in double spike-sample mix');
if isempty(errorratio)
	ylabel('Error in \alpha (1SD)');
else
	ylabel(['Error in ' rawdata.isolabel{errorratio(1)} '/' rawdata.isolabel{errorratio(2)} ' (1SD)']);
end
title([rawdata.isolabel{isoinv(1)} ', ' rawdata.isolabel{isoinv(2)} ', ' rawdata.isolabel{isoinv(3)} ', ' rawdata.isolabel{isoinv(4)} ' inversion' ])