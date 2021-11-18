function errorcurve2d(element,type,isospike,isoinv,errorratio,alpha,beta,resolution,threshold,ncontour,plottype,varargin)
%ERRORCURVE2D    A 2D contour plot of error as a function of double spike composition and double spike-sample proportions
%  ERRORCURVE2D(element,type,isospike,isoinv,errorratio,alpha,beta,resolution,threshold,ncontour,plottype,...)
%             element -- element used in double spike, e.g. 'Fe'
%             type -- type of spike, 'pure' or 'real'. Real spikes, such as those from
%                Oak Ridge National Labs, contain impurities (see 'maininput.csv'
%                or ISODATA.(element).rawspike) for the assumed compositions.
%                By default pure spikes are used.
%             isospike -- the isotopes used in the double spike e.g. [54 57].
%                By default the first two isotopes are chosen.
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default the first four isotopes are chosen.
%             errorratio -- by default, the error on the natural fractionation
%                factor (known as alpha) is given. Instead, the error on a
%                 particular ratio can be given by setting errorratio. e.g.
%                setting errorratio=[58 56] will give the error on 58Fe/56Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (natural and instrumental). Values of alpha and
%                beta can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%             resolution -- number of grid points in x and y. Default is 100.
%             threshold -- maximum contour to plot, relative to the minimum error.
%                Default is 0.25 i.e. 25% in excess of the minimum.
%             ncontour -- number of countours. Default is 25.
%             plottype -- by default, the error is plotted. By setting this to 'ppmperamu'
%                an estimate of the ppm per amu is plotted instead.
%             ... -- additional arguments are passed to contour command.
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%    errorcurve2d('Fe','pure',[57 58])
%
% See also errorestimate, contour
global ISODATA

% Set some default values
if isempty(ISODATA)
	dsstartup;
end
if (nargin<11) || isempty(plottype)
	plottype='default';
end
if (nargin<10) || isempty(ncontour)
	ncontour=25;
end
if (nargin<9) || isempty(threshold)
	threshold=0.25;
end
if (nargin<8) || isempty(resolution)
	resolution=100;
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
	isospike=[1 2];
end
if (nargin<2) || isempty(type)
	type='pure';
end
rawdata=ISODATA.(element);

% Convert isotope mass numbers to index numbers
errorratio=rawdata.isoindex(errorratio);
isoinv=rawdata.isoindex(isoinv);
isospike=rawdata.isoindex(isospike);

if strcmp(type,'pure')
	spike1=zeros(1,rawdata.nisos);
	spike2=zeros(1,rawdata.nisos);
	spike1(isospike(1))=1;
	spike2(isospike(2))=1;
else
	spike1=rawdata.rawspike(isospike(1),:);
	spike2=rawdata.rawspike(isospike(2),:);
end

prop=linspace(0.001,0.999,resolution);
spikeprop=linspace(0.001,0.999,resolution);

[iv jv]=meshgrid(1:length(prop),1:length(spikeprop));
[errvals ppmperamuvals]=arrayfun(@(i,j) errorestimate(element,prop(i),(spikeprop(j).*spike1 +(1-spikeprop(j)).*spike2),isoinv,errorratio,alpha,beta),iv,jv);
[optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalspike(element,type,isospike,isoinv,errorratio,alpha,beta);

%pcolor(prop,spikeprop,errvals);
%shading flat;
%caxis([opterr (1+threshold)*opterr]);

if strcmp(plottype,'ppmperamu')
	C=contour(prop,spikeprop,ppmperamuvals,linspace(optppmperamu,(1+threshold)*optppmperamu,ncontour+1),varargin{:});
else
	C=contour(prop,spikeprop,errvals,linspace(opterr,(1+threshold)*opterr,ncontour+1),varargin{:});
end

xlim([0 1]);
ylim([0 1]);
xlabel('proportion of double spike in double spike-sample mix');

if strcmp(type,'pure')
	ylabel(['proportion of ' rawdata.isolabel{isospike(1)} ' in ' rawdata.isolabel{isospike(1)} '-' rawdata.isolabel{isospike(2)} ' double spike']);
else
%	ylabel('proportion of first rawspike in double spike');
	ylabel(['proportion of ' rawdata.rawspikelabel{isospike(1)} ' in ' rawdata.rawspikelabel{isospike(1)} '-' rawdata.rawspikelabel{isospike(2)} ' double spike']);
end

invisostring=[rawdata.isolabel{isoinv(1)} ', ' rawdata.isolabel{isoinv(2)} ', ' rawdata.isolabel{isoinv(3)} ', ' rawdata.isolabel{isoinv(4)} ' inversion'];
if isempty(errorratio)
	title(['Error in \alpha (1SD) with ' invisostring]);
else
	title(['Error in ' rawdata.isolabel{errorratio(1)} '/' rawdata.isolabel{errorratio(2)} ' (1SD) with ' invisostring]);
end

hold on;
plot(optprop(1),optspikeprop(isospike(1)),'rx');
hold off;

