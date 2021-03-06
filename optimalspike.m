function [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalspike(element,type,isospike,isoinv,errorratio,alpha,beta)
%OPTIMALSPIKE    Find the optimal double spike composition and double spike-sample mixture proportions
% [optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]
%  =OPTIMALSPIKE(element,type,isospike,isoinv,errorratio,alpha,beta)
%             element -- element used in double spike, e.g. 'Fe'
%                This is the only mandatory argument.
%             type -- type of spike, 'pure' or 'real'. Real spikes, such as those from
%                Oak Ridge National Labs, contain impurities (see 'data/maininput.csv'
%                or ISODATA.(element).rawspike) for their assumed compositions.
%                By default pure spikes are used.
%             isospike -- the isotopes used in the double spike e.g. [54 57].
%                By default all choices of 2 isotopes are tried.
%             isoinv -- the isotopes used in the inversion, e.g. [54 56 57 58].
%                By default all choices of 4 isotopes are tried.
%             errorratio -- by default, the optimal double spike is chosen as that which
%                minimises the error on the natural fractionation factor (known as
%                alpha). Instead, the optimiser can be told to minimise the
%                error on a particular ratio by setting errorratio. e.g.
%                setting errorratio=[58 56] will minimise the error on 58Fe/56Fe.
%             alpha, beta -- there is a small dependance of the error on the fractionation
%                factors (instrumental and natural, or alpha and beta). Values of alpha and
%                beta can be set here if desired, although the effect on the optimal spikes
%                is slight unless the fractionations are very large. Default is zero.
%
% All the outputs are provided as matrices. Each column represents an isotope
% (see ISODATA.(element).isonum for the isotope numbers) e.g. for Fe the columns
% correspond to the isotopes 54Fe, 56Fe, 57Fe, 58Fe. The rows represent the
% different combinations of double spikes and isotopes being tried, in order of
% error: The first row is the best double spike, and the last row is the worst.
%     optspike -- the proportions of each isotope in the optimal double spike.
%     optprop -- the optimal proportion of spike in the double spike-sample mix.
%     opterr -- the error in the fractionation factor (or ratio if specified)
%               for the optimal spike.
%     optisoinv -- the 4 isotopes used in the inversion.
%     optspikeprop -- the proportion of each raw spike in the optimal double spike.
%     optppmperamu -- an alternative expression of the error in terms of ppm per amu.
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%   [optspike,optprop,opterr]=optimalspike('Fe')
%
% See also dsstartup
global ISODATA

% Have some default arguments
if isempty(ISODATA)
	dsstartup;
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
	isoinv=[];
end
if (nargin<3) || isempty(isospike)
	isospike=[];
end
if (nargin<2) || isempty(type)
	type='pure';
end

% Convert isotope mass numbers to index numbers
errorratio=ISODATA.(element).isoindex(errorratio);

if strcmp(type,'pure')
	[optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalpurespike(element,beta,alpha,errorratio,isospike,isoinv);
else
	[optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalrealspike(element,beta,alpha,errorratio,isospike,isoinv);
end
