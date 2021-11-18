function out=dsinversion(element,measured,spike,isoinv,standard)
%DSINVERSION    Do the double spike inversion for a given set of measurements
%    out=DSINVERSION(element,measured,spike,isoinv,standard)
%             element -- element used in double spike, e.g. 'Fe'
%             measured -- a matrix of beam intensities. Columns correspond to the
%                 different isotopes e.g. for Fe, first column is 54Fe, second is 56Fe,
%                 third is 57Fe, fourth is 58Fe. The matrix should have the same number
%                 of columns as there are isotopes available.
%             spike -- a composition vector for the spike. e.g. [0 0 0.5 0.5] is a 50-50
%                 mix of 57Fe and 58Fe.
%             isoinv -- the four isotopes to use in the inversion, e.g [54 56 57 58]. This
%                 defaults to the first four isotopes.
%             standard -- standard composition or unspiked run data. Defaults to the 
%                 value in ISODATA.(element).standard if not specified.
%
% This routine performs the double spike inversion on measured data to return the 
% "true" composition of the sample. Output is returned as a structure with the
% following fields
%             alpha -- the inferred natural fractionations
%             beta -- the inferred instrumental fractionations
%             prop -- the inferred proportions of spike to sample
%             sample -- the inferred compositions of the sample
%             mixture -- the inferred compositions of the mixture
%
% Example
%   out=dsinversion('Fe',measured,[0 0 0.5 0.5],[54 56 57 58]);
%   plot(out.alpha);
%
% See also dsstartup

global ISODATA

% Set some defaults
if isempty(ISODATA)
	dsstartup;
end
if (nargin<5) || (isempty(standard))
	standard=ISODATA.(element).standard;
end
if (nargin<4) || (isempty(isoinv))
	isoinv=[1 2 3 4];
end
rawdata=ISODATA.(element);
isoinv=rawdata.isoindex(isoinv);

% Avoid division by zero errors for small values
if any(spike(1,isoinv)<0.001)
	[m,ix]=max(spike(1,isoinv));  % find largest spike of those used in inversion
	deno=isoinv(ix);
	nume=isoinv(isoinv~=deno);
	isoinv=[deno nume];               % rearrange isoinv so we denominator with largest isotope 
end

% Duplicate so all vectors are same length
nmeasured=max([size(spike,1) size(standard,1) size(measured,1)]);
if (size(standard,1)==1)
	standard=repmat(standard,[nmeasured 1]);
end
if (size(spike,1)==1)
	spike=repmat(spike,[nmeasured 1]);
end
if (size(measured,1)==1)
	measured=repmat(measured,[nmeasured 1]);
end

% Take ratios based on the isotopes we are inverting
in=calcratioeddata(element,isoinv);
in.AT=spike(:,in.Ani)./repmat(spike(:,in.di),[1 in.nratios]);
in.An=standard(:,in.Ani)./repmat(standard(:,in.di),[1 in.nratios]);
in.Am=measured(:,in.Ani)./repmat(measured(:,in.di),[1 in.nratios]);

% Do the fractionation correction
z=fractionationcorrection(in.AP,in.An,in.AT,in.Am,in.srat);

% Create a structure of outputs for things we're interested in
lambda=z(:,1);
out.alpha=z(:,2);
out.beta=z(:,3);

% Calculate sample and mixture proportion, and proportion by mole
AM=zeros(size(in.Am));
AN=zeros(size(in.Am));
for i=1:nmeasured
	AM(i,:)=in.Am(i,:).*exp(-in.AP*out.beta(i));
	AN(i,:)=in.An(i,:).*exp(-in.AP*out.alpha(i));
	prop=ratioproptorealprop([lambda(i) (1-lambda(i))],[in.AT(i,:);AN(i,:)]);
	out.prop(i,1)=prop(1);
end
out.sample=zeros(size(measured));
out.mixture=zeros(size(measured));

out.sample(:,in.Ani)=AN;
out.sample(:,in.di)=1;
out.mixture(:,in.Ani)=AM;
out.mixture(:,in.di)=1;
out.sample=out.sample./repmat(sum(out.sample,2),[1 rawdata.nisos]);
out.mixture=out.mixture./repmat(sum(out.mixture,2),[1 rawdata.nisos]);