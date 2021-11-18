function cocktail(type,filename,elements)
% COCKTAIL   Generate double spike cocktail lists
%  COCKTAIL(type,filename,elements)
%             type -- type of spike, 'pure' or 'real'. Real spikes, such as those from
%                Oak Ridge National Labs, contain impurities (see 'maininput.csv'
%                or ISODATA.(element).rawspike) for their assumed compositions.
%                By default pure spikes are used.
%             filename -- file to store output (a CSV file).  Default is either 'cookbook_pure.csv'
%                or 'cookbook_real.csv' depending on the type.
%             elements -- which elements to include in the cookbook. Specify as a cell array
%                e.g. {'Ca','Fe'}. Default is all possible elements.
%
% This generates an exhaustive list of all possible double spikes for specified elements
% sorted in order of error.
%
% Note that a number of parameters are specified in the global variable ISODATA.
%
% Example
%   cocktail('real')
%
% See also dsstartup
global ISODATA

% default argument
if isempty(ISODATA)
	dsstartup;
end
if (nargin<1) || isempty(type)
	type='pure';
end
if (nargin<2) || isempty(filename)
	filename=['cookbook_' type '.csv'];
end
if (nargin<3)
	elements=fieldnames(ISODATA);
end
disp(['Writing to ' filename]);

title=['A double spike cocktail list: ' type ' spikes'];

fwritecell(filename,'%s','w',{title});
fwritecell(filename,'%s','a',{''});

for i=1:length(elements)
	element=elements{i};
	disp(element);
	in=ISODATA.(element);
	if (isequal(type,'pure')) || (isequal(type,'real') && in.nspikes>1)
		[optspike,optprop,opterr,optisoinv,optspikeprop,optppmperamu]=optimalspike(element,type);
		optisoinv=in.isoindex(optisoinv);
		optisonams=[{in.isoname{optisoinv(:,1)}}' {in.isoname{optisoinv(:,2)}}' {in.isoname{optisoinv(:,3)}}' {in.isoname{optisoinv(:,4)}}'];

		% write output to file
		isohead=strcat(repmat({'iso'},4,1),num2str((1:4)'))';

		if isequal(type,'pure')
			output=[optisonams num2cell([optspike optprop (1-optprop) opterr optppmperamu])];
			header=[isohead in.isoname {'spike'} {'sample'} {'error'} {'ppmperamu'}];
		
			fwritecell(filename,[repmat('%s,',1,4) repmat('%s,',1,in.nisos) '%s,%s,%s,%s'],'a',header);
			fwritecell(filename,[repmat('%s,',1,4) repmat('%f,',1,in.nisos) '%f,%f,%f,%f'],'a',output);
			fwritecell(filename,'%s','a',{''});
		else
			spikehead=strcat(repmat({'spike'},in.nspikes,1),num2str((1:in.nspikes)'))';
			output=[optisonams num2cell([optspike optspikeprop optprop (1-optprop) opterr optppmperamu])];
			header=[isohead in.isoname spikehead {'spike'} {'sample'} {'error'} {'ppmperamu'}];
		
			fwritecell(filename,[repmat('%s,',1,4) repmat('%s,',1,in.nisos) repmat('%s,',1,in.nspikes) '%s,%s,%s,%s'],'a',header);
			fwritecell(filename,[repmat('%s,',1,4) repmat('%f,',1,in.nisos) repmat('%f,',1,in.nspikes) '%f,%f,%f,%f'],'a',output);
			fwritecell(filename,'%s','a',{''});
		end
	end
end
disp(['Output written to ' filename]);