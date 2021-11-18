function data=loadrawdata()
% main routine for loading in the data about the different isotopic systems, such as standard compositions,
% atomic masses, and raw spike compositions

filename='maininput.csv';

fid=fopen(filename);
fgetl(fid);  % ignore header line

j=1;
data=[];
while 1
	tline = fgetl(fid);
	% bail out when file is finished
	if ~ischar(tline)
		break;
	end
	% split line based on commas
%	splitted=regexp(tline,',','split');  % This doesn't work on MATLAB's before 2007b...
	splitted=csvsplit(tline);

	if ~isempty(splitted{1})
		el=regexprep(splitted{1},'"',''); % remove quotes
		element{j}=el;
		j=j+1;
		i=1;
	end
	numbers=str2double({splitted{2:end}});
	data.(el).isonum(i)=numbers(1);
	data.(el).mass(i)=numbers(2);
	data.(el).standard(i)=numbers(3);
	data.(el).rawspike(:,i)=numbers(4:end);
	data.(el).isoname{i}=[el num2str(data.(el).isonum(i))];
	data.(el).isolabel{i}=['^{' num2str(data.(el).isonum(i)) '}' el];
	i=i+1;
end

% Renormalise compositions, add a few useful pieces of information
for j=1:length(element)
	el=element{j};
	data.(el).element=el;
	data.(el).standard=data.(el).standard/sum(data.(el).standard);
	nisos=length(data.(el).isonum);
	missingraws=any(isnan(data.(el).rawspike),2);
	data.(el).rawspike=data.(el).rawspike(~missingraws,:);
	data.(el).rawspike=data.(el).rawspike./repmat(sum(data.(el).rawspike,2),[1 nisos]);
	data.(el).isoindex(1,data.(el).isonum)=1:nisos;
	data.(el).isoindex(1,1:nisos)=1:nisos;
	data.(el).revisoindex(1,data.(el).isonum)=data.(el).isonum;
	data.(el).revisoindex(1,1:nisos)=data.(el).isonum;
	data.(el).nisos=length(data.(el).isonum);          % the number of isotopes
	data.(el).nratios=data.(el).nisos-1;               % the number of ratios 
	data.(el).nspikes=length(data.(el).rawspike(:,1)); % the number of spikes
	for k=1:data.(el).nspikes
		data.(el).rawspikelabel{k}=['spike ' num2str(k)];
	end
end

fclose(fid);

