function splitted=csvsplit(line)
% This is a replacement for the simple command 
%     regexp(tline,',','split')
% which doesn't seem to work on versions of MATLAB before 2007b
% Cumbersome code -- please improve!

splitted2=regexp(line,'([^,]*),?()|,','match');
splitted=splitted2;
for m=1:length(splitted2)
	spit=regexp(splitted2{m},'[^,]*','match');
	if isempty(spit)
		splitted{m}='';
	else
		splitted{m}=spit{1};
	end
end
if isequal(splitted2{end},',')
	splitted{length(splitted2)+1}='';
end
