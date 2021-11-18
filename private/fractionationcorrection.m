function z=fractionationcorrection(AP,An,AT,Am,srat)
% does fractionation correction for a vector of measured values
% srat specifies which 3 ratios you'd like to use

% Select appropriate ratios
P=AP(:,srat);                    % log of ratio of atomic masses
n=An(:,srat);                    % ratio standard
T=AT(:,srat);                    % ratio spike
m=Am(:,srat);                    % ratio measured

nmeasured=max([size(m,1) size(T,1) size(n,1)]);
% Duplicate so all matrices are same length
if (size(n,1)==1)
	n=repmat(n,[nmeasured 1]);
end
if (size(T,1)==1)
	T=repmat(T,[nmeasured 1]);
end
if (size(m,1)==1)
	m=repmat(m,[nmeasured 1]);
end

% dscorrection only works with a single value, here run over array
options = optimset('Display','off','Jacobian','on','TolFun',1e-15);
z=zeros(nmeasured,3);
for i=1:nmeasured
	z(i,:)=dscorrection(P,n(i,:),T(i,:),m(i,:),options);
end