function newdatacov=changedenomcov(data,datacov,olddi,newdi)
% change denominator of covariance matrix for given set of ratios 

nisos=length(data)+1;
oldni=[1:(olddi-1) (olddi+1):nisos];
dataplus=[data(1:(olddi-1)) 1 data((olddi):end)];
newni=[1:(newdi-1) (newdi+1):nisos];

datacovplus=zeros(nisos,nisos);
datacovplus(oldni,oldni)=datacov;
	
A=eye(nisos)./dataplus(newdi);
A(:,newdi)=A(:,newdi)-dataplus'./(dataplus(newdi).^2);

newdatacovplus=A*datacovplus*(A');
newdatacov=newdatacovplus(newni,newni);