function newdata=changedenom(data,olddi,newdi)
% change denominator for given set of ratios
nisos=length(data(1,:))+1;
ndata=length(data(:,1));
dataplus=[data(:,1:(olddi-1)) ones([ndata 1]) data(:,(olddi):end)];
denomplus=repmat(dataplus(:,newdi),[1 nisos]);
newdataplus=dataplus./denomplus;
newni=[1:(newdi-1) (newdi+1):nisos];
newdata=newdataplus(:,newni);