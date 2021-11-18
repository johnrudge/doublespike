function p=spikeprop(lambda,X1,X2)
% convert proportion in ratio space to proportion by moles
% X = lambda*X1 + (1-lambda)*X2 in ratio land

s1=(1+sum(X1));
s2=(1+sum(X2));
p1=lambda.*s1./(((1-lambda).*s2)+(lambda.*s1));

prop=ratioproptorealprop([lambda' (1-lambda)'],[X1;X2]);
p2=prop(:,1)';

p=p1;

