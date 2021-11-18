function z=dscorrection(P,n,T,m,options)
% Routine for double spike fractionation correction
% takes as input ratios: 
%       P -- log of ratio of atomic masses
%       n -- ratio of standard/ unspiked run
%       T -- ratio of spike
%       m -- ratio of measured
% outputs enriched ratio proportion (lambda), natural fractionation (alpha),
% and instrumental fractionation (beta) as a vector z=(lambda, (1-lambda)*alpha, beta)

% start by solving the linear problem
b=(m-n)';
A=[(T-n)' (-n.*P)' (m.*P)'];
y0=A\b;

% by starting at the linear solution, solve the non-linear problem
[y,fval,exitflag,output]=fsolve(@(y)F(y,P,n,T,m),y0,options);
z=y;
z(2)=y(2)/(1-y(1));  % z(2) is alpha

function [fval,Jac]=F(y,P,n,T,m)
% The nonlinear equations to solve
lambda=y(1);
alpha=y(2)/(1-lambda);
beta=y(3);

N=n.*exp(-alpha.*P);
M=m.*exp(-beta.*P);
fval=(lambda.*T)+((1-lambda).*N)-M;

% The Jacobian of the nonlinear equations -- can speed up root finding, but is not required
dfdlambdaprime=T - (N.*(1+alpha.*P));
dfdu=-N.*P;
dfdbeta=M.*P;
Jac=[dfdlambdaprime' dfdu' dfdbeta'];
