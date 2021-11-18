function [Vz, VAN, VAM]=fcerrorpropagation(z,AP,An,AT,Am,VAn,VAT,VAm,srat)
% linear error propagation for the fractionation correction

lambda=z(1);
alpha=z(2);
beta=z(3);

AM=Am.*exp(-AP*beta);
AN=An.*exp(-AP*alpha);

% Select appropriate ratios
P=AP(srat);                    % log of ratio of atomic masses
N=AN(srat);                    % ratio sample
T=AT(srat);                    % ratio spike
M=AM(srat);                    % ratio mixture
VT=VAT(srat,srat);             % covariance matrix of spike ratios
Vm=VAm(srat,srat);             % covariance matrix of measured ratios
Vn=VAn(srat,srat);             % covariance matrix of standard

% calculate various Jacobian matrices
dfdlambda=T - (N.*(1+alpha.*P));
dfdu=-N.*P;
dfdbeta=M.*P;
dfdy=[dfdlambda' dfdu' dfdbeta'];

dfdT=lambda.*eye(3);
dfdm=-diag(exp(-beta*P));
dfdn=(1-lambda)*diag(exp(-alpha*P));

% matrix to convert from (lambda, (1-lambda)alpha,beta) to (lambda,alpha,beta) 
K=[1 0 0
   (alpha/(1-lambda)) (1/(1-lambda)) 0
   0 0 1];

dzdT=-K*(dfdy\dfdT);
dzdm=-K*(dfdy\dfdm);
dzdn=-K*(dfdy\dfdn);

% Covariance matix for (lambda,beta,alpha)
Vz=dzdT*VT*(dzdT')+ dzdm*Vm*(dzdm') + dzdn*Vn*(dzdn'); 

% full matrices for all ratios
nratios=length(An);
dzdAT=zeros(3,nratios);
dzdAn=zeros(3,nratios);
dzdAm=zeros(3,nratios);
dzdAT(1:3,srat)=dzdT;
dzdAn(1:3,srat)=dzdn;
dzdAm(1:3,srat)=dzdm;

% Covariance matrix of sample
dalphadAT=dzdAT(2,:);
dalphadAn=dzdAn(2,:);
dalphadAm=dzdAm(2,:);

dANdAT=-(((AN.*AP)')*dalphadAT);
dANdAn=diag(exp(-AP.*alpha)) - (((AN.*AP)')*dalphadAn);
dANdAm=-(((AN.*AP)')*dalphadAm);

VAN=dANdAn*VAn*(dANdAn') + dANdAT*VAT*(dANdAT') + dANdAm*VAm*(dANdAm');

% Covariance matrix of mixture
dbetadAT=dzdAT(3,:);
dbetadAn=dzdAn(3,:);
dbetadAm=dzdAm(3,:);

dAMdAT=-(((AM.*AP)')*dbetadAT);
dAMdAn=-(((AM.*AP)')*dbetadAn);
dAMdAm=diag(exp(-beta.*AP)) -(((AM.*AP)')*dbetadAm);

VAM=dAMdAn*VAn*(dAMdAn') + dAMdAT*VAT*(dAMdAT') + dAMdAm*VAm*(dAMdAm');
