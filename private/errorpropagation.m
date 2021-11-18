function [Vz,VAN,VAM]=errorpropagation(z,AP,An,AT,Am,AVn,AVT,AVm,srat)
% does error propagation for the fractionation correction
% srat specifies which subset (3 ratios) you'd like to use

nratios=length(An);
nmeasured=length(Am(:,1));
Vz=zeros(nmeasured,3,3);
VAN=zeros(nmeasured,nratios,nratios);
VAM=zeros(nmeasured,nratios,nratios);
for i=1:nmeasured
	[Vz(i,:,:), VAN(i,:,:), VAM(i,:,:)]=fcerrorpropagation(z(i,:),AP,An,AT,Am(i,:),AVn,AVT,squeeze(AVm(i,:,:)),srat);
end