function prop=ratioproptorealprop(ratprop,ratios)
% convert a proportion in ratio space to one per mole

nratios=length(ratios(:,1));
nratprops=length(ratprop(:,1));

invpropdenom=(sum([ones(nratios,1) ratios],2))';
newinvpropdenom=repmat(invpropdenom,[nratprops 1]);

prop=(ratprop.*newinvpropdenom);
sprop=sum(prop,2);
sprop=repmat(sprop,[1 nratios]);

prop=prop./sprop;
