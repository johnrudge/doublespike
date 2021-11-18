function prop=realproptoratioprop(realprop,ratios)
% convert a proportion in moles to a proportion in ratio space

nratios=length(ratios(:,1));
nrealprops=length(realprop(:,1));

invpropdenom=1./(sum([ones(nratios,1) ratios],2))';
newinvpropdenom=repmat(invpropdenom,[nrealprops 1]);

prop=(realprop.*newinvpropdenom);
sprop=sum(prop,2);
sprop=repmat(sprop,[1 nratios]);

prop=prop./sprop;
