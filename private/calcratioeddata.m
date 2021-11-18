function out=calcratioeddata(element,isoinv)
% Calculates the relevant isotopic ratios given the four elements used in the inversion
% Assumes first of the four isotopes in isoinv is the denominator
global ISODATA
rawdata=ISODATA.(element);

if (nargin<2) || isempty(isoinv)
	isoinv=[1 2 3 4];
end

out.isoinv=rawdata.revisoindex(isoinv);
out.d=out.isoinv(1);
out.n=out.isoinv(2:end);
%out.An=rawdata.isonum(rawdata.isonum~=out.d);
out.di=find(rawdata.isonum==out.d);
out.ni(1)=find(rawdata.isonum==out.n(1));
out.ni(2)=find(rawdata.isonum==out.n(2));
out.ni(3)=find(rawdata.isonum==out.n(3));
out.Ani=find(rawdata.isonum~=out.d);
out.srat(1)=find(out.Ani==out.ni(1));
out.srat(2)=find(out.Ani==out.ni(2));
out.srat(3)=find(out.Ani==out.ni(3));
out.isoinv=rawdata.isoindex(out.isoinv);
out.nratios=rawdata.nisos-1;

out.An=rawdata.standard(out.Ani)./rawdata.standard(out.di);
out.AP=log(rawdata.mass(out.Ani)./rawdata.mass(out.di));
out.AS=rawdata.rawspike(:,out.Ani)./repmat(rawdata.rawspike(:,out.di),[1 out.nratios]);
