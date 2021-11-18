%SETERRORMODEL   Sets the coefficients of the error model in the global variable ISODATA
%    SETERRORMODEL(intensity,deltat,R,T,radiogenicisos)
%             intensity -- mean total beam intensity (volts). Default is 10 V.
%             deltat -- integration time (seconds). Default is 8 s.
%             R -- resistance (ohms). Default is 10^11 ohms.
%             T -- temperature (Kelvin). Default is 300 K.
%             radiogenicisos -- which isotopes use the radiogenic error model
%                  (errors on standard). Default is {'Pb','Sr','Hf','Os','Nd'}
%
% This function generates the coefficients of the ISODATA.(element).errormodel
% By default, only the measurements of the mixture have errors. For radiogenic
% isotopes, the standard (i.e. unspiked run) also has errors.
%
% The error model for the beams takes the form 
%     sigma_^2 = a + b * mu + c * mu^2
% Each beam is assumed independent.
%
% Example
%    seterrormodel(15,4);   % set 15 V total beam with 4 second integrations
%
% See also dsstartup
function seterrormodel(intensity,deltat,R,T,radiogenicisos,type)
global ISODATA

% Fundamental constants
elementarycharge=1.60217646e-19;  % Coulombs
k=1.3806504e-23;                  % Boltzman's constant (m^2 kg s^-2 K^-1)

if (nargin<6)
        type = 'fixed-mixture'; % by default, fix the intensity of beams for the mixture
end

if (nargin<5)
	radiogenicisos={'Pb','Sr','Hf','Os','Nd'};  % by default, use a different error model for these, as requires two runs.
end

% Mass spec properties
if (nargin<4) || isempty(T)
	T=300;            % Temperature (Kelvin)
end
if (nargin<3) || isempty(R)
	R=1e11;           % resistance (ohms)
end
if (nargin<2) || isempty(deltat)
	deltat=8;         % integration time (seconds)
end
if (nargin<1) || isempty(intensity)
	intensity=10;     % default mean total intensity, 10 volts
end

a=4*k*T*R/(deltat);             % Johnson-Nyquist noise in volts
b=elementarycharge*R/(deltat);  % Counting statistics prefactor

els=fieldnames(ISODATA);
for i=1:length(els)
	element=els{i};
	nisos=ISODATA.(element).nisos;

	% by default assume Johnson noise and counting statistics
        errormodel.measured.type='fixed-mixture';
	errormodel.measured.intensity=intensity;
	errormodel.measured.a=a.*ones(1,nisos);
	errormodel.measured.b=b.*ones(1,nisos);
	errormodel.measured.c=0.*ones(1,nisos);

	% by default, assume no noise on the spike
        errormodel.spike.type='fixed-mixture';
	errormodel.spike.intensity=intensity;
	errormodel.spike.a=0.*ones(1,nisos);
	errormodel.spike.b=0.*ones(1,nisos);
	errormodel.spike.c=0.*ones(1,nisos);

	% by default, assume no noise on standard unless it is radiogenic
        errormodel.standard.type='fixed-mixture';
	errormodel.standard.intensity=intensity;
	if isempty(intersect(radiogenicisos,{element}))
		errormodel.standard.a=0.*ones(1,nisos);
		errormodel.standard.b=0.*ones(1,nisos);
	else
		errormodel.standard.a=a.*ones(1,nisos);
		errormodel.standard.b=b.*ones(1,nisos);
	end
	errormodel.standard.c=0.*ones(1,nisos);

	ISODATA.(element).errormodel=errormodel;
end

