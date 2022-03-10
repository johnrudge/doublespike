%SETERRORMODEL   Sets the coefficients of the error model in the global variable ISODATA
%    SETERRORMODEL(element,intensity,deltat,R,T,radiogenicisos,type,R_reference)
%             element -- element used in double spike, e.g. 'Fe'
%             intensity -- mean total beam intensity (volts). Default is 10 V.
%             deltat -- integration time (seconds). Default is 8 s.
%             R -- resistance (ohms). Default is 10^11 ohms.
%                If a float then resistors are identical for all isotopes
%                If an array then individual resistances for each isotope
%             T -- temperature (Kelvin). Default is 300 K.
%             radiogenicisos -- which isotopes use the radiogenic error model
%                  (errors on standard). Default is {'Pb','Sr','Hf','Os','Nd'}
%             type -- If 'fixed-total' then total beam intensity of mixture is fixed,
%                If 'fixed-sample' then beam current from the sample is fixed.
%             R_reference -- reference resistance used for describing beam intensity.
%                The total beam current is intensity/R_reference
%                e.g. the 10 V default beam corresponds to 100 pA with R_reference=1e11 Ohms.
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
%    seterrormodel('Fe',15,4);   % set 15 V total beam with 4 second integrations
%
% See also dsstartup
function seterrormodel(element,intensity,deltat,R,T,radiogenicisos,type)
global ISODATA

% Fundamental constants
elementarycharge=1.60217646e-19;  % Coulombs
k=1.3806504e-23;                  % Boltzman's constant (m^2 kg s^-2 K^-1)

if (nargin<8)
        R_reference = 1e11; % default resistance used for describing beam intensity
end

if (nargin<7)
        type = 'fixed-total'; % by default, fix the total intensity of beams for the mixture
end

if (nargin<6)
	radiogenicisos={'Pb','Sr','Hf','Os','Nd'};  % by default, use a different error model for these, as requires two runs.
end

% Mass spec properties
if (nargin<5) || isempty(T)
	T=300;            % Temperature (Kelvin)
end
if (nargin<4) || isempty(R)
	R=1e11;           % resistance (ohms)
end
if (nargin<3) || isempty(deltat)
	deltat=8;         % integration time (seconds)
end
if (nargin<2) || isempty(intensity)
	intensity=10;     % default mean total intensity, 10 volts
end

nisos=ISODATA.(element).nisos;

if isscalar(R)
    R=R*ones(1,nisos);
end

a=4*k*T*R_reference^2./(deltat*R); % Johnson-Nyquist noise
b=elementarycharge*R_reference/(deltat).*ones(1,nisos); % Counting statistics
c=zeros(1,nisos);

% by default assume Johnson noise and counting statistics
errormodel.measured.type='fixed-total';
errormodel.measured.intensity=intensity;
errormodel.measured.a=a;
errormodel.measured.b=b;
errormodel.measured.c=c;

% by default, assume no noise on the spike
errormodel.spike.type='fixed-total';
errormodel.spike.intensity=intensity;
errormodel.spike.a=0.*ones(1,nisos);
errormodel.spike.b=0.*ones(1,nisos);
errormodel.spike.c=0.*ones(1,nisos);

% by default, assume no noise on standard unless it is radiogenic
errormodel.standard.type='fixed-total';
errormodel.standard.intensity=intensity;
if isempty(intersect(radiogenicisos,{element}))
	errormodel.standard.a=0.*ones(1,nisos);
	errormodel.standard.b=0.*ones(1,nisos);
	errormodel.standard.c=0.*ones(1,nisos);
else
	errormodel.standard.a=a;
	errormodel.standard.b=b;
	errormodel.standard.c=c;
end

ISODATA.(element).errormodel=errormodel;
