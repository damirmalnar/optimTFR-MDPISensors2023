function [varargout] = makeAFv2(s,varargin)
%MAKEAFV2 Make signal Ambiguity Function in Doppler-lag (nu,tau) domain. 
%   Detailed explanation goes here TODO
% 
% Inputs:
%   required:
%       s - input signal, real or analytic
%   optional:
%       Nf - lag-to-frequency FFT length, default 2^nextpow2(length(s))
%       ti - time instants for which to calculate AF, default 1:length(s) 
% Outputs:
%   required:
%       -
%   optional:
%       CAF - signal Symmetrical Ambiguity Function (AF), NuxNf array
%       WDF - signal Wigner-Ville Distribution (WVD), NfxNt array
%       TCF - signal Instantaneous Autocorrelation Function (IAF), NtxNf array
%       Nf - lag-to-frequency FFT length, supplied or default
%       Nt - IAF and WVD time support length, length(ti)
%       Nu - AF Doppler support length, 2^nextpow2(2*length(ti))
%       Ns - signal length
%       z - signal analytic associate, Nx1 complex valued vector
%
% 
% 2023-06-19 vMDPISensors
%   (C) Damir Malnar 2023. Supplementary materials to:
%   Citation: Malnar, D.; Vrankic, M. Optimising Time-Frequency Distributions: A Surface Metrology Approach. Sensors 2023, 1, 0. https://doi.org/

%% Input check & Defaults -----------------------------------------------%%

    % Number of inputs must be >=minargs and <=maxargs.
    narginchk(1, 4)
%     fprintf('Received 1 required, %d optional inputs.\n\n',...
%         size(varargin, 2))
    % Signal length - as input
    Ns = numel(s);
    % Set defaults for optional inputs: freq. bins, time instants        
    optargs = {2^ceil(log2(Ns)) 1:Ns};
    % Put defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:length(varargin)) = varargin;
    
    % Place optional args in memorable variable names
    [Nf, t] = optargs{:};
    clear optargs
    
    % TCF and WVF time support length
    Nt = numel(t);
%% z(n) - Analytic signal
    z = hilbert(real(s(:)));

%% TCF(tau,t) - Temporal Correlation Function
    halfNf = round(Nf/2);
    TCF = zeros(Nf,Nt,'double');
    for iT = 1:Nt
        ti = t(iT);
        tauMax = min([ti-1,Ns-ti,halfNf-1]);
        tau = -tauMax:tauMax;
        iL = rem(Nf+tau,Nf)+1;
        TCF(iL,iT) = z(ti+tau).*conj(z(ti-tau));
        if (ti<=Ns-halfNf) && (ti>=halfNf+1)
            TCF(halfNf+1,iT) = 0.5*(z(ti+halfNf)*conj(z(ti-halfNf)) + ...
                                    z(ti-halfNf)*conj(z(ti+halfNf)));
        end
    end
    WDF = real(fft(TCF,[],1));
    Nu = 2*Nt;

%% CAF(tau,nu) - Complex Ambiguity Function

% CAF = fft(TCF,[],2);%original
    CAF = fft(ifft(WDF,[],1),Nu,2);% now CAF is 'symmetric'

%% Optional output ------------------------------------------------------%%

    varargout = {CAF,WDF,TCF,Nf,Nt,Nu,Ns,z};

end %[EOF makeAF.m]
