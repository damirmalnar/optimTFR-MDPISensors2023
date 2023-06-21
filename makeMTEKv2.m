function MTEK = makeMTEKv2(p,Tau,Nu)
%MAKEMTEKV2 Make Multiform Tiltable Exponential Kernel in Doppler-Lag (tau,nu) domain.
%   Detailed explanation goes here TODO
%
% Inputs:
%   required:
%       p   - kernel parameters: tau0, nu0, r, beta, lambda
%       Tau - lag-to-frequency FFT length
%       Nu  - tfr time support length
%
%   optional:
%       -
% Outputs:
%   required:
%       MTEK - Nfx(2*Nt) array matching the signal ambiguity domain size
%
%   optional:
%       -
% Note:
%   uses MATLAB2020b automatic vector expansion. Use bsxfun() for older versions.
%
%
% 2023-06-19 vMDPISensors
%   (C) Damir Malnar 2023. Supplementary materials to:
%   Citation: Malnar, D.; Vrankic, M. Optimising Time-Frequency Distributions: A Surface Metrology Approach. Sensors 2023, 1, 0. https://doi.org/
    
%% Input check & Defaults -----------------------------------------------%%


%% Kernel preparation ---------------------------------------------------%%
    tau = [0:2:Tau-1, -Tau:2:-1].';
     nu = [0:0.5/Nu:0.5-0.5/Nu, -0.5:0.5/Nu:-0.5/Nu];
    tau = tau/p(1);
     nu = nu/p(2);
%% MTEK(tau,nu) - matrix of Multiform Tiltable Exponential kernel -------%%
    if p(4)<1.5
        MTEK = tau.*nu;
    else
        MTEK = abs(tau).*abs(nu);
    end

    tau = tau.*tau;
     nu = nu.*nu;
     
    % alpha
%     tau = tau.*(nu.^p(6));
%      nu = nu.*(tau.^p(6));

    MTEK = 2*p(3)*MTEK;
    
    MTEK = MTEK + tau;
    MTEK = MTEK + nu;
    
    MTEK = MTEK.*MTEK;
    MTEK = MTEK.^p(5);
%     MTEK = MTEK.^(2*p(5));% Attention! Complex results. Never like this - root of negative!
    MTEK = -pi*MTEK;
    MTEK = exp(MTEK);

%% Optional output ------------------------------------------------------%%

end %[EOF makeMTEK.m]