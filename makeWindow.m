function [ w, M ] = makeWindow( window, p )
%MAKEWINDOW Makes symmetric window function for given window name, length 
%and optional parameters. Handles odd and even lengths.
%   Detailed explanation goes here
%
%   TODO: add additional classic windows, add ultraspeherical window
%
% Inputs:
%   required:
%       window - window name, any of the supported windows:
%                rectangular, hann, hamming, blackman,
%                parametric blackman (pblackman), blackman-nuttall
%                (blacknuttall), blackman-harris (blackhariss), gauss,
%                approximate confined gauss (gaussac), kaiser, hann-poisson
%                (hannpoisson), planck-bessel (planckbessel)
%       p - vector of window parameters: window length (N), p1, p2,...
%   optional:
%        - 
% Outputs:
%   required:
%       w - window, Nx1 column vector
%   optional:
%       M - window half
% 
% 2023-06-19 vMDPISensors
%   (C) Damir Malnar 2023. Supplementary materials to:
%   Citation: Malnar, D.; Vrankic, M. Optimising Time-Frequency Distributions: A Surface Metrology Approach. Sensors 2023, 1, 0. https://doi.org/

%% Input check & Defaults -----------------------------------------------%%

    % Parse the inputs
    window = lower(window);
    % Cast to enforce Precision Rules on window parameters
    p = signal.internal.sigcasttofloat(p,...
        'double','makeWindow','P','allownumeric');
    if p(1) <= 0
        error('N must be greater than zero');
    end
    if odd(p(1)) == 1 % trivial window
        w = 1;
        return;
    end
    
%% Setup window calculation
       
%     N = p(1);                %window length
    N = odd(p(1));                %window length
    M = fix(0.5*(N+1));      %window half
    oddN = rem(N,2);
    w = zeros(1,N,'double');
    xi = 0:M-1;
    
%% Calcutate symmetric window -------------------------------------------%%

    switch window
        case 'rectangle' % Q.C. Pass
            w = ones(1,N,'double');
        case 'hann' % Q.C. Pass
            xi = xi/(N-1);
            w(1:M) = 0.5 - 0.5*cos(2*pi*xi);
            w(M+1:N) = w(M-oddN:-1:1);
        case 'hamming' % Q.C. Pass
            xi = xi/(N-1);
            w(1:M) = 0.54 - 0.46*cos(2*pi*xi);
            w(M+1:N) = w(M-oddN:-1:1);
        case 'blackman' % Q.C. Pass
            xi = xi/(N-1);
            w(1:M) = 0.42 - 0.5*cos(2*pi*xi) + 0.08*cos(4*pi*xi);
            w(1) = 0;
            w(M+1:N) = w(M-oddN:-1:1);
        case 'pblackman' % parametric blackman % 0<=p(2)<=1 % Q.C. Pass 
            xi = xi/(N-1);
            w(1:M) = 0.5*(1-p(2)) - 0.5*cos(2*pi*xi) + 0.5*p(2)*cos(4*pi*xi);
            w(1) = 0;
            w(M+1:N) = w(M-oddN:-1:1);       
        case 'blacknuttall' % Q.C. Pass
            xi = xi/(N-1);
            w(1:M) = 0.3635819 - 0.4891775*cos(2*pi*xi) +...
                     0.1365995*cos(4*pi*xi) - 0.0106411*cos(6*pi*xi);
            w(M+1:N) = w(M-oddN:-1:1);
        case 'blackharris' % Q.C. Pass
            xi = xi/(N-1);
            w(1:M) = 0.35875 - 0.48829*cos(2*pi*xi) +...
                     0.14128*cos(4*pi*xi) - 0.01168*cos(6*pi*xi);
            w(M+1:N) = w(M-oddN:-1:1);
        case 'gauss' % p(2)=(N-1)/(2*sigma_t), sigma_t=0.X*N % Q.C. Pass
            xi = xi/xi(M);
            w(M+1-oddN:N) = exp(-0.5*((p(2)*xi).^2));
            w(1:M-oddN) = w(N:-1:M+1);
        case 'gaussac' % (N-1)/(2*0.1*N)<=p(2)<=(N-1)/(2*0.145*N) %Q.C. Pass
            xi = xi/xi(M);
            G = @(x) exp(-0.5*((p(2)*x).^2));
            w(M+1-oddN:N) = G(xi) - G(-(2/(M-1)))*(G(xi(M:-1:1)+1));
            w(1:M-oddN) = w(N:-1:M+1);
        case 'kaiser' % p(2)>=0 % Q.C. Pass 
            xi = 4*((xi + 0.5*(1-oddN))/(N-1)).^2;
            bes = abs(besseli(0,p(2)));
            w(M+1-oddN:N) = abs(besseli(0,p(2)*sqrt(1-xi))/bes);
            w(1:M-oddN) = w(N:-1:M+1);
        case 'hannpoisson' % p(2)>=0  % Q.C. Pass
            xi = xi/(N-1);
            w(1:M) = (0.5 - 0.5*cos(2*pi*xi)).*exp(-p(2)*(1-2*xi));
            w(M+1:N) = w(M-oddN:-1:1);                
        case 'planckbessel' % 0<p(2)<=1, p(3)>=0 % Q.C. Pass 
            ix = ceil(p(2)*M)-1;
            xiTmp = 2*(xi(1:ix)/(N-1))-1;
            w(1:ix) = 1./(exp(2*p(2)*(1./(1+xiTmp) + 1./(1-2*p(2)+xiTmp)))+1);
            w(ix+1:M) = 1;
            w(M+1:N) = w(M-oddN:-1:1);
            xi = 4*((xi + 0.5*(1-oddN))/(N-1)).^2;
            bes = abs(besseli(0,p(3)));
            w(M+1-oddN:N) = w(M+1-oddN:N).*abs(besseli(0,p(3)*sqrt(1-xi))/bes);
            w(1:M-oddN) = w(N:-1:M+1);
        otherwise
            error('Unknown window name!');
    end
    
    % Normalise on max just in case
    w = w/w(M);

end %[EOF makeWindow.m]

