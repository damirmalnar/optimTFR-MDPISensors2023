function f = objectFun_cpu(p)
%OBJECTFUN_CPU Objective function that implements the proposed TFD performance measure
%   Detailed explanation goes here TODO
%
% Inputs:
%   p - 1xD vector of decision parameters
% Outputs:
%   f - 1x1 tfr performance measure
%
% Externals:
%  tfrsp - Spectrogram, requires tftb toolbox available from: https://tftb.nongnu.org/
%
%
% 2023-06-19 vMDPISensors
%   (C) Damir Malnar 2023. Supplementary materials to:
%   Citation: Malnar, D.; Vrankic, M. Optimising Time-Frequency Distributions: A Surface Metrology Approach. Sensors 2023, 1, 0. https://doi.org/

%% Common data

    global CAF Nt Nf truAvg z WDF

%% Make TFR
if 1
    % MTEK(tau,nu) kernel
    TFR = makeMTEKv2(p,Nf,Nt);
    TFR = TFR.*CAF;
    TFR = ifft(TFR,[],2);
    % smoothed TFD(f,t);
    TFR = real(fft(TFR(:,1:Nt),[],1));
end

% SPWVD
if 0
    g = makeWindow('hamming',p(1)).';
    h = makeWindow('hamming',p(3)).';
%     g = makeWindow('kaiser',p(1:2)).';
%     h = makeWindow('kaiser',p(3:4)).';
    TFR = tfrspwv([zeros(Nt/2,1);z;zeros(Nt/2,1)],1:2*Nt,Nf,g,h,0);%
    TFR = TFR(:,Nt/2+1:end-Nt/2);
end

% Spectrogram
if 0
    h = makeWindow('kaiser',p).';
%     h = makeWindow('hamming',p).';
    TFR = tfrsp(z,1:Nt,2*Nf,h,0); 
    TFR = TFR(1:Nf,:)/2;
end
%% ------------------------------------------------------------------------
if 1 % Proposed measure
    H = TFR(:);
    hBar = mean(H);
    ix = H>0;
    Vm = sum(H(ix));
    
    H = H-hBar;
    % Deal with roundoff-noise
    ix = abs(H)<=1e-13;
    H(ix) = 0.0;
    
    ix = H>=0;
%     Smr2 = sum(ix); % just for plot
    Vmp2 = sum(H(ix));
    c1 = abs(min(H));
    ix = H>=c1;
    Smr1 = sum(ix);    
    Vmp1 = sum(H(ix));

    p5 = abs(1.0-hBar/truAvg);
    p4 = 1.0-Vmp2/Vm;
    p3 = 1.0-Vmp1/Vm;
%     p2 = Smr2/numel(H); % just for plot
    p1 = Smr1/numel(H);
    
    f = p1 + p3 + p4 + p5;
    f = 0.25*f;
end

if 0 % Stankovic energy concentration measure
    ECM  = ((sum(realsqrt(abs(TFR/sum(TFR,'all'))),'all'))^2)/numel(TFR); % Stankovic
    f = ECM;
end
if 0 % volume normalised Renyi entropy
    RE3DV = log2(sum((TFR/sum(abs(TFR),'all')).^3,'all'))/(1-3); % Renyi a=3
%     RE3DV = log2(sum((TFR/sum(abs(TFR),'all')).^7,'all'))/(1-7); % Renyi a=7
    f = RE3DV;
end
end %[EOF objectFun_cpu.m]
