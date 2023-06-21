function [tfr, varargout] = optimTFR_cpu(s, nf, t, opts)
% Optimise time-frequency representation for a given signal
%
% Inputs:
%   required:
%       s - input signal, real or analytic
%   optional:
%       Nf - number of frequency bins, max number of lags
%       ti - time instants for which to calculate tfr
%       Np - population size used in optimisation, default Np=5
% Outputs:
%   required:
%       tfr     - optimised time-frequency representation, NfxNt array
%   optional:
%       z       - analytic associate of input signal
%       oPar    - optimised kernel parameters
%       qm      - tfr performance measure
%       kernel  - nu-tau domain zero-phase NfxNt kernel array
%       mask    - mask of wanted correlations, NfxNt logical array
%       outcome - evolution progress, objective function space
%       outputs - evolution progress, best parameter evolution
% Example:
%  load testSignals;
%  tfr = optimTFR_cpu(sigLFMX);
%  figure; mesh(tfr); xlabel('t'); ylabel('f'); colormap(jet);
%
% Externals:
%  nextname - uses Next Available Filename by Stephen Cobeldick available from https://github.com/mykappa/nextname
%
%
% 2023-06-19 vMDPISensors
%   (C) Damir Malnar 2023. Supplementary materials to:
%   Citation: Malnar, D.; Vrankic, M. Optimising Time-Frequency Distributions: A Surface Metrology Approach. Sensors 2023, 1, 0. https://doi.org/

%% ------------------ Input check & Defaults ----------------------------%%

    arguments
        s  {mustBeVector(s), mustBeNumeric(s), mustBeFinite(s), mustBeNonsparse(s), mustBeNonmissing(s)}
        nf (1,1) {mustBePositive(nf), mustBeInteger(nf)} = 2^ceil(log2(numel(s))) 
        t  {mustBeVector(t), mustBeInteger(t), mustBeNonsparse(t), mustBeNonmissing(t), mustBeInSignalRange(t,s)} = 1:numel(s)
        opts.PopSize  {mustBeInteger(opts.PopSize), mustBeGreaterThanOrEqual(opts.PopSize,4)} = 5
        opts.TestName {mustBeNonzeroLengthText(opts.TestName)} = [datestr(now, 'yyyy-mm-dd') '-testData']
        opts.FileSfx  {mustBeNonzeroLengthText(opts.FileSfx)} = '-1'
        opts.FileExt  {mustBeNonzeroLengthText(opts.FileExt), mustBeMember(opts.FileExt,['.mat','.txt','.xls'])} = '.mat'
        opts.TestSave {mustBeInteger(opts.TestSave), mustBeMember(opts.TestSave,[0 1])} = 0
        opts.StatName {mustBeNonzeroLengthText(opts.StatName)} = [datestr(now, 'yyyy-mm-dd') '-optiData']
        opts.StatSave {mustBeInteger(opts.StatSave), mustBeMember(opts.StatSave,[0 1])} = 0
    end
    
    global Ns Nf
    % Signal length
    Ns = numel(s);
    % lag-to-FFT bins
    Nf = nf;
    % Population size
%     Np = opts.PopSize;

%% ------------------ Signal preprocess ---------------------------------%%

    % Make signal Ambiguity function and return relevant parameters
    global CAF WDF Nt Nu z
    
    [CAF,WDF,~,Nf,Nt,Nu,~,z] = makeAFv2(s,Nf,t);
        
    % Signal energy
    global sigNrg
    sigNrg = sum(abs(z).^2);

    % Signal instantaneous energy from WVD
    global sigPwr
    sigPwr = sum(WDF,1);

    % Signal PSD from WVD
    global sigPSD
    sigPSD = sum(WDF,2);    
    
    global truAvg
    truAvg = sigNrg/Ns;

%% ------------------ DE engine setup -----------------------------------%%

    %SPWVD bound
if 0
    LU = [ 1, 0,  1, 0;
          Nt, 8, Nf, 8]; 
    lu = LU;
end

    %Spectrogram bounds
if 0
    LU = [ 1, 0;
          Nt, 8]; 
    lu = LU;
end

    % MTEK parameters bounds, no alpha
if 1 
    LU = [ eps(1), eps(1), -1.5, 1.0, 0.5;
               Nf,    0.5, +1.5, 2.0, 1.0]; 
    lu = [LU(1,1), LU(1,2), LU(1,3),  LU(1,4), eps(1);
          LU(2,1), LU(2,2), LU(2,1),  LU(2,4), LU(2,1)];
%         flintmax, flintmax, flintmax, LU(2,4), flintmax];
end    

%% ------------------ DE optimiser --------------------------------------%%

    % Return optimised parameters,optimisation statistics and final fitness
    [oPar, genStat] = deEngineS(LU,lu,opts.PopSize);

%% ------------------ Results -------------------------------------------%%

   % Make kernel using optimised parameters 
   kernel = makeMTEKv2(oPar, Nf, Nt);

   % Generate tfr using optimised kernel
   tfr = real(fft(ifft(kernel.*CAF,[],2),[],1));
   tfr = tfr(:,1:Nt);

%% ------------------ Formulate output & Print results ------------------%%
    
    varNames = {'Generations','OptimTime','AvgGenTime','BestFit','BestParameters'};
    optimStat = table(numel(genStat.No),...
                      sum(genStat.Wasted),...
                      mean(genStat.Wasted),...
                      genStat.BestFit(end),...
                      genStat.BestTarget(end,:),...
                      'VariableNames',varNames);
    % Optional outputs
    varargout{1} = optimStat;
    varargout{2} = kernel;
    varargout{3} = CAF;
    varargout{4} = WDF;
    varargout{5} = z;
    varargout{6} = t;
    varargout{7} = genStat;
%     varargout{8} = Nt;
%     varargout{9} = Ns;
    
    % Optionally save optimisation data
    if (isfield(opts,"TestSave") && (opts.TestSave == 1))
        fileName = nextname(opts.TestName,opts.FileSfx,opts.FileExt,true);
        save(fileName,'optimStat','kernel','CAF','WDF','z','t','genStat','Nt','Nf','Ns','oPar','tfr');
    end
    % Dangerous code here
    if (isfield(opts,"StatSave") && (opts.StatSave == 1))
        fileName = [opts.StatName opts.FileSfx opts.FileExt];        
        if isfile(fileName)
            fileObj = matfile(fileName,'Writable',true);
            fileObj.optimData = [fileObj.optimData;optimStat];                               
        else
            fileName = nextname(opts.StatName,opts.FileSfx,opts.FileExt,true);
            fileObj = matfile(fileName,'Writable',true);
            fileObj.optimData = optimStat;           
        end
    end    

    % Report optimised parameter values
%     printOptimisedParams(oPar, qm);
end % end optimTFR_cpu

%% ------------------ Helper functions ----------------------------------%%

function printOptimisedParams( params, qm )
% Prints optimised parameters formated

    fprintf('\nOptimised MTEK parameters:');
    fprintf('\n  tau0: %-11.9f',params(1));
    fprintf('\n   nu0: %-11.9f',params(2));
    fprintf('\n     r: %-11.9f',params(3));
    fprintf('\n  beta: %-11.9f',round(params(4)));
    fprintf('\n gamma: %-11.9f',1/round(params(4)));
    fprintf('\nlambda: %-11.9f',params(5));
    fprintf('\n\nproducing tfr having:');
    fprintf('\n     dispersion of wanted correlations: %-11.9f',qm(1));
    fprintf('\nand');
    fprintf('\nconcentration of unwanted correlations: %-11.9f',qm(2));
    fprintf('\nwith');
    fprintf('\n          agregated performace measure: %-11.9f',sum(qm,2));

end %end printOptimisedParams

% Custom validation function
function mustBeInSignalRange(a,b)
    % Test for range
    lowerBound = 1;
    upperBound = numel(b);
    if (a(1)<lowerBound && a(end)>upperBound)
        eid = 'SignalRange:notInSignalRange';
        msg = 'Time indices must be in signal length range.';
        throwAsCaller(MException(eid,msg))
    end
end
