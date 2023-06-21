function [parBest, varargout] = deEngineS(LU,lu,varargin)
% Differential Evolution optimiser according to Damir's strategy.
% Implemented as async evolution i.e. better trial replaces its target and
% immediately contributes in generating another trial during the same
% generation. Selection is by principle of dominance.
% Initial population sampled from Latin hypersphere giving a more even
% coverage of the solution space thus avoiding initial bias.
%
% parBest - best parameter set returned
% Parameters are handled directly in absolute values, no need for
% denormalisation before evaluation.
% LU - lower and upper parameter initialisation bounds. LU can not have
% -Inf/Inf. LU size is 2 x number of parameters.
% lu - lower and upper parameter bounds in "<" or ">" format.
% Unbounded parameter are indicateda as -Inf/Inf. lu size equals LU size.
% Np - population size
%
%
% 2023-06-19 vMDPISensors
%   (C) Damir Malnar 2023. Supplementary materials to:
%   Citation: Malnar, D.; Vrankic, M. Optimising Time-Frequency Distributions: A Surface Metrology Approach. Sensors 2023, 1, 0. https://doi.org/

%% ------------------ Input check & Defaults ----------------------------%%

    % Number of inputs must be >=minargs and <=maxargs.
    narginchk(2, 3)

    % Set defaults for optional inputs: population size
    optargs = {5};

    optargs(1:length(varargin)) = varargin;

    Np = optargs{:};
    clear optargs

%% ------------------ Setup ---------------------------------------------%%

    % Dimension of the problem
    [~, D] = size(LU);
    % Initialisation/Denormalisation bounds
        % Provided outside
    % Relative bounds
        % Provided outside
    % Population size
        % Provided outside

    % Auxiliary variables
%     outcome = [];
%     outputs = [];
    genMax = uint32(100000);
    genNum = uint32(0);

    % Setup RNG
    rng('shuffle');

%% ------------------ Implementation ------------------------------------%%

    % Overall header
    printOptimisationHeader();
    % Header for progress reporting
    printProgressHeader();

    % Initial population, using LHS to more evenly cover the
    % search space. Rows are individual members, columns are parameters.
    Target = lhsdesign(Np, D, 'criterion', 'maximin');
    Target = bsxfun(@times, LU(2,:)-LU(1,:), Target);
    Target = bsxfun(@plus, LU(1,:), Target);

    % Bulk evaluate initial target population on CPU
    tSG = tic;
    valTarget = cell2mat(arrayfun(@(ii) ...
        objectFun_cpu(Target(ii,:)),...
        (1:Np)','UniformOutput',false));
    tEG = toc(tSG);

    % Sort to obtain the values and indices of the best.
    [valBest, ixBest] = sort(valTarget,1,'ascend');
    
    % Objective function variance
    VarianceFit = var(valBest);
    
    % Progress tracking table
    varNames = {'No','Wasted','AverageFit','VarianceFit','WorstFit','BestFit','AverageTarget','DeviationTarget','WorstTarget','BestTarget'};
    generations = table(genNum,tEG,mean(valBest),VarianceFit,valBest(end),...
                        valBest(1),mean(Target),std(Target),...
                        Target(ixBest(end),:),Target(ixBest(1),:),...
                        'VariableNames',varNames);
                   
    % Progress report for 0th generation
    printProgressReport(genNum,valBest(1,1),mean(valBest(:,1)),valBest(end,1),tEG,Target(ixBest(1),:));

    % Search space - max and min difference vector
    deltaMinMax = mink(Target,2);
    % Minimum current difference vector
    deltaMinMax(2,:) = deltaMinMax(2,:)-deltaMinMax(1,:);
    % Maximum current difference vector
    deltaMinMax(1,:) = maxk(Target,1)-deltaMinMax(1,:);
    
    % Set all vector elemetns not feasible
    O = zeros(1,D,'logical');

    % Increment generation counter
    genNum = genNum+1;

    % Start evolution timer
    tSE = tic;
    % Start DE and run until NOT stopping conditions
    while (VarianceFit>1e-9) && (genNum<genMax)

        % Start generation timer
        tSG = tic;
        for ix = 1:Np
         
            % Async population model - local selection
            Trial = Target(ix,:);
            % mutation only
            O = ~O;
            while any(O)
                % Create step
                delta = sign(2*rand(1,sum(O)) - 1);
                delta = delta.*(deltaMinMax(2,O)+(deltaMinMax(1,O)-deltaMinMax(2,O)).*rand(1,sum(O)));
                % Create Trial solution
                Trial(O) = Target(ix,O) + delta;
                % Check if feasible
                O = Trial<lu(1,:) | Trial>lu(2,:);
            end
            % Evaluate Trial
            valTrial = objectFun_cpu(Trial);
            % Selection with the current Target
            if valTrial < valTarget(ix)
                % Replacement
                Target(ix,:) = Trial;
                valTarget(ix) = valTrial;
                % Update search space
                deltaMinMax = mink(Target,2);
                deltaMinMax(2,:) = deltaMinMax(2,:)-deltaMinMax(1,:);
                deltaMinMax(1,:) = maxk(Target,1)-deltaMinMax(1,:);
                % Track progess - not part of the optimisation algorithm
                [valBest, ixBest] = sort(valTarget,1,'ascend');
                % Objective function variance
                VarianceFit = var(valBest);
            end

        end % End current generation creation
        % Stop generation timer
        tEG = toc(tSG);

        % Update progress tracking table
        generationNow = table(genNum,tEG,mean(valBest),VarianceFit,valBest(end),...
                           valBest(1),mean(Target),std(Target),...
                           Target(ixBest(end),:),Target(ixBest(1),:),...
                           'VariableNames',varNames);        
        generations = [generations;generationNow];
        
        % Report progress for current generation
        printProgressReport(genNum,valBest(1,1),mean(valBest(:,1)),valBest(end,1),tEG,Target(ixBest(1),:));

        % Increment generation counter
        genNum = genNum+1;

    end % end DE
    % Stop evolution timer
    tEE = toc(tSE);

%% ------------------ Results------------------------------------------- %%

    % Best parameters returned
    parBest = Target(ixBest(1),:);

    % Optimisation progress, optional
    varargout{1} = generations;
%     varargout{2} = [];
%     varargout{3} = [];

    % Report optimisation ending
    printOptimisationEnding(tEE, Np, valBest(1));

end % end deEngine

%% Helper functions

function printProgressReport( generation, best, average, worst, tElapsed, bestTarget )
% Prints optimisation progress

    fprintf('\n%11u%3s%11.9f%3s%11.9f%3s%11.9f%3s%11.9f%3s%11.4f %11.4f %11.4f %11.4f %11.4f',...
        generation,' | ',...
        best,' | ',...
        average,' | ',...
        worst,' | ',...
        tElapsed,' | ',...
        bestTarget);

end % end printProgress

function printProgressHeader( )
% Prints progress header

    fprintf('\n%11s%3s%11s%3s%11s%3s%11s%3s%11s%3s%11s',...
        ' Generation',' | ',...
        '    Best   ',' | ',...
        '  Average  ',' | ',...
        '   Worst   ',' | ',...
        'Time       ',' | ',...
        'Parameters ');

end % end printProgressHeader

function printOptimisationHeader( )
% Prints optimisation header

    fprintf('\n--- DE Optimiser v13.0 by Damir - Serial Evolution ---\n');
    fprintf('\t\tPlease be patient...we have a lot to do...');
    fprintf('\nNOW OPTIMIZING...This may take a while...');
    fprintf('\nCoffee anyone?...\n');

end % end printOptimisationHeader

function printOptimisationEnding( tElapsed, Np, valEnd  )
% Prints optimisation ending

    fprintf('\n\nOPTIMIZED as best I could.');
    fprintf('\nUsing %u population members it took me %f seconds.', Np, tElapsed);
    fprintf('\nReached final fitness of: %11.9f\n\n', valEnd);

end % end printOptimisationEnding
