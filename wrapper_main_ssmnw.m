%WRAPPER_MAIN_SSMNW Runs function main_ssmnw.m
%
%   @INPUT
%   None. Hard coded script.
%   But do check that the file name and path for the paramter matrix is
%   correct. And check that the saveDir path is correct.
%
%   @OUTPUT
%   Saves the paramStruct to the created directory saveDir. 
%   Otherwise, the function simply runs main_ssmnw, which in turn generates outputs.
%
%   @AUTHORS
%   Alva Curtsdotter, Post doc @ BrosiLab, Dep of Environmental Sciences,
%   Emory University, Atlanta, Georgia, USA. Code initiated 2018-06-05.

% -------------------------------------------------------------------------

fprintf('\n\n>>> WRAPPER STARTS AT %s <<<\n', datestr(now))
wrapperStart = tic;

% Simulation parameters
paramStructName = 'paramStruct_SSMNW.mat';
load(paramStructName)
% numScen = size(params.data,1); %Number of scenarios; number corresponding to rows in parameter matrix. If set to 0, default values defined in msf_setup_ssmnw are used.
numScen =15; 

% Create directory to save results
saveDir = 'C:\Users\Caleb\Desktop\Research\Nestedness Matrices\';
timeStamp = datestr(now,'yyyymmmddTHHMM'); % '2018Jul13T0105';%
saveDir = [saveDir, 'testresults_', timeStamp, '\'];
if 7 ~= exist(saveDir, 'dir')
    mkdir(saveDir)
end


% Run main_ssmnw.
if numScen == 0
    fprintf('\n\n>>> RUNNING main_ssmnw WITH DEFAULT VALUES <<<\n')
    tic
    p = 0;
    main_ssmnw( p, saveDir )
    duration = toc;
    fprintf('\n\n>>> DURATION: %d <<<\n', duration)
else
    %save([saveDir, paramStructName], 'params')
    save([saveDir, paramStructName], 'params')
    fprintf('\n\n>>> RUNNING main_ssmnw, %i SCENARIOS & %i REPLICATES <<<\n', numScen, params.data(1,params.trials))
%     tic
    for p = 1:numScen%783:numScen%10%
        scenarioStart = tic;
        %for rep = 1:1 %Reps are run from within main.
        main_ssmnw( p, saveDir) % main_ssmnw( p, saveDir, 5)
        %        end
        duration = toc(scenarioStart);
        fprintf('\n\n>>> %i REPLICATES FOR SCENARIO %i TOOK %d <<<\n', params.data(1,params.trials), p, duration)
    end
end
% duration = toc;
% fprintf('\n\n>>> DURATION: %d <<<\n', duration)

fprintf('\n\n>>> WRAPPER ENDS AT %s <<<\n', datestr(now))
toc(wrapperStart)
