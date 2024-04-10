
% This function generates a network of gratings with certain contrast levels.
% For a set of gratings forming an annulus, the contrast level is different
% from that of other gratings. The network would be presented for 25 different
% mixtures of contrast levels inside the annulust and network spacings to
% measure the detection threshold.
sca;
clear all;
clc
clear PsychHID; % Force new enumeration of devices.
clear KbCheck; % Clear persistent cache of keyboard devices.

% AssertOpenGL;

PsychDefaultSetup(2);
rand('seed', sum(100 * clock));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experimental Parameters

pix_per_deg=28.31; %34;  %% our screens have 34pix per degree
width=6.2*pix_per_deg; %% the figure is 6.2 degrees square
height=6.2*pix_per_deg;
KbName('UnifyKeyNames');
Key1 = KbName('RightArrow'); % yes key
Key2 = KbName('LeftArrow'); % no key
spaceKey = KbName('space');
escKey = KbName('ESCAPE');
corrkey = [Key1, Key2]; % key for 101
textsize = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% experiment variables

nBlocks = 30;
nTrialsPerBlock = 25;
interTrialInterval =  0.5; %sec
fix_window=2; % degrees, this much deviation is allowed from the fixation point
Nsample = nBlocks; %nBlocks; % number of samples(stimuli) for each condition
eye_tracking = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir1 = 'D:\Maryam\Experiment\input\stimuli\';
dir2 = 'D:\Maryam\Experiment\input\variables\';
% dir1 = 'P:\FSE_MACSBIO\maryam.karimian\Arnold tongue project\experiments&simulations\6_June_2019\input\stimuli\';
% dir2 = 'P:\FSE_MACSBIO\maryam.karimian\Arnold tongue project\experiments&simulations\6_June_2019\input\variables\';

qall   = linspace(1,1.5,5);
pall   = linspace(0,0.99,5)+0.01;
spacingCoef = zeros(nTrialsPerBlock,1);
InContrange = zeros(nTrialsPerBlock,1);
j = (0:length(pall)*length(qall)-1);
source = [floor(j/length(qall)); mod(j,length(pall))]+1;
for j = 1:length(pall)*length(qall)
    InContRange(j) = pall(source(2,j));
    spacingCoef(j) = qall(source(1,j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Login prompt and open file for writing out data

prompt = {'Subject Name', 'Session Number', 'Subject Number', 'age', 'gender'};
defaults = {'MaryamKarimian', '1', '1', '29', 'f'};
answer = inputdlg(prompt, 'MaryamKarimian', 2, defaults);
[subName, sessionid, subid, subage, gender] = deal(answer{:}); % all input variables are strings
out_dir = ['output\psychophysics\subject' num2str(subid) '_session' num2str(sessionid)];
if ~exist(out_dir,'dir'), mkdir(out_dir); end
outputname = fullfile(out_dir, [subName gender subid]);

if exist([outputname '.xls'],'file') % check to avoid overiding an existing file
    fileproblem = input('That file already exists! Append a .x (1), overwrite (2), or break (3/default)?');
    if isempty(fileproblem) || fileproblem==3
        return;
    elseif fileproblem==1
        outputname = [outputname '.x'];
    end
end
outfile = fopen([outputname '.xls'],'w'); % open a file for writing data out
fprintf(outfile, 'subid\t subage\t gender\t session\t blockNumber\t trialNumber\t ContrastSDev\t Spacing\t Response\t CorrectAns\t ReactionTime\t \n');

outfile_gaze = fopen([outputname '_gaze.xls'],'w');
fprintf(outfile_gaze, 'blockNumber\t trialNumber\t timeStamp\t gazeAngle\t gazeX\t gazeY\t \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experimental Instructions

text = sprintf(['In this experiment you are asked to FIXATE on the small dot in the center.\n' ...
    '\n Then 25 networks of gratings will be presented.'...
    '\n There is a rectangle with a more homogeneous contrast difference' ...
    '\n inside its elements, which appears in a random location of the network. \n' ...
    '\n\n Press "-->" if the rectangle is HORIZONTAL\n' ...
    ' Press "<--" if it is VERTIVAL \n\n' ...
    'Press any key to start the experiment.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block Loop


screen_setup;
%----------------------------------------------------------------------
%                       set up Eyelink
%----------------------------------------------------------------------

if eye_tracking
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    
    el=EyelinkInitDefaults(mainwin);
    
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    dummymode=0;
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    
    % open file to record data to
    Eyelink('openfile', ['data.edf']);
    
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % do a final check of calibration using driftcorrection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  EyelinkDoDriftCorrection(el);
    
    WaitSecs(0.1);
    Eyelink('StartRecording');
    
    eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
    if eye_used == el.BINOCULAR; % if both eyes are tracked
        eye_used = el.LEFT_EYE; % use left eye
    end
end

%% Experiment loop
index     = 0;

start_time = GetSecs;


for a = 1:nBlocks
    
    filename=[dir2,'orientation.mat'];
    load(filename,'vertical');

%     if mod(a,2)==1 && a>1
%         screen_setup;
%     end
    if a==1
        WaitSecs(1);
        Screen('TextSize', mainwin, textsize);
        DrawFormattedText(mainwin, text, 'center', 'center', textcolor,[],[],[],2);
        Screen(mainwin, 'Flip');
        
        KbWait();
    end
    %%%%%%%%%%%%%%%%%%%%%%
    
    trialorder = [1:nTrialsPerBlock];
    trialorder = trialorder(randperm(length(trialorder)));
    if (GetSecs-start_time)>(10*60)
        DrawFormattedText(mainwin, ['Take A Short Break'], 'center', yCenter-(pix_per_deg*2), black);
        DrawFormattedText(mainwin, [num2str(a-1) 'rounds are done'], 'center', 'center', black);
        DrawFormattedText(mainwin, ['Press Any Key To Restart'], 'center', yCenter+(pix_per_deg*2), black);
        Screen('Flip', mainwin);
        if eye_tracking
            Eyelink('StopRecording');
        end
        KbStrokeWait
        if eye_tracking
            EyelinkDoTrackerSetup(el);
            WaitSecs(0.1);
            Eyelink('StartRecording');
            WaitSecs(0.1);
        end
        start_time=GetSecs;
    end
    
    b=1; plotting_index = 1; count = 1; %current trial number to check
    
    orientation = vertical(:,a);
    
    while b <= length(trialorder) %&& reversal_point <= max_rps
        
        %%%%%% grey screen
        Screen('FillRect', mainwin, bgcolor); Screen(mainwin, 'Flip');
        WaitSecs(0.1); % 100 ms
        % Flip again to sync us to the vertical retrace at the same time as drawing our fixation point
        timeStart = GetSecs;
        while GetSecs <= timeStart +  1 %&& ~KbCheck
            dotColor = [0 1 1];
            draw_dot;
            Screen('Flip', mainwin);
        end
        if eye_tracking
            wait_fixation=1;
            while wait_fixation
                if Eyelink( 'NewFloatSampleAvailable') > 0
                    % get the sample in the form of an event structure
                    evt = Eyelink( 'NewestFloatSample');
                    if eye_used ~= -1 % do we know which eye to use yet?
                        % if we do, get current gaze position from sample
                        if (abs(evt.gx(eye_used+1)-xCenter)/pix_per_deg)<fix_window  & (abs(evt.gy(eye_used+1)-yCenter)/pix_per_deg)<fix_window
                            wait_fixation=0;
                        end
                    end
                end
                dotColor = [0 1 1];
                draw_dot;
                Screen('Flip', mainwin);
            end
        end
        
        timeStart = GetSecs;
        
        count = count+1;
        keyIsDown=0;
        rt=0;
        response = nan;
        keypressed = nan;
        
        break_fixation=0;

        while GetSecs <= timeStart +  1 && ~KbCheck
            
            Project_Stimulus;
            
            if eye_tracking
                if Eyelink( 'NewFloatSampleAvailable') > 0
                    % get the sample in the form of an event structure
                    evt = Eyelink( 'NewestFloatSample');
                    if eye_used ~= -1 % do we know which eye to use yet?
                        %%%%% if the eye position is beyond the fixation window stop the trial and play a error noise
                        if (abs(evt.gx(eye_used+1)-xCenter)/pix_per_deg)>fix_window  & (abs(evt.gy(eye_used+1)-yCenter)/pix_per_deg)>fix_window
                            break_fixation=1; %% use the break fixation flag to stop other loops running
                            beep
                            break
                        end
                    end
                end
            end
        end 
        
        if ~break_fixation
            
            timeStartx = GetSecs;
            while GetSecs <= timeStartx +  0.5 %%/ifi
                
                if eye_tracking
                    if Eyelink( 'NewFloatSampleAvailable') > 0
                        % get the sample in the form of an event structure
                        evt = Eyelink( 'NewestFloatSample');
                        if eye_used ~= -1 % do we know which eye to use yet?
                            % if we do, get current gaze position from sample
                            timeStamp = GetSecs-timeStartx; 
                            xGaze = evt.gx(eye_used+1);
                            yGaze = evt.gy(eye_used+1);
                            gazeAngle = ((evt.gx(eye_used+1)-xCenter)^2 + ...
                                (evt.gy(eye_used+1)-yCenter)^2)/pix_per_deg;
                            
                            index = index + 1;
                            
                            fprintf(outfile_gaze, '%d\t %d\t %d\t %8.3f\t %8.3f\t %8.2f\t \n', a, trialorder(b), ...
                                    timeStamp, gazeAngle, xGaze, yGaze);
                        end
                    end
                end
                
                % =========== get response
                
                [keyIsDown, secs, keyCode] = KbCheck;
                
                if keyIsDown
                    nKeys = sum(keyCode);
                    if nKeys==1
                        if keyCode(Key1)||keyCode(Key2)
                            rt = 1000.*(GetSecs-timeStart); %% I think here I should have used timeStartx-timestart! because i have also accounted for the systems' processing time.
                            keypressed=find(keyCode);
                            Screen('Flip', mainwin);
                            break % for outer while
                        elseif keyCode(escKey)
                            % fclose(outfile);
                            if eye_tracking
                                Eyelink('StopRecording');
                                Eyelink('closefile');
                            end
                            ShowCursor;
                            Screen('CloseAll');
                            sca;
                            return
                        end
                        keyIsDown=0;
                        keyCode=0;
                    end
                else
                    response = 1000;
                end
                
                % ======================================
            end
            
            if isnan(keypressed)
                trialorder(end+1)    = trialorder(b);
                InContRange(end+1)   = InContRange(trialorder(b));
                spacingCoef(end+1)   = spacingCoef(trialorder(b));
                orientation(end+1)   = orientation(trialorder(b));

            else
                
                if keypressed==corrkey(1)
                    response = 1;
                elseif keypressed==corrkey(2)
                    response = 0;
                end
                
                if orientation(trialorder(b))==1
                    corrAns = 0;
                elseif orientation(trialorder(b))==0
                    corrAns = 1;
                end
                
                if response == corrAns
                    dotColor = [0 1 0];
                else
                    dotColor = [1 0 0];
                end
                
                draw_dot;
                Screen('Flip', mainwin);
                WaitSecs(0.5);

                % writing data out
                fprintf(outfile, '%s\t %s\t %s\t %s\t %d\t %d\t %8.3f\t %8.3f\t %d\t %d\t %8.2f\t \n', subid, subage,...,
                     gender, sessionid, a, trialorder(b), InContRange(trialorder(b)), spacingCoef(trialorder(b)), response, corrAns, rt);
                %**********************************************************
                
                % store variable for plotting
                curr_cont(a,plotting_index) = InContRange(trialorder(b));
                curr_spc(a,plotting_index) = spacingCoef(trialorder(b));
                
                plotting_index = plotting_index+1;
            end
            
        else
            trialorder(end+1)  = trialorder(b);
            InContRange(end+1) = InContRange(trialorder(b));
            spacingCoef(end+1) = spacingCoef(trialorder(b));
            orientation(end+1)   = orientation(trialorder(b));

            
        end % fixation
        
        Screen('FillRect', mainwin ,bgcolor);    Screen(mainwin, 'Flip');
        WaitSecs(interTrialInterval);
        
        b = b+1;
    end % end of trial loop
    
    % looking for end
    if a == nBlocks
        DrawFormattedText(mainwin, 'End of experiment. Press any key to exit.', 'center', 'center',textcolor);
        Screen(mainwin, 'Flip');
        KbWait([], 3); %wait for keystroke
        
    else
        DrawFormattedText(mainwin, ['End of the block ', num2str(a), '. Press any key to continue.'], 'center', 'center',textcolor);
        Screen(mainwin, 'Flip');
        KbWait(); %wait for keystroke
        
    end
%     if mod(a,2)==0 
%         Screen('CloseAll'); 
%     end
end % end of block loop

Screen('CloseAll');

fprintf('\n\n\n\n\nFINISHED!...\n\n');

% time_elapsed = toc;
if eye_tracking
    Eyelink('StopRecording');    Eyelink('closefile');

    %% gaze plots
%     figure;
%     [meanX, meanY] = gazePlots(xGaze, yGaze, nBlocks, nTrialsPerBlock); %length(curr_cont));
%     savefig([outputname '_gazePlots.fig']);
%     fprintf(outfile_gaze, '\n %s\t %dt \n', 'Mean X', meanX');
%     fprintf(outfile_gaze, '\n %s\t %d\t \n', 'Mean Y', meanY);
%     fclose(outfile_gaze);
end


