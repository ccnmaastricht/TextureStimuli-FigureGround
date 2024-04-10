
% while GetSecs <= timeStart +  1 && ~KbCheck
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', mainwin);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(screenRect);
    
    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', mainwin, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
%     filename = [dir,'texture_cont_',num2str(trialorder(b)),'.png']
%     Stimuli = imread(filename);
%     imageTexture = Screen('MakeTexture', mainwin, Stimuli);
% f = 1;
    
%     filename=[dir,'texture_trial',num2str(trialorder(b)),'_block',num2str(f(b)),'.mat'];
    filename=[dir1,'texture_trial',num2str(trialorder(b)),'_block',num2str(a),'.mat'];
    load(filename,'imageArray')
    imageTexture = Screen('MakeTexture', mainwin , imageArray);%% put the images on a screen
    
    Screen('DrawTexture', mainwin, imageTexture, [], [], 0);
%     Screen('DrawTexture', mainwin, imageTexture, [], [xCenter-(screenRect(3)/2)...
%         yCenter-(screenRect(4)/2) xCenter+(screenRect(3)/2) yCenter+(screenRect(4)/2)], 0);
    
    dotColor = [0 1 1];
    draw_dot;
    
    % Flip to the screen
    Screen('Flip', mainwin);

% end