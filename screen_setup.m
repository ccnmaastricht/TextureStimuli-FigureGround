
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Screen Parameters

% == colors of screen == %
screenNumber=max(Screen('Screens'));

% Find the color values which correspond to white and black.
white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);

% Round gray to integral number, to avoid roundoff artifacts with some
% graphics cards:
gray=round((white+black)/2);

% This makes sure that on floating point framebuffers we still get a
% well defined gray.
if gray == white
    gray=white / 2;
end

bgcolor = gray; 
textcolor = white;
% [mainwin, screenRect] = PsychImaging('OpenWindow', screenNumber, gray);

% == Open a double buffered fullscreen window with a gray background == %
% [mainwin, screenRect]=Screen('OpenWindow',screenNumber, bgcolor);
[mainwin, screenRect] = PsychImaging('OpenWindow', screenNumber, bgcolor);
[xCenter, yCenter] = RectCenter(screenRect);
center = [screenRect(3)/2 screenRect(4)/2];
HideCursor(mainwin);
% Screen(mainwin, 'Flip');
