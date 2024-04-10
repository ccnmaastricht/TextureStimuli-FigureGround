

[xCenter, yCenter] = RectCenter(screenRect);

% dotColor = [0 1 0];

dotSizePix = 10;

Screen('DrawDots', mainwin, [xCenter, yCenter], dotSizePix, dotColor, [], 2);
