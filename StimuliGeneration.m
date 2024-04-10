% This function generates 25 network of gratings with 5 different contrast
% level and 5 different spacing.
close all;
clearvars;
clear all;
sca;

dir1 = 'P:\FSE_MACSBIO\maryam.karimian\Arnold tongue project\experiments&simulations\6_June_2019\input\contrast\';
dir2 = 'P:\FSE_MACSBIO\maryam.karimian\Arnold tongue project\experiments&simulations\6_June_2019\input\spacing\';
dir3 = 'P:\FSE_MACSBIO\maryam.karimian\Arnold tongue project\experiments&simulations\6_June_2019\input\stimuli\';
dir4 = 'P:\FSE_MACSBIO\maryam.karimian\Arnold tongue project\experiments&simulations\6_June_2019\input\variables\';

% Luminannce correction
x     = linspace(0,1,15)';
y1234 = [.32 1.17 4.05 10 20.4 33 47 56 78.3 100 119 142 164 178 191;...
    .33 1.61 5.96 9.55 17.1 27.8 41.5 60.6 79.5 100 118 134 155 169 182;...
    .33 1.56 5.52 12 21.8 30.9 40.9 51.9 68.3 87.5 114 137 156 172 183;...
    .33 1.39 4.98 11.4 20.9 31 43 58.6 78.7 96.5 115 133 153 167 177];
y     = mean(y1234,1)';
p = polyfit(x,y,3);
RGB2Scr = @(x) p(1)*(x.^3)+p(2)*(x.^2)+p(3)*x+p(4);
p = polyfit(y,x,3);
Scr2RGB = @(x) p(1)*(x.^3)+p(2)*(x.^2)+p(3)*x+p(4);
black = max(0,Scr2RGB(min(y)));
white = min(1,Scr2RGB(max(y)));
grey  = Scr2RGB((max(y)-min(y))/2);
inc = white - grey;


rand('seed', sum(100 * clock));
% annulusDimPix = 30;
CircleSizemm = 7;
frequency = 2;
offset    = pi;
InContRange = 0.01;
nblocks = 30; 
prelearningSteps   = 1;
% Number of gratings in the length and width of the biggest possible rectangle.
Length = 4; %Only even numbers. note that (StimLength/2)<=(N/2)-barLoc
Width  = 2; %Only even numbers. note that (StimWidth/2)<=(N/2)-barLoc
courserAreaRmm = 70;
dev        = .3;
qall       = linspace(1,1.5,5);
pall       = linspace(0,0.99,5)+0.01;
vertical   = ones(numel(pall)*numel(qall),nblocks)*2;

for l = 1:1%nblocks
    k = 0;
    for q = 1:5
        spacingCoef = qall(q);
        for p = 1:5
            InContRange = pall(p);
            k = k + 1;
            
            PsychDefaultSetup(2);
            
            Screen('Preference', 'SkipSyncTests', 1);
            
            screens = Screen('Screens');
            
            screenNumber = max(screens);
            
            % %             white = Scr2RGB(WhiteIndex(screenNumber));
            % %             black = Scr2RGB(BlackIndex(screenNumber));
            % %             grey  = white/2;
            %             white = WhiteIndex(screenNumber);
            %             black = BlackIndex(screenNumber);
            %             x = linspace(0,1,256);
            %             gr = (max(RGB2Scr(x)) - min(RGB2Scr(x))) /2;
            %             [~,idx] = min(abs(RGB2Scr(x) - gr));
            %             grey = x(idx);
            %             inc = white - grey;
            
            [mainwin, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
            %             [screenXpixels, screenYpixels] = Screen('WindowSize', mainwin);
            %             [screenXmm, screenYmm] = Screen('DisplaySize', mainwin);
            screenXpixels = 1280;
            screenYpixels = 1024;
            screenXmm = 452;
            screenYmm = 361;
            ifi = Screen('GetFlipInterval', mainwin);
            [xCenter, yCenter] = RectCenter(windowRect);
            %             xCenter = 640;
            %             yCenter = 512;
            Screen('BlendFunction', mainwin, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
            PixelsPermm = screenXpixels/screenXmm;
            annulusDimPix = CircleSizemm*PixelsPermm;
            courserAreaRpixels = courserAreaRmm*PixelsPermm;
            
            % Define the Grating
            r         = linspace(-1,1,2*annulusDimPix+1);
            [x, y]    = meshgrid(r);
            R         = sqrt(x.^2 + y.^2);
            
            gridDim = spacingCoef * annulusDimPix;
            X = [xCenter];
            Xfill = 0;
            i=1;
            while Xfill<(screenXpixels-annulusDimPix)/2
                X = [X,X(i)+gridDim];
                Xfill = Xfill + gridDim;
                i = i+1;
            end
            DeviationFromCenterX = i;
            Xfill = 0;
            while Xfill<(screenXpixels-annulusDimPix)/2
                X = [X(1)-gridDim,X];
                Xfill = Xfill + gridDim;
                i = i+1;
            end
            NX = i;
            
            Y = [yCenter];
            Yfill = 0;
            j=1;
            while Yfill<(screenYpixels-annulusDimPix)/2
                Y = [Y,Y(j)+gridDim];
                Yfill = Yfill + gridDim;
                j = j+1;
            end
            DeviationFromCenterY = j;
            Yfill = 0;
            while Yfill<(screenYpixels-annulusDimPix)/2
                Y = [Y(1)-gridDim,Y];
                Yfill = Yfill + gridDim;
                j = j+1;
            end
            NY = j;
            X = repmat(X,j,1);
            Y = repmat(Y',1,i);
            
            Xshift = (rand(size(X))*2-1)*(gridDim-annulusDimPix)/2;
            Yshift = (rand(size(Y))*2-1)*(gridDim-annulusDimPix)/2;
            
            X = X + Xshift; %X = X*PixelsPermm;
            Y = Y + Yshift; %Y = Y*PixelsPermm;
            Xpos = X(:);   Ypos = Y(:);
            
            MaxgridDim      = qall(end)*annulusDimPix;
            MaxStimLength = Length*MaxgridDim + annulusDimPix/2;
            MaxStimWidth  = Width*MaxgridDim + annulusDimPix/2;
            MaxStimLengthmm = MaxStimLength/PixelsPermm;
            MaxStimWidthmm  = MaxStimWidth/PixelsPermm;
            StimLength = round((MaxStimLength-annulusDimPix/2)/gridDim);
            StimWidth  = round((MaxStimWidth -annulusDimPix/2)/gridDim);
            StimLengthmm = (StimLength*gridDim +annulusDimPix/2)/PixelsPermm;
            StimWidthmm  = (StimWidth*gridDim + annulusDimPix/2)/PixelsPermm;
            prop    = (2*Width+1)/(2*Length+1);
            newprop = (2*StimWidth+1)/(2*StimLength+1);
            
            radius = abs((Xpos-xCenter)-(Ypos-yCenter).*1i);
            % for Simona (sessions 1-8: 2d quadrant)
%             G = radius>=courserAreaRpixels-(10*PixelsPermm) & radius<=courserAreaRpixels+(10*PixelsPermm) & Ypos<yCenter-MaxStimLength & Xpos<xCenter-MaxStimLength;
            % for Simona (session 9: 3d quadrant)
%             G = radius>=courserAreaRpixels-(10*PixelsPermm) & radius<=courserAreaRpixels+(10*PixelsPermm) & Ypos>yCenter+MaxStimLength & Xpos<xCenter-MaxStimLength;
            % for the last session (2d quadrant)
%             G = radius>=courserAreaRpixels-(10*PixelsPermm) & radius<=courserAreaRpixels+(10*PixelsPermm) & Ypos<yCenter-MaxStimLength & Xpos<xCenter-MaxStimLength;
            % for the first 8 sessions (4th quadrant)
            G = radius>=courserAreaRpixels-(10*PixelsPermm) & radius<=courserAreaRpixels+(10*PixelsPermm) & Ypos>yCenter+MaxStimLength & Xpos>xCenter+MaxStimLength;
            B = find(G==1);
            rectCent = B(floor(rand*nnz(G))+1);
            yRectCent = floor(rectCent/size(X,1))+1;
            xRectCent = mod(rectCent,size(X,1));
            M = zeros(size(X));
            if rand>=0.5
                vertical(k,l) = 1;
            else
                vertical(k,l) = 0;
            end
            
            switch vertical(k,l)
                case 0
                    
                    if newprop < prop + 0.05 && newprop > prop - 0.05
                        M((xRectCent-StimWidth):(xRectCent+StimWidth),(yRectCent-StimLength):(yRectCent+StimLength)) = 1;
                    else
                        
                        while newprop >= prop + 0.05
                            ii = 1; jj = 1;
                            if StimLengthmm >= MaxStimLengthmm
                                M((xRectCent-StimWidth):(xRectCent+StimWidth-ii),(yRectCent-StimLength):(yRectCent+StimLength)) = 1;
                                num = 2*StimWidth+1 - ii; denum = 2*StimLength+1;
                                newprop = num/denum;
                                ii = ii + 1;
                            elseif StimLengthmm < MaxStimLengthmm
                                M((xRectCent-StimWidth):(xRectCent+StimWidth),(yRectCent-StimLength-jj):(yRectCent+StimLength)) = 1;
                                num = 2*StimWidth+1; denum = 2*StimLength+1 + jj;
                                newprop = num/denum;
                                jj = jj + 1;
                                StimLengthmm = StimLengthmm + CircleSizemm/2;
                            end
                        end
                        
                        while newprop <= prop - 0.05
                            ii = 1; jj = 1;
                            if StimLengthmm >= MaxStimLengthmm
                                M((xRectCent-StimWidth):(xRectCent+StimWidth),(yRectCent-StimLength):(yRectCent+StimLength-ii)) = 1;
                                num = 2*StimWidth+1; denum = 2*StimLength+1 - ii;
                                newprop = num/denum;
                                ii = ii + 1;
                                StimLengthmm = StimLengthmm - CircleSizemm/2;
                            elseif StimLengthmm < MaxStimLengthmm
                                M((xRectCent-StimWidth-jj):(xRectCent+StimWidth),(yRectCent-StimLength):(yRectCent+StimLength)) = 1;
                                num = 2*StimWidth+1 + jj; denum = 2*StimLength+1;
                                newprop = num/denum;
                                jj = jj + 1;
                            end
                        end
                        
                    end
                    
                otherwise
                    
                    if newprop <= prop + 0.05 && newprop >= prop - 0.05
                        M((xRectCent-StimLength):(xRectCent+StimLength),(yRectCent-StimWidth):(yRectCent+StimWidth)) = 1;
                    else
                        
                        while newprop >= prop + 0.05
                            ii = 1; jj = 1;
                            if StimLengthmm >= MaxStimLengthmm
                                M((xRectCent-StimLength):(xRectCent+StimLength),(yRectCent-StimWidth):(yRectCent+StimWidth-ii)) = 1;
                                num = 2*StimWidth+1 - ii; denum = 2*StimLength+1;
                                newprop = num/denum;
                                ii = ii + 1;
                            elseif StimLengthmm < MaxStimLengthmm
                                M((xRectCent-StimLength-jj):(xRectCent+StimLength),(yRectCent-StimWidth):(yRectCent+StimWidth)) = 1;
                                num = 2*StimWidth+1; denum = 2*StimLength+1 + jj;
                                newprop = num/denum;
                                jj = jj + 1;
                                StimLengthmm = StimLengthmm + CircleSizemm/2;
                            end
                        end
                        
                        while newprop <= prop - 0.05
                            ii = 1; jj = 1;
                            if StimLengthmm >= MaxStimLengthmm
                                M((xRectCent-StimLength):(xRectCent+StimLength-ii),(yRectCent-StimWidth):(yRectCent+StimWidth)) = 1;
                                num = 2*StimWidth+1; denum = 2*StimLength+1 - ii;
                                newprop = num/denum;
                                ii = ii + 1;
                                StimLengthmm = StimLengthmm - CircleSizemm/2;
                            elseif StimLengthmm < MaxStimLengthmm
                                M((xRectCent-StimLength):(xRectCent+StimLength),(yRectCent-StimWidth-jj):(yRectCent+StimWidth)) = 1;
                                num = 2*StimWidth+1 + jj; denum = 2*StimLength+1;
                                newprop = num/denum;
                                jj = jj + 1;
                            end
                        end
                        
                    end
                    
            end
            
            Ngratings = numel(X);
            contrast = rand(size(X));%*0.2 + (0.8);
            %                 contrast  = contrast(:);
            %             U = contrast <= .5;
            %             contrast(U) = rand(nnz(U),1)* (0.5-dev);
            %             contrast(~U)= rand(nnz(~U),1)*(0.5-dev)+(0.5+dev);
            contrast(M==1) = rand(nnz(M),1) * InContRange + .5 - InContRange/2;
            contrast = min(max(contrast,0),1);
            centerX  = size(contrast,1)-DeviationFromCenterX+1;
            centerY  = size(contrast,2)-DeviationFromCenterY+1;
            contrast = reshape(contrast,[size(X,1),size(X,2)]);
            contrast(size(X,1)-DeviationFromCenterY+1,size(X,2)-DeviationFromCenterX+1) = 0;
            contrast  = contrast(:);
            
            % Making a circle mask for each annulus
            [x, y] = meshgrid(-annulusDimPix:annulusDimPix, -annulusDimPix:annulusDimPix);
            H = sqrt(x.^2 + y.^2)<annulusDimPix;
            [s1, s2] = size(H);
            mask1 = ones(s1, s2, 2) * grey;
            mask1(:, :, 2) = white * (1 - H);
            
            % Make the destination rectangles for all the Annuluses in the array
            baseRect = [0 0 annulusDimPix annulusDimPix];
            allRects = nan(4, Ngratings);
            for b = 1:Ngratings
                allRects(:, b) = CenterRectOnPointd(baseRect, Xpos(b), Ypos(b));
            end
            
            HideCursor;
            
            % Define each annulus texture and merge it with the mask.
            for b = 1 : Ngratings
                
                annulus =  grey + inc * contrast(b) * cos(2 * pi * R * frequency + offset);
                
                % Make our radial gabor (Annulus) texure into a screen texture for drawing
                Texture = Screen('MakeTexture', mainwin, annulus);
                Screen('DrawTexture', mainwin, Texture, [], allRects(:,b));
                
                % Make our mask1 into a screen texture for drawing
                mask1tex = Screen('MakeTexture', mainwin, mask1);
                Screen('DrawTextures', mainwin, mask1tex, [], allRects(:,b));
                
            end
            
            
            % Flip to the screen
            Screen('Flip', mainwin);
            
            % Get images array to save
            imageArray=Screen('GetImage', mainwin);
            
            %Wait for a key press
            %             KbStrokeWait;
            WaitSecs(0.5);
            
            % Clear the screen
            sca
            
            % Write the image array on an image format
            
            fname1 = [dir3,'texture_trial',num2str(k),'_block',num2str(l),'.mat'];
            save(fname1,'imageArray','-v7.3');
%             fname2 = [dir3,'texture_trial',num2str(k),'_block',num2str(l),'.png'];
%             imwrite(imageArray,fname2);
            
            %             if (p==5 && pall(5)==1)
            %                vertical(k,l) = 5;
            %             end
            
            % Save contrast and Spacings for the computational work
            fname3 = [dir2,'Pos4Spacing_trial',num2str(k),'_block',num2str(l),'_p',num2str(p),'_q',num2str(q),'.mat'];
            save(fname3,'X','Y');
            fname4 = [dir1,'Contrast_trial',num2str(k),'_block',num2str(l),'_p',num2str(p),'_q',num2str(q),'.mat'];
            save(fname4,'contrast');
            fname5 = [dir4,'Exper_Variables_trial','_block',num2str(l),'_p',num2str(p),'_q',num2str(q),'.mat'];
            save(fname5,'Ngratings','PixelsPermm','NX','NY','M','vertical');
            fname6 = [dir4,'orientation.mat'];
            save(fname6,'vertical');
            
        end
    end
end