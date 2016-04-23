% Cluster Fix algorithm for Ephiz Recordings
% Seth König, January, 2015
% Runs the same as previous versions of Cluster Fix, but sets eye data
% points outside of desired boundaries (i.e. outside of the image) to NaNs.
% Code here is  generically organized since everyone may have slightly
% different ways of organizing their data and/or different tasks require
% different organization. Code assumes eye data was collected at 1000 Hz and
% is already aligned to spike times.
%
% Cluster Fix requires "findgaps.m" to run! It should be in the same folder
% you found this script. 
%
% Cluster Fix should run automatically. Cluster Fix is not that fast and it
% will take approximatley as long as it takes to collect the data. If you 
% experience errors or think parameters should be adjusted please look at 
% the comments first. If you still experience errors please contact me. 
%
% sample calibrated eye data. You should start with calibrated eye data. The
% data can be in pixels or in dva but make sure the desired boundaries are 
% correctly used in the for-loop below. Eye data is in the cell array eyedat
% with each cell corresponding to a different trial; in the 1st row of
% each cell is the horizontal (x) eye trace and in the 2nd row is the
% vertical (y) eye trace. Here, eye data is orginally in cortex dva.  
load('SampleEyeData')

%%%---Pre-Process eye data data---%%%
% define the size of the ROI/image. These will be used as the boundaries in
% the next section. Make sure they match up appropriately. For example,
% this code is using pixel coordinates where [1,1] represents the lower
% left corner of the screen and [800,600] represents the upper right corner
% of the screen. 
imageX = 800; %horizontal size of the image 
imageY = 600; % vertical size of the image

% For-loop replaces eye data points outside of the desired boundaries/ROI/Image
% with NaNs. By default I use a buffer of 24 pixels (~ 1 dva) to 48 pixels (~ 2 dva).
% 0 pixels may be too restrictive if calibration is an issue and/or the monkey
% is fixating an item on the border of the image.
for t = 1:length(eyedat);%iterate by trial
    x = eyedat{t}(1,:); %horizontal eye data;
    y = eyedat{t}(2,:); %vertical eye data;
    
    %convert eye data from dva coordinates to pixel coordinates
    x = 24*x; % ~24 pixels/dva 
    y = 24*y; % ~24 pixels/dva 
    x = x+imageX/2; % re-center
    y = y+imageY/2; % re-center
    
    %replace eye data too far outside ROI/image/Boundary with NaNs
    %too far left
    y(x < -24) = NaN; 
    x(x < -24) = NaN;
    %too far right
    y(x > imageX+24) = NaN;
    x(x > imageX+24)= NaN;
    %too far below 
    x(y < -24) = NaN;
    y(y < -24) = NaN;
    %too far above
    x(y > imageY+24) = NaN;
    y(y > imageY+24) = NaN;
    
    %store readjusted data
    eyedat{t}(1,:) = x;
    eyedat{t}(2,:) = y;
end
clear x y

%%%---Run Cluster Fix---%%%
% use Cluster Fix to Detect fixations and saccades in XY eye data
% since timing is very important NaNs must remain. Trials with NaNs or with eye
% data lasting less than 200 ms (approximately less than 1 saccade and fixation)
% are ignored, and then trials are processed in chunks seperated by NaNs. Chunks
% are then finally recombined together into their original format. For more
% details refer the the description in "ClusterFix_Plexon.m". 
fixationstats = cell(1,length(eyedat)); %structure array containing important data
% fixationstats.fixationtimes: indexes when fixation is occuring, row 1 is
% start of fixation/row 2 end of fixation, coloumn by ordinal fixation #

% fixationstats.fixations: mean position of fixation, row 1 is mean
% horizontal position/row 2 mean vertical position, columb by ordinal fixation #

% fixationstats.saccadetimes: indexes of saccade periods, row 1 when
% saccade starts/row 2 when saccade ends, column by ordinal saccade #

% fixationstats.XY = [x;y]; virtually same as eye data. May have more NaNs
% if segments of eye data within boundaries are shorter than 200 ms. 
for t = 1:length(eyedat);%runs 1 trail at a time. 
    if ~isempty(eyedat{t})
        % Attempt not to break up a trial into more segments than necessary 
        % becauuse it may be detrimental to clustering. 
        fixationstats{t} = ClusterFix_Plexon(eyedat{t});
    end
end
clear t 
%%
%%%---Visually Verify that Cluster Fix is working "properly"---%%%
% fyi plots raw eyetrace
% Artifact of indexing: these are segments either next to the NaNs or
% inbetween a fixation and a saccade that are an artifact of the way the
% code plots. If there is more than 1 data point in a row that is black then
% there is a concern.  
for t = 1:10%length(fixationstats) %change if you want more plots.
    fixations = fixationstats{t}.fixations; %mean fixation location
    fixationtimes = fixationstats{t}.fixationtimes; %indexes of fixation times
    saccadetimes = fixationstats{t}.saccadetimes; %indexes of saccadetimes
    x = fixationstats{t}.XY(1,:); %horizontal eye trace
    y = fixationstats{t}.XY(2,:); % vertical eye trace
    
    figure
    
    %2D plot with x vs time
    subplot(2,2,1)
    hold on
    p(1)= plot(x,'k');
    for f = 1:size(fixationtimes,2)
        p(2) = plot(fixationtimes(1,f):fixationtimes(2,f),...
            x(fixationtimes(1,f):fixationtimes(2,f)),'r');
    end
    for s = 1:size(saccadetimes,2)
        p(3) = plot(saccadetimes(1,s):saccadetimes(2,s),...
            x(saccadetimes(1,s):saccadetimes(2,s)),'g');
    end
    hold off
    box off
    xlabel('Time (ms)')
    ylabel('Horizontal eye position (pixels)')
    legend(p,'Artifact of indexing','fixations','saccades')
    
    %2D plot with x vs time
    subplot(2,2,3)
    hold on
    p(1) = plot(y,'k');
    for f = 1:size(fixationtimes,2)
        p(2) = plot(fixationtimes(1,f):fixationtimes(2,f),...
            y(fixationtimes(1,f):fixationtimes(2,f)),'r');
    end
    for s = 1:size(saccadetimes,2)
        p(3) = plot(saccadetimes(1,s):saccadetimes(2,s),...
            y(saccadetimes(1,s):saccadetimes(2,s)),'g');
    end
    hold off
    box off
    xlabel('Time (ms)')
    ylabel('Horizontal eye position (pixels)')
    legend(p,'Artifact of indexing','fixations','saccades')
    
    %2D plot of x vs y
    subplot(2,2,[2 4])
    p(1) = plot(x,y,'k');
    hold on
    for f = 1:size(fixationtimes,2)
        p(2) = plot(x(fixationtimes(1,f):fixationtimes(2,f)),...
            y(fixationtimes(1,f):fixationtimes(2,f)),'r');
    end
    for s = 1:size(saccadetimes,2)
        p(3) = plot(x(saccadetimes(1,s):saccadetimes(2,s)),...
            y(saccadetimes(1,s):saccadetimes(2,s)),'g');
    end
    for f = 1:size(fixations,2)
        p(4) = plot(fixations(1,f),fixations(2,f),'*k');
    end
    hold off
    axis equal
    xlabel('Horizontal eye position (pixels)')
    ylabel('Vertical eye position (pixels)')
    legend(p,'Artifact of indexing','fixations','saccades',...
        'mean fixatoin location')
end
