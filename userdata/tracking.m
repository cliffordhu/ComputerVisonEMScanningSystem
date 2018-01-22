%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Name : Red Object Detection and Tracking                        %
% Author       : Arindam Bose                                             %
% Version      : 1.05                                                     %
% Description  : How to detect and track red objects in Live Video        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
redThresh = 0.1; % Threshold for red detection
% vidDevice = imaq.VideoDevice('winvideo', 1, 'YUY2_640x480', ... % Acquire input video stream
%                     'ROI', [1 1 640 480], ...
%                     'ReturnedColorSpace', 'rgb');

%p.ip=[str ':6060'];
% resuliton 32 full size. 
% resultion 16 half size 8: 320x240
p.cam=ipcam(['http://141.121.239.62:6060/videostream.cgi?user=admin&pwd=P0o9i8u7&resolution=32&rate=3']);
p.cam.Timeout=60;
msgbox('connected!')
vidDevice=p.cam
%vidInfo = imaqhwinfo(vidDevice); % Acquire input video property
hblob = vision.BlobAnalysis('AreaOutputPort', false, ... % Set blob analysis handling
                                'CentroidOutputPort', true, ... 
                                'BoundingBoxOutputPort', true', ...
                                'MinimumBlobArea', 400, ...
                                'MaximumBlobArea', 16000, ...
                                'MaximumCount', 1);
hshapeinsRedBox = vision.ShapeInserter('BorderColor', 'Custom', ... % Set Red box handling
                                        'CustomBorderColor', [1 0 0], ...
                                        'Fill', true, ...
                                        'FillColor', 'Custom', ...
                                        'CustomFillColor', [1 0 0], ...
                                        'Opacity', 0.4);
htextins = vision.TextInserter('Text', 'Number of Red Object: %2d', ... % Set text for number of blobs
                                    'Location',  [7 2], ...
                                    'Color', [1 0 0], ... // red color
                                    'FontSize', 12);
htextinsCent = vision.TextInserter('Text', '+      X:%4d, Y:%4d', ... % set text for centroid
                                    'LocationSource', 'Input port', ...
                                    'Color', [1 1 0], ... // yellow color
                                    'FontSize', 14);
% hVideoIn = vision.VideoPlayer('Name', 'Final Video', ... % Output video player
%                                 'Position', [100 100 vidInfo.MaxWidth+20 vidInfo.MaxHeight+30]);
nFrame = 0; % Frame number initialization
 rgbFrame = snapshot(vidDevice); % Acquire single frame
    h1= imshow(rgbFrame);
   
%% Processing Loop
while(nFrame < 2000)
    rgbFrame = snapshot(vidDevice); % Acquire single frame
    rgbFrame = flipdim(rgbFrame,2); % obtain the mirror image for displaying
    diffFrame = imsubtract(rgbFrame(:,:,1), rgb2gray(rgbFrame)); % Get red component of the image
    diffFrame = medfilt2(diffFrame, [3 3]); % Filter out the noise by using median filter
    binFrame = im2bw(diffFrame, redThresh); % Convert the image into binary image with the red objects as white
    [centroid, bbox] = step(hblob, binFrame); % Get the centroids and bounding boxes of the blobs
    centroid = uint16(centroid); % Convert the centroids into Integer for further steps 
    rgbFrame(1:20,1:165,:) = 0; % put a black region on the output stream
    vidIn = step(hshapeinsRedBox, rgbFrame, bbox); % Instert the red box
    for object = 1:1:length(bbox(:,1)) % Write the corresponding centroids
        centX = centroid(object,1); centY = centroid(object,2);
        vidIn = step(htextinsCent, vidIn, [centX centY], [centX-6 centY-9]); 
    end
    vidIn = step(htextins, vidIn, uint8(length(bbox(:,1)))); % Count the number of blobs
 %   step(hVideoIn, vidIn); % Output video stream

   h1.CData=vidIn ;
  nFrame = nFrame+1;
end
%% Clearing Memory
% release(hVideoIn); % Release all memory and buffer used
% release(vidDevice);
% % clear all;
% clc;