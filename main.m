function main
global g
% setup working folder and main window figure  interface. 
g.f = figure('OuterPosition',[1 1 1200 800]);
g.f.Name='NearField Scan by Computer Vision  KEYSIGHT COS HTC'
f=uimenu('Label','New');
uimenu(f,'Label','Near Field Scan by Computer Vision');
folderName1='./objs';
folderName2='./userdata';
folderName3='./lib';
g.ifgrid=1;

g.fname.data='./userdata/*.mat';
g.fname.watchlist='./userdata/*.mat';
g.fname.order='./userdata/*.mat';

%if ~isdeployed
addpath(folderName1,folderName2,folderName3);
%end

% Greate GUI interface

b=uix.HBoxFlex( 'Parent', g.f);
g.p1 = uix.VBox( 'Parent', b, 'Padding', 5 );
g.p2 = uix.TabPanel( 'Parent', b, 'Padding', 5 );
set( b, 'Widths', [-1 -8], 'Spacing', 5 );
% setup button on left panel
g.b1=uicontrol( 'Parent', g.p1, 'String', 'Configuration','Callback', @fb1 );
g.b3=uicontrol( 'Parent', g.p1, 'String', 'Measurement','Callback', @fb3 );        
g.b6=uicontrol( 'Parent', g.p1, 'String', 'Data Visualization','Callback', @fb6 );
g.tx11=uicontrol( 'Style','text','Parent', g.p1, 'String', '' );
g.tx12=uicontrol( 'Style','text','Parent', g.p1, 'String', '' );
set( g.p1, 'Heights', 30*ones(1,5) );
% setup UITab on right panel. 
g.p2c.t1 = uix.Panel( 'Parent', g.p2, 'Title', 'Configuration', 'Padding', 5 );
g.p2c.t2 = uix.Panel( 'Parent', g.p2, 'Title', 'Data Collection', 'Padding', 5 );
g.p2c.t3 = uix.Panel( 'Parent', g.p2, 'Title', 'Data Visualization', 'Padding', 5 );
g.p2.TabTitles(1)={'Config'};
g.p2.TabTitles(2)={'Collection'};
g.p2.TabTitles(3)={'Visualization'};
%% setup tab 1  Connection to Camera and SA
g.p2c.t1c=uix.VBoxFlex( 'Parent', g.p2c.t1);
uicontrol( 'Style','text','Parent', g.p2c.t1c, 'String', 'IP Camera IP address' );
g.b1c.rx1=uicontrol('Style','edit','Parent', g.p2c.t1c, 'String', 'http://10.112.89.123:6060/videostream.cgi?user=admin&pwd=P0o9i8u7&resolution=32&rate=3');
uicontrol( 'Style','text','Parent', g.p2c.t1c, 'String', 'Spectrum Analyzer address' );
g.b1c.rx2=uicontrol('Style','edit','Parent', g.p2c.t1c, 'String', 'TCPIP0::localhost::inst0::INSTR');
g.b1=uicontrol( 'Parent', g.p2c.t1c, 'String', 'Connect','Callback', @fconnect );
set( g.p2c.t1c, 'Heights', 30*ones(1,5) );
%% setup tab 2  Data Collection Interface Setup
g.p2c.t2c=uix.VBoxFlex( 'Parent', g.p2c.t2);
g.p2c.t2d.c=uicontainer('Parent',g.p2c.t2c);
g.ax0=axes('Parent',g.p2c.t2d.c);
rgbFrame=imshow('testsetup.jpg');
uistack(g.ax0,'bottom');
set(g.ax0,'handlevisibility','off','visible','off')
g.ax1=axes('Parent',g.p2c.t2d.c);     
g.p2c.t2d.btc = uix.HButtonBox( 'Parent', g.p2c.t2c, 'Padding', 5 );
set( g.p2c.t2c, 'Heights', [-1 30] );
g.p2c.t2d.bt1=uicontrol( 'Parent',g.p2c.t2d.btc, 'String', 'preview' ,'UserData',[1 2 1],'Callback',@bt_measure);
g.p2c.t2d.bt2=uicontrol( 'Parent',g.p2c.t2d.btc, 'String', 'setup mesh' ,'UserData',[1 2 2],'Callback',@bt_measure);
g.p2c.t2d.bt3=uicontrol( 'Parent',g.p2c.t2d.btc, 'String', 'Start Scan' ,'UserData',[1 2 3],'Callback',@bt_measure);
g.p2c.t2d.bt4=uicontrol( 'Parent',g.p2c.t2d.btc, 'String', 'Stop scan' ,'UserData',[1 2 4],'Callback',@bt_measure);
g.p2c.t2d.bt5=uicontrol( 'Parent',g.p2c.t2d.btc, 'String', 'Resume' ,'UserData',[1 2 5],'Callback',@bt_measure);
g.p2c.t2d.bt6=uicontrol( 'Parent',g.p2c.t2d.btc, 'String', 'Finish' ,'UserData',[1 2 6],'Callback',@bt_measure);

%% setup tab 3 Post Measurement Data Visualization
g.p2c.t3c=uix.VBoxFlex( 'Parent', g.p2c.t3);
g.p2c.t3d.c=uicontainer('Parent',g.p2c.t3c);
g.p2c.t3d.c1=uicontainer('Parent',g.p2c.t3c);
g.ax20=axes('Parent',g.p2c.t3d.c);
rgbFrame=imshow('testsetup.jpg');
uistack(g.ax20,'bottom');
set(g.ax20,'handlevisibility','off','visible','off')
g.ax3=axes('Parent',g.p2c.t3d.c);
g.ax4=axes('Parent',g.p2c.t3d.c1); 
g.p2c.t3d.btc = uix.HButtonBox( 'Parent', g.p2c.t3c, 'Padding', 5 );
set( g.p2c.t3c, 'Heights', [-5 -2 30] );
g.p2c.t3d.bt1=uicontrol( 'Parent',g.p2c.t3d.btc, 'String', 'load data' ,'UserData',[1 2 1],'Callback',@bt_analysis);
g.p2c.t3d.bt2=uicontrol( 'Parent',g.p2c.t3d.btc, 'String', 'maximum hold' ,'UserData',[1 2 2],'Callback',@bt_analysis);
g.p2c.t3d.bt3=uicontrol( 'Parent',g.p2c.t3d.btc, 'String', 'spectrum' ,'UserData',[1 2 3],'Callback',@bt_analysis);
g.p2c.t3d.bt4=uicontrol( 'Parent',g.p2c.t3d.btc, 'String', 'transparency level' ,'UserData',[1 2 4],'Callback',@bt_analysis);
g.p2c.t3d.bt5=uicontrol( 'Parent',g.p2c.t3d.btc, 'String', 'Grid' ,'UserData',[1 2 5],'Callback',@bt_analysis);
g.p2c.t3d.bt6=uicontrol( 'Parent',g.p2c.t3d.btc, 'String', 'Reserved' ,'UserData',[1 2 6],'Callback',@bt_analysis);
             
end

function fb1(~,h) % Click the connection button, bring the UI tab 1 to active
global g

g.p2.Selection=1;

end
function fb3(~,h)% Click the Data collection button, bring the UI tab 2 to active
global g

g.p2.Selection=2;

end

function fb6(h,~)% Click the Visualization button, bring the UI tab 3 to active
global g

g.p2.Selection=3;

end

function fconnect(~,h) 

global g
global p
global d
%%connect to camera % Make connection to Camera Need to load the IP Camera Driver from matlab support package installer. I used the FosCam IP Camera

p.cam=ipcam(g.b1c.rx1.String);
p.cam.Timeout=60;
%% connect to SA
% remove all previous connection on USB
a=instrfindall;
for i=1:length(a)
    fclose(a(i))
end
% saObj is the OOB of SA stored under ./objs folder.
% setup connection to SA and initilze the SA. You need to modify the code
% to work with your SA model

p.sa=saObj;
p.sa.setaddress(g.b1c.rx2.String);
p.sa.ini;
p.sa.configRE;
p.sa.setsweep2(1)
pause(1)
p.sa.takedata(1);
p.sa.takefreq;
d.freq=p.sa.freq;
msgbox('connected!')
g.p2.Selection=2;

end

function bt_measure(h,~)
% when Click the button on UITab 2 . 
global p
global g
global d
switch h.UserData(3)
   
    case 1 % when click the preview button
         h=preview(p.cam);
         
    case 2 % when click the meash button. 
         % get customer input for the mesh size
         prompt={'Enter the mesh size Veritcal direction','Enter the mesh size Horizontal direction'};
         name='Input for Mesh function';
         numlines=1;
         defaultanswer={'11','22'};
         answer=inputdlg(prompt,name,numlines,defaultanswer);
        
         mm=str2num(answer{1});nn=str2num(answer{2});
         d.mm=mm;d.nn=nn;
        % snapshot the DUT 
        g.ax0.Children.CData=snapshot(p.cam);
        axes(g.ax1)
        % overlay the meshgrid onto the DUT
        [g.h g.txt]=overlaygrid_setup(d.mm,d.nn);
   
        
    case 3 % click measurement. start the measurement
         %% Tracking
          % initilize the data holder. 
           d.data(1:d.mm*d.nn)={zeros(length(d.freq),1)};
          % start measurement. 
          tracking;
    case 4  % toggle the pause flag
        g.ifrun=0;
    case 5  % resume the measurement. 
        tracking;
    case 6  % at the end of data collection save the data
        g.ifrun=0;
        % get the info where to save the data
        [filename, pathname] = uiputfile('./userdata/*.mat', 'Pick an file name to save the scanning data');
        if isequal(filename,0) || isequal(pathname,0)
           disp('User pressed cancel')
        else
           disp(['User selected ', fullfile(pathname, filename)])
        end
        fname=fullfile(pathname, filename);
        % save the measured matrix into fname
        save(fname, 'd');
        % save the snapshot of DUT to the same fname but end with .jpg.
        rgbFrame=snapshot(p.cam);
        imwrite(rgbFrame,[fname(1:end-4) '.jpg']);
end

end

function tracking 
% functin of tracking the probe location and save trace data. this is the
% main funcion. 
global g
global d
global p

                redThresh = 0.2; % Threshold for red detection
                % vidDevice = imaq.VideoDevice('winvideo', 1, 'YUY2_640x480', ... % Acquire input video stream
                %                     'ROI', [1 1 640 480], ...
                %                     'ReturnedColorSpace', 'rgb');
                % 
                % %p.ip=[str ':6060'];
                % % resuliton 32 full size. 
                % % resultion 16 half size 8: 320x240
                % p.cam=ipcam(['http://141.121.239.62:6060/videostream.cgi?user=admin&pwd=P0o9i8u7&resolution=32&rate=3']);
                % p.cam.Timeout=60;
                % msgbox('connected!')
                % vidDevice=p.cam
                %vidInfo = imaqhwinfo(vidDevice); % Acquire input video property
                gridsize=640*480/d.mm/d.nn;
                hblob = vision.BlobAnalysis('AreaOutputPort', false, ... % Set blob analysis handling
                                                'CentroidOutputPort', true, ... 
                                                'BoundingBoxOutputPort', true', ...
                                                'MinimumBlobArea', floor(gridsize/2), ...
                                                'MaximumBlobArea', floor(gridsize*2), ...
                                                'MaximumCount', 1);
                hshapeinsRedBox = vision.ShapeInserter('BorderColor', 'Custom', ... % Set Red box handling
                                                        'CustomBorderColor', [1 0 0], ...
                                                        'Fill', true, ...
                                                        'FillColor', 'Custom', ...
                                                        'CustomFillColor', [1 0 0], ...
                                                        'Opacity', 0.6);
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
                g.ifrun=1;

                %% Processing Loop stop after 400 frame or toggle the pause button
                
               % while(nFrame < 40000 && g.ifrun)
                    while(g.ifrun)
                    rgbFrame = snapshot(p.cam); % Acquire single frame
                    %rgbFrame = flipdim(rgbFrame,2); % obtain the mirror image for displaying
                    diffFrame = imsubtract(rgbFrame(:,:,1), rgb2gray(rgbFrame)); % Get red component of the image
                    diffFrame = medfilt2(diffFrame, [3 3]); % Filter out the noise by using median filter
                    binFrame = im2bw(diffFrame, redThresh); % Convert the image into binary image with the red objects as white
                   [centroid, bbox] = step(hblob, binFrame); % Get the centroids and bounding boxes of the blobs
                    centroid = uint16(centroid); % Convert the centroids into Integer for further steps 
                    vidIn = step(hshapeinsRedBox, rgbFrame, bbox); % Instert the red box
                    for object = 1:1:length(bbox(:,1)) % Write the corresponding centroids
                        centX = centroid(object,1); centY = centroid(object,2);
                        centX1 = ceil((centX)/d.dx); centY1 = ceil((centY)/d.dy);
                        vidIn = step(htextinsCent, vidIn, [centX centY], [centX-6 centY-9]); 
                        title(gca, ['(X,Y) = (', num2str(centX1), ', ',num2str(centY1), ')']);
                        if centX1>0 && centX1<=d.nn && centY1>0 && centY1<=d.mm
                            ind=sub2ind(size(d.M),centY1,centX1);
                            % update the mask
                            d.M(ind)=1;
                            mask=d.M*0.5;
                            g.h.AlphaData=mask;
                            %take the trace out of SA
                            tmp=[d.data{ind} p.sa.takedata(1)];
                            d.data(ind)={max(tmp,[],2)};
                            
                            % update the reading
                            d.MData(ind)=max(d.data{ind});
                            g.txt(ind).String=num2str(d.MData(ind)) ;
                            %update the plot
                            g.h.CData(ind)=d.MData(ind);

                        end


                    end
                    vidIn = step(htextins, vidIn, uint8(length(bbox(:,1)))); % Count the number of blobs
                 %   step(hVideoIn, vidIn); % Output video stream
                    g.ax0.Children.CData=vidIn ; % plot the updated frame to the axis on GUI
                  nFrame = nFrame+1; % update the frame number.
                end

                
end

function [h txt]=overlaygrid_setup(m,n)
global g
global d
d.M = zeros(m,n);
d.MData = zeros(m,n);
[r c] = size(d.M);
X=640; Y=480;
d.dx=X/(c);
d.dy=Y/(r);
hold(g.ax1,'on');
%# text location and labels
[xloc yloc] = meshgrid((1:c)*d.dx,(1:r)*d.dy);
g.xloc = xloc(:); g.yloc = yloc(:);
str = strtrim(cellstr( num2str(d.MData(:),'%.3g') ));
xticklabels = cellstr( num2str((1:c)','X%d') );
yticklabels = cellstr( num2str((1:r)','Y%d') );

%# plot colored cells
h = imagesc((1:c)*d.dx, (1:r)*d.dy, d.MData);
%mask = M>0.9;               %# or any other mask
mask = d.M*0.5;               %# or any other mask

h.AlphaData=mask;

colormap(summer)            %# colormap([0 1 0])
set(gca, 'Box','on', 'XAxisLocation','top', 'YDir','reverse', ...
    'XLim',([0 c]+0.5)*d.dx, 'YLim',([0 r]+0.5)*d.dy, 'TickLength',[0 0], ...
    'XTick',(1:c)*d.dx, 'YTick',(1:r)*d.dy, ...
    'XTickLabel',xticklabels, 'YTickLabel',yticklabels, ...
    'LineWidth',2, 'Color','none', ...
    'FontWeight','bold', 'FontSize',8, 'DataAspectRatio',[1 1 1]);

%# plot grid
xv1 = repmat(((2:c)-0.5)*d.dx, ([2 1])); xv1(end+1,:) = NaN;
xv2 = repmat([0.5*d.dx;(c+0.5)*d.dx;NaN], [1 r-1]);
yv1 = repmat([0.5*d.dy;(r+0.5)*d.dy;NaN], [1 c-1]);
yv2 = repmat(((2:r)-0.5)*d.dy, [2 1]); yv2(end+1,:) = NaN;
line([xv1(:);xv2(:)], [yv1(:);yv2(:)], 'Color','k', 'HandleVisibility','off')

%# plot text
str = strtrim(cellstr( num2str(d.MData(:),'%.3g') ));
txt=text(g.xloc, g.yloc, str, 'FontSize',8, 'HorizontalAlignment','center');

end

function bt_analysis(h,~)
global g
global d
switch h.UserData(3)
      case 1
        
         [filename, pathname] = uigetfile('./userdata/*.mat', 'Pick a MATLAB code file');
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
    else
       disp(['User selected ', fullfile(pathname, filename)])
    end
    fname=fullfile(pathname, filename);
     load(fname);
     rgbframe1=imread([fname(1:end-4) '.jpg']);
     tmp=size(rgbframe1);
       g.sx=640;g.sy=floor(640/tmp(2)*tmp(1));
     rgbframe=imresize(rgbframe1,[ g.sy g.sx ]);
     g.ax20.XLim=[0 640];
     g.ax20.YLim=[0 640/tmp(2)*tmp(1)];
     g.ax20.Children.CData=rgbframe;
     axes(g.ax3);
     

     [g.h1 g.txt1]=overlaygrid_setup(d.mm,d.nn);
     g.h1.AlphaData=ones(d.mm,d.nn)*0.5;
      d.A=cell2mat(d.data);
        maxhold=max(d.A,[],2);
        dcm_obj = datacursormode(gcf);
        axes(g.ax4);
        plot(g.ax4,d.freq/1e6,maxhold);
        grid on
        xlim=[d.freq(1) d.freq(end)];
        xlabel('Frequency(MHz)');ylabel('Amplitude(dBuV)')
        set(dcm_obj,'UpdateFcn',@NewCallback);
        g.title1=title('Maximum hold value')
        
  
    case 2
        set (gcf, 'WindowButtonMotionFcn','');
        d.A=cell2mat(d.data);
        maxhold=max(d.A,[],2);
        dcm_obj = datacursormode(gcf);
        axes(g.ax4);
        plot(g.ax4,d.freq/1e6,maxhold);
        grid on
        xlim=[d.freq(1) d.freq(end)];
        xlabel('Frequency(MHz)');ylabel('Amplitude(dBuV)')
        set(dcm_obj,'UpdateFcn',@NewCallback);
        g.title1=title('Maximum hold value')
    case 3
        axes(g.ax3);
        set (gcf, 'WindowButtonMotionFcn', @mouseMove);
    
    case 4
         prompt={'Enter the transparent ratio for field overlay(0 totally transparent, 1 totally blocked)'};
   name='Input for Peaks function';
   numlines=1;
   defaultanswer={'0.5'};
 
   answer=inputdlg(prompt,name,numlines,defaultanswer);
      g.h1.AlphaData=ones(size(g.h1.CData))*str2num(answer{1});
    case 5
   togglegrid;
end

end

function output_txt = NewCallback(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
global g
global d

pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

I=find(d.freq>=pos(1)*1e6);
ind=I(1);
B=d.A(ind,:);
C=reshape(B,[d.mm d.nn]);
g.h1.CData=C;

 g.h1.XData=[0:g.sx/(d.nn-1):g.sx];
 g.h1.YData=[0:g.sy/(d.mm-1):g.sy];
smoothheatmap;
          
          
g.title1.String=['Selected freq=' num2str(floor(d.freq(ind)/1e6)) 'MHz'];
%update the reading.
 
if    g.ifgrid
 for i=1:d.mm
     for k=1:d.nn
         ind=sub2ind([d.mm d.nn],i,k);
         g.txt1(ind).String=num2str(C(i,k),4);
        
     end
 end
 end
 
end

function smoothheatmap
global g
global d

        data=g.h1.CData;
        dx=g.h1.XData(2);dy=g.h1.YData(2);
         sx= g.h1.XData(end);sy= g.h1.YData(end);
        [X,Y] = meshgrid([1:d.nn]*dx,[1:d.mm]*dy);
        [X2,Y2]=meshgrid([1:0.2:d.nn]*dx,[1:0.2:d.mm]*dy);
         outData = interp2(X, Y, data, X2, Y2, 'linear');
         tmp=size(outData);
         dx=sx/(tmp(2)-1);dy=sy/(tmp(1)-1);
         outDataX=[0:dx:sx]; outDataX(end)=sx+2.5*dx;
         outDataY=[0:dy:sy];outDataY(end)=sy+2.5*dy;
        
         g.h1.XData=outDataX;
         g.h1.YData=outDataY;
         g.h1.CData=outData;
         g.h1.AlphaData=ones(size(g.h1.CData))*g.h1.AlphaData(1,1);
end
function togglegrid
global g
 if g.ifgrid
            g.ifgrid=0;
            for i=1:length( g.txt1)
            g.txt1(i).Text.Visible='off';
            g.line1.Visible='off';
            end
            
        else
            g.ifgrid=1;
             for i=1:length( g.txt1)
            g.txt1(i).Text.Visible='on';
            g.line1.Visible='on';
            end
 end
        
        
end


function mouseMove (object, eventdata)
global g
global d
axes(g.ax3);
C = get (g.ax3, 'CurrentPoint');
b1=ceil(C(1,1)/d.dx);
b2=ceil(C(1,2)/d.dy);
[mm n]=size(d.M);
g.title1.String=['(X,Y) = (', num2str(b1), ', ',num2str(b2), ')'];
    if b1>0 && b1<=d.nn && b2>0 && b2<=d.mm
        
        g.ax4.Children(end).YData=d.data{sub2ind(size(d.M),b2,b1)}';
       
    end
end
