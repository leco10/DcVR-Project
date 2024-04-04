function varargout = Tail_Tracker_GUI(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Tail_Tracker_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @Tail_Tracker_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function Tail_Tracker_GUI_OpeningFcn(hObject, ~, handles, varargin)
% pathnames
clc; % clear the command window
disp([datestr(datetime('now')) ': started Tail Tracker GUI']);
handles.pathname = 'C:\Users\jlab\Desktop\Leo\TailTracker\';
cd(handles.pathname);

% create camera streaming object
imaqreset;
handles.Vid = videoinput('gentl','1','Mono8');
handles.VidSource = getselectedsource(handles.Vid);
handles.VidSource.BinningHorizontal = 2;
handles.VidSource.BinningHorizontalMode = 'Sum';
handles.VidSource.BinningVertical = 2;
handles.VidSource.BinningVerticalMode = 'Sum';
handles.VidSource.GainAuto = 'Off';
handles.VidSource.Gain = 1;
vid_rez=handles.Vid.VideoResolution;

% load or create config and state file
handles.GenericAcquisitionROI = [0 0 vid_rez];
if exist([handles.pathname 'Tail_Tracker_config.mat'],'file') == 2
    load([handles.pathname 'Tail_Tracker_config.mat'],'config');
else
    % default config
    config.Acquisition.CameraFPS = 20;
    config.Acquisition.CameraExposureTime = 1;
    config.Acquisition.AcquisitionROI = handles.GenericAcquisitionROI;
    config.Tracking.Point1 = [10,10] ;
    config.Tracking.Point2 = [100,100];
    config.Tracking.nSegments = 4;
    config.Tracking.Blurring = 5;
    config.Tracking.TrackingAngle = deg2rad(60);
    config.Estimator.vigor_window = 0.25; % [s]
    config.Estimator.reset_window = 0.05; % [s]
    config.Estimator.vigor_thresh = 0.25; % [degrees], forward bout threshold
    config.Estimator.forward_coef = 2; % [cm/(s*degrees)], transformation from vigor into forward velocity
    config.Estimator.turning_coef = 18.5; % [1/s], transformation from smoothed tail into angular velocity
    config.Estimator.v_ang_kernel_tau = 2; % time constant of inertia kernel
    config.Estimator.inertia_kernel_tau = 1.5; % time constant of inertia kernel
    config.Estimator.inertia_kernel_decay = 0.01; % decay to 1/100 after tau
    config.Estimator.tail_thresh = 0.1; % angular velocity threshold
    config.Estimator.tail_smoothing = 0.5; % [a.u.], low-pass-filtering of the tail to reduce wiggling
end
handles.config = config; clear config;
% create line and add listener function
handles.hImage = image(handles.cam_image, zeros(flip(vid_rez),'uint8'));
hold(handles.cam_image,'on');
handles.line = images.roi.Line(handles.cam_image,'Position',[handles.config.Tracking.Point1(1),handles.config.Tracking.Point1(2);handles.config.Tracking.Point2(1),handles.config.Tracking.Point2(2)]);
addlistener(handles.line,'ROIMoved',@line_moved);

% starting GUI state and variables
handles.run = false;
handles.use_video = false;
% current data
handles.Data.CurrentData.t = datetime('now');
handles.Data.CurrentData.tail = 0;
handles.Data.CurrentData.vigor = 0;
handles.Data.CurrentData.velocity = 0;
handles.Data.CurrentData.angular_velocity = 0;
handles.Data.CurrentData.time_since_bout_onset = 0;
handles.Data.CurrentData.tail_derivative = 0;
handles.Data.CurrentData.swim_forw = false; 
handles.Data.CurrentData.swim_ang = false; 
handles.Data.CurrentData.tail_LPF = 0;

% buffers
handles.Data.BufferData.t = [];
handles.Data.BufferData.tail = [];
handles.Data.BufferData.vigor = [];
handles.Data.BufferData.velocity = [];
handles.Data.BufferData.angular_velocity = [];
handles.Data.BufferData.tail_derivative = [];

% decay kernels
handles.Data.ang_vel_dampening_y = exp((0.01:0.01:1)*-log(1/handles.config.Estimator.inertia_kernel_decay)/(handles.config.Estimator.v_ang_kernel_tau)); % 0.55: mean bout duration, 100: ang_vel dampening 0.01 at end of mean bout
handles.Data.inertia_kernel = exp((0.01:0.01:1)*-log(1/handles.config.Estimator.inertia_kernel_decay)/(handles.config.Estimator.inertia_kernel_tau));
handles.Data.v_forward_inertia_list = [];

if isequal(handles.GenericAcquisitionROI, handles.config.Acquisition.AcquisitionROI)
    handles.SetAcquisitionROIButton.String = 'Draw acquisition ROI';
else
    handles.SetAcquisitionROIButton.String = 'Reset acquisition ROI';
end
handles.pi2 = pi/2;
handles.two_pi = 2*pi;

% create graphical objects
disableDefaultInteractivity(handles.cam_image);
handles.semicircles_obj = plot(handles.cam_image,0,0,'color','w','linewidth',1);
handles.tail_point_obj = plot(handles.cam_image,0,0,'marker','.','color','g','LineStyle','-','MarkerSize',20);
handles.tail_axis.XLim = [-10 0];
axes(handles.tail_axis)
disableDefaultInteractivity(handles.tail_axis);
set(handles.tail_axis,'xlim',[-10,0],'color','none','xcolor','white');
hold(handles.tail_axis,'on');
xlabel('time [s]');
ylabel('tail sum [rad]');
handles.tail_plot = plot(handles.tail_axis,0,0,'color','g');
handles.turn_line_obj = line([-10 0],[handles.config.Estimator.tail_thresh handles.config.Estimator.tail_thresh],'color','g','linestyle',':');
yyaxis right;
set(handles.tail_axis,'ycolor','r');
hold(handles.tail_axis,'on');
ylabel(handles.tail_axis,'vigor [rad]');
handles.vigor_plot = plot(handles.tail_axis,0,0,'color','r');
handles.swim_line_obj = line([-10 0],[handles.config.Estimator.vigor_thresh handles.config.Estimator.vigor_thresh],'color','r','linestyle',':');
yyaxis left;
set(handles.tail_axis,'ycolor','g');

%velocity axis
axes(handles.vel_axis)
handles.vel_axis.XLim = [-10 0];
disableDefaultInteractivity(handles.vel_axis);
set(handles.vel_axis,'xlim',[-10,0],'color','none','xcolor','white');
hold(handles.vel_axis,'on');
xlabel(handles.vel_axis,'time [s]');
handles.vel_plot = plot(handles.vel_axis,0,0,'color','g');
set(handles.vel_axis,'ycolor','g');
ylabel(handles.vel_axis,'velocity [cm/s]');
yyaxis right;
set(handles.vel_axis,'ycolor','r');
hold(handles.vel_axis,'on');
ylabel(handles.vel_axis,'angular velocity [°/s]');
handles.ang_vel_plot = plot(handles.vel_axis,0,0,'color','r');
yyaxis left;
set(handles.vel_axis,'ycolor','g');

% set video acquisition parameters
handles.VidSource.ExposureAuto = 'Off';
handles.VidSource.AcquisitionFrameRateEnable = 'True';
set_and_save_config(handles);
CameraFPS_Callback(handles.CameraFPS, [], handles);
CameraExposureTime_Callback(handles.CameraExposureTime, [], handles);
handles = get_handles(handles);

% change video ROI
handles.Vid.ROIPosition = handles.config.Acquisition.AcquisitionROI;
handles.sz = flip(handles.config.Acquisition.AcquisitionROI(3:4));
temp = struct;
temp.Source = handles.line;
temp.CurrentPosition = [handles.config.Tracking.Point1;handles.config.Tracking.Point2];
update_handles(handles);
line_moved([],temp);
handles = get_handles(handles);
handles.run = true;
update_handles(handles);
save_config(handles);
handles = get_handles(handles);

% start streaming video
setappdata(handles.hImage,'UpdatePreviewWindowFcn',@GUI_update_fcn);
preview(handles.Vid, handles.hImage);

% update guidata and set output
handles.output = hObject;
update_handles(handles);
guidata(hObject, handles);

%% --- Main function
function GUI_update_fcn(~,event,himage)
% get handles
tic
handles = getappdata(himage,'HandleToHandles');
if handles.run
    % timing
    t0 = handles.Data.CurrentData.t;
    handles.Data.CurrentData.t = datetime(event.Timestamp,'InputFormat','HH:mm:ss.SSS');
    handles.Data.CurrentData.dt = seconds(handles.Data.CurrentData.t - t0);
    handles.ActualFPS.String = round(1/handles.Data.CurrentData.dt);
    % get the image
    if handles.use_video
        if handles.video.t >= handles.video.dur-1
            handles.video.t = 0;
        end
        handles.video.t = handles.video.t + handles.Data.CurrentData.dt;
        id = ceil(handles.video.FPS*handles.video.t);
        A = handles.video.data(:,:,id);
    else
        A = event.Data;
    end

    % blur the image
    if handles.config.Tracking.Blurring > 0
        A = uint8(conv2(A,handles.filter,'same'));
    end

    % pre allocate vectors for tail coordinates
    r = [handles.config.Tracking.Point1(2); nan(handles.config.Tracking.nSegments,1)];
    c = [handles.config.Tracking.Point1(1); nan(handles.config.Tracking.nSegments,1)];
    alphas = [0; nan(handles.config.Tracking.nSegments,1)];

    for i=2:handles.config.Tracking.nSegments+1
        % find circle around the previous points
        px_ids = sub2ind(handles.sz,r(i-1),c(i-1))+handles.offsets;
        [r0,c0] = ind2sub(handles.sz,px_ids);
        theta = cart2pol(c0-c(i-1),r0-r(i-1));
        theta = theta-handles.pi2;
        theta(theta<-pi) = theta(theta<-pi)+handles.two_pi;

        % keep pixels within a segment
        keep_these = theta-alphas(i-1)>=-handles.config.Tracking.TrackingAngle & theta-alphas(i-1)<=handles.config.Tracking.TrackingAngle;
        px_ids = px_ids(keep_these);

        % find the darkest pixel
        [~,darkest_pixel] = min(A(px_ids));
        darkest_pixel = px_ids(darkest_pixel);
        [r(i),c(i)] = ind2sub(handles.sz, darkest_pixel);
        if i==2
            alphas(i) = 0;
        else
            alphas(i) = cart2pol(c(i)-c(i-1),r(i)-r(i-1))-handles.pi2;
            if alphas(i)<-pi
                alphas(i) = alphas(i) + handles.two_pi;
            end
            % comment this out
            [rr,cc] = ind2sub(handles.sz,px_ids);
            handles.semicircles_obj(i-1).XData = cc;
            handles.semicircles_obj(i-1).YData = rr;
        end
    end
    r(1) = [];
    c(1) = [];
    handles.Data.CurrentData.tail = sum(alphas);
    
    [Data] = tail_tracker_estimator(handles.Data,handles.config.Estimator);
    handles.Data = Data;
    
    % display things
    himage.CData = A;
    handles.tail_point_obj.XData = c;
    handles.tail_point_obj.YData = r;
    handles.tail_plot.XData = handles.Data.BufferData.t;
    handles.tail_plot.YData = handles.Data.BufferData.tail;
    handles.vigor_plot.XData = handles.Data.BufferData.t;
    handles.vigor_plot.YData = handles.Data.BufferData.vigor;
    handles.vel_plot.XData = handles.Data.BufferData.t;
    handles.vel_plot.YData = handles.Data.BufferData.velocity;
    handles.ang_vel_plot.XData = handles.Data.BufferData.t;
    handles.ang_vel_plot.YData = handles.Data.BufferData.angular_velocity;

    % update handles %
    update_handles(handles);

    % how much time did the whole iteration take %
    handles.ComputationFPS.String = round(1/toc);
end

% --- Set and save config
function set_and_save_config(handles)
handles.CameraFPS.String = handles.config.Acquisition.CameraFPS;
handles.CameraExposureTime.String = handles.config.Acquisition.CameraExposureTime;
handles.nSegments.String = handles.config.Tracking.nSegments;
handles.Blurring.String = handles.config.Tracking.Blurring;
handles.TrackingAngle.String = rad2deg(handles.config.Tracking.TrackingAngle);
handles.VigorWindow.String = handles.config.Estimator.vigor_window;
handles.VigorWindow_s = handles.config.Estimator.vigor_window;
handles.SwimThresh.String = handles.config.Estimator.vigor_thresh;
handles.vel_coefficient.String = handles.config.Estimator.forward_coef;
handles.ang_vel_coefficient.String = handles.config.Estimator.turning_coef;
handles.decay_duration.String = handles.config.Estimator.inertia_kernel_tau;
handles.decay_coefficient.String = handles.config.Estimator.inertia_kernel_decay;
handles.tail_treshold_edit.String = handles.config.Estimator.tail_thresh;
handles.reset_window.String = handles.config.Estimator.reset_window;
handles.smoothing.String = handles.config.Estimator.tail_smoothing;
handles.v_ang_kernel_duration.String = handles.config.Estimator.v_ang_kernel_tau;

handles.filter = ones(handles.config.Tracking.Blurring,handles.config.Tracking.Blurring)/(handles.config.Tracking.Blurring^2);
save_config(handles);
function save_config(handles)
update_handles(handles);
config = handles.config;
save([handles.pathname 'Tail_Tracker_config.mat'],'config');

% --- Tracking
function line_moved(~,evtData)
handles = getappdata(findobj(evtData.Source.Parent,'type','image'),'HandleToHandles');
currentPosition = round(evtData.CurrentPosition);
handles.config.Tracking.Point1 = currentPosition(1,:);
handles.config.Tracking.Point2 = currentPosition(2,:);
% compute the segment length
nSegments = handles.config.Tracking.nSegments;
dTip = sqrt((currentPosition(1,1)-currentPosition(2,1))^2+(currentPosition(1,2)-currentPosition(2,2))^2);

[columnsInImage, rowsInImage] = meshgrid(1:handles.sz(2), 1:handles.sz(1));
circlePixels = (rowsInImage - currentPosition(1,2)).^2 + (columnsInImage - currentPosition(1,1)).^2 <= (dTip/nSegments)^2;
circlePixels = edge(circlePixels);
circlePixels_ids = find(circlePixels);
cent = sub2ind(handles.sz,currentPosition(1,2),currentPosition(1,1));
handles.offsets = circlePixels_ids-cent;

delete(handles.semicircles_obj);
for i=1:nSegments
    handles.semicircles_obj(i) = plot(handles.cam_image,0,0,'color','w','linestyle','none','marker','.','MarkerSize',1);
end

set_and_save_config(handles);

function nSegments_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Tracking.nSegments = round(str2double(hObject.String));
temp = struct;
temp.Source = handles.line;
temp.CurrentPosition = handles.line.Position;
update_handles(handles);
line_moved([],temp);
handles = get_handles(handles);
set_and_save_config(handles);

function Blurring_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Tracking.Blurring = round(str2double(hObject.String));
set_and_save_config(handles);

function TrackingAngle_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Tracking.TrackingAngle = deg2rad(str2double(hObject.String));
temp = struct;
temp.Source = handles.line;
temp.CurrentPosition = handles.line.Position;
update_handles(handles);
line_moved([],temp);
handles = get_handles(handles);
set_and_save_config(handles);


% --- Estimator
function VigorWindow_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.vigor_window = str2double(hObject.String);
set_and_save_config(handles);

function SwimThresh_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.vigor_thresh = str2double(hObject.String);
handles.swim_line_obj.YData = [handles.config.Estimator.vigor_thresh handles.config.Estimator.vigor_thresh];
set_and_save_config(handles);

function tail_treshold_edit_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.tail_thresh = str2double(hObject.String);
handles.turn_line_obj.YData = [handles.config.Estimator.tail_thresh handles.config.Estimator.tail_thresh];
set_and_save_config(handles);

function decay_duration_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.inertia_kernel_tau = str2double(hObject.String);
set_and_save_config(handles);

function vel_coefficient_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.forward_coef = str2double(hObject.String);
set_and_save_config(handles);

function decay_coefficient_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.inertia_kernel_decay = str2double(hObject.String);
set_and_save_config(handles);

function ang_vel_coefficient_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.turning_coef = str2double(hObject.String);
set_and_save_config(handles);

function reset_window_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.reset_window = str2double(hObject.String);
set_and_save_config(handles);

function smoothing_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.tail_smoothing = str2double(hObject.String);
set_and_save_config(handles);

function v_ang_kernel_duration_Callback(hObject, ~, handles)
handles = get_handles(handles);
handles.config.Estimator.v_ang_kernel_tau = str2double(hObject.String);
set_and_save_config(handles);

% --- Acquisition
function CameraFPS_Callback(hObject, ~, handles)
handles = get_handles(handles);
fps = str2double(hObject.String);
handles.VidSource.AcquisitionFrameRate = fps;
fps = handles.VidSource.AcquisitionFrameRate;
handles.VidSource.AcquisitionFrameRate = fps;
handles.config.Acquisition.CameraFPS = fps;
set_and_save_config(handles);
function CameraExposureTime_Callback(hObject, ~, handles)
handles = get_handles(handles);
ExposureTime = str2double(hObject.String);
ExposureTime = min(ExposureTime,1/75*1000);
handles.VidSource.ExposureTime = ExposureTime*1000;
ExposureTime = handles.VidSource.ExposureTime/1000;
handles.VidSource.ExposureTime = ExposureTime*1000;
handles.config.Acquisition.CameraExposureTime = ExposureTime;
set_and_save_config(handles);
function SetAcquisitionROIButton_Callback(hObject, ~, handles)
handles = get_handles(handles);
if strcmp(hObject.String, 'Draw acquisition ROI')
    hObject.String = 'Set acquisition ROI';
    handles.aq_ROI_obj = images.roi.Rectangle(handles.cam_image,'Position',[1,1,100,100]);
    update_handles(handles);
else
    switch hObject.String
        case 'Reset acquisition ROI'
            hObject.String = 'Draw acquisition ROI';
            handles.config.Acquisition.AcquisitionROI = handles.GenericAcquisitionROI;
        case 'Set acquisition ROI'
            hObject.String = 'Reset acquisition ROI';
            aq_ROI = round(handles.aq_ROI_obj.Position);
            if rem(aq_ROI(3),2)~=0
                aq_ROI(3) = aq_ROI(3) - 1;
            end
            if rem(aq_ROI(4),2)~=0
                aq_ROI(4) = aq_ROI(4) - 1;
            end
            delete(handles.aq_ROI_obj);
            aq_ROI(1:2) = aq_ROI(1:2)-1;
            handles.config.Acquisition.AcquisitionROI = aq_ROI;
    end
    update_handles(handles);
    update_acquisition_ROI(handles)
end
function update_acquisition_ROI(handles)
handles = get_handles(handles);
handles.run = false;
% change video ROI
handles.Vid.ROIPosition = handles.config.Acquisition.AcquisitionROI;
handles.config.Acquisition.AcquisitionROI = handles.Vid.ROIPosition;
handles.line.Position = [10,10;100,100];
handles.sz = flip(handles.config.Acquisition.AcquisitionROI(3:4));
update_handles(handles);
temp = struct;
temp.Source = handles.line;
temp.CurrentPosition = handles.line.Position;
update_handles(handles);
line_moved([],temp);
handles = get_handles(handles);
handles.run = true;
update_handles(handles);
save_config(handles);

% --- Load video
function LoadVideoButton_Callback(~, ~, handles)
handles = get_handles(handles);
handles.run = false;
update_handles(handles);
handles.use_video = true;
[filename,pathname] = uigetfile(handles.pathname);
v = VideoReader([pathname, filename]);
handles.video.nFrames = v.NumFrames;
handles.video.FPS = v.FrameRate;
handles.video.t = 0;
handles.video.dur = handles.video.nFrames/handles.video.FPS;
handles.video.data = zeros(v.Height,v.Width,handles.video.nFrames,'uint8');
handles.line.Position = [10,10;100,100];
handles.sz = [v.Height,v.Width];
update_handles(handles);
temp = struct;
temp.Source = handles.line;
temp.CurrentPosition = handles.line.Position;
update_handles(handles);
line_moved([],temp);
handles = get_handles(handles);
nbytes = fprintf('Loading video frames: 0 of %d', handles.video.nFrames);
for i=1:handles.video.nFrames
    this_frame = read(v,i);
    handles.video.data(:,:,i) = uint8(this_frame(:,:,1));
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('Loading video frames: %d of %d\n',i,handles.video.nFrames);
end
handles.run = true;
update_handles(handles);

% --- Update and retrieve handles
function update_handles(handles)
setappdata(handles.hImage,'HandleToHandles',handles);
function [handles] = get_handles(handles)
handles = getappdata(handles.hImage,'HandleToHandles');

function nSegments_CreateFcn(hObject, eventdata, handles)
function Blurring_CreateFcn(hObject, eventdata, handles)
function TrackingAngle_CreateFcn(hObject, eventdata, handles)
function VigorWindow_CreateFcn(hObject, eventdata, handles)
function SwimThresh_CreateFcn(hObject, eventdata, handles)
function CameraFPS_CreateFcn(hObject, eventdata, handles)
function CameraExposureTime_CreateFcn(hObject, eventdata, handles)
function vel_coefficient_CreateFcn(hObject, eventdata, handles)
function decay_duration_CreateFcn(hObject, eventdata, handles)
function decay_coefficient_CreateFcn(hObject, eventdata, handles)
function ang_vel_coefficient_CreateFcn(hObject, eventdata, handles)
function tail_treshold_edit_CreateFcn(hObject, eventdata, handles)
function reset_window_CreateFcn(hObject, eventdata, handles)
function smoothing_CreateFcn(hObject, eventdata, handles)
function v_ang_kernel_duration_CreateFcn(hObject, eventdata, handles)


function varargout = Tail_Tracker_GUI_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
