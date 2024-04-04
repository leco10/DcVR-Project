function Tail_Tracker_script(exp_time)
disp([datestr(datetime('now')) ': started Tail Tracker script']);

exp_id = "Exp08_2D_OMR";
fish_id = input("Enter Fish ID: ","s");

Data.pathname = 'C:\Users\jlab\Desktop\Leo\TailTracker\';
cd(Data.pathname);

load([Data.pathname 'Tail_Tracker_config.mat'],'config');
Data.filter = ones(config.Tracking.Blurring,config.Tracking.Blurring)/(config.Tracking.Blurring^2);
% create camera streaming object
imaqreset;
Vid = videoinput('gentl','1','Mono8');
VidSource = getselectedsource(Vid);
VidSource.BinningHorizontal = 2;
VidSource.BinningHorizontalMode = 'Sum';
VidSource.BinningVertical = 2;
VidSource.BinningVerticalMode = 'Sum';
VidSource.GainAuto = 'Off';
VidSource.Gain = 1;
VidSource.ExposureAuto = 'Off';
VidSource.AcquisitionFrameRateEnable = 'True';
VidSource.AcquisitionFrameRate = 75;%config.Acquisition.CameraFPS; %use max fps by default
VidSource.ExposureTime = config.Acquisition.CameraExposureTime*1000;
Vid.ROIPosition = config.Acquisition.AcquisitionROI;

% some variables

%toggle save video
Data.save_video = false;

% current data
Data.CurrentData.t = 0;
Data.CurrentData.tail = 0;
Data.CurrentData.vigor = 0;
Data.CurrentData.velocity = 0;
Data.CurrentData.angular_velocity = 0;
Data.CurrentData.time_since_bout_onset = 0;
Data.CurrentData.tail_derivative = 0;
Data.CurrentData.swim_forw = false; 
Data.CurrentData.swim_ang = false; 
Data.CurrentData.tail_LPF = 0;

% buffers
Data.BufferData.t = [];
Data.BufferData.tail = [];
Data.BufferData.vigor = [];
Data.BufferData.velocity = [];
Data.BufferData.angular_velocity = [];
Data.BufferData.tail_derivative = [];

% decay kernels
%Data.ang_vel_dampening_y = exp((0.01:0.01:1)*-log(1/config.Estimator.inertia_kernel_decay)/(config.Estimator.v_ang_kernel_tau)); % 0.55: mean bout duration, 100: ang_vel dampening 0.01 at end of mean bout
c = [ones(1,20) zeros(1,80)];
d = sgolayfilt(c,1,11);
Data.ang_vel_dampening_y = sgolayfilt(d,1,7);
Data.inertia_kernel = exp((0.01:0.01:1)*-log(1/config.Estimator.inertia_kernel_decay)/(config.Estimator.inertia_kernel_tau));
Data.v_forward_inertia_list = [];
% saving
Data.SaveData.Time = single(0);
Data.SaveData.Tail = single(0);
Data.SaveData.Vigor = single(0);
Data.SaveData.Tail_LPF = single(0);
Data.SaveData.Tail_derivative = single(0);
Data.SaveData.V_forward = single(0);
Data.SaveData.V_angular = single(0);

Data.pi2 = pi/2;
Data.two_pi = 2*pi;
Data.exp_duration = exp_time;
Data.filename = 'vigor_data.txt';
% offsets for tail segments
Data.sz = flip(config.Acquisition.AcquisitionROI(3:4));
point1 = config.Tracking.Point1;
point2 = config.Tracking.Point2;

% compute the segment length
dTip = sqrt((point1(1)-point2(1))^2+(point1(2)-point2(2))^2);

[columnsInImage, rowsInImage] = meshgrid(1:Data.sz(2), 1:Data.sz(1));
circlePixels = (rowsInImage - point1(2)).^2 + (columnsInImage - point1(1)).^2 <= (dTip/config.Tracking.nSegments)^2;
circlePixels = edge(circlePixels);
circlePixels = find(circlePixels);
Data.offsets = circlePixels-sub2ind(Data.sz,point1(2),point1(1));
Data.config = config;

% create graphical objects
% figure with tail plot
figure('color','k','Position',[100 100 300 300]);
tail_axis = axes('color','k','xcolor','w'); hold on;
subplot(1,2,1, tail_axis)
disableDefaultInteractivity(gca);
set(gca,'xlim',[-10 0]);
xlabel('time [s]');
yyaxis left;
set(gca,'ycolor','g');
ylabel('tail sum [°]');
Data.tail_plot_obj = plot(0,0,'g');
line([-10 0],[Data.config.Estimator.tail_thresh Data.config.Estimator.tail_thresh],'color','g','linestyle',':');
yyaxis right;
set(gca,'ycolor','r');
ylabel('vigor [°]');
Data.vigor_plot_obj = plot(0,0,'color','r');
line([-10 0],[Data.config.Estimator.vigor_thresh Data.config.Estimator.vigor_thresh],'color','r','linestyle',':');

vel_axis = axes('color','k','xcolor','w'); hold on;
subplot(1,2,2, vel_axis)
vel_axis.XLim = [-10 0];
disableDefaultInteractivity(vel_axis);
set(vel_axis,'xlim',[-10,0],'color','none','xcolor','white');
hold(vel_axis,'on');
xlabel(vel_axis,'time [s]');
Data.vel_plot_obj = plot(vel_axis,0,0,'color','g');
set(vel_axis,'ycolor','g');
ylabel(vel_axis,'velocity [cm/s]');
yyaxis right;
set(vel_axis,'ycolor','r');
hold(vel_axis,'on');
ylabel(vel_axis,'angular velocity [°/s]');
Data.ang_vel_plot_obj = plot(vel_axis,0,0,'color','r');
yyaxis left;
set(vel_axis,'ycolor','g');

% figure with camera image
figure('color','k','Position',[300 100 500 500]);
axes;
Data.im_obj = imshow(getsnapshot(Vid));
hold on;
for in=1:config.Tracking.nSegments
    Data.semicircles_obj(in) = plot(0,0,'color','w','linewidth',1);
end
Data.tail_point_obj = plot(0,0,'marker','.','color','g','LineStyle','-','MarkerSize',20);
Data.fps_txt_obj = text(0,0,'FPS: 0 Hz','horizontalalignment','left','verticalalignment','top','fontsize',10,'FontWeight','bold','color','w');
Data.time_txt_obj = text(0,15,'Time: 0 s','horizontalalignment','left','verticalalignment','top','fontsize',10,'FontWeight','bold','color','w');
disableDefaultInteractivity(gca);
Data.display_timer = 0;
Data.display_rate = 10; % Hz
Data.display_period = 1/Data.display_rate;

% start video recording
if Data.save_video
    vid_filename = [char(datetime,'yyMMddHHmmss') '_fish_1_video.mat']; %change fish number!!
    Data.VideoWriter = VideoWriter(append('E:\',vid_filename),'MPEG-4');
    Data.VideoWriter.FrameRate = 75;
    Data.VideoWriter.Quality = 70;
    open(Data.VideoWriter);
end
pause(5)
% start streaming video
Vid.FramesPerTrigger = inf;
Vid.TimerPeriod = 1/VidSource.AcquisitionFrameRate;
Vid.TimerFcn = @update_fcn;
clear circlePixels columnsInImage rowsInImage config dTip i point1 point2;
Data.start_time = datetime('now');
Vid.UserData = Data;
clear Data;
start(Vid);
    function obj = update_fcn(obj,event)
        Data = obj.UserData;

        % timing
        t0 = Data.CurrentData.t;
        Data.CurrentData.t = seconds(datetime(event.Data.AbsTime,'InputFormat','HH:mm:ss.SSS')-Data.start_time);
        Data.CurrentData.dt = Data.CurrentData.t - t0;


        % get the image
        A = getdata(obj,1);

        % blur the image
        if Data.config.Tracking.Blurring > 0
            A = uint8(conv2(A,Data.filter,'same'));
        end

        % pre allocate vectors for tail coordinates
        r = [Data.config.Tracking.Point1(2); nan(Data.config.Tracking.nSegments,1)];
        c = [Data.config.Tracking.Point1(1); nan(Data.config.Tracking.nSegments,1)];
        alphas = [0; nan(Data.config.Tracking.nSegments,1)];

        for i=2:Data.config.Tracking.nSegments+1
            % find circle around the previous points
            px_ids = sub2ind(Data.sz,r(i-1),c(i-1))+Data.offsets;
            [r0,c0] = ind2sub(Data.sz,px_ids);
            theta = cart2pol(c0-c(i-1),r0-r(i-1));
            theta = theta-Data.pi2;
            theta(theta<-pi) = theta(theta<-pi)+Data.two_pi;

            % keep pixels within a segment
            keep_these = theta-alphas(i-1)>=-Data.config.Tracking.TrackingAngle & theta-alphas(i-1)<=Data.config.Tracking.TrackingAngle;
            px_ids = px_ids(keep_these);

            % find the darkest pixel
            [~,darkest_pixel] = min(A(px_ids));
            darkest_pixel = px_ids(darkest_pixel);
            [r(i),c(i)] = ind2sub(Data.sz, darkest_pixel);
            if i==2
                alphas(i) = 0;
            else
                alphas(i) = cart2pol(c(i)-c(i-1),r(i)-r(i-1))-Data.pi2;
                if alphas(i)<-pi
                    alphas(i) = alphas(i) + Data.two_pi;
                end
                % comment this out
                [rr,cc] = ind2sub(Data.sz,px_ids);
                Data.semicircles_obj(i-1).XData = cc;
                Data.semicircles_obj(i-1).YData = rr;
            end
        end
        r(1) = [];
        c(1) = [];
        Data.CurrentData.tail = sum(alphas);
        
        % call estimator function
        [Data] = tail_tracker_estimator(Data,Data.config.Estimator);

        % disply things
        Data.display_timer = Data.display_timer + Data.CurrentData.dt;
        if Data.display_timer >= Data.display_period
            Data.display_timer = 0;
            Data.tail_plot_obj.XData = Data.BufferData.t;
            Data.tail_plot_obj.YData = Data.BufferData.tail;
            Data.vigor_plot_obj.XData = Data.BufferData.t;
            Data.vigor_plot_obj.YData = Data.BufferData.vigor;
            Data.vel_plot_obj.XData = Data.BufferData.t;
            Data.vel_plot_obj.YData = Data.BufferData.velocity;
            Data.ang_vel_plot_obj.XData = Data.BufferData.t;
            Data.ang_vel_plot_obj.YData = Data.BufferData.angular_velocity;
            Data.im_obj.CData = A;
            Data.fps_txt_obj.String = ['FPS: ' num2str(round(1/Data.CurrentData.dt)) ' Hz'];
            Data.time_txt_obj.String = ['Time: ' num2str(round(Data.CurrentData.t)) ' s'];
            Data.tail_point_obj.XData = c;
            Data.tail_point_obj.YData = r;
            fileID = fopen(Data.filename,'w');
            fprintf(fileID,'%6.3f,%6.3f,%5.2f',[Data.CurrentData.v_forward,Data.CurrentData.v_angular,Data.CurrentData.t]); %velocity,ang velo, time
            fclose(fileID);
        end
        %save video
        if Data.save_video
            nDroppedFrames = round(Data.CurrentData.dt*75)-1; % max frame rate = 75
            writeVideo(Data.VideoWriter,getframe(gcf));
            if nDroppedFrames>0
                for ind=1:nDroppedFrames
                    writeVideo(Data.VideoWriter,getframe(gcf));
                end
            end
        end
      
       % save variables
       Data.SaveData.Time = [Data.SaveData.Time single(Data.CurrentData.t)];
       Data.SaveData.Tail = [Data.SaveData.Tail single(Data.CurrentData.tail)];
       Data.SaveData.Vigor = [Data.SaveData.Vigor single(Data.CurrentData.vigor)];
       Data.SaveData.Tail_LPF =[Data.SaveData.Tail_LPF single(Data.CurrentData.tail_LPF)];
       Data.SaveData.Tail_derivative =[Data.SaveData.Tail_derivative single(Data.CurrentData.tail_derivative)];
       Data.SaveData.V_forward =[Data.SaveData.V_forward single(Data.CurrentData.v_forward)];
       Data.SaveData.V_angular =[Data.SaveData.V_angular single(Data.CurrentData.v_angular)];

        % flush the data
        flushdata(obj,'all');

        % save the data after exp_duration
        if Data.CurrentData.t > Data.exp_duration %exp_duration is given by python script
            local_save_data (Data);
            close all
            stop(obj);
            delete(obj)
            clear obj
        end

        % this should be the last line
        obj.UserData = Data;
    end

    function local_save_data(Data)
        alt_path = Data.pathname;
        filename = [char(datetime,'yyMMddHHmmss') '_data.mat'];
        config = Data.config;
        StartTime = Data.start_time;
        Data = Data.SaveData;
        try
            save(fullfile('E:\leos data', exp_id,fish_id, filename),'Data','config','StartTime');
        catch
            save([alt_path filename],'Data','config','StartTime');
        end
    end
end