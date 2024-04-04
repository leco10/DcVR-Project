function [Data] = tail_tracker_estimator(Data,par)

% add time to decay kernel
Data.CurrentData.time_since_bout_onset= Data.CurrentData.time_since_bout_onset + Data.CurrentData.dt; % in s

% smooth the tail and compute tail derivative
tail_LPF0 = Data.CurrentData.tail_LPF;
Data.CurrentData.tail_LPF = Data.CurrentData.tail_LPF*par.tail_smoothing + Data.CurrentData.tail*(1-par.tail_smoothing);
Data.CurrentData.tail_derivative = Data.CurrentData.tail_LPF - tail_LPF0;

% buffer the tail for computing the vigor
Data.BufferData.t = [Data.BufferData.t Data.CurrentData.dt];
Data.BufferData.t = Data.BufferData.t - Data.CurrentData.dt;
Data.BufferData.tail = [Data.BufferData.tail Data.CurrentData.tail];
Data.BufferData.tail_derivative = [Data.BufferData.tail_derivative Data.CurrentData.tail_derivative];

% compute and buffer the vigor
tf_vigor = Data.BufferData.t>-par.vigor_window;
Data.CurrentData.vigor = std(Data.BufferData.tail(tf_vigor));
Data.BufferData.vigor = [Data.BufferData.vigor Data.CurrentData.vigor];

% decide if we are swimming forward and compute forward velocity
last_swim_forw = Data.CurrentData.swim_forw;
Data.CurrentData.swim_forw = Data.CurrentData.vigor > par.vigor_thresh;
if Data.CurrentData.swim_forw
    Data.CurrentData.v_forward = Data.CurrentData.vigor*par.forward_coef;
else
    %calculate inertia
    if last_swim_forw && ~Data.CurrentData.swim_forw
        Data.v_forward_inertia_list = Data.BufferData.vigor(end-1)*par.forward_coef*Data.inertia_kernel;
    end
    if ~isempty(Data.v_forward_inertia_list)
        Data.CurrentData.v_forward = Data.v_forward_inertia_list(1);
        Data.v_forward_inertia_list(1) = [];
    else
        Data.CurrentData.v_forward = 0;
    end
end

% decide is we are turning
if ~Data.CurrentData.swim_ang && abs(Data.CurrentData.tail_derivative)>par.tail_thresh
    Data.CurrentData.swim_ang = true;
    Data.CurrentData.time_since_bout_onset = 0;
end
tf_not_turning = Data.BufferData.t>-par.reset_window;
if Data.CurrentData.swim_ang
    %check if bout ended
    if all(abs(Data.BufferData.tail_derivative(tf_not_turning))<par.tail_thresh)
        Data.CurrentData.swim_ang = false;
        Data.CurrentData.v_angular = 0;
    else
        if sign(Data.CurrentData.tail_derivative) == sign(Data.CurrentData.tail)
            n_10ms = min(100,round(Data.CurrentData.time_since_bout_onset*100)+1);
            Data.CurrentData.v_angular = Data.CurrentData.tail_derivative*-par.turning_coef*Data.ang_vel_dampening_y(n_10ms);
        else
            Data.CurrentData.v_angular = 0;
        end
    end
else
    Data.CurrentData.v_angular = 0;
end

% buffer velocities for plotting
Data.BufferData.velocity = [Data.BufferData.velocity Data.CurrentData.v_forward];
Data.BufferData.angular_velocity = [Data.BufferData.angular_velocity Data.CurrentData.v_angular];

% remove old values from the 10 s buffer
tf_remove = Data.BufferData.t<-10;
Data.BufferData.t(tf_remove) = [];
Data.BufferData.tail(tf_remove) = [];
Data.BufferData.tail_derivative(tf_remove) = [];
Data.BufferData.vigor(tf_remove) = [];
Data.BufferData.velocity(tf_remove) = [];
Data.BufferData.angular_velocity(tf_remove) = [];

% % save things
% Data.SaveData.Time = [Data.SaveData.Time single(Data.CurrentData.t)];
% Data.SaveData.Tail = [Data.SaveData.Tail single(Data.CurrentData.tail)];
% Data.SaveData.Vigor = [Data.SaveData.Vigor single(Data.CurrentData.vigor)];
% Data.SaveData.Tail_LPF =[Data.SaveData.Tail_LPF single(Data.CurrentData.tail_LPF)];
% Data.SaveData.Tail_derivative =[Data.SaveData.Tail_derivative single(Data.CurrentData.tail_derivative)];
% Data.SaveData.V_forward =[Data.SaveData.V_forward single(Data.CurrentData.v_forward)];
% Data.SaveData.V_angular =[Data.SaveData.V_angular single(Data.CurrentData.v_angular)];


