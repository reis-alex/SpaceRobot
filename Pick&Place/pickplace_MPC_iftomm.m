addpath(genpath([pwd '\Toolbox_SMS_sym']));
clear all, clc
close all
import casadi.*

% load robot through URDFs, generate model
robot_path_pre = [pwd '\sms_simple_iftomm.urdf'];
[satelite,~] = urdf2robot_flex_visu(robot_path_pre);
[robot_pre_sim]= robot_slx_format(satelite);
satelite = SPART_casadi(robot_path_pre);

% auxiliary skew symmetric matrix
SkewSym = @(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

% for quaternion integration and conversions
quaternion = quaternion();


%% System dimensions
xsat_ini = vertcat(zeros(6,1), ...                  % null initial displacements
                    zeros(3,1), ... % first 3: q_wheels,
                    20*pi/180,30*pi/180,50*pi/180, ... % last 3: q_arm
                    zeros(length(satelite.idx.velocities),1));                               % null velocities

% used as reference, must be commented
% xsat_ini = vertcat(zeros(6,1), ...                  % null initial displacements
%                     zeros(3,1), ... % first 3: q_wheels,
%                     -40*pi/180,150*pi/180,20*pi/180, ... % last 3: q_arm
%                     zeros(length(satelite.idx.velocities),1)); 

% initialize satelite's state vector
R0s = eye(3);
q0s = quaternion.rotm_to_quat(R0s); %initial quaternion, obtained directly from R0

Hs = satelite.dynamics.H;
Cs = satelite.dynamics.C;

link_positions(:,:,1) = full(satelite.kinematics.rL(R0s,xsat_ini(1:6+satelite.robot.n_q)));
pos_t = full(link_positions(:,end,1));

% since r0 = 0, the reference computed here is w.r.t. the body frame

plot3(link_positions(1,4:end,1),link_positions(2,4:end,1),link_positions(3,4:end,1),'-og')
hold on
plot3(link_positions(1,end,1),link_positions(2,end,1),link_positions(3,end,1),'-xr','LineWidth',10)
plot3(link_positions(1,4,1),link_positions(2,4,1),link_positions(3,4,1),'or','LineWidth',10)
axis auto

%% Auxiliary functions
H_sym = satelite.dynamics.H(satelite.R0,satelite.state_vars);
C_sym = satelite.dynamics.C(satelite.R0,satelite.state_vars);

n0 = 6;
nq = satelite.robot.n_q;
nr = 3; 

H00 = H_sym(1:n0, 1:n0);
H0q = H_sym(1:n0, n0+1 : n0+nq);
Hqq = H_sym(n0+1 : n0+nq, n0+1 : n0+nq);

moment_satelite = casadi.Function('h_sat',{satelite.R0,satelite.state_vars},{H00*satelite.state_vars(satelite.idx.velocities(1:6))});
moment_wheels = casadi.Function('h_wheel',{satelite.R0,satelite.state_vars},{H0q(:,1:3)*satelite.state_vars(satelite.idx.velocities(7:9))});
moment_arm = casadi.Function('h_arm',{satelite.R0,satelite.state_vars},{H0q(:,4:end)*satelite.state_vars(satelite.idx.velocities(10:12))});

hsat = H0q(:,4:end)*satelite.state_vars(satelite.idx.velocities(10:12)) + H0q(:,1:3)*satelite.state_vars(satelite.idx.velocities(7:9));
hwheel = H0q(:,1:3)*satelite.state_vars(satelite.idx.velocities(7:9));
harm = H0q(:,4:end)*satelite.state_vars(satelite.idx.velocities(10:12));

% R0 function
angles = satelite.state_vars(1:3);
yaw = angles(3);
pitch = angles(2);
roll = angles(1);
c1 = cos(yaw);
s1 = sin(yaw);
c2 = cos(pitch);
s2 =  sin(pitch);
c3 = cos(roll);
s3 =  sin(roll);
R = [c1*c2, c1*s2*s3 - s1*c3, c1*s2*c3 + s1*s3;
     s1*c2, s1*s2*s3 + c1*c3, s1*s2*c3 - c1*s3;
     -s2,   c2*s3,             c2*c3];
Rf = casadi.Function('R0',{satelite.state_vars(1:3)},{R});

EE_idx = find(strcmp({satelite.robot.links.name},'Link_EE'));
J = satelite.Jacob(EE_idx);
Jm = J.Jm;
Jm = casadi.substitute(Jm, satelite.R0,  R);
mu = casadi.Function('manipulability',{satelite.state_vars([1:12 19:end])},{det(Jm*Jm'+(1e-6)*eye(6))});

%% MPC
tau = casadi.SX.sym('torque_arm',satelite.robot.n_q,1);

% Define MPC problem
opt.model.states   = [satelite.state_vars([1:12 19:end])];
f_mpc = vertcat(satelite.dynamics.ddX(satelite.R0,satelite.state_vars,tau)); 

% texp is the prediction model. The lines below serve as simplifications
texp = [satelite.state_vars(satelite.idx.velocities);f_mpc(7:end);];
texp = casadi.substitute(texp, satelite.state_vars(13:18),  -H00\H0q*satelite.state_vars(19:end));
texp = casadi.substitute(texp, satelite.R0,  R);

% acceleration function
accel =  casadi.substitute(f_mpc,satelite.state_vars(13:18),  -H00\H0q*satelite.state_vars(19:end));
accel = casadi.substitute(accel, satelite.R0,  R);
accel = casadi.Function('accel',{satelite.state_vars([1:12 19:end]), tau},{accel});

% opt is used in building the MPC
opt.model.function = casadi.Function('fmpc',{vertcat(satelite.state_vars([1:12 19:end])),tau},{texp});
opt.model.controls = tau;
opt.continuous_model.integration = 'euler';

opt.dt          = 0.1;
opt.n_controls  = satelite.robot.n_q;          
opt.n_states    = length(opt.model.states);
opt.N           = 4;

% Define parameters to the MPC
opt.parameters.name{1} = 'target';
opt.parameters.name{2} = 'itm_target';
opt.parameters.dim = vertcat([3,1], [3,1]);

% Cost functions

% define the function computing the position of the end effector
ee_fun = satelite.kinematics.rL(R,vertcat(satelite.state_vars(satelite.idx.positions,1)));
ee_fun = casadi.Function('ee',{opt.model.states(1:6+satelite.robot.n_q)},{R'*ee_fun(:,end)});

opt.costs.stage.parameters = opt.parameters.name(2);
opt.costs.stage.function   = @(x,u,varargin) 1e3*sum((ee_fun(x(1:6+satelite.robot.n_q))-varargin{:}(1:3)).^2) + u'*1e3*u;% + (x(1:3)'*1e7*x(1:3) + omega0(x)'*1e6*omega0(x) + u(1:3)'*1e3*u(1:3));%(hsatf(x)'*hsatf(x))*1e4;

opt.costs.general.parameters = opt.parameters.name(1:2);
opt.costs.general.function   = @(x,u,varargin) 1e5*(varargin{end}-varargin{end-1})'*(varargin{end}-varargin{end-1});

% Constraints
opt.constraints.states.upper  = vertcat( 50*pi/180*ones(3,1),  inf*ones(3,1), inf*ones(3,1),  inf*ones(satelite.robot.n_q-3,1),  3600*2*pi/60*ones(3,1),  0.09*ones(satelite.robot.n_q-3,1));
opt.constraints.states.lower  = vertcat(-50*pi/180*ones(3,1), -inf*ones(3,1), -inf*ones(3,1), -inf*ones(satelite.robot.n_q-3,1), -3600*2*pi/60*ones(3,1), -0.09*ones(satelite.robot.n_q-3,1));
opt.constraints.control.upper = vertcat(0.175*ones(3,1),50*ones(3,1));
opt.constraints.control.lower = -opt.constraints.control.upper;

opt.constraints.general.parameters  = opt.parameters.name(2);
opt.constraints.general.function{1} = @(x,varargin) (ee_fun(x(1:6+satelite.robot.n_q))-varargin{:}(1:3));
opt.constraints.general.elements{1} = 'end';
opt.constraints.general.type{1} = 'equality';

opt.constraints.parameters.name  = opt.parameters.name(2);
opt.constraints.parameters.upper =  vertcat(inf*ones(3,1));
opt.constraints.parameters.lower = -vertcat(inf*ones(3,1));

% Define inputs to optimization
opt.input.vector = opt.parameters.name(1);
opt.solver = 'ipopt';

% build MPC
[solver_mpc,args_mpc] = build_mpc(opt);


%% simulation loop
% simulation parameters

ine_wheels = {satelite.robot.links.inertia};
ine_wheels = ine_wheels{2};

mom_arm = 0.01*ones(6,1);
mom_wheels = 0.01*ones(6,2);
mpc_x0   = zeros(1,length(args_mpc.vars{1}));

dq0 = casadi.Function('dq0', {satelite.R0,satelite.state_vars([1:12 19:end])}, {-H00\H0q*satelite.state_vars(19:end)});
xsat(:,1) = vertcat(xsat_ini([1:12 19:end],1));
feval = zeros(31,1);
ref = [-0.1384; -0.4628; -0.7569];

%%
    
for k = 1:2000
    k

    link_positions(:,:,k) = full(R0s(:,:,k)'*satelite.kinematics.rL(R0s(:,:,k),vertcat(xsat(1:6+satelite.robot.n_q,k))));
    mpc_input = vertcat(xsat(:,k),ref);

    sol = solver_mpc('x0', mpc_x0, 'lbx', args_mpc.lbx, 'ubx', args_mpc.ubx, ...
        'lbg', args_mpc.lbg, 'ubg', args_mpc.ubg, 'p', mpc_input);
    mpc_x0        = full(sol.x);
    torque(:,k) = full(sol.x(args_mpc.vars{3}));
    feval = full(opt.model.function(xsat(:,k),torque(:,k)));
    xsat(:,k+1) = xsat(:,k) + opt.dt*feval; %vertcat(xsat(satelite.idx.velocities,k),full(feval));
    
    check_feas(solver_mpc.stats())


    tt(:,k) = full(dq0(R0s(:,:,k),xsat(1:18,k)));
    % update satelite quaternion, R0 according to next omega
    [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(tt(1:3,k),q0s(:,k),opt.dt);
    Rtest = Rf(xsat(1:3,k));
    xsat_full(:,k) = vertcat(xsat(1:12,k),tt(:,k),xsat(13:end,k));
    mom_arm(:,k) = full(moment_arm(R0s(:,:,k),xsat_full(:,k)));
    mom_sate(:,k) = full(moment_satelite(R0s(:,:,k),xsat_full(:,k)));
    mom_wheels(:,k) = full(moment_wheels(R0s(:,:,k),xsat_full(:,k)));
    itmt(:,k) = full(sol.x(end-2:end));

    manip(k) = full(mu(xsat(:,k)));
    accelerations(:,k) = full(accel(xsat(:,k),torque(:,k)));
    cost(k) = full(sol.f);
end

%%
% save('states','xsat')
% save('accelerations','accelerations')
% save('tau','torque')
% ee_pos = squeeze(link_positions(:,end,:));
% save('EE','ee_pos');
% save('inter_ref','itmt')
% save('mom_arm','mom_arm')
% save('mom_wheel','mom_wheels')
% save('mom_sat','mom_sate')
% save('mom_arm','mom_arm')
% save('mom_wheel','mom_wheels')
% save('mom_sat','mom_sate')

%%
close all

figure
plot(itmt(1,:),'-.r')
hold on
plot(ref(1,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--r')
plot(squeeze(link_positions(1,end,:)),'-r')
plot(itmt(2,:),'-.g')
plot(ref(2,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--g')
plot(squeeze(link_positions(2,end,:)),'-g')
plot(itmt(3,:),'-.b')
plot(ref(3,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--b')
plot(squeeze(link_positions(3,end,:)),'-b')
xlabel('Time [s]')
ylabel('Position [m]')
grid on
title('Position of end effector')

figure
plot(xsat(1,:),'r')
hold on
plot(xsat(2,:),'g')
plot(xsat(3,:),'b')
xlabel('Time [s]')
ylabel('Position [rad]')
grid on
title('Angular position of the base')

figure
plot(xsat(13,:),'r')
hold on
plot(xsat(14,:),'g')
plot(xsat(15,:),'b')
xlabel('Time [s]')
ylabel('Velocity [rad/s]')
grid on
title('Angular velocity of the reaction wheels')

figure
plot(xsat(16,:),'r')
hold on
plot(xsat(17,:),'g')
plot(xsat(18,:),'b')
xlabel('Time [s]')
ylabel('Velocity [rad/s]')
grid on
title('Angular velocity of the manipulator joints')

figure
plot(torque(1,:),'r')
hold on
plot(torque(2,:),'g')
plot(torque(3,:),'b')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Reaction wheels torques')

figure
stairs(torque(4,:),'r')
hold on
stairs(torque(5,:),'g')
stairs(torque(6,:),'b')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Manipulator joint torques')

% figure
% plot(xsat(23,:),'r')
% hold on
% plot(xsat(24,:),'g')
% plot(xsat(25,:),'b')
% title('Integral of h_{sat}')


figure
subplot(131)
plot(mom_sate(1,:),'b'); hold on; plot(mom_arm(1,:),'-g'),plot(mom_wheels(1,:),'--r')
subplot(132)
plot(mom_sate(2,:),'b'); hold on; plot(mom_arm(2,:),'-g'),plot(mom_wheels(2,:),'--r')
subplot(133)
plot(mom_sate(3,:),'b'); hold on; plot(mom_arm(3,:),'-g'),plot(mom_wheels(3,:),'--r')
legend('mom sat','mom bras', 'mom roues')

figure
plot(mom_wheels(1,:),'-r')
hold on
plot(mom_wheels(2,:),'-g')
plot(mom_wheels(3,:),'-b')

%%
figure
hold on
plot3(link_positions(1,7,1),link_positions(2,7,1),link_positions(3,7,1),'or','LineWidth',10)
hold on
for i = 1:k
    plot3(link_positions(1,7:end,i),link_positions(2,7:end,i),link_positions(3,7:end,i),'-og')
    plot3(link_positions(1,end,i),link_positions(2,end,i),link_positions(3,end,i),'-xr','LineWidth',10)
    axis auto
    view([-37.5 30])
    drawnow
    grid
end

%%
function [] = check_feas(stats)
        if ~stats.success
            disp('Solver failed.');
            disp(['Return status: ', stats.return_status]);

            if contains(stats.return_status, 'Infeasible')
                error('Problem is infeasible!');
            end
        end
end