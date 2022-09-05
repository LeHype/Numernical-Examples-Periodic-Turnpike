
import casadi.*
slacks = MX;
% DGL creation here 
OsciDGL= @(x,u,d) [x(2);
                   
                u;
                cost_energy(x,u,d)
    ];

x0 = [0 0]';

% Iteration to get multiple trajectories
for i = 1:5
% Initial state
x0 = [-0.04 0 0];
% Time horizon and increments
tspan = 35-5*i;
dt = 0.05;
NumInc = floor(tspan/dt);
time = [0:dt:NumInc*dt];

% Box constraints on the ocp 
x_box = [-Inf Inf; -Inf, Inf; -Inf Inf];
u_box = [-Inf Inf ];
% Create the OCP from the ODE 
[ocp, x, u,varout{1:6}] = ode2ocp(OsciDGL, 3, 1, NumInc, dt, x_box=x_box, u_box=u_box,x0='param',foh=true,nd=1);

ocp.set_value(u(1),0); %% first value for u has to be zero. 
x0_p = varout{2};
xf_p = varout{3};
% set x0
ocp.set_value(x0_p,x0);
% set cost function
ocp.minimize(x(3,end));
ocp.set_value(varout{4},time);

sol = ocp.solve()


xxx = ocp.value(x);
uuu = ocp.value(u);


    figure(1)
%     subplot(2,1,1)

    scatter(time,xxx(1,:),30,repmat(tspan,NumInc+1,1),'filled')
    xlabel('Time')
    ylabel('Position')
%     legend('x0 = 0.1 , T = 200','x0 = -0.1, T = 150','xf = 0.1 , T = 150')
    hold on
%     subplot(2,1,2)
%     plot(time,xxx(2,:))
%     legend('x0 = 0.1 , T = 200','x0 = -0.1, T = 150','xf = 0.1 , T = 150')
%     hold on
%     xlabel('Time')
%     ylabel('Input u')
 
end  
 colormap(gca,'turbo')
 cb = colorbar; % creates the colorbar on side
 cb.Label.String = 'Time Horizon';
% use caxis to change range on colorbar
caxis([15,34]); % low end is 10, high end is 30

function [cost] = cost_energy(x,u,t)
    cost = 0.5*(x(1)-cos(2*pi*t))^2+(x(2)-sin(2*pi*t))^2+u^2;
    
end

