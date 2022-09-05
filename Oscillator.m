% U= [u_1 u_2]       
% X = [x1 x2]

import casadi.*
xi = 2;         %speed constant
X_norm = 1;     %Radius
w_norm = 3;     %resonance frequency
slacks = MX;


% DGL creation here 
OsciDGL = @(x,u) [xi*(2*X_norm^2-norm(x)^2) -w_norm; w_norm xi*(2*X_norm^2-norm(x)^2)]*x+u;



%%
for i = 1:16
%[DISREGARD]Used for storing the results     
varmap = strcat('x',num2str(i));
varmap1 = strcat('u',num2str(i));
varmap2 = strcat('t',num2str(i));
varmap3 = strcat('V',num2str(i));
% Time horizon and increments
tspan = 18;
dt = 0.05;
NumInc = floor(tspan/dt);;
time = [0:dt:NumInc*dt];
%Tolerance for goal position 
Epsilon = 0.1; 
% Goal and start position
Goal = [-0.1 -2+i]; 
Start = [-3 0.01];
% Box constraints on the ocp
x_box = [-Inf Inf; -Inf, Inf];
u_box = [-Inf Inf ; -Inf, Inf];
% Create the OCP from the ODE 
[ocp, x, u,varout{1:6}] = ode2ocp(OsciDGL, 2, 2, NumInc, dt, x_box=x_box, u_box=u_box,x0='param',foh=true,verbose=1);
ocp.solver('ipopt')
ocp.set_value(u(:,1),[0 0]); %% first value for u has to be zero. 

x0_p = varout{2};
% set x0
ocp.set_value(x0_p,Start);
% set cost function
ocp.minimize(sum(u(1,:).^2+u(2,:).^2));
%set terminal constraint
ocp.subject_to(x(1,end)<=Goal(1)+Epsilon);
ocp.subject_to(x(1,end)>=Goal(1)-Epsilon);
% % % 
ocp.subject_to(x(2,end)<=Goal(2)+Epsilon);
ocp.subject_to(x(2,end)>=Goal(2)-Epsilon);

%set minimal distance constraint
ocp.subject_to((x(1,:).^2+x(2,:).^2)>=0.2);

% For iterative solutions try to set the initial condition
% to the one from the previous iteration

% IMPORTANT!! You need a good initialization for good results
% System is very prone to local minima.R un the ocp
% without goal constraints and take the result as an initial guess
try
ocp.set_initial(x(:,2:end),xxx(:,2:end));
disp('ok')
catch 
    ocp.set_initial(x(:,2:end),repmat(Start,NumInc,1)')
disp('nok')
end
%[DISREGARD] This is just for plotting. It is very inelegant
% so probably best to leave alone. 
sol = ocp.solve()
xxx = ocp.value(x);
uuu = ocp.value(u);
Coordinates= zeros(2,length(xxx));
Velocity= zeros(1,length(xxx));
for j = 1:length(xxx)
        [xx,yy]= cart2pol((xxx(1,j)),(xxx(2,j)));
        Coordinates(1,j)=xx;
        Coordinates(2,j)=yy;
end
V = diff(Coordinates(1,:))/dt;
V = [V V(end)];  %% Add one value for plotting purposes
Values.(varmap) =xxx;
Values.(varmap1)=uuu;
Values.(varmap2)=time;
Values.(varmap3)=V;
AngleFull = [0:0.001:pi*2];
[Atmo1 Atmo2] = arrayfun(@(x) pol2cart(x,sqrt(2)*X_norm),AngleFull);
figure(4)
scatter(0,0,10,'blue','filled')
hold on 
scatter(Goal(1),Goal(2),50,'+')
% plot([Goal(1)*(1-Epsilon) Goal(1)*(1-Epsilon) ],[Goal(2)*(1-Epsilon)  Goal(2)*(1+Epsilon)],'r',LineWidth=5)
% plot([Goal(1)*(1+Epsilon)  Goal(1)*(1+Epsilon) ],[Goal(2)*(1-Epsilon)  Goal(2)*(1+Epsilon)],'r',LineWidth=5)
% plot([Goal(1)*(1-Epsilon)  Goal(1)*(1+Epsilon) ],[Goal(2)*(1-Epsilon)  Goal(2)*(1-Epsilon)],'r',LineWidth=5)
% plot([Goal(1)*(1-Epsilon)  Goal(1)*(1+Epsilon) ],[Goal(2)*(1+Epsilon)  Goal(2)*(1+Epsilon)],'r',LineWidth=5)
rectangle('Position',[Goal(1)-Epsilon Goal(2)-Epsilon (2*Epsilon) (2*Epsilon)])
scatter(Atmo1,Atmo2,10,'black')
scatter(Start(1),Start(2),500,'blued','filled')
axis equal
scatter(xxx(1,end),xxx(2,end),500,'blue^','filled')
scatter(xxx(1,:),xxx(2,:),50,V,'filled')
 colormap(gca,'turbo')
 cb = colorbar; % creates the colorbar on side
 cb.Label.String = 'Velocity ';
% use caxis to change range on colorbar
caxis([3-0.5,3+0.5]); 
end
%%
xxx = ocp.value(x);
uuu = ocp.value(u);


    figure(1)
    subplot(2,1,1)

    scatter(time,xxx(1,:))
    xlabel('Time')
    ylabel('Position')
    legend('x0 = 0.1 , T = 200','x0 = -0.1, T = 150','xf = 0.1 , T = 150')
    hold on
    subplot(2,1,2)
    plot(time,xxx(2,:))
    legend('x0 = 0.1 , T = 200','x0 = -0.1, T = 150','xf = 0.1 , T = 150')
    hold on
    xlabel('Time')
    ylabel('Input u')
 
    
    %%
    CCC = zeros(length(xxx),1);
    for i = 1:length(xxx)
        CCC(i) = cost_energy(xxx(:,i),uuu(i),time(i));
    end


function x_end = integrator_step_disturbed(x0, u, dt, odefun, d)
% calculate one integration step with step size dt
import casadi.*

x0_rk = x0;
k = casadi.MX( size( x0, 1 ), 4 );

if size(u,2) == 1
    k(:,1) = odefun(x0_rk(:,end)                  , u, d);
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u, d);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u, d);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u, d);
else
    k(:,1) = odefun(x0_rk(:,end)                  , u(:,1),       d(:,1));
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u(:,2),       d(:,2));
end

x_end  = x0_rk + dt / 6 * k * [1 2 2 1]';

end