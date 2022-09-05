% X= [Vs V_phi]
% U= [u_s u_phi]
         
% X = [s s_dot phi phi_dot]
m= [1.523598269*10^13 1];     %Mass of planet/spacecrafe
Gamma = 6.6743*10^(-11);      %Gravitational constant
import casadi.*

slacks = MX;
% Createion of DGL
KeplerDGL = @(x,u) [x(2);
                     x(1)*x(4)^2-(Gamma*m(1)*m(2))/(m(2)*x(1)^2)+(1/m(2))*u(1);
                     x(4);
                     -(2/x(1))*x(4)*x(2)+(1/(m(2)*x(1)^2))*u(2)];









%Iteration to create multiple trajectories
for i = 1:20

 
%[DISREGARD]Used for storing the results
varmap = strcat('x',num2str(i));
varmap1 = strcat('u',num2str(i));
varmap2 = strcat('t',num2str(i));
varmap3 = strcat('Coord',num2str(i));
%time horizon is decreased incrementally
tspan = 58+i*2;
%timestep
dt = 0.2;
NumInc = floor(tspan/dt);;
time = [0:dt:NumInc*dt];
%Starting Position
Start = [-100 -40];
% StartVelocities = The spacecraft is on an approaching trajectory
StartVelocities = [-5 0];

Goal = [-100 100];  %% in cartesian coordinates only defines distace from earth

% Some trigonometrics to converto to polar coordinates
[StartPolar(1) StartPolar(2)] = cart2pol(Start(1),Start(2));
[GoalPolar(1) GoalPolar(2)] = cart2pol(Goal(1),Goal(2)) ;
if GoalPolar(1)<0
    GoalPolar(1)=GoalPolar(1)+2*pi;
end

% define tolerance to Goal zone
Length_epsilon = 12;
Angle_epsilon = 2*pi/500;                 %% How far of can the final position be

% Define desired exit angle
GoalPolar(1)= pi/10*13;
x0 = [StartPolar(2) StartVelocities(1) StartPolar(1) StartVelocities(2)];

MinDist = 6; %% Minimal Distance to Earth

% initial guess has to be nonzero because of singularity
x_initial = repmat(x0,NumInc,1)';
 
x_box = [-Inf Inf; -Inf, Inf; -Inf Inf; -Inf, Inf];
% Limit rocket output 
u_box = [-5 5 ; -5 5];
% Create OCP from ODE 
[ocp, x, u,varout{1:6}] = ode2ocp(KeplerDGL, 4, 2, NumInc, dt, x_box=x_box, u_box=u_box,x0='param',foh=true,verbose=1)

ocp.solver('ipopt')
ocp.set_value(u(1:2,1),[0 0]); %% first value for u has to be zero. 

x0_p = varout{2};
xf_p = varout{3};
% set x0
ocp.set_value(x0_p,x0);

%Deviation from Stable Orbit
DevStOr = (sqrt(Gamma*m(1)./(x(1,:).^3))-x(4,:)).^2;
%Minimizer
costfun = (u(1,:).^2+u(2,:).^2);
ocp.minimize(sum(costfun));

%Box constraints for final position
ocp.subject_to (x(1,end)>= GoalPolar(2)-Length_epsilon);
ocp.subject_to (x(1,end)<= GoalPolar(2)+Length_epsilon);
% % % % %  
ocp.subject_to (x(3,end)>= GoalPolar(1)-Angle_epsilon);
ocp.subject_to (x(3,end)<= GoalPolar(1)+Angle_epsilon);
% Minimal Distance to Earth
ocp.subject_to(x(1,:) >= MinDist);
ocp.subject_to(x(1,:) <= 200);
% I want the spacecraft to leave at a certain angle. 
% So i set the angular velocity in the final increment to 0
ocp.subject_to(x(4,end) >=0.0);
ocp.subject_to(x(4,end) <=0.01);

ocp.set_initial(x(:,2:end),x_initial(:,:));




sol = ocp.solve()

%[DISREGARD] This is just for plotting. It is very inelegant
% so probably best to leave alone. 


figure(i)
% xxx = ocp.debug.value(x);

xxx = ocp.value(x);
uuu = ocp.value(u);
Values.(varmap) =xxx;
Values.(varmap1)=uuu;
Values.(varmap3)=Coordinates;
Values.(varmap2)=time;
% uuu = ocp.debug.valmaue(u);
Coordinates= zeros(2,length(xxx));
Velocity= zeros(1,length(xxx));
for j = 1:length(xxx)
        [xx,yy]= pol2cart((xxx(3,j)),(xxx(1,j)));
        Coordinates(1,j)=xx;
        Coordinates(2,j)=yy;
        
        [vx,vy] = pol2cart((xxx(4,j)),(xxx(2,j)));
        Velocity(j) = sqrt(vx^2+vy^2);


   

end
Angle = [GoalPolar(1)-Angle_epsilon:0.001:GoalPolar+Angle_epsilon];
[arc11 arc12] = arrayfun(@(x) pol2cart(x,GoalPolar(2)-2*Length_epsilon),Angle);
[arc21 arc22] = arrayfun(@(x) pol2cart(x,GoalPolar(2)+2*Length_epsilon),Angle);

AngleFull = [0:0.001:pi*2];
[Atmo1 Atmo2] = arrayfun(@(x) pol2cart(x,MinDist),AngleFull);


scatter(0,0,100,'blue','filled')
hold on 
scatter(arc11,arc12,10,'r');
scatter(arc21,arc22,10,'r');
scatter(Goal(1),Goal(2),50,'+')
plot([arc11(end) arc21(end)],[arc12(end) arc22(end)]','r',LineWidth=5)
plot([arc11(1) arc21(1)],[arc12(1) arc22(1)]','r',LineWidth=5)
scatter(Atmo1,Atmo2,10,'black')
scatter(Start(1),Start(2),500,Velocity(1),'d','filled')
axis equal
scatter(Coordinates(1,end),Coordinates(2,end),500,Velocity(end),'^','filled')
scatter(Coordinates(1,:),Coordinates(2,:),80,Velocity,'filled')
 colormap(gca,'turbo')
 cb = colorbar; % creates the colorbar on side
 cb.Label.String = 'Velocity [m/s]';
% use caxis to change range on colorbar
caxis([min(Velocity),max(Velocity)]); % low end is 10, high end is 30
xlabel('x')
ylabel('y')
set(gcf,'color','w');


    figure(444)
    clf
    subplot(2,1,1)
    plot(time,ocp.value(costfun))
    ylim([-max(ocp.value(costfun))/4 max(ocp.value(costfun))])
    xlabel('time')
    ylabel('$l(x,u)$',Interpreter='latex')
    subplot(2,1,2)
    plot(time,ocp.value(DevStOr))
    xlabel('time')
    ylabel('Distance from stable orbit',Interpreter='latex')
    ylim([-max(ocp.value(DevStOr))/4 max(ocp.value(DevStOr))])
    
    
end
