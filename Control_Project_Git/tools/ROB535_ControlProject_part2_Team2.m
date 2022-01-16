function [sol_2, FLAG_terminate] = ROB535_ControlProject_part2_Team2(TestTrack, Xobs_seen, curr_state)

curr_state = curr_state';

%% parameters
Nw=2;
f=0.01;
Iz=2667;
a=1.35;
b=1.45;
By=0.27;
Cy=1.2;
Dy=0.7;
Ey=-1.6;
Shy=0;
Svy=0;
m=1400;
g=9.806;

track = TestTrack;
l1 = size(track.cline,2);
x_orig = linspace(1,10,l1);
x_interp = linspace(1,10,l1*20); %previously X20
r1 = interp1(x_orig,track.cline(1,:),x_interp);
r2 = interp1(x_orig,track.cline(2,:),x_interp);
track.cline = [r1;r2];
r1 = interp1(x_orig,track.bl(1,:),x_interp);
r2 = interp1(x_orig,track.bl(2,:),x_interp);
track.bl = [r1;r2];
r1 = interp1(x_orig,track.br(1,:),x_interp);
r2 = interp1(x_orig,track.br(2,:),x_interp);
track.br = [r1;r2];

%number of states and inputs in dynamic model
nstates=6;
ninputs=2;

%input ranges (first row is wheel angle, second row forward force)
input_range=[-0.499,   0.499;...
            -4999,4999];

%load('project_traj_01_timestepfast_centerline.mat'); %loads Eric's fast trajectory
%load('project_traj_01_timestepfast.mat'); %loads Eric's fast trajectory
%load('project_traj_01_timestepslow.mat'); %loads Nikhil's slow trajectory
%load('project_traj_01_timestepfast_centerline2.mat'); %loads Meet's stable trajectory

persistent U_ref Y_ref

if isempty(U_ref)
    %load('ROB535_ControlProject_part1_Team2.mat', 'U', 'Y');
    %load('project_traj_01_timestep_centerline_ND.mat', 'U', 'Y');
    load('refs/RefTraj.mat', 'U', 'Y');
    U_ref = U'; %already 0.01 timestep
    Y_ref = Y'; %already 0.01 timestep
end

dt=0.01;  % time discretization
persistent T
if isempty(T)
    len1 = size(Y_ref,2);
    T = 0:dt:dt*(len1-1); %time span
end

persistent dfdx dfdu

if isempty(dfdx)
    %Setting up the dynamics and partials
    syms x u y v phi r Fx delta_f
    % x u y v phi r
    %slip angle functions in degrees
    a_f=rad2deg(delta_f-atan2(v+a*r,u));
    a_r=rad2deg(-atan2((v-b*r),u));

    %Nonlinear Tire Dynamics
    phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
    phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

    F_zf=b/(a+b)*m*g;
    F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

    F_zr=a/(a+b)*m*g;
    F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

    % F_total=sqrt((Nw*Fx)^2+(F_yr^2));
    % F_max=0.7*m*g;

    % if F_total>F_max
    %     
    %     Fx=F_max/F_total*Fx;
    %   
    %     F_yr=F_max/F_total*F_yr;
    % end

    %vehicle dynamics
    % x u y v phi r
     f1 = u*cos(phi)-v*sin(phi);
     f2 = (-f*m*g+Nw*Fx-F_yf*sin(delta_f))/m+v*r;
     f3 = u*sin(phi)+v*cos(phi);
     f4 = (F_yf*cos(delta_f)+F_yr)/m-u*r;
     f5 = r;
     f6 = (F_yf*a*cos(delta_f)-F_yr*b)/Iz;

    dfdx = [diff(f1,x)   diff(f1,u)   diff(f1,y)   diff(f1,v)   diff(f1,phi)   diff(f1,r)
            diff(f2,x)   diff(f2,u)   diff(f2,y)   diff(f2,v)   diff(f2,phi)   diff(f2,r)
            diff(f3,x)   diff(f3,u)   diff(f3,y)   diff(f3,v)   diff(f3,phi)   diff(f3,r)
            diff(f4,x)   diff(f4,u)   diff(f4,y)   diff(f4,v)   diff(f4,phi)   diff(f4,r)
            diff(f5,x)   diff(f5,u)   diff(f5,y)   diff(f5,v)   diff(f5,phi)   diff(f5,r)
            diff(f6,x)   diff(f6,u)   diff(f6,y)   diff(f6,v)   diff(f6,phi)   diff(f6,r)];

    dfdu = [diff(f1,delta_f) diff(f1,Fx)
            diff(f2,delta_f) diff(f2,Fx)
            diff(f3,delta_f) diff(f3,Fx)
            diff(f4,delta_f) diff(f4,Fx)
            diff(f5,delta_f) diff(f5,Fx)
            diff(f6,delta_f) diff(f6,Fx)];


    %these are the system linearized in discrete time about the reference
    %trajectory i.e. x(i+1)=A_i*x_i+B_i*u_i
    dfdx = matlabFunction(dfdx);
    dfdu = matlabFunction(dfdu);

end

A=@(i) eye(nstates) + dt*dfdx(U_ref(1,i),Y_ref(5,i),Y_ref(6,i),Y_ref(2,i),Y_ref(4,i));
B=@(i) dt*dfdu(U_ref(1,i),Y_ref(6,i),Y_ref(2,i),Y_ref(4,i));


%11 timesteps for 3 states, 10 timesteps for 2 inputs
npred=20;
Ndec=(npred+1)*nstates+ninputs*npred;
%decision variable will be z=[x_1...x_11;u_1...u_10] (x refers to state
%vector, u refers to input vector)

batch = 0.5/dt+1;

%final trajectory
Y=NaN(nstates,batch);

%applied inputs
U=NaN(ninputs,batch);

%input from QP
u_mpc=NaN(ninputs,batch);

%error in states (actual-reference)
eY=NaN(nstates,batch);

%set initial condition
Y(:,1)=curr_state;

%find start point
diffs2 = (Y_ref - curr_state).^2;
sum_diffs2 = diffs2(1,:) + diffs2(3,:);
[~, t_index] = min(sum_diffs2);

FLAG_terminate = 0;

%for i=t_index:min(t_index+batch-1, length(T)-1)
for local_i=1:batch
    diffs2 = (Y_ref(:,t_index:min(t_index+batch*2, size(Y_ref, 2))) - Y(:,local_i)).^2;
    sum_diffs2 = diffs2(1,:) + diffs2(3,:);
    [~, j] = min(sum_diffs2);
    %i = t_index + j - 1;
    i = local_i + t_index - 1;
    
    %shorten prediction horizon if we are at the end of trajectory
    npred_i=min([npred,length(T)-i-1]);
    
    %calculate error
    eY(:,local_i)=Y(:,local_i)-Y_ref(:,i);

    %generate equality constraints
    [Aeq,beq]=eq_cons(i,A,B,eY(:,local_i),npred_i,nstates,ninputs);
    
    %generate inequality constraints
    car = [Y(1,local_i); Y(3,local_i)];
    [A_con,b_con] = find_constraints(car,track,Xobs_seen);
    
    %generate boundary constraints
    [Lb,Ub]=bound_cons(i,U_ref,npred_i,input_range,nstates,ninputs);
    
    %cost for states [x, longi vel, y, lat vel, heading, rot vel] 
    Q=[8, 0.4, 8, 0.4, 160, 2];
    
    %cost for inputs [steering, throttle]
    R=[150, 0];
    
    H=diag([repmat(Q,[1,npred_i+1]),repmat(R,[1,npred_i])]);
    
    func=zeros(nstates*(npred_i+1)+ninputs*npred_i,1);
    
    while 1
        [Aineq,bineq]=ineq_cons(A_con,b_con, npred_i,nstates,ninputs);
        horizon_ref = [reshape(Y_ref(:,i:i+npred_i),1,[]), zeros(1, npred_i*ninputs)];
        bineq = bineq - Aineq*horizon_ref';
        [x,fval,exitflag,output] = quadprog(H,func,Aineq,bineq,Aeq,beq,Lb,Ub);
        if exitflag == -2
            n = size(A_con, 1);
            disp(['reducing constraints to ' num2str(n-1)]);
            A_con = A_con(1:n-1,:);
            b_con = b_con(1:n-1,:);
            if n <= 2
                FLAG_terminate = 1;
                break;
            end
        else
            break;
        end
    end
    
    
    if FLAG_terminate == 1
        break;
    end
    
    if size(x) < nstates+ninputs
        break
    end
    
    %get linearized input
    u_mpc(:,local_i)=x(nstates*(npred_i+1)+1:nstates*(npred_i+1)+ninputs);
    
    %get input
    U(:,local_i)=u_mpc(:,local_i)+U_ref(:,i);
    
    
    %simulate model
    [~,ztemp]=ode45(@(t,z)kinematic_bike_dynamics(t,z,U(:,local_i),0,Nw,f,Iz,a,b,By,Cy,Dy,Ey,Shy,Svy,m,g),[0 dt],Y(:,local_i));
    
    %store final state
    Y(:,local_i+1)=ztemp(end,:)';
    disp(['Done loop ' num2str(i) ' out of ' num2str(length(T)-1)])
end

sol_2 = U';
t_index = t_index + batch;
if t_index >= length(T)
    FLAG_terminate = 1;
end

if FLAG_terminate == 1
    clear ROB535_ControlProject_part2_Team2;
end


end












