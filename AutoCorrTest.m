%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Stochastic Differential Equation Model  %%%%%%%%
%%%%%%%%         Author: Ryan Houlihan            %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
format long;

%%%%%%%% Values that need to be set %%%%%%%%%%%%%%%%%%%%%%
N = 10000;             % number of springs
N_runs = 1;         % Number of runs to estimate auto-corr over
h = 0.0005;           % time step
T = 10;               % maximum time

alpha = 0.8;          % value for alpha
zeta = 35;

m=1; x=2; v=3; s_len=4; k=5; sigma=0.0; % index's

objects=zeros(N,5);   % prepare place to store locations:
                      % 1 big and (Num_objts-1) small

objects(1,m)=10;    % Large Mass
objects(1,k)=0;       % Large spring constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NS= T/h;             % number of steps to take
t=(0:h:T-h);          % t is the vector [0 1h 2h 3h ... Nh]

del_omega = zeta/N;      % Delta omega value (time/number of springs)
omega = zeros(1,N);   % frequency array
k_all = zeros(1,N);   % spring constant array
m_all = zeros(1,N);   % little mas array
R = zeros(2,NS);
R_exact=zeros(2,NS);

acc=zeros(N_runs,NS);
corr=zeros(NS,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Find Springs %%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N
    % compute all the small objects frequencies
    omega(i)= i*del_omega;
end
for i=1:N
    % compute all the small objects k's and mass's
    k_all(i) = (2*alpha*del_omega)/((alpha^2 + omega(i)^2)*pi);
    m_all(i) = k_all(i)/(omega(i)^2);
end

% Compute the exact memory kernel
% and the estimated memory kernel
for i=0:NS
    R(1,i+1)=alpha*i*h;
    R_exact(1,i+1)=alpha*i*h;
    for j=1:N
        R(2,i+1)=R(2,i+1)+k_all(j)*cos(omega(j)*i*h);
    end
    R_exact(2,i+1) = exp(-alpha*i*h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%

for r=1:N_runs
    % set all the small objects to their respective values
    % also set the inital position and velocity of each 
    % of the small objects

% initalize small object    
% find xi and eta for the simulation
xi = randn(N-1,1);
eta = randn(N-1,1);

    for j=2:N
        objects(j,m)=m_all(j-1);
        objects(j,s_len)=0;
        objects(j,k)=k_all(j-1);
        objects(j,x)=sqrt(1/objects(j,k))*xi(j-1,1);
        objects(j,v)=1/sqrt(objects(j,m))*eta(j-1,1);
    end
    
    for i=1:NS+1
        % compute all the forces on all the small objects
        % and update their position and velocities
        for j=2:N
            big_x=objects(1,x);
            lil_x=objects(j,x);
            lil_m=objects(j,m);
            lil_v=objects(j,v);
            
            % calculate movement and change in velocity of little mass
            rel_dist=lil_x-big_x;
            % calculte spring force
            f_spring=objects(j,k)*(objects(j,s_len)-rel_dist);
            
            % Compute velocity
            objects(j,v)=lil_v+h*f_spring/lil_m;
            
            % Compute position
            objects(j,x)=lil_x+lil_v*h+((h^2)/2)*f_spring/lil_m;
        end
        
        % sum up the the n+1 forces of each small mass on the main mass
        f_all=0;
        for j=2:(N)
            big_x=objects(1,x);
            big_v=objects(1,v);
            lil_x=objects(j,x);
            lil_v=objects(j,v);
            lil_m=objects(j,m);
            
            % calculate distance from large mass to small
            rel_dist=lil_x-big_x;
            
            % compute the addition of the small mass to ttl force
            f_all=f_all+objects(j,k)*rel_dist;
        end
        
        % Add contribution of big mass spring to ttl force
        f_all=((f_all)-objects(1,k)*big_x)/objects(1,m);
        
        % save acceleration
        acc(r,i)=f_all;
        
        % update the marked object's position and velocity
        objects(1,v)=big_v+h*f_all;
        objects(1,x)=big_x+big_v*h+((h^2)/2)*f_all;
    end
end

% Compute force auto-correlation
for i=1:NS
    % mean of F(0)
    m1=mean(acc(:,1));
    % mean of F(t)
    m_curr=mean(acc(:,i));
    for r=1:N_runs
        corr(i)=corr(i)+(acc(r,1)-m1)*(acc(r,i)-m_curr)/((acc(r,1)-m1)^2);
    end
    corr(i)=corr(i)/(N_runs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
% Plot R(t) estimate
plot(R(1,:), R(2,:), 'Color', [.1 .4 0], 'Linewidth', 2);
hold on;
% Plot R(t) exact
plot(R_exact(1,:), R_exact(2,:), '--', 'Color', [.5 .3 1], 'Linewidth', 2);
hold on;
% Plot force autocorrelation
plot((alpha*t),corr, 'Color', [1 0 .5], 'Linewidth', 2);

% Create xlabel, title, legend
xlabe = xlabel('$$\alpha t$$');
set(xlabe, 'Interpreter', 'latex', 'fontsize', 14);
title(sprintf('alpha = %0.2f, # Springs = %d, Time = %d, TS = %d, # Runs = %d, Zeta = %d, M = %d', ...
							 alpha, N, T, h, N_runs, zeta, objects(1,m)));
l=legend('$$R(t) = \sum_{j=1}^N k_j cos(\omega_j t)$$', '$$R(t) = e^{{-} \alpha |t|}$$', ...
							           '$$M^2<\ddot Q(0)\ddot Q(t)>$$',3);
set(l,'Interpreter', 'latex', 'fontsize', 14, 'Location', 'NorthEast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
