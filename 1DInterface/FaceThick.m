% compute the interface thickness and the interface energy
clear all; clc;

nx = 1000;	% x range
dx = 0.005;	% grid spacing of x
L = 0.01; 	% mobility
k = 20; 	% interface energy coefficient
u = 200; 	% chemical energy coefficient

% max dt to ensure numerical stablity
dt = 1 / (2*L) / (2*u + k/dx^2) ;

% total time steps (should be large enough to reach equilibrium)
total_step = 50000;

phi = zeros(nx,1);
grad_phi = zeros(nx,1);
lap_phi = zeros(nx,1);
x = zeros(nx,1);

%==============initial value of phi =======================
for i = 1:nx
    x(i) = i * dx;
    phi(i) = 0.5;
end
%==============setting the initial disturbance=============
f = 0.02;
phi(1) = phi(1)*(1+f);  phi(nx) = phi(nx)*(1-f);
phi(2) = phi(2)*(1+f);  phi(nx-1) = phi(nx-1)*(1-f);

%===============gradient&Laplace===========================
for times = 1:total_step
    
    for i=1:nx
        %Numan boundary condition
        ip = i+1; im = i-1;
        if im == 0
            im = 1;
        elseif ip == nx+1
            ip = nx;
        end
        
        %gradient
        grad_phi(i) = (phi(ip) - phi(im))/(2*dx);
        
        %laplacian
        lap_phi(i) = ( phi(ip)-2*phi(i)+phi(im) )/(dx^2);
    end
    
    %=============evolution================================
    for i=1:nx
        %time evolution
        phii=phi(i);
        term1 = u*(1-2*phii);   % driving force of chemical energy
        term2 = -k*lap_phi(i);  % driving force of gradient energy
        term = term1 + term2;
        phi(i) = phii - L*term*dt;% update phi
        
        % set the range for phi
        if phi(i) > 1
            phi(i) = 1;
        elseif phi(i) < 0
            phi(i) = 0;
        end
    end
    
    
    %visualization of the output disp(phi);
    if mod(times,300)==0
        figure(1);
        plot(x, phi, 'lineWidth', 2);
        axis([1*dx nx*dx -0.05 1.05]);
        
        if mod(times,500)==0
            times;
        end
    end
    
end

% compare numerical & analytical value of interface thicknss
am1 = 0;
for i=1:nx
    if (0.0001 < phi(i)) && (phi(i) < 0.9999)
        am1=am1+1;
    end
end
fprintf('the numerical face thickness is %d \n', am1 * dx);
fprintf('the analitical face thickness is %d \n\n', pi * sqrt(k / (2*u)));


ener = 0; % energy
for i = 1:nx
    term = u*phi(i)*(1-phi(i)) + k/2*(grad_phi(i))^2;
    ener = ener + term * dx;
end
fprintf('the numerical face energy is %d \n', ener);
fprintf('the analitical face energy is %d \n', pi / (4*sqrt(2)) * sqrt(k*u));
