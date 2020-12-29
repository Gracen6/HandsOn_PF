% generate the polycrystal using Steinbach's multi-phase field model
clear all; clc;

%=========== parameters for geometric model================
am = 0;
Q = 20;     % number of different phases
nx = 120;   % number of grid in x-direction
ny = 120;   % number of grid in y-direction
dx = 0.1;
dy = 0.1;
dt = 0.01;  % time step

%============ material parameters==========================
L = 0.02;   % mobility
k = 4.0;    % gradient coefficient
u = 64.0;   % energy barrier, chemical coefficient
u_p = 3*u;  % coeff. for penalty of 3 phases coexistence

times = 0;
phi = zeros(Q,nx,ny);  % order parameter

AMphi = 0;   % total amount of different phases
AMdphi = 0;  % total amount of the increments of different phases
AMva = 0;    % total amount of variational derivative
varia = zeros(Q,1); % variational derivative for each phase

grad_phi = zeros(Q,nx,ny,2);

lap_phi = zeros(Q,nx,ny);

x = zeros(nx);
gap = zeros(nx,ny);      % show the grain boundary
polyCrys = zeros(nx,ny); % show different grains as different colors
phase = zeros(Q,1);

% output the system's total energy versus time
file_eng = fopen('TotalEnergy.txt', 'w');
fprintf(file_eng,'time energy \n');
eng = zeros(nx, ny);     % energy of each element

%==========initialize parameters for phase field===========
for i = 1:nx
    for j = 1:ny
        for p = 1:Q
            phi(p,i,j) = 1 / Q;
        end
    end
    x(i)=i;
end
%===========set the initial perturbatiom===================
f = 0.01;   %relative proportion for perturbation
for i=1:nx
    for j=1:ny
        for p=1:Q
            phi(p,i,j) = phi(p,i,j) - f + 2*f*rand(1,1);
        end
        % total phases
        am = 0;
        for p = 1:Q
            am = am + phi(p,i,j);
        end
        % normalization
        for p = 1:Q
            phi(p,i,j) = phi(p,i,j)/am;
        end
        % verify
        AMphi = 0;
        for p = 1:Q
            AMphi = AMphi + phi(p,i,j);
        end
        if abs(AMphi-1) > 0.0001
            AMphi
            error('error, total phi does not equal 1');
        end
    end
end

%==========gradient&Laplace================================
for times = 1:12000 %20000
    
    for i = 1:nx
        for j = 1:ny
            % Periodic boundary condition
            ip = i+1;  im = i-1;  jp = j+1;  jm = j-1;
            if (im == 0)
                im = nx;
            elseif (ip == nx+1)
                ip = 1;
            end
            if jm == 0
                jm = ny;
            elseif jp == ny+1
                jp = 1;
            end
            
            %gradient
            for p=1:Q
                grad_phi(p,i,j,1) = (phi(p,ip,j) - phi(p,im,j))/(2*dx);
                grad_phi(p,i,j,2) = (phi(p,i,jp) - phi(p,i,jm))/(2*dy);
            end
            
            %laplacian
            for p = 1:Q
                lap_phi(p,i,j) = ( phi(p,ip,j)-2*phi(p,i,j)...
                    +phi(p,im,j) )/(dx^2) + ( phi(p,i,jp)...
                    -2*phi(p,i,j)+phi(p,i,jm) )/(dy^2);
            end
            
        end
    end
    
    %========evolution=====================================
    for i = 1:nx
        for j = 1:ny
            
            AMphi = 0;
            for p = 1:Q
                AMphi = AMphi + abs( phi(p,i,j) );
            end
            
            % compute variational derivative
            AMva = 0;
            for p = 1:Q
                
                term2 = 0; % the penalty term
                for m = 1:Q
                    if m ~= p
                        for n = (m+1):Q
                            if n ~= p
                                term2 = term2 + abs(phi(m,i,j)) * abs(phi(n,i,j));
                            end
                        end
                    end
                end
                term2 = term2 * sign(phi(p,i,j));
                term2 = term2 * u_p;
                phii = phi(p,i,j);
                
                % variational derivative of chemical potential, 
                varia(p) = u * sign(phii) * (AMphi-abs(phii)) ...
                    + term2 - k * lap_phi(p,i,j); % variation derivative of penalty term, & gradient energy, 
                AMva = AMva + varia(p);                
            end            
            % update order parameter
            for p = 1:Q
                term = varia(p) - AMva / Q;
                phi(p,i,j) = phi(p,i,j) - L*term*dt;
            end            
            % total phases
            am = 0;
            for p = 1:Q
                am = am + phi(p,i,j);
            end
            % normalization
            for p = 1:Q
                phi(p,i,j) = phi(p,i,j) / am;
            end
            
            % verify
            AMphi = 0;
            for p = 1 : Q
                AMphi = AMphi + phi(p,i,j);
            end
            if abs(AMphi-1) > 0.001
                AMphi
                error('error, total phi does not equal 1');
            end
            
            
            % =================== the energy ==============
            eng(i, j) = 0;
            for p1 = 1:Q
                for p2 = (p1 + 1):Q
                    eng(i, j) = eng(i, j) + u * abs(phi(p1, i, j) ...
                        * phi(p2, i, j)); % chemical potential term
                end
            end
            for p1 = 1:Q
                for p2 = (p1 + 1):Q
                    for p3 = (p2 + 1):Q
                        eng(i, j) = eng(i, j) + u_p * abs(phi(p1, i, j)...
                            * phi(p2, i, j) * phi(p3, i, j)); 
							% penalty term of 3 phases coexistence
                    end
                end
            end
            for p = 1:Q
                eng(i, j) = eng(i, j) + k / 2 * ((grad_phi(p,i,j,1))^2 ...
                    + (grad_phi(p,i,j,2))^2); % gradient energy term
            end
            
            % ============= grain boundary ================
            gap(i, j) = max(phi(:,i,j)) - min(phi(:,i,j));
            
            % === different grains in different colors ====
            [aa, bb] = max(phi(:,i,j));
            polyCrys(i, j) = bb; % show different grains
            
        end
    end
    
    
    %visualize output disp(phi);
    if mod(times,5) == 0
        times
        % output energy
        fprintf(file_eng, '%d    %d \n', times, sum(sum(eng)));
        % ==========plot phase 1 along x at y = 6==========
        figure(1);  plot(x,phi(1,:,ny/2)); axis([1 nx -0.01 1.01]);
        % ==========plot the grain boundary================
        figure(2);
        colormap('jet(64)'); pcolor(gap);shading flat;
        if mod(times,8) == 0  % save figures
            picname = ['GrainBoundary_' num2str(times) '.png'];
            saveas(gcf, picname);
        end
        % ==========plot different grains(phases)==========
        figure(3);
        colormap('jet(64)'); pcolor(polyCrys);shading flat;
        if mod(times,8) == 0  % save figures
            picname = ['PolyCrys_' num2str(times) '.png'];
            saveas(gcf, picname);
        end
        
    end
    
end

fclose(file_eng); % close the energy file
