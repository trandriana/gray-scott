% Gray-Scott Model using Spectral Fourier Method
clc; clear; close all;

% Parameters
Nx = 256;            % Number of grid points in x
Ny = Nx;            % Number of grid points in y
T = 800000;            % Number of time steps
N = 10;
dt = T/T/N ;       % Time step size

% Reaction-Diffusion Parameters
Du = 2.e-5;           % Diffusion coefficient for u
Dv = Du/2;           % Diffusion coefficient for v
F = 0.078;           % Feed rate
k = 0.061;           % Kill rate

a1 =  -F; a2 = 0; g1 = F;
b1 = 0; b2 = -(F + k); g2 = 0;

% Wavenumbers for Fourier transform (2D grid)
A = 2 ;
x = linspace(-A, A-2*A/Nx, Nx)'; dx = x(2) - x(1);
mu = 1i*pi/A*(-Nx/2:Nx/2-1)';
[Kx, Ky] = ndgrid(mu, mu);
lambda = (Kx.^2)+ (Ky.^2);
lambda = fftshift(lambda);

% Initialize u and v
u = ones(Nx, Ny);  
v = zeros(Nx, Ny);  

% Random perturbations in the center of the grid
r = round(Nx/20);
% u(mid-10:mid+10, mid-10:mid+10) = 0.30; %+ 0.02*rand(21, 21);
% v(mid-10:mid+10, mid-10:mid+10) = 0.15;% + 0.02*rand(21, 21);
inside_the_small_box = Nx/2-r:Nx/2+r;
nodes_in_the_small_box = length(inside_the_small_box);
rng("default")
u(inside_the_small_box, inside_the_small_box) = 0.5  + rand(nodes_in_the_small_box, nodes_in_the_small_box);
v(inside_the_small_box, inside_the_small_box) = 0.25 + rand(nodes_in_the_small_box, nodes_in_the_small_box);
% Precompute the Fourier transforms of initial u and v

% clear mov;
movieName = 'Worms'

mov = VideoWriter(['movies\_',movieName]);
open(mov);
%
my_colormap = jet;
pcolor(x, x, real(u));
set(gca, 'clim', [0 1]);
my_ticks = linspace(x(1), x(end), 5)';
% xticks(my_ticks); yticks(my_ticks);
axis off
xlim([x(1), x(end)])
ylim([x(1), x(end)])
colormap(my_colormap);
colorbar;
shading("interp")
title(['Number of Modes = ', '(', int2str(Nx), ', ', int2str(Ny), '), ', 't = ', num2str(0)]);

input('Type Enter to continue!')
im = getframe(gcf);
writeVideo(mov, frame2im(im));
drawnow;
%
% print ('-f1', '-r600', '-dpdf',strcat(['GrayScott_sample' '_IC']))

disp('Computations ongoing.')
% tau = zeros(1,T);

% Time stepping loop
for t = 1:T
    % disp([int2str(t), '/', int2str(T)])
    % Reaction terms
    uv = u.*v;
    uv2 = uv.*v;

    yinf = norm(uv(:),'inf');
    tau = min(dt, dt/(2*yinf));
    
    u = u + tau * (b1*v - uv2 + g1);
    v = v + tau * (a2*u + uv2 + g2);

    u_hat = squeeze(fftn(u));
    v_hat = squeeze(fftn(v));

    u_hat = u_hat.*exp(tau * (Du*lambda + a1));
    v_hat = v_hat.*exp(tau * (Dv*lambda + b2));
    
    u = ifftn(u_hat);
    v = ifftn(v_hat);
    
    % Display the result every 100 time steps
    if mod(t, 500 ) == 0
        pcolor(x, x, real(u)); 
        set(gca, 'clim', [0 1]);
        % xticks(my_ticks); yticks(my_ticks);
        axis off

        xlim([x(1), x(end)])
        ylim([x(1), x(end)])
        colormap(my_colormap);
        colorbar;
        shading("interp")
        title(['Number of Modes = ', '(', int2str(Nx), ', ', int2str(Ny), '), ', 't = ', num2str(t*dt)]);
        im = getframe(gcf);
        writeVideo(mov, frame2im(im));
        drawnow;
    end
end

close(mov);
print ('-f1', '-r600', '-djpeg', strcat([movieName,'_']))