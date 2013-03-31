%-------------------------------------------------%
% Solve 3D incompresible N-S equations            %
% in vorticity_y - laplacian v:                   %
%                                                 %
% d psi/dt = h_psi + 1/Re Lap(psi)                %
% d o_y/dt = h_oy  + 1/Re Lap(o_y)                %
% div(v)   = 0                                    %
% psi      = Lap(v)                               %
% oy       = du/dz - dw/dx                        %
%                                                 %
% Spatial discretization: Fourier                 %
% Temporal discretization: 3 order Runge-kutta    %
% Low storage: van der Houwen (2R)                %
%                                                 %
% ()_hat = Fourier coefficient                    %
%                                                 %
% adrian Dic-2012                                 %
%-------------------------------------------------%

function vortices3D

clear all

global poisson poisson2D kx ky kz Lap dealias Re u v w u_hat v_hat w_hat ndim

%-----------domain-----------%
nsave = 10;
nmax  = 2;  % max steps
Re    = 10000;  % viscosity
T     = 10000;  % total time
CFL   = 0.5;    % CFL
Lx    = 1*pi;   % x length
Ly    = 1*pi;   % y length
Lz    = 1*pi;   % z length
nx    = 128;    % x-modes
ny    = 128;    % y-modes
nz    = 128;    % z-modes
%----------------------------%
dx    = Lx/(nx-1);
dy    = Ly/(ny-1);
dz    = Lz/(nz-1);
ndim  = [nx ny nz];

%---------display------------%
disp(['Reynolds number: '     ,num2str(Re)])
disp(['x length Lx: '         ,num2str(Lx)])
disp(['y length Ly: '         ,num2str(Ly)])
disp(['y length Lz: '         ,num2str(Lz)])
disp(['CFL number: '          ,num2str(CFL)])
disp(['number of nx modes: '  ,int2str(nx)])
disp(['number of ny modes: '  ,int2str(ny)])
disp(['number of nz modes: '  ,int2str(nz)])
disp(['total time simulated: ',num2str(T)])

%-----------mesh-------------%
%xx      = linspace(-Lx/2,Lx/2,nx);
%yy      = linspace(-Ly/2,Ly/2,ny);
%zz      = linspace(-Lz/2,Lz/2,nz);
%[x y z] = meshgrid(xx,yy,zz); clear xx yy zz

%--initialize low storage RK--%
b   = [0 1/6 1/3 1/3 1/6];
a   = [0 1/2 1/2 1];
dtv = CFL*min([dx dy dz])^2*Re/pi^2;

%-----initialize Fourier------%
[kx ky kz] = meshgrid( mod((1:nx)-ceil(nx/2+1),nx)-floor(nx/2) , ...
                       mod((1:ny)-ceil(ny/2+1),ny)-floor(ny/2) , ...
                       mod((1:nz)-ceil(nz/2+1),nz)-floor(nz/2) );

% Cutting frequencies 2/3 rule: from n/3 to n/2-1
dealias = ( abs(kx)<2/3*(nx/2) & abs(ky)<2/3*(ny/2) & abs(kz)<2/3*(nz/2) );  

kx   =   2*pi*kx/Lx;
ky   =   2*pi*ky/Ly;
kz   =   2*pi*kz/Lz;
Lap  = -(kx.^2 + ky.^2 + kz.^2);

% singular point psi 0-mode undetermined
poisson          = Lap; 
poisson(1,1,1)   = 1; 
kzs              = kz; 
kzs(kzs==0)      = 1; 
poisson2D        = kx.^2+kz.^2;
poisson2D(poisson2D==0) = 1;

max(max(max(poisson2D)))

%------initial conditions----%
u_hat   = (rand(nx,ny,nz)-0.5)+1i*(rand(nx,ny,nz)-0.5); u_hat(-Lap> (nx/30)^2 ) = 0;
v_hat   = (rand(nx,ny,nz)-0.5)+1i*(rand(nx,ny,nz)-0.5); v_hat(-Lap> (nx/30)^2 ) = 0;
u       = 5e4*real(ifftn(u_hat,ndim));
v       = 5e4*real(ifftn(v_hat,ndim));
u_hat   = fftn(u,ndim);
v_hat   = fftn(v,ndim); 
w_hat   = -( kx.*u_hat + ky.*v_hat )./kzs; clear kzs 
w       = real( ifftn(w_hat,ndim) );

oy_hat  = 1i*kz.*u_hat - 1i*kx.*w_hat;
psi_hat = Lap.*v_hat;


disp('!------------------------------!')
disp('      starting computation      ')
disp('!------------------------------!')

%------temporal loop---------%
t      = 0; 
ii     = 0; 
j      = 0;
S1_oy  = oy_hat; 
S1_psi = psi_hat; 

%------%
% ii      = 121;
% load(['data/vortices3D.RK.',int2str(ii),'.mat'])
% disp(['data/vortices3D.RK.',int2str(ii),'.mat'])
% u_hat   = fftn(u,ndim);
% v_hat   = fftn(v,ndim); 
% w_hat   = fftn(w,ndim); 
% oy_hat  = 1i*kz.*u_hat - 1i*kx.*w_hat;
% psi_hat = Lap.*v_hat;
% j       = 10*ii;
% S1_oy   = oy_hat; 
% S1_psi  = psi_hat; 
%------%

tic
while t<T 

    j = j+1; disp(j)
    if j>nmax,break,end    

    % compute time step
    umax = max( abs(u(:)) );
    vmax = max( abs(v(:)) );
    wmax = max( abs(w(:)) );
    dtx  = CFL*dx/umax/pi;
    dty  = CFL*dy/vmax/pi; 
    dtz  = CFL*dz/wmax/pi;
    dt   = min( [dtv dtx dty dtz] ); 

    max(max(max(oy_hat)))
    % RK substeps
    for i=1:4
        S1_oy   =  oy_hat + (a(i)-b(i))*dt*S1_oy;
        S1_psi  = psi_hat + (a(i)-b(i))*dt*S1_psi;
        [S1_psi S1_oy] = Fw(S1_psi,S1_oy);
        oy_hat  = oy_hat  + b(i+1)*dt*S1_oy;
        psi_hat = psi_hat + b(i+1)*dt*S1_psi;
        max(max(max(oy_hat)));dt
    end
    
    
    
    % Explicit Euler
    %[S1_psi S1_oy] = Fw(psi_hat,oy_hat);
    %psi_hat        = psi_hat + dt*S1_psi;
    %oy_hat         = oy_hat  + dt*S1_oy; 

    t = t + dt;
           
    % write to disk
    if mod(j,nsave)==0,
        toc 
        disp(['total time simulated: ',num2str(t)] )
        disp(['dt: '                  ,num2str(dt)])
        v_hat  = psi_hat./poisson;
        u_hat  = ( -1i*oy_hat.*kz - kx.*ky.*v_hat )./poisson2D;
        w_hat  = (  1i*oy_hat.*kx - ky.*kz.*v_hat )./poisson2D;
        u      = real( ifftn(u_hat,ndim) );
        v      = real( ifftn(v_hat,ndim) );
        w      = real( ifftn(w_hat,ndim) );
        E      = (sum(u(:).^2) + sum(v(:).^2) + sum(w(:).^2))/(Lx*Ly*Lz);
        compute_dissp;
        disspt = sum(dissp(:))/(Lx*Ly*Lz);
        disp(['total energy     : ',num2str(E)])
        disp(['total dissipation: ',num2str(disspt)])
        ii    = ii+1;
        save(['./data/vortices3D.RK.',int2str(ii),'.mat'],'u','v','w','Lx','Ly','Lz','Re','t','E','disspt');
        tic
    end
    
%    drawnow
%    pcolor(squeeze(v(1,:,:)))%,shading flat,caxis([-1.5 1.5])
    
end

end
