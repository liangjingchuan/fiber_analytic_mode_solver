%function Field=fiber_field(xx,coeffs,orderi,mesh_X,mesh_Y)
%This function calculate the field distribution of fiber analytically
%function [Field,N_indx]=fiber_field(fiber_geom,n_eff,coeffs,orderi)
function Field=fiber_field(fiber_geom,n_eff,coeffs,orderi)

% %%% parameters for test
% fiber_geom.n_core=1.45;
% fiber_geom.n_cladding=1.2;
% fiber_geom.core_width=1e-6;%radius
% fiber_geom.lambda=1.55e-6;
% fiber_geom.num_grids=201;
% fiber_geom.mesh_grids=linspace(-3*fiber_geom.core_width,3*fiber_geom.num_grids,fiber_geom.num_grids); 
% fiber_geom.num_region=2;
% n_eff=1.372786791039351;
% coeffs=[0.098506725303696 + 0.000000000000000i; -1.209159534986718e-17 - 3.164610092071529e-04i; 0.995131199635683 + 0.000000000000000i; 0.000000000000000 - 0.003196941353591i];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    e0=(1.0E-9)/(36*pi);
    u0=(1.0E-7)*4*pi;
    c=1/sqrt(e0*u0);

    k0=2*pi/fiber_geom.lambda;
    omega=2*pi*c/fiber_geom.lambda;
    n1=fiber_geom.n_core;
    n2=fiber_geom.n_cladding;
    ep1=e0*n1^2;
    ep2=e0*n2^2;
    a=fiber_geom.core_width;
    mesh_X=fiber_geom.mesh_grids;
    mesh_Y=fiber_geom.mesh_grids;
    beta=n_eff*k0;
    A=coeffs(1);
    B=coeffs(2);
    C=coeffs(3);
    D=coeffs(4);
    m=orderi;
    
    [ax,ay]=meshgrid(mesh_X,mesh_Y);
    ax_size=size(ax);
    int_size=ax_size(1)*ax_size(2);
    ax=reshape(ax,int_size,1);
    ay=reshape(ay,int_size,1);
    a_rho=sqrt(ax.^2+ay.^2);
    a_phi=angle(ax+ay*1j);
    F_tmp=zeros(int_size,6);
    % N_indx=zeros(int_size,1);
    % int_coor=[ax ay];
    
    ind_core=find(a_rho<=a);
    p2=(n1^2-n_eff^2)*k0^2;
    u=sqrt(p2);
    Jm=besselj(m,u*a_rho(ind_core));
    if m~=0
        Jm_prime=(besselj(m-1,u*a_rho(ind_core))-besselj(m+1,u*a_rho(ind_core)))/2.0;
    else
        Jm_prime=-besselj(1,u*a_rho(ind_core));
    end

    % N_indx(ind_core)=n1;
    F_tmp(ind_core,1)=-1j/p2*(A*beta*u*Jm_prime+B*1j*m*u0*omega./a_rho(ind_core).*Jm);
    F_tmp(ind_core,2)=-1j/p2*(A*1j*m*beta./a_rho(ind_core).*Jm-B*omega*u0*u*Jm_prime);
    F_tmp(ind_core,3)=A*Jm;
    F_tmp(ind_core,4)=-1j/p2*(B*beta*u*Jm_prime-A*1j*m*omega*ep1./a_rho(ind_core).*Jm);
    F_tmp(ind_core,5)=-1j/p2*(B*1j*m*beta./a_rho(ind_core).*Jm+A*omega*ep1*u*Jm_prime);
    F_tmp(ind_core,6)=B*Jm;
    
    ind_cladding=find(a_rho>a);
    q2=(n2^2-n_eff^2)*k0^2;
    w=sqrt(-q2);
    Km=besselk(m,w*a_rho(ind_cladding));
    if m~=0
        Km_prime=-(besselk(m-1,w*a_rho(ind_cladding))+besselk(m+1,w*a_rho(ind_cladding)))/2.0;
    else
        Km_prime=-besselk(1,w*a_rho(ind_cladding));
    end

    % N_indx(ind_cladding)=n2;
    F_tmp(ind_cladding,1)=-1j/q2*(C*beta*w*Km_prime+D*1j*m*u0*omega./a_rho(ind_cladding).*Km);
    F_tmp(ind_cladding,2)=-1j/q2*(C*1j*m*beta./a_rho(ind_cladding).*Km-D*omega*u0*w*Km_prime);
    F_tmp(ind_cladding,3)=C*Km;
    F_tmp(ind_cladding,4)=-1j/q2*(D*beta*w*Km_prime-C*1j*m*omega*ep2./a_rho(ind_cladding).*Km);
    F_tmp(ind_cladding,5)=-1j/q2*(D*1j*m*beta./a_rho(ind_cladding).*Km+C*omega*ep2*w*Km_prime);
    F_tmp(ind_cladding,6)=D*Km;
    
    for count=1:6
        F_tmp(:,count)=F_tmp(:,count).*exp(1j*m*a_phi);    
    end
    
    tmp1=(F_tmp(:,1).*cos(a_phi))-(F_tmp(:,2).*sin(a_phi)); % Ex <= E_rho+E_phi
    tmp2=(F_tmp(:,1).*sin(a_phi))+(F_tmp(:,2).*cos(a_phi)); % Ey <= E_rho+E_phi
    tmp3=(F_tmp(:,4).*cos(a_phi))-(F_tmp(:,5).*sin(a_phi)); % Hx <= H_rho+H_phi
    tmp4=(F_tmp(:,4).*sin(a_phi))+(F_tmp(:,5).*cos(a_phi)); % Hy <= H_rho+H_phi
    F_tmp(:,1)=tmp1; % Ex
    F_tmp(:,2)=tmp2; % Ey
    F_tmp(:,4)=tmp3; % Hx
    F_tmp(:,5)=tmp4; % Hy

    Field=zeros(ax_size(1),ax_size(2),6);
    for count=1:6
        Field(:,:,count)=reshape(F_tmp(:,count),ax_size(1),ax_size(2));
    end
    % N_indx=reshape(N_indx,ax_size(1),ax_size(2));
    
    
    
%     for regioni=1:LAST_REGION
%         ri_ind=nodesinregion(int_coor,regioni);
%         n_i=REGION_INDX(regioni);
%         if n_i>NEFF
%             kc1=sqrt(n_i^2-NEFF^2)*k0;
%             J_prime=(besselj(orderi-1,kc1*ar(ri_ind))-besselj(orderi+1,kc1*ar(ri_ind)))/2.0;
%             F_tmp(ri_ind,3)=coeffs(1)*(besselj(orderi,kc1*ar(ri_ind))); %Ez
%             F_tmp(ri_ind,6)=coeffs(2)*(besselj(orderi,kc1*ar(ri_ind))); %Hz
%             F_tmp(ri_ind,1)=(-j*NEFF*k0*coeffs(1)/kc1)*J_prime; %Erho
%             F_tmp(ri_ind,1)=F_tmp(ri_ind,1)+(k0*orderi*coeffs(2)/kc1^2)*(besselj(orderi,kc1*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,2)=(orderi*NEFF*k0*coeffs(1)/kc1^2)*(besselj(orderi,kc1*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,2)=F_tmp(ri_ind,2)+(j*k0*coeffs(2)/kc1)*J_prime;            
%             F_tmp(ri_ind,4)=-(j*NEFF*k0*coeffs(2)/kc1)*J_prime;
%             F_tmp(ri_ind,4)=F_tmp(ri_ind,4)-(k0*n_i^2*orderi*coeffs(1)/kc1^2)*(besselj(orderi,kc1*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,5)=(orderi*NEFF*k0*coeffs(2)/kc1^2)*(besselj(orderi,kc1*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,5)=F_tmp(ri_ind,5)-(j*n_i^2*k0*coeffs(1)/kc1)*J_prime;            
%         else
%             kc2=sqrt(NEFF^2-n_i^2)*k0;
%             K_prime=-(besselk(orderi-1,kc2*ar(ri_ind))+besselk(orderi+1,kc2*ar(ri_ind)))/2.0;
%             F_tmp(ri_ind,3)=coeffs(3)*(besselk(orderi,kc2*ar(ri_ind)));
%             F_tmp(ri_ind,6)=coeffs(4)*(besselk(orderi,kc2*ar(ri_ind)));
%             F_tmp(ri_ind,1)=(j*NEFF*k0*coeffs(3)/kc2)*K_prime;
%             F_tmp(ri_ind,1)=F_tmp(ri_ind,1)-(k0*orderi*coeffs(4)/kc2^2)*(besselk(orderi,kc2*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,2)=(-orderi*NEFF*k0*coeffs(3)/kc2^2)*(besselk(orderi,kc2*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,2)=F_tmp(ri_ind,2)-(j*k0*coeffs(4)/kc2)*K_prime;
%             F_tmp(ri_ind,4)=(j*NEFF*k0*coeffs(4)/kc2)*K_prime;
%             F_tmp(ri_ind,4)=F_tmp(ri_ind,4)+(n_i^2*k0*orderi*coeffs(3)/kc2^2)*(besselk(orderi,kc2*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,5)=(-orderi*NEFF*k0*coeffs(4)/kc2^2)*(besselk(orderi,kc2*ar(ri_ind))./ar(ri_ind));
%             F_tmp(ri_ind,5)=F_tmp(ri_ind,5)+(j*n_i^2*k0*coeffs(3)/kc2)*K_prime;
%         end
%     end;
%     
%     for count=1:6
%         F_tmp(:,count)=F_tmp(:,count).*exp(j*orderi*a_phi);    
%     end   
%    % 
%     tmp1=(F_tmp(:,1).*cos(a_phi))-(F_tmp(:,2).*sin(a_phi));
%     tmp2=(F_tmp(:,1).*sin(a_phi))+(F_tmp(:,2).*cos(a_phi));
%     tmp3=(F_tmp(:,4).*cos(a_phi))-(F_tmp(:,5).*sin(a_phi));
%     tmp4=(F_tmp(:,4).*sin(a_phi))+(F_tmp(:,5).*cos(a_phi));
%     F_tmp(:,1)=tmp1;
%     F_tmp(:,2)=tmp2;
%     F_tmp(:,4)=tmp3;
%     F_tmp(:,5)=tmp4;
%     
%     Field=zeros(ax_size(1),ax_size(2),6);
%     for count=1:6
%         Field(:,:,count)=reshape(F_tmp(:,count),ax_size(1),ax_size(2));
%     end
