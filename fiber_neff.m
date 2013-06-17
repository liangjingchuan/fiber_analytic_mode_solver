%function y=fiber_neff(x,orderi,fiber_para)
%This function find the analytic solution of Neff in fiber
function y=fiber_neff(x,orderi,fiber_geom)

NEFF=x;
n1=fiber_geom.n_core;
n2=fiber_geom.n_cladding;
lambda=fiber_geom.lambda;
core_a=fiber_geom.core_width;

k0=2*pi/lambda;
k1=n1*k0;
k2=n2*k0;
kc1=sqrt(n1^2-NEFF.^2)*k0; % u
kc2=sqrt(NEFF.^2-n2^2)*k0; % w

if orderi~=0
    J=besselj(orderi,kc1*core_a);
    K=besselk(orderi,kc2*core_a);
    J_prime=(besselj(orderi-1,kc1*core_a)-besselj(orderi+1,kc1*core_a))/2.0;
    K_prime=-(besselk(orderi-1,kc2*core_a)+besselk(orderi+1,kc2*core_a))/2.0;
    J_math=(1./(kc1)).*(J_prime./J);
    K_math=(1./(kc2)).*(K_prime./K);
    y=(J_math+K_math).*(k1^2.*J_math+k2^2.*K_math)-(NEFF.*k0.*orderi./core_a).^2.*(1./kc1.^2+1./kc2.^2).^2;
else
    J0=besselj(0,kc1*core_a);
    K0=besselk(0,kc2*core_a);
    J1=besselj(1,kc1*core_a);
    K1=besselk(1,kc2*core_a);
    y=(J1./(kc1.*J0)+K1./(kc2.*K0)).*(k1^2.*J1./(kc1.*J0)+k2^2.*K1./(kc2.*K0));
end


    
    
    
    
%     %%% Tao %%%
%     M=zeros(4,4);
%     M(1,1)=J;
%     M(1,3)=-K;
%     M(2,1)=((k0*orderi*NEFF)/(kc1^2*core_a))*J;
%     M(2,2)=(1j*omega*mu/kc1)*J_prime;
%     M(2,3)=((k0*orderi*NEFF)/(kc2^2*core_a))*K;
%     M(2,4)=(1j*omega*mu/kc2)*K_prime;
%     M(3,2)=M(1,1);
%     M(3,4)=M(1,3);
%     M(4,1)=-(1j*omega*ep1/kc1)*J_prime;
%     M(4,2)=M(2,1);
%     M(4,3)=-(1j*omega*ep2/kc2)*K_prime;
%     M(4,4)=M(2,3);
%     
%     ratio=-M(1,3)/M(1,1);
%     M2=[M(1,:);M(3,:);M(2,:);M(4,:)];
%     M3=M2(3:4,3:4)+ratio*M2(3:4,1:2);
%     X_prime=(1/kc2)*(K_prime/K);
%     Y_prime=(1/kc1)*(J_prime/J);
%     rho=(n1/n2)^2;
%     char_fac=(1-rho)^2*X_prime^2/4.0;
%     char_fac=char_fac+orderi^2*((NEFF/n2)*((1/(kc1*core_a)^2)+(1/(kc2*core_a)^2)))^2;
%     char_sqrt=Y_prime+(1+rho)*X_prime/2.0+sqrt(char_fac);
%     
%     y=log(1+10*abs(det(M)))-1;
%     y2=abs(char_sqrt);
%     y3=abs(det(M3));
