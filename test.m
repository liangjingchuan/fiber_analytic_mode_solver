% test
close all;
clear all;

% n_temp=(1:1:search_num-1);
% figure;
% plot(n_temp,search_product);

% % bessel function
% x=linspace(0,10,100);
% y0=besselj(0,x);
% y1=besselj(1,x);
% y2=besselj(2,x);
% y3=besselj(3,x);
% figure; hold on;
% plot(x,y0,'k');
% plot(x,y1,'r');
% plot(x,y2,'g');
% plot(x,y3,'b');
% hold off;
% grid on;
% figure;
% plot(x,y0+y1);

%%%
% %%% test
% x=linspace(n2,n1,100);
% NEFF=x;
% core_a=a;
% omega=2*pi*c/lambda;
% k0=2*pi/lambda;
% k1=n1*k0;
% k2=n2*k0;
% kc1=sqrt(n1^2-NEFF.^2)*k0; % u
% kc2=sqrt(NEFF.^2-n2^2)*k0; % w
% orderi=1;
% J=besselj(orderi,kc1*core_a);
% K=besselk(orderi,kc2*core_a);
% J_prime=(besselj(orderi-1,kc1.*core_a)-besselj(orderi+1,kc1.*core_a))/2.0;
% K_prime=-(besselk(orderi-1,kc2.*core_a)+besselk(orderi+1,kc2.*core_a))/2.0;
% J_math=(1./(kc1.*core_a)).*(J_prime./J);
% K_math=(1./(kc2.*core_a)).*(K_prime./K);
% Y=(J_math+K_math).*(k1^2.*J_math+k2^2.*K_math)-(NEFF.*k0.*orderi./core_a).^2.*(1./kc1.^2+1./kc2.^2).^2;
% figure;
% plot(NEFF,Y);


% x1=linspace(0,a,100);
% x2=linspace(a,5*a,400);
% rho=[x1,x2(2:end)];
% bessel_1=vec_coeff(1)*besselj(m,u*x1);
% bessel_3=vec_coeff(3)*besselk(m,w*x2);
% Ez=[bessel_1,bessel_3(2:end)];
% bessel_2=vec_coeff(2)*besselj(m,u*x1);
% bessel_4=vec_coeff(4)*besselk(m,w*x2);
% Hz=[bessel_2,bessel_4(2:end)];
% figure;
% plot(rho,Ez);
% title('Ez');
% figure;
% plot(rho,abs(Hz));
% title('Hz');

ntype=2;
a=0;
b=0;

switch ntype
    case 1
        a=2;
    case 2
        a=3;
    otherwise
         disp('Unknown method.')
end

a



