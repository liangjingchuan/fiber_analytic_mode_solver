% fiber mode

close all;
clear all;

%%% parameters input
fiber_para.n_core=1.45;
fiber_para.n_cladding=1.0;
fiber_para.core_width=1e-6; % radius
fiber_para.lambda=1.55e-6;
fiber_para.num_grids=401;
fiber_para.mesh_grids=linspace(-5*fiber_para.core_width,5*fiber_para.core_width,fiber_para.num_grids); 
fiber_para.V=(2*pi*fiber_para.core_width/fiber_para.lambda)*sqrt(fiber_para.n_core^2-fiber_para.n_cladding^2); % single mode when V<2.405
fiber_para.M=round(fiber_para.V^2/2); % number of guided modes
fiber_para.num_region=2;

%%% output
% plot_type='3D'; % 2D,3D
result={};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e0=(1.0E-9)/(36*pi);
u0=(1.0E-7)*4*pi;
Z0=sqrt(u0/e0); % impedance
c=1/sqrt(e0*u0); % 3e8 speed of light in free space

a=fiber_para.core_width;
lambda=fiber_para.lambda;
omega=2*pi*c/lambda;
k0=2*pi/lambda;
n1=fiber_para.n_core;
n2=fiber_para.n_cladding;
mu=u0;
ep1=e0*n1^2;
ep2=e0*n2^2;

num_bessel=fiber_para.M;

% i_orderm=1; % nu or m
for i_orderm=0:num_bessel
    
    %%% find the effective refractive index
    n_min=n2;
    n_max=n1;
    %n_bottom=fminbnd(@fiber_neff,n_min,n_max,optimset('TolX',1e-10),orderi,fiber_para);
    %n_search=(n_min+n_max)/2;
    search_num=1000;
    search_array=linspace(n_max,n_min,search_num);
    search_value=fiber_neff(search_array,i_orderm,fiber_para);
    search_product=search_value(1:end-1).*search_value(2:end);
    search_indx=find(search_product<=0);
    
    n_array=search_array(search_indx);
    num_zero=length(search_indx);
    sub_order=ones(4,1); %% suborder ot bessel function
    
    % i_nindx=1;
    for i_ordern=1:num_zero
        
        n_search=[search_array(search_indx(i_ordern)) search_array(search_indx(i_ordern)+1)];
        n_eff=fzero(@(x) fiber_neff(x,i_orderm,fiber_para),n_search);
        
        %%% calculate the coefficients (A,B,C,D)
        u=sqrt(n1^2-n_eff^2)*k0; % u
        w=sqrt(n_eff^2-n2^2)*k0; % w
        
        m=i_orderm;
        
        Ja=besselj(m,u*a);
        Ka=besselk(m,w*a);
        if m~=0
            Ja_prime=(besselj(m-1,u*a)-besselj(m+1,u*a))/2.0;
            Ka_prime=-(besselk(m-1,w*a)+besselk(m+1,w*a))/2.0;
        else
            Ja_prime=-besselj(1,u*a);
            Ka_prime=-besselk(1,w*a);
        end
        
        M=zeros(4,4);
        M(1,1)=Ja;
        M(1,3)=-Ka;
        M(2,1)=((n_eff*k0*m)/(a*u^2))*Ja;
        M(2,2)=(1j*omega*mu/u)*Ja_prime;
        M(2,3)=((n_eff*k0*m)/(a*w^2))*Ka;
        M(2,4)=(1j*omega*mu/w)*Ka_prime;
        M(3,2)=M(1,1);
        M(3,4)=M(1,3);
        M(4,1)=-(1j*omega*ep1/u)*Ja_prime;
        M(4,2)=M(2,1);
        M(4,3)=-(1j*omega*ep2/w)*Ka_prime;
        M(4,4)=M(2,3);
        vec_coeff=null(M);
        
        Field=fiber_field(fiber_para,n_eff,vec_coeff,i_orderm);
        
        %%% results
        result(end+1).Neff=n_eff;
        result(end).Field=Field;
        result(end).order_m=i_orderm;
        result(end).ABCD=vec_coeff;
        
        %%% HE or EH
        if i_orderm~=0
            if vec_coeff(1)>vec_coeff(2)
                result(end).type='HE';
            else
                result(end).type='EH';
            end
        else
            if vec_coeff(1)==0
                result(end).type='TE';
            else
                result(end).type='TM';
            end
        end % if
        %%%
        switch result(end).type
            case 'TE'
                result(end).order_n=sub_order(1);
                sub_order(1)=sub_order(1)+1;
            case 'TM'
                result(end).order_n=sub_order(2);
                sub_order(2)=sub_order(2)+1;
            case 'HE'
                result(end).order_n=sub_order(3);
                sub_order(3)=sub_order(3)+1;
            case 'EH'
                result(end).order_n=sub_order(4);
                sub_order(4)=sub_order(4)+1;
            otherwise
                disp('Unknown modes.')
        end % switch
    end % i_ordern
end % i_orderm

temp_cell=struct2cell(result);
temp_mat=cell2mat(temp_cell(1,:));

[array_Neff, indx_Neff]=sort(temp_mat,'descend');
analytic_result=result(indx_Neff);

%%% plot field
figid=1;
plot_type='2D'; % 2D,3D
ModeNo=1;
plot_result=analytic_result; % specify the result to be ploted
Field=analytic_result(ModeNo).Field;
m_order=analytic_result(ModeNo).order_m;
n_order=analytic_result(ModeNo).order_n;
type=analytic_result(ModeNo).type;
Neff=analytic_result(ModeNo).Neff;

figid=plotfields(fiber_para,plot_result,figid,plot_type);
%figid=plotfields2(fiber_para,Field,figid);

rescale=1e6;
mesh_X=fiber_para.mesh_grids*rescale;
mesh_Y=fiber_para.mesh_grids*rescale;
x_min=min(mesh_X);
x_max=max(mesh_X);
y_min=min(mesh_Y);
y_max=max(mesh_Y);
%%% wave guide
radius=fiber_para.core_width*rescale;
center_X=0;
center_Y=0;
num_point=600;
a_angle=linspace(0,2*pi,num_point);
a_y=radius.*sin(a_angle)+center_Y;
a_x=radius.*cos(a_angle)+center_X;
Linewidth=1;
Intensity=0.5.*(Field(:,:,1).*conj(Field(:,:,5))-Field(:,:,2).*conj(Field(:,:,4)));

figure(figid);figid=figid+1;
hold on;
imagesc(mesh_X,mesh_Y,real(Intensity));  title(['Intensity of ',type ,'_',num2str(m_order),'_',num2str(n_order),'   n_e=',num2str(Neff)]);
axis([x_min x_max y_min y_max]);
xlabel('X [\mum]');ylabel('Y [\mum]');
colorbar;
plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
hold off;

% figure(figid);figid=figid+1;
% hold on;
% imagesc(mesh_X,mesh_Y,N_distribution); title('refractive index (n)');
% axis([x_min x_max y_min y_max]);
% xlabel('X [\mum]');ylabel('Y [\mum]');
% colorbar;
% plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
% hold off;

% %%% test plot
xxx=fiber_para.mesh_grids;
yyy=diag(Field(:,:,1));
% yyy=Field(round(fiber_para.num_grids/2),:,1);
figure;
plot(xxx,real(yyy));
title(['Cross view: Ez','   n_e=',num2str(Neff)]);


