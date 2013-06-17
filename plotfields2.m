%$Id: plotfields.m,v 1.13 2006/09/08 11:53:25 taolu Exp $
%$Revision: 1.13 $
%$Author: taolu $
%$Date: 2006/09/08 11:53:25 $
%function figid=plotfields(mesh_X,mesh_Y,Field,Phi,figid,levels,file_label)
function figid=plotfields2(fiber_geom,Field,figid)

% %%% parameters
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
rescale=1e6;
mesh_X=fiber_geom.mesh_grids*rescale;
mesh_Y=fiber_geom.mesh_grids*rescale;
x_min=min(mesh_X);
x_max=max(mesh_X);
y_min=min(mesh_Y);
y_max=max(mesh_Y);

%%% wave guide
radius=fiber_geom.core_width*rescale;
center_X=0;
center_Y=0;
num_point=600;
a_angle=linspace(0,2*pi,num_point);
a_x=radius.*cos(a_angle)+center_X;
a_y=radius.*sin(a_angle)+center_Y;
Linewidth=1;


% field_str={'E_x';'E_y';'E_z';'H_x';'H_y';'H_z'};
% z_min=min(min(real(Field)));
% z_max=max(max(real(Field)))+eps;

% a_max=max(max(abs(real(Field))));
% b_max=max(max(abs(Field)));
% all_max=[reshape(a_max,6,1) reshape(b_max,6,1)];
% z_max=max(max(real(Field)))+eps;
% z_min=min(min(real(Field)));
% z_max=reshape(z_max,6,1)+eps;
% z_min=reshape(z_min,6,1);

figure(figid);clf;
figid=figid+1;
%subplot(231);
hold on;
imagesc(mesh_X,mesh_Y,real(Field(:,:,1))); 
title('E_x');
axis([x_min x_max y_min y_max]);
xlabel('X [\mum]');ylabel('Y [\mum]');
colorbar;
%contour(mesh_X,mesh_Y,N_indx);
plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
hold off;

figure(figid);clf;
figid=figid+1;
%subplot(232);
hold on;
imagesc(mesh_X,mesh_Y,real(Field(:,:,2))); title('E_y');
axis([x_min x_max y_min y_max]);
xlabel('X [\mum]');ylabel('Y [\mum]');
colorbar;
plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
hold off;

figure(figid);clf;
figid=figid+1;
%subplot(233);
hold on;
imagesc(mesh_X,mesh_Y,real(Field(:,:,3))); title('E_z');
axis([x_min x_max y_min y_max]);
xlabel('X [\mum]');ylabel('Y [\mum]');
colorbar;
plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
hold off;

figure(figid);clf;
figid=figid+1;
%subplot(234);
hold on;
imagesc(mesh_X,mesh_Y,real(Field(:,:,4))); title('H_x');
axis([x_min x_max y_min y_max]);
xlabel('X [\mum]');ylabel('Y [\mum]');
colorbar;
plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
hold off;


figure(figid);clf;
figid=figid+1;
%subplot(235);
hold on;
imagesc(mesh_X,mesh_Y,real(Field(:,:,5))); title('H_y');
axis([x_min x_max y_min y_max]);
xlabel('X [\mum]');ylabel('Y [\mum]');
colorbar;
plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
hold off;

figure(figid);clf;
figid=figid+1;
%subplot(236);
hold on;
imagesc(mesh_X,mesh_Y,real(Field(:,:,6))); title('H_z');
axis([x_min x_max y_min y_max]);
xlabel('X [\mum]');ylabel('Y [\mum]');
colorbar;
plot(a_x,a_y,'k','LineWidth',Linewidth); % plot waveguide
hold off;

% if exist('file_label')
%     if ~isempty(file_label)
%         saveas(gcf,[file_label int2str(figid-1) '_Et.fig'], 'fig');
%     end
% end
% 
% figid=figid+1;
% I=(Field(:,:,1).*conj(Field(:,:,5)))-(Field(:,:,2).*conj(Field(:,:,4)));
% figure(figid);clf;
%        subplot(211);mesh(mesh_X,mesh_Y,real(I));
%       subplot(212);
% contour(mesh_X,mesh_Y,real(I),levels); plot_waveguide;
% title(file_label);

% otherwise
%     disp('In plotfields.m not implemented');
% end

