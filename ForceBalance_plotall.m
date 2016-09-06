ss = 1; %sumbsample time
tbegin=find(time_all>2.38e4,1);
tend=find(time_all>2.42e4,1);
sigbegin=1;
sigend=parameters.grid.n_elements;
%%
%calculate terms of force balance
h_g = -(parameters.rho_w/parameters.rho)*Base(xg_all,parameters);
% GL_long = 2.*parameters.B_Glen.*h_g;
% GL_driving = 0.5.*parameters.rho.*...
%     (1-(parameters.rho/parameters.rho_w)).*parameters.g.*(h_g.^2);

B_Glen = zeros(size(h_all));
for q=1:length(time_all)
    B_Glen(:,q) = mean(set_B_Glen(squeeze(T_all(:,:,q)),parameters),2);
end

dudx = diff(squeeze(u_full_all(:,1,:)))./(repmat(xg_all,parameters.grid.n_nodes-1,1).*...
    repmat(diff(parameters.grid.sigma_node),1,length(xg_all)));
stress_long = diff(2.*B_Glen.*h_all.*...
    (abs(dudx).^((1/parameters.n_Glen)-1)).*dudx)./...
    (repmat(xg_all,parameters.grid.n_nodes-2,1).*...
    repmat(diff(parameters.grid.sigma_node(1:end-1)),1,length(xg_all)));

u_mean = (squeeze(u_full_all(1:end-1,1,:)+u_full_all(2:end,1,:)))./2;
stress_taub = parameters.a_till.*exp(parameters.b_till.*e_all).*...
    u_mean./sqrt(u_mean.^2+parameters.u_eps.^2);
% u_mean(w_all==parameters.w_min) = 0;

stress_taul = parameters.Dupont_G.*h_all.*...
    (abs(u_mean).^((1/parameters.n_Glen)-1)).*u_mean;

bed_elev = Base(repmat(xg_all,parameters.grid.n_elements,1).*...
    repmat(parameters.grid.sigma_element,1,length(xg_all)),parameters);
h_mean = (h_all(1:end-1,:)+h_all(2:end,:))./2;
stress_driving = parameters.rho.*parameters.g.*h_mean.*diff(h_all+bed_elev)./...
    (repmat(xg_all,parameters.grid.n_elements-1,1).*...
    repmat(diff(parameters.grid.sigma_element),1,length(xg_all)));

%calculate terms of MW budget

HF_friction = u_mean.*stress_taub;

HF_geo = parameters.G.*ones(size(stress_taub));

dh_basal = h_all.*(parameters.grid.eta_element(2)-parameters.grid.eta_element(1));
HF_conduction = parameters.k_i.*(squeeze(T_all(:,2,:)-T_all(:,1,:))./dh_basal);

% ResidualHeating = MWHF_all-(HF_friction+HF_geo+HF_conduction);

%%
% figure(3)
% subplot(2,1,1)
% plot(time_all(tbegin:ss:tend),GL_driving(tbegin:ss:tend),'k.');hold on
% xlabel('time (years)','fontsize',16);
% ylabel('GL Driving','fontsize',16);

figure(4);set(4,'units','normalized','position',[0 0.1 1 0.75]);
subplot(2,2,1)
pcolor(time_all(tbegin:ss:tend),parameters.grid.sigma_element(1:end-1),stress_long(:,tbegin:ss:tend)./1000);hold on
shading('flat')
colorbar;
xlabel('time (years)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
set(gca,'YDir','Reverse')
title('Longitudinal Stress (kPa)','fontsize',16);
caxis([0 20])

subplot(2,2,2)
pcolor(time_all(tbegin:ss:tend),parameters.grid.sigma_element,stress_taub(:,tbegin:ss:tend)./1000);hold on
shading('flat')
colorbar;
set(gca,'YDir','Reverse')
xlabel('time (years)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Basal Shear Stress (kPa)','fontsize',16);
caxis([0 20])

subplot(2,2,3)
pcolor(time_all(tbegin:ss:tend),parameters.grid.sigma_element,stress_taul(:,tbegin:ss:tend)./1000);hold on
shading('flat')
colorbar;
set(gca,'YDir','Reverse')
xlabel('time (years)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Lateral Shear Stress (kPa)','fontsize',16);
caxis([0 20])

subplot(2,2,4)
pcolor(time_all(tbegin:ss:tend),parameters.grid.sigma_element(1:end-1),-stress_driving(:,tbegin:ss:tend)./1000);hold on
shading('flat')
colorbar;
set(gca,'YDir','Reverse')
xlabel('time (years)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Driving Stress (kPa)','fontsize',16);
caxis([0 20])
%%
figure(5)
subplot(3,1,1)
pcolor(time_all(tbegin:ss:tend)./1000,parameters.grid.sigma_element(sigbegin:sigend-1),HF_geo((sigbegin:sigend-1),tbegin:ss:tend));hold on
shading('flat')
colorbar;
xlabel('time (kyr)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
set(gca,'YDir','Reverse')
title('Geothermal Heat Flux (W/m^2)','fontsize',16);
caxis([0 1])
set(gca,'fontsize',20)

subplot(3,1,2)
pcolor(time_all(tbegin:ss:tend)./1000,parameters.grid.sigma_element(sigbegin:sigend-1),HF_friction((sigbegin:sigend-1),tbegin:ss:tend));hold on
shading('flat')
colorbar;
set(gca,'YDir','Reverse')
xlabel('time (kyr)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Frictional Heat Flux (W/m^2)','fontsize',16);
caxis([0 1])
set(gca,'fontsize',20)

subplot(3,1,3)
pcolor(time_all(tbegin:ss:tend)./1000,parameters.grid.sigma_element(sigbegin:sigend-1),-HF_conduction((sigbegin:sigend-1),tbegin:ss:tend));hold on
shading('flat')
colorbar;
set(gca,'YDir','Reverse')
xlabel('time (kyr)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Conductive Heat Flux (W/m^2)','fontsize',16);
caxis([0 1])
set(gca,'fontsize',20)
%%
figure(6)
contourf(time_all(tbegin:ss:tend),parameters.grid.sigma_element(sigbegin:sigend-1),...
    HF_geo((sigbegin:sigend-1),tbegin:ss:tend)+...
    HF_friction((sigbegin:sigend-1),tbegin:ss:tend)+...
    HF_conduction((sigbegin:sigend-1),tbegin:ss:tend),linspace(-0.1,0.1,11));hold on
colorbar;
set(gca,'YDir','Reverse')
xlabel('time (years)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Net Heat Flux (W/m^2)','fontsize',16);
caxis([-0.1 0.1])
colormap(redblue)
set(gca,'fontsize',20)

% figure(8)
% subplot(2,1,1)
% pcolor(time_all(tbegin:ss:tend),parameters.grid.sigma_element(sigbegin:sigend-1),...
%     dh_basal((sigbegin:sigend-1),tbegin:ss:tend));hold on
% shading('flat')
% colorbar;
% set(gca,'YDir','Reverse')
% xlabel('time (years)','fontsize',16);
% ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
% title('Basal Ice thickness (m)','fontsize',16);
% 
% subplot(2,1,2)
% pcolor(time_all(tbegin:ss:tend),parameters.grid.sigma_element(sigbegin:sigend-1),...
%     squeeze(T_all((sigbegin:sigend-1),2,tbegin:ss:tend)-T_all((sigbegin:sigend-1),1,tbegin:ss:tend)));hold on
% shading('flat')
% colorbar;
% set(gca,'YDir','Reverse')
% xlabel('time (years)','fontsize',16);
% ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
% title('Basal Temperature Gradient (K)','fontsize',16);



%%
pt = 30;
tbegin_zoom=15200;
tend_zoom=16000;
% tbegin_zoom = 11500;
% tend_zoom = 13000;

figure(10)
subplot(3,1,1)
plotyy(time_all(tbegin_zoom:ss:tend_zoom),h_all(pt,tbegin_zoom:ss:tend_zoom),time_all(tbegin_zoom:ss:tend_zoom),parameters.year.*u_mean(pt,tbegin_zoom:ss:tend_zoom));hold on
ylabel('Thickness','fontsize',16)
% ylim([195 260])
% xlim([time_all(tbegin_zoom) time_all(tend_zoom)])

% subplot(5,1,2)
% plot(,'r','linewidth',2)
% % legend('Thickness Up','Velocity Up','Thickness Down','Velocity Down')
% ylabel('GL Velocity','fontsize',16)

subplot(3,1,2)
plotyy(time_all(tbegin_zoom:ss:tend_zoom),w_all(pt,tbegin_zoom:ss:tend_zoom),time_all(tbegin_zoom:ss:tend_zoom),-squeeze(T_all(pt,2,tbegin_zoom:ss:tend_zoom)-T_all(pt,1,tbegin_zoom:ss:tend_zoom)));hold on
ylabel('Till Water Content','fontsize',16)
ylim([0.4 1.0])
% xlim([time_all(tbegin_zoom) time_all(tend_zoom)])

% subplot(3,1,4)
% plot(time_all(tbegin_zoom:ss:tend_zoom),-squeeze(T_all(pt,2,tbegin_zoom:ss:tend_zoom)-T_all(pt,1,tbegin_zoom:ss:tend_zoom)),'r','linewidth',2);
% ylabel('GL Basal Temp Gradient (K)','fontsize',16)

subplot(3,1,3)
plot(time_all(tbegin_zoom:ss:tend_zoom),HF_geo(pt,tbegin_zoom:ss:tend_zoom),'b');hold on
plot(time_all(tbegin_zoom:ss:tend_zoom),HF_friction(pt,tbegin_zoom:ss:tend_zoom),'r')
plot(time_all(tbegin_zoom:ss:tend_zoom),-HF_conduction(pt,tbegin_zoom:ss:tend_zoom),'m')
% plot(time_all(tbegin_zoom:ss:tend_zoom),HF_geo(pt,tbegin_zoom:ss:tend_zoom)+HF_friction(pt,tbegin_zoom:ss:tend_zoom)+HF_conduction(pt,tbegin_zoom:ss:tend_zoom),'k','linewidth',2);
% plot(time_all(tbegin_zoom:ss:tend_zoom),ResidualHeating(pt,tbegin_zoom:ss:tend_zoom),'Color',[.7 .7 .7]);hold on
% plot(time_all(tbegin_zoom:ss:tend_zoom),MWHF_all(pt,tbegin_zoom:ss:tend_zoom),'k','linewidth',2);hold on
% legend('Geothermal','Frictional','Diffusive','Residual Heating','Net')
% plot(time_all(tbegin_zoom:ss:tend_zoom),zeros(size(tbegin_zoom:ss:tend_zoom)),'k--');hold on
xlim([time_all(tbegin_zoom) time_all(tend_zoom)])
ylim([-0.15 0.15])
%%
pt = 20;
pt2 = 40;
tbegin=15200;
tend=16000;

figure(6);
subplot(3,1,1)
[AX,H1,H2] = plotyy(time_all(tbegin:ss:tend),h_all(pt,tbegin:ss:tend),time_all(tbegin:ss:tend),parameters.year.*u_mean(pt,tbegin:ss:tend));hold on
set(H1,'linewidth',2);set(H2,'linewidth',2)
set(get(AX(1),'Ylabel'),'String','Thickness (m)','Fontsize',16) 
set(get(AX(2),'Ylabel'),'String','Velocity (m/yr)','Fontsize',16)
title('~3 km upstream of GL','Fontsize',16)
[AX,H1,H2] = plotyy(time_all(tbegin:ss:tend),h_all(pt2,tbegin:ss:tend),time_all(tbegin:ss:tend),parameters.year.*u_mean(pt2,tbegin:ss:tend));hold on
set(H1,'linewidth',2,'linestyle','--');set(H2,'linewidth',2,'linestyle','--')
set(get(AX(1),'Ylabel'),'String','Thickness (m)','Fontsize',16) 
set(get(AX(2),'Ylabel'),'String','Velocity (m/yr)','Fontsize',16)
title('Upstream point solid, downstream point dashed','Fontsize',16)

subplot(3,1,2)
plot(time_all(tbegin:ss:tend),stress_long(pt,tbegin:ss:tend),'b','linewidth',2);hold on
plot(time_all(tbegin:ss:tend),stress_taub(pt,tbegin:ss:tend),'r','linewidth',2)
plot(time_all(tbegin:ss:tend),stress_taul(pt,tbegin:ss:tend),'m','linewidth',2)
plot(time_all(tbegin:ss:tend),-stress_driving(pt,tbegin:ss:tend),'k','linewidth',2)
plot(time_all(tbegin:ss:tend),stress_long(pt,tbegin:ss:tend)+stress_taub(pt,tbegin:ss:tend)+stress_taul(pt,tbegin:ss:tend)+stress_driving(pt,tbegin:ss:tend),'Color',[.7 .7 .7],'linewidth',2)
% plot(time_all(tbegin:ss:tend),stress_long(pt2,tbegin:ss:tend),'b--','linewidth',2);hold on
% plot(time_all(tbegin:ss:tend),stress_taub(pt2,tbegin:ss:tend),'r--','linewidth',2)
% plot(time_all(tbegin:ss:tend),stress_taul(pt2,tbegin:ss:tend),'m--','linewidth',2)
% plot(time_all(tbegin:ss:tend),-stress_driving(pt2,tbegin:ss:tend),'k--','linewidth',2)
% plot(time_all(tbegin:ss:tend),stress_long(pt,tbegin:ss:tend)+stress_taub(pt,tbegin:ss:tend)+stress_taul(pt,tbegin:ss:tend)+stress_driving(pt,tbegin:ss:tend),'Color',[.7 .7 .7],'linestyle','--','linewidth',2)
legend('Long Up','Basal Up','Lateral Up','Driving Up','Total Up');%...
%     ,'Long Down','Basal Down','Lateral Down','Driving Down','Total Down')
subplot(3,1,3)
plot(time_all(tbegin:ss:tend),HF_geo(pt,tbegin:ss:tend),'b','linewidth',2);hold on
plot(time_all(tbegin:ss:tend),HF_friction(pt,tbegin:ss:tend),'r','linewidth',2)
plot(time_all(tbegin:ss:tend),-HF_conduction(pt,tbegin:ss:tend),'m','linewidth',2)
plot(time_all(tbegin:ss:tend),HF_geo(pt2,tbegin:ss:tend),'b--','linewidth',2);hold on
plot(time_all(tbegin:ss:tend),HF_friction(pt2,tbegin:ss:tend),'r--','linewidth',2)
plot(time_all(tbegin:ss:tend),-HF_conduction(pt2,tbegin:ss:tend),'m--','linewidth',2)
legend('Geothermal Up','Frictional Up','Conductive Up','Geothermal Down','Frictional Down','Conductive Down')

% subplot(2,1,1)
% [AX,H1,H2] = plotyy(time_all(tbegin:ss:tend)-1e4,stress_long(pt,tbegin:ss:tend),time_all(tbegin:ss:tend)-1e4,HF_friction(pt,tbegin:ss:tend));hold on
% set(H1,'linewidth',2);set(H2,'linewidth',2)
% set(get(AX(1),'Ylabel'),'String','Longitudinal Stress (Pa)','Fontsize',16) 
% set(get(AX(2),'Ylabel'),'String','Frictional Heat Flux (W/m^2)','Fontsize',16)
% title('~3 km upstream of GL','Fontsize',16)
% subplot(2,1,2)
% [AX,H1,H2] = plotyy(time_all(tbegin:ss:tend)-1e4,stress_long(pt2,tbegin:ss:tend),time_all(tbegin:ss:tend)-1e4,HF_friction(pt2,tbegin:ss:tend));hold on
% set(H1,'linewidth',2);set(H2,'linewidth',2)
% set(get(AX(1),'Ylabel'),'String','Longitudinal Stress (Pa)','Fontsize',16) 
% set(get(AX(2),'Ylabel'),'String','Frictional Heat Flux (W/m^2)','Fontsize',16)
% title('~23 km upstream of GL','Fontsize',16)

% figure(7)
% subplot(3,1,1)
% plot(time_all(tbegin:ss:tend),h_all(pt,tbegin:ss:tend),'b','linewidth',2);hold on
% plot(time_all(tbegin:ss:tend),parameters.year.*u_mean(pt,tbegin:ss:tend),'r','linewidth',2)
% plot(time_all(tbegin:ss:tend),h_all(pt2,tbegin:ss:tend),'b--','linewidth',2);hold on
% plot(time_all(tbegin:ss:tend),parameters.year.*u_mean(pt2,tbegin:ss:tend),'r--','linewidth',2)
% legend('Thickness Up','Velocity Up','Thickness Down','Velocity Down')
% 
% subplot(3,1,2)
% plot(time_all(tbegin:ss:tend),w_all(pt,tbegin:ss:tend),'b','linewidth',2);hold on
% plot(time_all(tbegin:ss:tend),-squeeze(T_all(pt,2,tbegin:ss:tend)-T_all(pt,1,tbegin:ss:tend)),'r','linewidth',2);
% plot(time_all(tbegin:ss:tend),w_all(pt2,tbegin:ss:tend),'b--','linewidth',2);hold on
% plot(time_all(tbegin:ss:tend),-squeeze(T_all(pt2,2,tbegin:ss:tend)-T_all(pt2,1,tbegin:ss:tend)),'r--','linewidth',2);
% legend('Till Water Content Up','T1-T2 Up','Till Water Content Down','T1-T2 Down')
% 
% subplot(3,1,3)
% plot(time_all(tbegin:ss:tend),HF_geo(pt,tbegin:ss:tend),'b','linewidth',2);hold on
% plot(time_all(tbegin:ss:tend),HF_friction(pt,tbegin:ss:tend),'r','linewidth',2)
% plot(time_all(tbegin:ss:tend),-HF_conduction(pt,tbegin:ss:tend),'m','linewidth',2)
% plot(time_all(tbegin:ss:tend),HF_geo(pt2,tbegin:ss:tend),'b--','linewidth',2);hold on
% plot(time_all(tbegin:ss:tend),HF_friction(pt2,tbegin:ss:tend),'r--','linewidth',2)
% plot(time_all(tbegin:ss:tend),-HF_conduction(pt2,tbegin:ss:tend),'m--','linewidth',2)
% legend('Geothermal Up','Frictional Up','Conductive Up','Geothermal Down','Frictional Down','Conductive Down')