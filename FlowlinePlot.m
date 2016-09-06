normalPlots = 1;

if(normalPlots)

h_g=-(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
bed_elev = Base(x_g.*parameters.grid.sigma_node,parameters);

z_element=(parameters.grid.eta_node*(h'))';
z_element=z_element+repmat(Base(x_g.*parameters.grid.sigma_element,parameters),[1 parameters.grid.n2_nodes]);
x_element=repmat(x_g*parameters.grid.sigma_element,[1 parameters.grid.n2_nodes]);

T_s = parameters.T_s;
T_lin=zeros(parameters.grid.n_elements,parameters.grid.n2_nodes);
for i=1:parameters.grid.n_elements
    T_lin(i,:) = linspace(parameters.MP,T_s(i),parameters.grid.n2_nodes);
end
T_lin = (T_lin(:,1:parameters.grid.n2_nodes-1)+T_lin(:,2:parameters.grid.n2_nodes))./2;

%%

for i=1:parameters.grid.n_elements
    T_lin(i,:) = linspace(0,T_s(i),parameters.grid.n2_elements);
end

if(t~=1)
    set(p1,'color','k','linewidth',1)
    set(p2,'color','k','linewidth',1)
%     set(p2b,'color','k')
    set(p3b,'color','k')
    set(p3,'color','k','linewidth',1)
    set(p4,'color','k','linewidth',1)
    set(p6,'color','k','linewidth',1)
end
    
%%

figure(1);set(1,'units','normalized','position',[0 0.1 1 0.75]);
subplot(3,2,1)
% p1 = plot(parameters.grid.sigma_node,u_old*parameters.year,'r');hold on
p1 = plot(parameters.grid.sigma_node,u_old_mean*parameters.year,'r','linewidth',3);hold on
xlabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
ylim([min(u_old_mean*parameters.year) max([u_old_mean*parameters.year;10])]);
ylabel('Velocity (m/yr)','fontsize',16)

subplot(3,2,2)
p2 = plot(parameters.grid.sigma_element,MW_HF,'r','linewidth',3);hold on
% p2b = plot(parameters.grid.sigma_node,ExtraHF_warm,'r');
% p2c = plot(parameters.grid.sigma_node,ExtraHF_cold,'b');
xlabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
ylabel('Heat Flux (W/m^2)','fontsize',16)
ylim([-0.1 0.1]);
drawnow

subplot(3,2,3)
plot(x_g.*parameters.grid.sigma_node/1000,bed_elev,'k','linewidth',2);hold on
p3=plot([x_g.*parameters.grid.sigma_node;x_g]/1000,[bed_elev+[h;h_g];bed_elev(end)],'r','linewidth',3);hold on
p3b=plot(x_g/1000,bed_elev(end),'r.','markersize',8)
xlabel('x (km)','fontsize',16);
ylabel('Ice Height (m)','fontsize',16)

% subplot(3,2,3)
% % plot(time,parameters.buttress,'k.');hold on
% plot(time,x_g,'k.');hold on
% xlabel('Years','fontsize',16);
% ylabel('GL Vel (m/yr)','fontsize',16)

subplot(3,2,4)
p4 = plot(parameters.grid.sigma_element,parameters.taub,'r','linewidth',3);hold on
xlabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
ylabel('$\tau_b$','Interpreter','LaTeX','fontsize',16)

% subplot(3,2,3)
% plot(time,x_g/1e3,'k.');hold on
% xlabel('Years','fontsize',16);
% ylabel('GL Position (km)','fontsize',16)

% subplot(3,2,5)
% pcolor(x_element./1000,z_element,[T,T_s]);shading('flat');hold on
% plot(x_g.*parameters.grid.sigma_node/1000,bed_elev,'k','linewidth',2);
% plot([x_g.*parameters.grid.sigma_node;x_g]/1000,[bed_elev+[h;h_g];bed_elev(end)],'m','linewidth',2);
% xlabel('x (km)','fontsize',16);
% ylabel('z (m)','fontsize',16)
% caxis([min(parameters.T_s) 0]);
% colorbar
% hold off

subplot(3,2,6)
p6 = plot(parameters.grid.sigma_element,e.*h_till,'r','linewidth',3);hold on
xlabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
ylabel('Till Water Content (m)','fontsize',16)

% subplot(3,2,5)
% p5 = plot(parameters.grid.sigma_element,h_till,'r');hold on
% xlabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
% ylabel('Till thickness (m)','fontsize',16)
subplot(3,2,5)
plot(time,x_g./1e3,'b.');hold on
xlabel('Time (yr)','Interpreter','LaTeX','fontsize',16);
ylabel('GL Pos (km)','fontsize',16)

drawnow

else

[SIG,ETA] = meshgrid(parameters.grid.sigma_element,parameters.grid.eta_element);
TempDiag
figure(2);set(2,'units','normalized','position',[1 0.1 0.5 0.5]);
subplot(6,1,1);contourf(SIG',ETA',T-parameters.MP);colorbar;
subplot(6,1,2);contourf(SIG',ETA',-parameters.dtau.*(Uadv));colorbar;%caxis(parameters.dtau.*[-1e-9 1e-9])
subplot(6,1,3);contourf(SIG',ETA',-parameters.dtau.*(Vadv));colorbar;%caxis(parameters.dtau.*[-1e-9 1e-9])
subplot(6,1,4);contourf(SIG',ETA',-parameters.dtau.*Ustretch);colorbar;%caxis(parameters.dtau.*[-1e-9 1e-9])
subplot(6,1,5);contourf(SIG',ETA',-parameters.dtau.*Ustretch2);colorbar;%caxis(parameters.dtau.*[-1e-9 1e-9])
subplot(6,1,6);contourf(SIG',ETA',-parameters.dtau.*Fdiff);colorbar;

end