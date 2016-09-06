T_s = parameters.T_s;
T_lin=zeros(parameters.grid.n_elements,parameters.grid.n2_nodes);
for i=1:parameters.grid.n_elements
    T_lin(i,:) = linspace(0,T_s(i),parameters.grid.n2_nodes);
end
T_lin = (T_lin(:,1:parameters.grid.n2_nodes-1)+T_lin(:,2:parameters.grid.n2_nodes))./2;
T_init = squeeze(T_all(:,:,1));
   
ss = 1; %sumbsample time
% tbegin=8000;
% tend=1.9e4;
% find(time_all>2e4,1);
% tbegin = 9386;
% tend = 12566;
% tbegin=15260;
% tend=17930;
% tbegin=find(time_all>22.6e3,1);
% tend=length(time_all)-37;
tbegin=find(time_all>1.5e4,1);
tend=length(time_all);

sigbegin = 1;
sigend = parameters.grid.n_nodes;

%%
figure(1);set(1,'units','normalized','position',[0 0.1 0.8 0.75]);
subplot(2,2,1)
pcolor(time_all(tbegin:ss:tend)./1000,parameters.grid.sigma_element(sigbegin:sigend-1),h_all((sigbegin:sigend-1),tbegin:ss:tend));hold on
shading('flat')
colorbar;
xlabel('time (kyr)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
set(gca,'YDir','Reverse','fontsize',20)
title('Ice Thickness (m)','fontsize',16);

subplot(2,2,2)
pcolor(time_all(tbegin:ss:tend)./1000,parameters.grid.sigma_element(sigbegin:sigend-1),squeeze(mean(T_all((sigbegin:sigend-1),:,tbegin:ss:tend),2)));hold on
shading('flat')
colorbar;
% caxis([-10 -2])
set(gca,'YDir','Reverse','fontsize',20)
xlabel('time (kyr)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Mean Column Temperature (deg C)','fontsize',16)

subplot(2,2,3)
pcolor(time_all(tbegin:ss:tend)./1000,parameters.grid.sigma_node(sigbegin:sigend),parameters.year.*squeeze(u_full_all((sigbegin:sigend),1,tbegin:ss:tend)));hold on
shading('flat')
colorbar;
caxis([0 1000])
set(gca,'YDir','Reverse','fontsize',20)
xlabel('time (kyr)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Velocity (m/yr)','fontsize',16);

subplot(2,2,4)
pcolor(time_all(tbegin:ss:tend)./1000,parameters.grid.sigma_element(sigbegin:sigend-1),e_all((sigbegin:sigend-1),tbegin:ss:tend));hold on
shading('flat')
colorbar;
% caxis([0.4 0.8])
set(gca,'YDir','Reverse','fontsize',20)
xlabel('time (kyr)','fontsize',16);
ylabel('$\sigma$','Interpreter','LaTeX','fontsize',16);
title('Void Ratio','fontsize',16);

figure(2)
plot(time_all(tbegin:1:tend)./1000,xg_all(tbegin:1:tend)/1000,'k','linewidth',4);hold on
xlabel('time (kyr)','fontsize',16);
ylabel('GL Position (km)','fontsize',16);
xlim([time_all(tbegin)./1000 time_all(tend)./1000])
set(gca,'fontsize',20)

% figure(5)
% plot(tbegin:ss:tend,u_all(end,tbegin:ss:tend)/1000,'k.');hold on
% xlabel('time (years)','fontsize',16);
% ylabel('GL Position (km)','fontsize',16);

%%
close all
[SIG,ETA] = meshgrid(parameters.grid.sigma_element,parameters.grid.eta_element);
figure(3);
set(3,'units','normalized','position',[0.1 0.5 0.6 0.4]);

q=0;
timeset=0;

for t=tbegin:1:tend
        
    timeset=timeset+diff(time_all(t-1:t));
    if(timeset<10);
        continue
    else
        timeset=0;
    end
    q=q+1;
    

    x_g = xg_all(t);
%     T = squeeze(T_all(:,:,t));
    h = h_all(:,t);
    T_s = parameters.T_s;

    h_g=-(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
    bed_elev = Base(x_g.*parameters.grid.sigma_node,parameters);

    z_element=(parameters.grid.eta_node*(h'))';
    z_element=z_element+repmat(Base(x_g.*parameters.grid.sigma_element,parameters),[1 parameters.grid.n2_nodes]);
    x_element=repmat(x_g*parameters.grid.sigma_element,[1 parameters.grid.n2_nodes]);
    
%     subplot(2,1,1)
    
%     subplot(3,1,1:2)
%     pcolor(x_element./1000,z_element,[T,T_s]);shading('flat');hold on
    plot(x_g.*parameters.grid.sigma_node/1000,bed_elev,'k','linewidth',5);hold on
    plot([x_g.*parameters.grid.sigma_node;x_g]/1000,[bed_elev+[h;h_g];bed_elev(end)],'k','linewidth',2);


%     plot(xg_all(t).*parameters.grid.sigma_element/1000,bed_elev,'k','linewidth',2);hold on
%     plot([xg_all(t).*parameters.grid.sigma_element]/1000,[bed_elev+h_all(:,t)],'m','linewidth',2);hold on
%     pcolor(x_element./1000,z_element,[squeeze(T_all(:,:,t)),parameters.T_s]);
    shading('flat')
    xlabel('x (km)','fontsize',26);
    ylabel('z (m)','fontsize',26)
    caxis([min(parameters.T_s) 0]);
    set(gca,'fontsize',26,'linewidth',2)
%     caxis([-4 0]);
    axis([0 max(xg_all)/1000 min(Base(xg_all,parameters)) parameters.icedivide+max(max(h_all))])
    colormap('redblue')
    h=colorbar;
    set(h,'fontsize',26);
    title(['Time =',num2str(time_all(t)),' years'],'fontsize',32)
    drawnow
    
%     set(3,'units','pixels','position',[0 0 1600 961])
%     saveSameSize(3, 'format', 'png', 'renderer', 'opengl','file',['/Users/alexanderrobel/ResearchDataFiles/FlowlineExpOutput/AnimationTest/FL_frame_',sprintf('%3.3d',q),'.png']);


    clf
end