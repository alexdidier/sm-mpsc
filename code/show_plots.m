%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linear SM-MPSC Plots
%
% Alexandre Didier
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure(1)
hold on
xlabel('$p_{k,1}$','Interpreter','latex')
ylabel('$p_{k,2}$','Interpreter','latex')
zlabel('$p_{k,3}$','Interpreter','latex')
k=10;
p2=S_feas{k+1}.plot('Color','blue','alpha',0.2);
p3=S_feas{k+2}.plot('Color','gray','alpha',0.4);
p1=S_feas{1}.plot('Color','red','alpha',1);
legend([p1,p2,p3],{'$X_{\textup{feas}}(\Theta_0)$','$X_{\textup{feas}}(\Theta_{150})$','$X_{\textup{feas,f}}(\Theta_{150})$'},'Interpreter','latex','FontSize',11)
hold off

%%
figure(2)
colour=[0.6350 0.0780 0.1840];

subplot(2,1, [2])
hold on
u_diff_lin_norm = [];
for i=1:length(u_diff_lin)
    u_diff_lin_norm(i) = norm(u_diff_lin(:,i));
end
u_diff_norm = [];
for i=1:length(u_diff)
    u_diff_norm(i) = norm(u_diff(:,i));
end
bar(plot_t_lin, u_diff_norm,'b', 'FaceAlpha', 0.5)
bar(plot_t_lin, u_diff_lin_norm,'r', 'FaceAlpha', 0.5)
ax=gca;
ax.XAxis.FontSize=10;
ax.YAxis.FontSize=10;
xlabel('Time [s]','FontSize',11,'Interpreter','latex')
axis([0 time 0 7.5])
hold off

i=1;
subplot(2,1,[1])
hold on
p4=plot(plot_t_lin,plot_x_lin(2*i-1,:),'r','LineWidth',2);


p3=plot(plot_t, plot_x(2*i-1,:),'b','LineWidth',2);
p1=plot(plot_t,x_ref_PGSD(2*i-1,:),'k--');
p2=plot(plot_t,2.3*ones(size(plot_t)),'k', 'LineWidth', 1.5 );
plot(plot_t,-2.3*ones(size(plot_t)),'k', 'LineWidth', 1.5 );

axis([0 time -3 3])
ax=gca;
ax.XAxis.FontSize=10;
ax.YAxis.FontSize=10;
legend([p1,p2,p3, p4],{'Reference','Constraints','Adaptive MPSC','MPSC in [6]'},'Interpreter','latex','FontSize',16); 
set(ax.XAxis,'visible','off' );
hold off

%%
figure(3)
fontsize=20;

hold on
p1=Polyhedron(Omega.A, plot_Omega(:,1)).plot('Color','gray','alpha',0.1);
p2=Polyhedron(Omega.A, plot_Omega(:,5)).plot('Color','gray','alpha',0.3);
p3=Polyhedron(Omega.A, plot_Omega(:,11)).plot('Color','gray','alpha',0.5);
p4=Polyhedron(Omega.A, plot_Omega(:,151)).plot('Color','gray','alpha',0.7);
for i=[1 5 11 151]
p6=plot((plot_Omega(1,i)-plot_Omega(2,i))/2,(plot_Omega(3,i)-plot_Omega(4,i))/2,'b.','MarkerSize',10);
end
p5=plot(k_true(1),k_true(2),'r*','MarkerSize',10,'linewidth',0.3);
legend([p5,p6,p1,p2,p3,p4],{'$\theta^*$','$\bar{\theta}_k$','$\Theta_0$','$\Theta_{4}$','$\Theta_{10}$','$\Theta_{150}$'},'Interpreter','latex','FontSize',fontsize,'Location','nw')

axis([0.02 0.28 0.02 0.28])
set(gca,'FontSize',fontsize-8)
xlabel('$\theta_1$','Interpreter','latex','FontSize',fontsize);
ylabel('$\theta_2$','Interpreter','latex','FontSize',fontsize);
set(gcf,'position',[0,0,600,400])
hold off
