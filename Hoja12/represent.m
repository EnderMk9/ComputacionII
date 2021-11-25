% figure(1);
% plot(t,x,t,dx,t,ddx);
% xlabel('time (s)'); ylabel('x (m)')
% legend('$x$','$\dot{x}$','$\ddot{x}$','Interpreter','latex');
% figure(2);
% plot(t,log10(x),t,log10(dx),t,log10(ddx));
% xlabel('time (s)'); ylabel('log_{10}(x (m))')
% legend('$log_{10}(x)$','$log_{10}(\dot{x})$','$log_{10}(\ddot{x})$','Interpreter','latex');
figure(1);
plot(t,log10(ddiff));
xlabel('time (s)'); ylabel('log_{10}(error)')
legend('$\log_{10}(\Delta \dot{x})$','Interpreter','latex');
figure(2);
plot(t,log10(dddiff));
xlabel('time (s)'); ylabel('log_{10}(error)')
legend('$\log_{10}(\Delta \ddot{x})$','Interpreter','latex');
