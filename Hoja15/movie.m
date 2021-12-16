% fig = figure(4)
% anim = plot([0,0],[5,y2(1)+y1(1)],'-k',[0,y1(1)+y2(1)],[0,y2(1)],'-k',0,y1(1),'*',0,-10+y1(1)+y2(1),'*');
% grid on; xlim([-0.5 0.5]); ylim([-20 5]);
% xlabel('x (m)'), ylabel('y (m)');
% text = annotation('textbox', [0.1, 0.8, 0.3, 0.03], 'String', strcat('t = ',32,num2str(0),' s'), ...
%     'FontSize', 15, 'FitBoxToText','off','Interpreter','latex', 'BackgroundColor', 'w');
% grid off;
% tic;
% for i = 2:numel(t)
%     set(anim(3), 'XData',0, 'YData',y1(i));
%     set(anim(4), 'XData',0, 'YData',-10+y1(i)+y2(i));
%     set(anim(1), 'XData',[0,0], 'YData',[5,y1(i)]);
%     set(anim(2), 'XData',[0,0], 'YData',[y1(i),-10+y1(i)+y2(i)]);
%     set(text, 'String', strcat('t = ',32,num2str(t(i),'%10.1f'),' s'));
%     drawnow;
%     frame = getframe(fig);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     imwrite(imind,cm,strcat('movie/',num2str(i-1),'.png'));
% end
% toc;