load('dat.mat');
plot(xs0,it0,'*',xs1,it1,'*',xs2,it2,'*',-4.5641,1,'+k',3.432,1,'+k',-0.25907,1,'+k');
legend('root 0', 'root 1', 'root 2');
xlabel("x"); ylabel("iterations")