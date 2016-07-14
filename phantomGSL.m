

% % Trial 1: Force sweeps (first 2 may be bad)
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid/volume.mat')
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid/volume_forcedata.mat')
% 
% for i = 1:size(volume,3)
%     gsl1(i) = mean2(volume(10:60,:,i));
% end
% 
% figure(1)
% subplot(3,2,1); hold on;
% plot(volumeForces(:,2),'k','LineWidth',2)
% plot(30*gsl1,'b','LineWidth',2)
% legend('Force','Mean')
% title('1')
% 
% r = corrcoef(volumeForces(:,2),gsl1);
% fprintf('\n1: %.2f',r(2,1))


% Trial 2: Force sweeps
load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid2/volume.mat')
load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid2/volume_forcedata.mat')
load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid2/xy.mat')

for i = 1:size(volume,3)
    gsl2(i) = mean2(volume(y(i):y(i)+50,x(i)-50:x(i)+50,i));
end

subplot(3,2,3); hold on;
plot(volumeForces(:,2),'k','LineWidth',2)
plot(30*gsl2,'b','LineWidth',2)
legend('Force','Mean')
title('2')

r = corrcoef(volumeForces(:,2),gsl2);
fprintf('\n2: %.2f',r(2,1))

figure(2); subplot(2,1,1);
plot(volumeForces(:,1),gsl2,'k','LineWidth',2);
ylabel('GSL')
axis([0 20 .2 .4]);
subplot(2,1,2)
plot(volumeForces(:,1),volumeForces(:,2),'k','LineWidth',2);
xlabel('Time (s)')
ylabel('Force (N)')
axis([0 20 0 15]);
suptitle('Phantom experiment')
% 
% % Trial 3: Force sweeps
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid3/volume.mat')
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid3/volume_forcedata.mat')
% 
% for i = 1:size(volume,3)
%     gsl3(i) = mean2(volume(10:60,:,i));
% end
% 
% figure(1);
% subplot(3,2,5); hold on;
% plot(volumeForces(:,2),'k','LineWidth',2)
% plot(30*gsl3,'b','LineWidth',2)
% legend('Force','Mean')
% title('3')
% 
% r = corrcoef(volumeForces(:,2),gsl3);
% fprintf('\n3: %.2f\n',r(2,1))
% 
% 
% % Trial 4: Indirectly applied force, with paper plate
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid4/volume.mat')
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid4/volume_forcedata.mat')
% 
% for i = 1:size(volume,3)
%     gsl4(i) = mean2(volume(10:60,:,i));
% end
% 
% subplot(3,2,2); hold on;
% plot(volumeForces(:,2),'k','LineWidth',2)
% plot(30*gsl4,'b','LineWidth',2)
% legend('Force','Mean')
% title('Plate - disregard data')
% 
% 
% % Trial 5: Indirectly applied force
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid5/volume.mat')
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid5/volume_forcedata.mat')
% 
% for i = 1:size(volume,3)
%     gsl5(i) = mean2(volume(10:60,:,i));
% end
% 
% subplot(3,2,4); hold on;
% plot(volumeForces(:,2),'k','LineWidth',2)
% plot(30*gsl5,'b','LineWidth',2)
% legend('Force','Mean')
% title('5')
% 
% 
% % Trial 6: Indirectly applied force
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid6/volume.mat')
% load('/home/sisir/code/data/DMD_Complete_20Aug2015/PhantomCompressions_18May2016/deltoid6/volume_forcedata.mat')
% 
% for i = 1:size(volume,3)
%     gsl6(i) = mean2(volume(10:60,:,i));
% end
% 
% subplot(3,2,6); hold on;
% plot(volumeForces(:,2),'k','LineWidth',2)
% plot(30*gsl6,'b','LineWidth',2)
% legend('Force','Mean')
% title('6')
% 
% 
