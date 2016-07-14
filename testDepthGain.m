function [slope_raw, slope_bg, slope_profile,sweep_marked,depth_marked,GSLs_profile,top,bot,profile] = testDepthGain

load('/home/sisir/code/forceGSL/BM_background/BMgain.mat')
load('/home/sisir/code/forceGSL/directories/GSL_crosssection.mat')


GSLs = [];
bgsweep = [];

for h = 1:18
    bgsweep(:,:,h) = BMgain/max(max(BMgain));
end

for j = 1:length(d)
%     [GSLs(j,:),sweep_marked,GSLs_depth(j,:),depth_marked] = analyzeGSLs(d(j),bgsweep);
    [slope_raw(j), slope_bg(j), slope_profile(j),sweep_marked,depth_marked,GSLs_profile,top,bot,profile] = analyzeGSLs(d(j),bgsweep);
end



corrcoef(slope_raw,slope_bg)
corrcoef(slope_raw,slope_profile)
corrcoef(slope_bg,slope_profile)















% 
% figure; plot(GSLs(:,:)')
% figure; plot(GSLs_depth(:,:)')
% figure; plot(GSLs_adj(:,:)')
% mean(GSLs(:,18)-GSLs(:,1));
% mean(GSLs_depth(:,18)-GSLs_depth(:,1));
% mean(GSLs_adj(:,18)-GSLs_adj(:,1));


% age = [];
% th0 = [];
% mw6 = [];
% hasmw6 = [1 2 4 6 7 8 9 10 19 22 24 26 27];



% for k = 1:57
%     
%     age(end+1) = d(k).age;
%     th0(end+1) = d(k).th0;
%     mw6(end+1) = d(k).mw6;
    
%     plotScatter(5,d(k).age,(GSLs(k,18)-GSLs(k,1))/GSLs(k,1));
%     plotScatter(6,d(k).th0,(GSLs(k,18)-GSLs(k,1))/GSLs(k,1));
%     plotScatter(7,d(k).mw6,(GSLs(k,18)-GSLs(k,1))/GSLs(k,1));
%     plotScatter(8,d(k).age,(GSLs_adj(k,18)-GSLs_adj(k,1))/GSLs_adj(k,1));
%     plotScatter(9,d(k).th0,(GSLs_adj(k,18)-GSLs_adj(k,1))/GSLs_adj(k,1));
%     plotScatter(10,d(k).mw6,(GSLs_adj(k,18)-GSLs_adj(k,1))/GSLs_adj(k,1));
    
% end

% figure; plot(age(1:30),GSLs(1:30,1),'ro')
% hold on; plot(age(31:57),GSLs(31:57,1),'bx')
% 
% figure; plot(age(1:30),GSLs_adj(1:30,1),'ro')
% hold on; plot(age(31:57),GSLs_adj(31:57,1),'bx')
% 
% corr(age(1:30)',GSLs(1:30,1))
% corr(age(31:57)',GSLs(31:57,1))
% 
% corr(age(1:30)',GSLs_adj(1:30,1))
% corr(age(31:57)',GSLs_adj(31:57,1))

% corr(GSLs(31:57,1),age')
% figure; plot(age,GSLs(31:57,1),'o')
% 
% corr(GSLs(31:57,18),age')
% figure; plot(age,GSLs(31:57,18),'o')
% 
% corr(GSLs_adj(31:57,1),age')
% figure; plot(age,GSLs(31:57,1),'o')
% 
% corr(GSLs_adj(31:57,18),age')
% figure; plot(age,GSLs(31:57,18),'o')



end




% function [GSLs,sweep_marked,GSLs_depth,depth_marked] = analyzeGSLs(entry,bgsweep)
function [p1, p2, p3,sweep_marked,depth_marked,GSLs_profile,top,bot,profile] = analyzeGSLs(entry,bgsweep)

sweep = entry.sweep;
yc = entry.ellipse(:,1);
xc = entry.ellipse(:,2);
a = entry.ellipse(:,3);
b = entry.ellipse(:,4);
alpha = entry.ellipse(:,5);
interMus = entry.interMus;
subQ = entry.subQ;
GSLs = [];
sweep_marked = [];
top = [];
bot = [];

for i = 1:18
    [GSLs(i),sweep_marked(:,:,i),top(i),bot(i)] = getGSL(sweep(:,:,i),xc(i),yc(i),a(i),b(i),alpha(i),interMus(i),subQ(i));
    [GSLs_depth(i),depth_marked(:,:,i),~,~] = getGSL(bgsweep(:,:,i),xc(i),yc(i),a(i),b(i),alpha(i),interMus(i),subQ(i));
%     [GSLs_adj(i),adj_marked(:,:,i),~,~] = getGSL(sweep(:,:,i)./bgsweep(:,:,i),xc(i),yc(i),a(i),b(i),alpha(i),interMus(i),subQ(i));
end

for i = 1:411
    profile(i) = mean2(sweep(i,:,:));
end
for i = 1:18
    GSLs_profile(i) = median(profile(top(i):bot(i)));
end
 
% figure; hold on
% plot([1.5:.5:10],GSLs,'k')
% plot([1.5:.5:10],GSLs_depth,'g')
% plot([1.5:.5:10],GSLs_profile,'b')

p1 = polyfit([1.5:.5:10],GSLs,1);
p1 = p1(1);
p2 = polyfit([1.5:.5:10],GSLs_depth,1);
p2 = p2(1);
p3 = polyfit([1.5:.5:10],GSLs_profile,1);
p3 = p3(1);


end



function [gsl,img,top,useDepth] = getGSL(img,xc,yc,a,b,alpha,interMus,subQ)

npts = 800;
t = linspace(0, 2*pi, npts);
Q = [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];
X = Q * [a * cos(t); b * sin(t)] + repmat([yc;xc], 1, npts);

Y = round(X(1,:));
X = round(X(2,:));

Y = max(min(Y,interMus),subQ);
X = max(min(X,315),1);

for i = 1:npts
    img(Y(i),X(i)) = 0;
end

bottom = max(Y);
top = min(Y);
useDepth = round(top + (bottom-top)/3);

ROIcontents = [];

for yy = top:useDepth
    rowIndex = find(Y == yy);
    bookendLeft = min(X(rowIndex));
    bookendRight = max(X(rowIndex));
    
%     img(yy,bookendLeft:bookendRight) = 1;

    ROIcontents = [ROIcontents img(yy,bookendLeft:bookendRight)];
    
end

gsl = median(ROIcontents);

end



function plotScatter(f,patientData,observedData)

figure(f);
hold on;
plot(patientData,observedData,'o')

end