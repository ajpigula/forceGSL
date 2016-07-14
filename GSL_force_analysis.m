
function [sweep_tem,sweep_marked] = GSL_force_analysis

clear
warning off

load('/home/sisir/code/forceGSL/directories/GSL_crosssection.mat')

GSLs = [];

% for j = 1:length(d)
%     [GSLs(j,:),sweep_tem,sweep_marked] = analyzeGSLs(d(j)); 
% end

% absChange = GSLs - repmat(GSLs(:,1),1,18);
% percChange = (GSLs - repmat(GSLs(:,1),1,18))./repmat(GSLs(:,1),1,18);

GSLs = zeros(57,18);
absChange = zeros(57,18);
percChange = zeros(57,18);

for i = 1:57
    
    GSLs(i,:) = d(i).GSL;
    absChange(i,:) = d(i).absChange;
    percChange(i,:) = d(i).percChange;
        
end

plotSummary(1,GSLs(1:30,:))
title('DMD')
xlabel('Force (N)')
ylabel('GSL')
axis([0 12 .3 .8])

plotSummary(2,GSLs(31:57,:))
title('Healthy')
xlabel('Force (N)')
ylabel('GSL')
axis([0 12 .3 .8])

plotSummary(3,absChange(1:30,:))
title('DMD')
xlabel('Force (N)')
ylabel('GSL, change')

plotSummary(4,absChange(31:57,:))
title('Healthy')
xlabel('Force (N)')
ylabel('GSL, change')


% plotErrorBars(1,GSLs)
% suptitle('Raw data')
% plotErrorBars(2,absChange)
% suptitle('Absolute change')
% plotErrorBars(3,percChange)
% suptitle('Percent change')
% 
% figure(5);
% 
% for k = 31:57
%     plotScatter(4,d(k).age,GSLs(k,18));
%     plotScatter(7,d(k).age,absChange(k,18));
%     plotScatter(10,d(k).age,percChange(k,18));
%     
%     plotScatter(5,d(k).th0,GSLs(k,18));
%     plotScatter(8,d(k).th0,absChange(k,18));
%     plotScatter(11,d(k).th0,percChange(k,18));
%     
%     plotScatter(6,d(k).mw6,GSLs(k,18));
%     plotScatter(9,d(k).mw6,absChange(k,18));
%     plotScatter(12,d(k).mw6,percChange(k,18));
%     
% end

% figure(4); xlabel('Age'); ylabel('GSL'); title('DMD')
% figure(5); xlabel('Age'); ylabel('Abs change'); title('DMD')
% figure(6); xlabel('Age'); ylabel('Perc change'); title('DMD')
% 
% figure(7); xlabel('Th0'); ylabel('GSL'); title('DMD')
% figure(8); xlabel('Th0'); ylabel('Abs change'); title('DMD')
% figure(9); xlabel('Th0'); ylabel('Perc change'); title('DMD')
% 
% figure(10); xlabel('6MW'); ylabel('GSL'); title('DMD')
% figure(11); xlabel('6MW'); ylabel('Abs change'); title('DMD')
% figure(12); xlabel('6MW'); ylabel('Perc change'); title('DMD')
% 


end

function A = depthCorrection(sweep)

for k = 1:18
    for i = 1:411
        D = i*50/411;
        for j = 1:315
            A(i,j,k) = sqrt(sweep(i,j,k)*D);
        end
    end
end
    
end


function plotErrorBars(f,data)

figure(f); subplot(1,2,1); hold on;
plot(1.5:.5:10,data(1:30,:)','r');
plot(1.5:.5:10,data(31:57,:)','b');

subplot(1,2,2); hold on;
errorbar([1.5:.5:10]',mean(data(1:30,:),1)',std(data(1:30,:),1)','Color','r','LineWidth',2,'LineStyle','--','Marker','.');
errorbar([1.5:.5:10]',mean(data(31:57,:),1)',std(data(31:57,:),1)','Color','b','LineWidth',2,'LineStyle','--','Marker','.');

end


function plotSummary(f,data)

figure(f); hold on;

errorbar([1.5:.5:10]',mean(data,1)',std(data,1)','Color','k','LineWidth',2,'Marker','.');
plot([1.5:.5:10]',max(data,[],1),'k*')
plot([1.5:.5:10]',min(data,[],1),'k*')

end


function [GSLs,sweep_tem,sweep_marked] = analyzeGSLs(entry)

sweep = entry.sweep;

A = depthCorrection(sweep);

yc = entry.ellipse(:,1);
xc = entry.ellipse(:,2);
a = entry.ellipse(:,3);
b = entry.ellipse(:,4);
alpha = entry.ellipse(:,5);
interMus = entry.interMus;
subQ = entry.subQ;
GSLs = [];
sweep_marked = [];

for i = 1:18
%     [GSLs(i),sweep_marked(:,:,i)] = getGSL(sweep(:,:,i),xc(i),yc(i),a(i),b(i),alpha(i),interMus(i),subQ(i));
    [~,sweep_tem(:,:,i)] = getGSL(sweep(:,:,i),xc(i),yc(i),a(i),b(i),alpha(i),interMus(i),subQ(i));
    [GSLs(i),sweep_marked(:,:,i)] = getGSL(A(:,:,i),xc(i),yc(i),a(i),b(i),alpha(i),interMus(i),subQ(i));
end

implay(sweep_tem)
implay(sweep_marked/max(max(max(sweep_marked))))

end


function [gsl,img] = getGSL(img,xc,yc,a,b,alpha,interMus,subQ)

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

% figure(f);
subplot(3,3,f-3)
hold on;
plot(patientData,observedData,'o')

end

