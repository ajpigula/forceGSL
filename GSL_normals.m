
% function [atForces_thick,atForces_fixed,atForces_frac_20,atForces_frac_33,atForces_frac_50,percChanges_fixed,percChanges_20,percChanges_33,percChanges_50] = GSL_normals
function th0 = GSL_normals


clear
warning off

% load('/home/sisir/code/dmd/directories/normals_crosssection.mat')
% load('/home/sisir/code/dmd/directories/DMD_GSL_sample.mat')
load('/home/sisir/code/dmd/directories/GSL_crosssection.mat')

basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/';
runID1 = 'Sweeps_22Mar2016/';
runID2 = 'Compression_Jan2016_quads/';

figure(6); hold on;
plot([.3 .6],[.3 .6],'k--')
figure(7); hold on;
plot([.3 .6],[.3 .6],'k--')
figure(8); hold on;
plot([.3 .6],[.3 .6],'k--')

atForces_fixed = [];
atForces_frac_20 = [];
atForces_frac_33 = [];
atForces_frac_50 = [];
atForces_thick = [];
percChanges_fixed = [];
percChanges_20 = [];
percChanges_33 = [];
percChanges_50 = [];
age = [];
mw = [];
th0 = [];

for i = 1:length(d)
    
    folder = strcat(d(i).name,'_',d(i).visit);
    tree = dir(fullfile(basePath,runID1,folder));
    
    idx = [];
    for j = 3:length(tree)
        p = findstr(tree(j).name,'_');
        thisMuscle = tree(j).name(p(1)+1:p(end)-1);
        if length(d(i).muscle) > length(thisMuscle)
            continue
        end
        if strcmp(thisMuscle(1:length(d(i).muscle)),d(i).muscle)
            load(fullfile(basePath,runID1,folder,tree(j).name))
        end
    end
    
    load(fullfile( basePath,runID2, strcat(d(i).name,'_',d(i).muscle,'_',d(i).visit) ))
    
    thickness = interMus(:,1)-subQ(:,1);
    fprintf('Final thickness: %i\n',thickness(end))
    
    th0(end+1) = interpolateTh0(thickness,sweepForces(1:length(thickness),2));
    
    if thickness(end) > 45
        
        fprintf('%s\n',d(i).name)
        
        GSLs_ROI = placeROI_fixed(sweep,subQ,0);
        if strcmp(d(1).name,'01001')
            GSLs_ROI = placeROI_fixed(sweep,subQ,15);
        end
        if strcmp(d(1).name,'01008')
            GSLs_ROI = placeROI_fixed(sweep,subQ,5);
        end
        if strcmp(d(1).name,'01017')
            GSLs_ROI = placeROI_fixed(sweep,subQ,10);
        end
        change_ROI = 100*(GSLs_ROI-GSLs_ROI(1))/GSLs_ROI(1);
        
        GSLs_frac_20 = placeROI_frac(sweep,subQ,interMus,.2);
        change_frac_20 = 100*(GSLs_frac_20-GSLs_frac_20(1))/GSLs_frac_20(1);
        GSLs_frac_33 = placeROI_frac(sweep,subQ,interMus,.33);
        change_frac_33 = 100*(GSLs_frac_33-GSLs_frac_33(1))/GSLs_frac_33(1);
        GSLs_frac_50 = placeROI_frac(sweep,subQ,interMus,.5);
        change_frac_50 = 100*(GSLs_frac_50-GSLs_frac_50(1))/GSLs_frac_50(1);
        
        atForces_fixed(end+1,:) = atForces(GSLs_ROI,sweepForces);
        atForces_frac_20(end+1,:) = atForces(GSLs_frac_20,sweepForces);
        atForces_frac_33(end+1,:) = atForces(GSLs_frac_33,sweepForces);
        atForces_frac_50(end+1,:) = atForces(GSLs_frac_50,sweepForces);
        
        atForces_thick(end+1,:) = atForces(thickness,sweepForces);
        
        percChanges_fixed(end+1,:) = atForces(change_ROI,sweepForces);
        percChanges_20(end+1,:) = atForces(change_frac_20,sweepForces);
        percChanges_33(end+1,:) = atForces(change_frac_33,sweepForces);
        percChanges_50(end+1,:) = atForces(change_frac_50,sweepForces);
                
      
        age(end+1) = d(i).age;
        mw(end+1) = d(i).mw6;
%         th0(end+1) = interpolateTh0(thickness,sweepForces(1:length(thickness),2));
        

        plotvForce(1,GSLs_ROI,sweepForces);
        plotvForce(2,change_ROI,sweepForces);
        
        plotvForce(3,.2*thickness,sweepForces);
        plotvForce(4,.33*thickness,sweepForces);
        plotvForce(5,.5*thickness,sweepForces);
        
        plotROIvfrac(6,GSLs_ROI,GSLs_frac_20)
        plotROIvfrac(7,GSLs_ROI,GSLs_frac_33)
        plotROIvfrac(8,GSLs_ROI,GSLs_frac_50)
        
        plotvForce(9,GSLs_frac_20,sweepForces);
        plotvForce(10,GSLs_frac_33,sweepForces);
        plotvForce(11,GSLs_frac_50,sweepForces);
        
        plotvForce(12,change_frac_20,sweepForces);
        plotvForce(13,change_frac_33,sweepForces);
        plotvForce(14,change_frac_50,sweepForces);
        
       

    else
        fprintf('Muscle too thin\n\n')
    end
    
end

fprintf('\n')

figure(1); 
errorbar([1.5:.5:10],mean(atForces_fixed,1),std(atForces_fixed,1)/sqrt(19),'ko-','LineWidth',2)
title('Fixed-size ROI')
xlabel('Applied force (N)')
ylabel('GSL')
set(gca,'fontsize', 12)


figure(2); 
errorbar([1.5:.5:10],mean(percChanges_fixed,1),std(percChanges_fixed,1)/sqrt(19),'ko-','LineWidth',2)
title('Fixed ROI, % change')
xlabel('Applied force (N)')
ylabel('% change')
set(gca,'fontsize', 12)


figure(3); title('ROI height, d = 20%')
xlabel('Applied force (N)')
ylabel('ROI height (pixels)')
set(gca,'fontsize', 12)


figure(4); title('ROI height, d = 33%')
xlabel('Applied force (N)')
ylabel('ROI height (pixels)')
set(gca,'fontsize', 12)


figure(5); title('ROI height, d = 50%')
xlabel('Applied force (N)')
ylabel('ROI height (pixels)')
set(gca,'fontsize', 12)


figure(6); title('GSL: fixed-size ROI vs d = 20%')
xlabel('Fixed')
ylabel('d = 20')
set(gca,'fontsize', 12)


figure(7); title('GSL: fixed-size ROI vs d = 33%')
xlabel('Fixed')
ylabel('d = 33')
set(gca,'fontsize', 12)


figure(8); title('GSL: fixed-size ROI vs d = 50%')
xlabel('Fixed')
ylabel('d = 50')
set(gca,'fontsize', 12)

figure(9); title('d = 20')
errorbar([1.5:.5:10],mean(atForces_frac_20,1),std(atForces_frac_20,1)/sqrt(19),'ko-','LineWidth',2)
xlabel('Applied force (N)')
ylabel('GSL')
set(gca,'fontsize', 12)

figure(10); 
errorbar([1.5:.5:10],mean(atForces_frac_33,1),std(atForces_frac_33,1)/sqrt(19),'ko-','LineWidth',2)
title('d = 33')
xlabel('Applied force (N)')
ylabel('GSL')
set(gca,'fontsize', 12)

figure(11); 
errorbar([1.5:.5:10],mean(atForces_frac_50,1),std(atForces_frac_50,1)/sqrt(19),'ko-','LineWidth',2)
title('d = 50')
xlabel('Applied force (N)')
ylabel('GSL')
set(gca,'fontsize', 12)

figure(12); 
errorbar([1.5:.5:10],mean(percChanges_20,1),std(percChanges_20,1)/sqrt(19),'ko-','LineWidth',2)
title('d = 20')
xlabel('Applied force (N)')
ylabel('% change')
set(gca,'fontsize', 12)

figure(13); 
errorbar([1.5:.5:10],mean(percChanges_33,1),std(percChanges_33,1)/sqrt(19),'ko-','LineWidth',2)
title('d = 33')
xlabel('Applied force (N)')
ylabel('% change')
set(gca,'fontsize', 12)

figure(14); 
errorbar([1.5:.5:10],mean(percChanges_50,1),std(percChanges_50,1)/sqrt(19),'ko-','LineWidth',2)
title('d = 50')
xlabel('Applied force (N)')
ylabel('% change')
set(gca,'fontsize', 12)


corrs(15, atForces_fixed, atForces_frac_20, atForces_frac_33, atForces_frac_50);


figure(16);
plot(age,percChanges_fixed(:,end),'o')
xlabel('Age')
ylabel('Total % change')
r = corr(age',percChanges_fixed(:,end));
fprintf('Correlation with age: %f\n',r)

% figure(17);
% plot(th0,percChanges_fixed(:,end),'o')
% xlabel('Original thickness (th0)')
% ylabel('Total % change')
% r = corr(th0',percChanges_fixed(:,end));
% fprintf('Correlation with th0: %f\n',r)

figure(18);
plot(mw(find(mw)),percChanges_fixed(find(mw),end),'o')
xlabel('6MWT')
ylabel('Total % change')
r = corr(mw(find(mw))',percChanges_fixed(find(mw),end));
fprintf('Correlation with 6MWT: %f\n',r)

end



function GSLs_ROI = placeROI_fixed(sweep,subQ,extra)

GSLs = [];

for i = 1:length(subQ)
        
    ROI = sweep(subQ(i,1)+2+extra:subQ(i,1)+42+extra, subQ(i,2)-50:subQ(i,2)+50,i);
    GSLs_ROI(i) = median(reshape(ROI,1,size(ROI,1)*size(ROI,2)));

end


end

function GSLs_frac = placeROI_frac(sweep,subQ,interMus,frac)

GSLs = [];

ROI_depth = interMus(:,1) - subQ(:,1);

for i = 1:length(subQ)
    
    fractional = sweep(subQ(i,1)+2:subQ(i,1)+2+frac*ROI_depth(i), subQ(i,2)-50:subQ(i,2)+50,i);
    GSLs_frac(i) = median(reshape(fractional,1,size(fractional,1)*size(fractional,2)));

%     sweep([subQ(i,1)+2 subQ(i,1)+2+round(frac*ROI_depth(i))], subQ(i,2)-50:subQ(i,2)+50,i) = 0;
%     sweep(subQ(i,1)+subQ(i,1)+2+round(frac*ROI_depth(i)), [subQ(i,2)-50 subQ(i,2)+50],i) = 0;
    
end


% implay(sweep(:,:,1:length(subQ)))
% figure; imshow(sweep(:,:,length(subQ)));

end

function percs = atForces(GSL,sweepForces)

forces = sweepForces(1:length(GSL),2);
percs = zeros(1,18);

for F = 1:18
    
    thisForce = 1+F/2;
    
    idx = find( abs(forces-thisForce) == min(abs(forces-thisForce)) );
    idx = idx(1);
    percs(F) = GSL(idx);

end

end


function corrs(f, fixed, frac20, frac33, frac50)

for i = 1:18
    
    r20(i) = corr(fixed(:,i),frac20(:,i));    
    r33(i) = corr(fixed(:,i),frac33(:,i));
    r50(i) = corr(fixed(:,i),frac50(:,i));
    
end

SE = 1/sqrt(19-3);
conf = 1.96; % z for 95% confidence interval

z20 = .5*(log(1+r20)-log(1-r20));
z20(2,:) = z20(1,:) + SE*conf;
z20(3,:) = z20(1,:) - SE*conf;
r20 = (exp(2*z20)-1) ./ (exp(2*z20)+1);

z33 = .5*(log(1+r33)-log(1-r33));
z33(2,:) = z33(1,:) + SE*conf;
z33(3,:) = z33(1,:) - SE*conf;
r33 = (exp(2*z33)-1) ./ (exp(2*z33)+1);

z50 = .5*(log(1+r50)-log(1-r50));
z50(2,:) = z50(1,:) + SE*conf;
z50(3,:) = z50(1,:) - SE*conf;
r50 = (exp(2*z50)-1) ./ (exp(2*z50)+1);


figure(f); subplot(1,3,1)
errorbar([1.5:.5:10],r20(1,:),r20(1,:)-r20(3,:),r20(2,:)-r20(1,:),'ko-')
title('20')
xlabel('Force (N)')
ylabel('Correlation coefficient')
axis([1 10 .75 1])

subplot(1,3,2)
errorbar([1.5:.5:10],r33(1,:),r33(1,:)-r33(3,:),r33(2,:)-r33(1,:),'kx-')
title('33')
xlabel('Force (N)')
axis([1 10 .75 1])

subplot(1,3,3)
errorbar([1.5:.5:10],r50(1,:),r50(1,:)-r50(3,:),r50(2,:)-r50(1,:),'kd-')
title('50')
axis([1 10 .75 1])

xlabel('Force (N)')
ylabel('Correlation coefficient')
suptitle('Correlation coefficients')

end


function plotvForce(f, GSLs, sweepForces)

figure(f); hold on;
plot(sweepForces(1:length(GSLs),2),GSLs)

end

function plotvStrain(f, GSLs, thickness, th0)

strain = (th0-thickness(1:length(GSLs)))/th0;

figure(f); hold on;
plot(strain,GSLs,'LineWidth',2)

end

function plotvThickness(f,GSLs,thickness)

figure(f); hold on;
plot(thickness(1:length(GSLs)),GSLs,'LineWidth',2)

end

function plotROIvfrac(f,GSLs_ROI,GSLs_frac)

figure(f); hold on;
plot(GSLs_ROI,GSLs_frac,'LineWidth',2)
xlabel('GSL from fixed-size ROI')
ylabel('GSL from upper third of muscle')

end






