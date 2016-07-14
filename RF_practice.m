
dindex = 42;
patient = '02016';
visitID = '130204/';
muscle = 'quadsB/';

load('/home/sisir/code/dmd/directories/GSL_crosssection.mat')
load('/home/sisir/code/dmd/RF_background/RFgain_adj.mat')

dynamic_range = 80;
toData = '/media/U/Output_Complete/';

f = fopen( strcat(toData,patient,'-',visitID,muscle,'force_data.lvm') );
C = textscan(f,'%f %f %f %f %f %f');
fclose(f);
forces = C{2};

numFrames = length(forces);
vid1 = zeros(1210,256,numFrames);
vid2 = zeros(1210,256,numFrames);
dropped = [];

for i = 1:numFrames

fid = fopen( strcat(toData,patient,'-',visitID,muscle,'US_Images/',patient,'_',num2str(i),'RF.dat') );
RFdata = fread(fid, inf, 'short');
fclose(fid);

RFdata = RFdata.';
if size(RFdata,2) ~= 464642
    dropped = [dropped i];
    continue
end

env = abs(hilbert(RFdata));
log_env = 20*log10(env/max(env));
log_env = 255/dynamic_range*(log_env+dynamic_range);
img3 = reshape(log_env(3:end), [], 256);
img3 = img3./RFgain_adj;

vid1(:,:,i) = img3([1:280 606:930 1211:1815],:)/255;
vid2(:,:,i) = img3([1:280+161 607+161:931+140 1211+140:1815],:)/255;

end


keepframes = setdiff(1:numFrames,dropped);
vid1 = vid1(:,:,keepframes);
vid2 = vid2(:,:,keepframes);
forces = forces(keepframes);
fmax = find(forces == max(forces));
fidx = find(round(forces(1:fmax))==2,1,'last');
fdiff = diff(forces(1:fidx));
fstart = find(fdiff < 0,1,'last');
forces = forces(fstart:fmax);
vid1 = vid1(:,:,fstart:fmax);
vid2 = vid2(:,:,fstart:fmax);


storedforces = d(dindex).forces;
idx = [];

for i = 1:18
    foundme = find( forces == storedforces(i) );
    
    if isempty(foundme)
        foundme = find( abs(forces-storedforces(i)) ==  min(abs(forces-storedforces(i))) );
    end
    
    idx = [idx foundme(1)];
end

vid1 = vid1(:,:,idx);
vid2 = vid2(:,:,idx);
forces = forces(idx);

ellipse_params = d(dindex).ellipse;
sweep = d(dindex).sweep;

playMeEllipse;

ymax = size(vid1,1);
xmax = size(vid1,2);

z = ellipse_params(:,[1 2]);
a = ellipse_params(:,3);
b = ellipse_params(:,4);
alpha = ellipse_params(:,5);

npts = 5000;
t = linspace(0, 2*pi, npts);

img1 = vid1;
img2 = vid2;

for k = 1:18
    
        % Rotation matrix
    Q = [cos(alpha(k)), -sin(alpha(k)); sin(alpha(k)) cos(alpha(k))];
    % Ellipse points
    X = Q * [a(k) * cos(t); b(k) * sin(t)] + repmat(z(k,:)', 1, npts);
    
    Y = round(X(1,:)*1210/220);
    X = round(X(2,:)*256/315);
    
    
    Y = max(min(Y,ymax),1);
    X = max(min(X,xmax),1);
    
    
    for m = 1:size(X,2)
        img1(Y(m),X(m),k) = 0;
        img2(Y(m),X(m),k) = 0;
    end
   
    bottom = max(Y);
    top = min(Y);
    useDepth = round(top + (bottom-top)/3);
    
    ROIcontents1 = [];
    ROIcontents2 = [];
    
    for yy = top:useDepth
        rowIndex = find(Y == yy);
        bookendLeft = min(X(rowIndex));
        bookendRight = max(X(rowIndex));

%         img(yy,bookendLeft:bookendRight,k) = 1;

        ROIcontents1 = [ROIcontents1 vid1(yy,bookendLeft:bookendRight,k)];
        ROIcontents2 = [ROIcontents2 vid2(yy,bookendLeft:bookendRight,k)];

    end

    gsl1(k) = median(ROIcontents1);
    gsl2(k) = median(ROIcontents2);
    
end

implay(img1);
implay(img1);

clear data log_env img3 ans c fs dynamic_range f fid i env foundme 
clear di keepframes img3 one storedforces patient muscle visitID k Q
clear numFrames toData RFdata C idx xmax ymax a alpha b bookendLeft bookendRight
clear ROIcontents clear npts t top useDepth X Y yy z rowIndex sweep bottom


% 
% 
% r(dindex).sweep = vid;
% r(dindex).forces = forces;
% r(dindex).gsl = gsl';
% r(dindex).done = 1;








