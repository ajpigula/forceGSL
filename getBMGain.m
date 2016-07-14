basePath = '/media/U/Output_Complete';

dynamic_range = 80;

d = dir(basePath);
imgs = [];


for i = 3:length(d)-42
    d2 = dir(fullfile(basePath,d(i).name));
    
    for j = 4:length(d2)
        d3 = dir(fullfile(basePath,d(i).name,d2(j).name,'US_Images','*_1BM.dat'));
        
        fprintf('\n%s',fullfile(basePath,d(i).name,d2(j).name,'US_Images',d3.name))
        
        fid = fopen(fullfile(basePath,d(i).name,d2(j).name,'US_Images',d3.name));
        
        BMData = fread(fid);
        fclose(fid);
        
        if isempty(BMData)
            continue;
        end

        BMData = reshape(BMData(1:(end-4)), 640,480);
        BMData = imrotate(BMData(95:409,20:430), 90);
        BMData = mat2gray(BMData, [0 255]);
        imgs(:,:,end+1) = BMData;
        
    end
    
    
end

BMgain = imgs(:,:,2:2043);
% The first image is all 0s, and the images after 2043 are from data
% intended to be discarded (indicated by an X at the beginning of the
% filename)
% 
% BMgain_adj = BMgain;
% for i = 1:110
%     BMgain_adj(i,:) = BMgain(111,:);
% end
% Because most images have the subcutaneous layer in the same place so it
% doesn't average out

% BMgain_adj = BMgain_adj/(mean2(BMgain_adj));
% Normalization

% fprintf('\n\n')
% figure; imshow(mean(imgs,3)/255);

