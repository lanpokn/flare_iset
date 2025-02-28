loop = 5000;
% Define a directory to save images
outputDir = './output_images/flare2';  
outputDir_aper = './output_images/aperture2';  
outputDir_light = './output_images/light2';
apertureSize = 4800;
imageSize = 1200;
%from 532 is new
for i = 726:loop
    %scene = sceneCreate('hdr image',...
     %   'dynamic range',5,...
      %  'patch shape','circle','npatches',5,'patch size',10);
    
    scene = sceneCreate('hdr','ncircles',1,'nlines',0,'radius',0.01,'circlecolors',{'random'},'imagesize',imageSize);
    %scene = sceneCreate('hdr','ncircles',1,'nlines',1);
    %scene = sceneSet(scene,'fov',1);
    scene = sceneSet(scene,'fov',1);
    % create oi
    [oi,wvf] = oiCreate('wvf');
    
    % This sets the sample density of the aperture image. We use a high
    % sampling density to show the flare clearly.1024*4
    %1000 4000 is good
    %800 6000 is good
    wvf = wvfSet(wvf,'npixels',apertureSize);
    
    
    
    % % Try change to dot and line scratches, as you like
    nsides = [0 randi([3,10]) randi([30,60])];
    img = cell(numel(nsides),1);
    sensor = [];
    sensor_l = [];
    

    
    % Create the directory if it doesn't exist
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    % Iterate over nsides to generate and save images
    % 到时候如果做单张仿真图片的对齐，只需要把下边涉及随机的一些关键函数，一律开一个非随机版本即可
    % ghost和光圈这些还有关，但就是没找到调整光源大小的方法，无敌了
    for ii = 1:numel(nsides)
        %scratch
        % 'dot mean',200, 'dot sd',100, 'dot opacity',0.6,'dot radius',15+50*rand(),...
        % 'line mean',200, 'line sd',100, 'line opacity',0.6,'linewidth',5+25*rand(),'segmentlength',3000+4500*rand());
        [aperture, params] = wvfAperture(wvf,'nsides',nsides(ii),...
            'dot mean',100, 'dot sd',100, 'dot opacity',0.6,'dot radius',5+25*rand(),...
            'line mean',300, 'line sd',150, 'line opacity',0.6,'linewidth',2+10*rand(),'segmentlength',3000+4500*rand());
        apertureFilename = fullfile(outputDir_aper, sprintf('image_%d_%d.png', i,nsides(ii)));
        imwrite(aperture, apertureFilename);  % Save as PNG (you can change the format)
    
        %无敌了，根本找不到合理的改zcoeff的接口，只能在wvfComputer里自己搞了
        oi = oiSet(oi,'optics wvf',wvf);
    
        %oi = oiSet(oi,'fnumber',1.5);
        %这个改完，晕的纹路真实了一些，但是依然没有调整光源大小的作用
        %能改晕的样子，也可以随机
        oi = oiSet(oi,'fnumber',1.5);
        
        %这个影响光源大小？好像不是，倒是让结果图片变成800x800的标准版了
        %oi = oiSet(oi,'focal length',4.38e-3,'m');
        oi = oiSet(oi,'focal length',2.38e-4,'m');
        %defocusD = 0; vertical_astigmatism = 0; oblique_astigmatism = 0;
        %oi = oiSet(oi,'zcoeffs',[defocusD,vertical_astigmatism,oblique_astigmatism]
        %     oi = oiCompute(wvf,scene);
        %这个和大小关系也不大，但是和颜色相关的稳定性关系大
        %oi = oiCompute(oi, scene,'crop',true,'pixel size',3e-6,'aperture',aperture);%3e-6
        oi = oiCompute(oi, scene,'crop',true,'pixel size',3e-6,'aperture',aperture);%3e-6
        % To see the flare in the OI, we must use hdr rendering mode.
        %this matters the size!
        %这个不给。或者给的很小，比如0.1，就是光源，但0.1可能过于小了.好像并不行，四边形时不行
        Illumi_level = 5+45*rand();
        oi = oiAdjustIlluminance(oi, Illumi_level);
        %oiWindow(oi,'render flag','hdr');
        % To visualize, we used the ip window image
        if isempty(sensor)
            % First time through, create a sensor using this function.
            [ip, sensor] = piRadiance2RGB(oi, 'etime', 1/10);
        else
            % Next time, use the sensor
            sensor = sensorCompute(sensor, oi);
            ip = ipCompute(ip, sensor);
        end
        %[ip, sensor] = piRadiance2RGB(oi, 'etime', 1/10);
    
        % Store the image into the img cell array
        img{ii} = ipGet(ip, 'srgb');
        
        % Save the image to the output directory
        imgFilename = fullfile(outputDir, sprintf('image_%d_%d.png', i,nsides(ii)));
        img{ii} = imresize(img{ii},[512,512]);
        imwrite(img{ii}, imgFilename); 

        % get a ideal pupil
        im = ones([apertureSize,apertureSize], 'single');
        centerPoint = [apertureSize/2 + 1, apertureSize/2+1];
        radius = (apertureSize - 1)/3;
        [X,Y] = meshgrid((1:apertureSize) - centerPoint(1),(1:apertureSize) - centerPoint(2));
        imRadius = sqrt(X.^2 + Y.^2);
        % ieNewGraphWin; imagesc(imRadius); colormap(gray); colorbar; axis image
        idx = (imRadius > radius);
        im(idx) = 0;
        aperture_ideal = im;

        oi = oiCompute(oi, scene,'crop',true,'pixel size',3e-6,'aperture',aperture_ideal);%3e-6
        oi = oiAdjustIlluminance(oi, Illumi_level/5.0);

        if isempty(sensor_l)
            % First time through, create a sensor using this function.
            [ip_l, sensor_l] = piRadiance2RGB(oi, 'etime', 1/10);
        else
            % Next time, use the sensor
            sensor_l = sensorCompute(sensor_l, oi);
            ip_l = ipCompute(ip_l, sensor_l);
        end
        %[ip, sensor] = piRadiance2RGB(oi, 'etime', 1/10);
    
        % Store the image into the img cell array
        img{ii} = ipGet(ip_l, 'srgb');
        img{ii} = imresize(img{ii},[512,512]);
        % Save the image to the output directory
        imgFilename = fullfile(outputDir_light, sprintf('image_%d_%d.png', i,nsides(ii)));

        % % 读取或获取三通道图像 img{ii} (已调整为 512x512)
        % I = 0.2989 * img{ii}(:,:,1) + 0.5870 * img{ii}(:,:,2) + 0.1140 * img{ii}(:,:,3); % 计算强度
        % 
        % max_intensity = max(I(:));  % 计算最大强度
        % threshold = 0.8 * max_intensity; % 计算阈值
        % 
        % mask = I >= threshold; % 生成掩码（保留强度大于阈值的部分）
        % 
        % % 扩展掩码到三通道
        % mask3 = cat(3, mask, mask, mask);
        % 
        % % 进行归零处理
        % img{ii} = img{ii} .* mask3;
        % 
        imwrite(img{ii}, imgFilename); 

        
        % Optionally, you can display the image
        %ipWindow(ip);
    end
end
