function [im, params] = wvfAperture(wvf, varargin)
% Create a synthetic aperture
%
% Synopsis
%   [im, params] = wvfAperture(wvf, varargin)
%
% Brief description:
%   Returns a pupil aperture function. A main application is creating
%   scattering flare. Returns an N x N monochromatic image emulating a
%   pupil aperture with small disks (dust) and lines (scratches).  The
%   disks and polylines have random size and opacity to an otherwise clear
%   aperture.
%
% Inputs
%   wvf: wavefront structure
%
% Optional key/val parameters
%   line mean - Number of lines
%   line sd
%   line opacity
%   line width
%   segment length
%
%   dot mean  - Number of dots
%   dot sd
%   dot opacity
%   dot radius
%
%   n sides - Number of aperture sides
%
% Output
%  im: A matrix of values in the [0, 1] range. 0 means completely
%     opaque and 1 means completely transparent. The returned matrix
%     is real-valued.  We ignore any phase shift that may be
%     introduced by the "dust" and "scratches".
%
% See also
%   piFlareApply, RandomDirtyAperture
%

% Examples:
%{
    wvf = wvfCreate;
    im = wvfAperture(wvf); % Default with some scratches and dust.
    ieNewGraphWin; imagesc(im); colormap(gray); axis image
%}
%{
    % Diffraction limited, circular (no scratches or dust)
    wvf = wvfCreate;
    im = wvfAperture(wvf,'dot mean',0,'dot sd',0,'line mean',0,'line sd',0); 
    ieNewGraphWin; imagesc(im); colormap(gray); axis image
%}
%{
    wvf = wvfCreate;
    im = wvfAperture(wvf,'segment length',100); % Default
    ieNewGraphWin; imagesc(im); colormap(gray); axis image
%}
%{
    wvf = wvfCreate;
    [im,params] = wvfAperture(wvf,'n sides',8); 
    ieNewGraphWin; imagesc(im); colormap(gray); axis image
%}

%% Inputs

varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('wvf',@isstruct);

p.addParameter('dotmean',20,@isnumeric);
p.addParameter('dotsd',5,@isnumeric);
p.addParameter('dotopacity',0.5,@isnumeric);
p.addParameter('dotradius',[],@isnumeric);

p.addParameter('linemean',20,@isnumeric);
p.addParameter('linesd',5,@isnumeric);
p.addParameter('lineopacity',0.5,@isnumeric);
p.addParameter('linewidth',2,@isnumeric);

p.addParameter('segmentlength',600,@isnumeric);

p.addParameter('texfile',[]);

p.addParameter('nsides',0, @(x)(isnumeric(x) && (x > 2 || x <= 0)));

p.parse(wvf,varargin{:});
dotMean     = p.Results.dotmean;
dotSD       = p.Results.dotsd;
dotOpacity  = p.Results.dotopacity;
dotRadius   = p.Results.dotradius;
lineMean       = p.Results.linemean;
lineSD         = p.Results.linesd;
lineOpacity    = p.Results.lineopacity;
lineWidth      = p.Results.linewidth;
segmentLength  = p.Results.segmentlength;
nSides         = p.Results.nsides;

texFile = p.Results.texfile;

% Adjust
imageSize = wvfGet(wvf, 'spatial samples');
im = ones([imageSize,imageSize], 'single');

if isempty(texFile)
    if isempty(dotRadius), dotRadius = round(imageSize/200); end
    %% Add dots (circles), simulating dust.

    % We should do more to control the random variable
    num_dots = randn(1,1)*dotSD + dotMean;
    num_dots = round(num_dots);

    % Make circles, but limit their size.
    % BW: Not sure why we don't just clip.
    max_radius = dotRadius*5;
    for i = 1:num_dots
        radius = max(round(dotRadius + rand*5),0);
        radius = min(radius,max_radius);
        circle_xyr = rand(1, 3, 'single') .* [imageSize, imageSize, radius];
        opacity = dotOpacity + (rand * 0.5);
        opacity = min(opacity,1);

        % Computer vision toolbox
        im = insertShape(im, 'FilledCircle', circle_xyr, 'Color', 'black', ...
            'Opacity', opacity);
    end

    % %% Add polylines, simulating scratches.
    % 
    num_lines = randn(1,1)*lineSD + lineMean;
    num_lines = round(num_lines);
    
    num_groups = randi([1, 7]);  % 生成1到5组
    if rand>0.5
        num_groups = 2;
    end
    lines_per_group = floor(num_lines / num_groups)-1;  % 每组的划痕数
    
    % Generate direction of each group (consistent within each group)
    directions = rand(1, num_groups) * 2 * pi;  % 随机生成 1 到 5 组的方向（角度）
    
    % For each group, generate scratches with the same direction
    for group_idx = 1:num_groups
        % Get number of lines for this group
        num_lines_in_group = lines_per_group;
        
        % Generate scratches in the same direction
        for i = 1:num_lines_in_group
            % Random length of the scratch segment
            % num_segments = 1;
            segment_length = rand * segmentLength;
    
            % Start position
            start_xy = rand(2, 1) * imageSize;
    
            % Calculate the direction vector for this group
            direction_vector = [cos(directions(group_idx)); sin(directions(group_idx))];
    
            % Generate points along the line (same direction)
            segments_xy = direction_vector * segment_length;
    
            % Vertices of the scratch line
            vertices_xy = cumsum([start_xy, segments_xy], 2);
            vertices_xy = reshape(vertices_xy, 1, []);
            
            % Width of the scratches
            width = round(max(1, lineWidth + randn * (lineWidth / 2)));
    
            % Set opacity of the line
            opacity = lineOpacity + (rand * 0.5);
            
            % Insert the scratch line into the image
            im = insertShape(im, 'Line', vertices_xy, 'LineWidth', width, 'Color', [opacity, opacity, opacity]);
        end
    end
    % 参数设置
    % num_lines = round(randn * lineSD + lineMean);  % 总线条数
    % num_groups = 2;  % 仅两组，互相垂直
    % lines_per_group = floor(num_lines / num_groups);  % 每组的线条数
    % 
    % % 创建空白图像
    % im = zeros(imageSize, imageSize, 3);
    % 
    % % 计算延申长度，确保线条超出边界 (2倍图像对角线长度)
    % extend_length = sqrt(2) * imageSize * 2;
    % 
    % % 随机生成第一个组的角度，然后计算互相垂直的角度
    % angle1 = rand * 180;  % 随机生成一个基准角度 (0~180度)
    % angle2 = mod(angle1 + 90, 180);  % 互相垂直的角度
    % 
    % % 设置两组方向向量
    % directions = [angle1, angle2];
    % dir_vectors = [cosd(directions); sind(directions)];  % 方向向量
    % 
    % % 等间距计算
    % spacing = imageSize / (lines_per_group + 1);
    % 
    % % 为每一组生成网格线
    % for group_idx = 1:num_groups
    %     direction_vector = dir_vectors(:, group_idx);
    %     perpendicular_vector = [-direction_vector(2); direction_vector(1)];  % 垂直方向向量
    % 
    %     % 生成当前组的线条
    %     for i = -lines_per_group:lines_per_group
    %         % 计算平行线的偏移位置
    %         offset = i * spacing;
    %         start_xy = imageSize / 2 + offset * perpendicular_vector;  % 起点沿垂直方向偏移
    % 
    %         % 计算无限延申的起点和终点
    %         extended_start = start_xy - extend_length * direction_vector;
    %         extended_end = start_xy + extend_length * direction_vector;
    % 
    %         % 转换为插入格式
    %         vertices_xy = [extended_start; extended_end];
    %         vertices_xy = reshape(vertices_xy, 1, []);
    % 
    %         % 设置宽度和透明度
    %         width = round(max(1, lineWidth + randn * (lineWidth / 2)));
    %         opacity = lineOpacity + (rand * 0.5);
    % 
    %         % 插入线条到图像
    %         im = insertShape(im, 'Line', vertices_xy, 'LineWidth', width, ...
    %                          'Color', [opacity, opacity, opacity]);
    %     end
    % end
else
    im = imread(texFile);
    try
        im = rgb2gray(im);
    catch
        % do nothing
    end
    im = double(im);
    im = im/max2(im);

    if mean2(im)<0.2
        im = 1-im;
    end
    im(im<0.95) = im(im<0.95)-rand(1);
    im(im<0)=0;
    im = imresize(im, [imageSize, imageSize]);
end
centerPoint = [imageSize/2 + 1, imageSize/2+1];
radius = (imageSize - 1)/3;

% Clip the image with a bounding polygon
if nSides > 0
    % create n-sided polygon
    pgon1 = nsidedpoly(nSides, 'Center', centerPoint, 'radius', radius);
    % create a binary image with the polygon
    pgonmask = poly2mask(floor(pgon1.Vertices(:,1)), floor(pgon1.Vertices(:,2)), imageSize, imageSize);
    im = im.*pgonmask;
    %clip with another polygon
    %add by hhq
    % centerPoint = [ 1, 1];
    % radius = (imageSize - 1)/2;
    % pgon1 = nsidedpoly(nSides, 'Center', centerPoint, 'radius', radius);
    % % create a binary image with the polygon
    % pgonmask = poly2mask(floor(pgon1.Vertices(:,1)), floor(pgon1.Vertices(:,2)), imageSize, imageSize);
    % im = im.*pgonmask;
end



% Color image to gray. In some cases, when there are no dots or scratches,
% im is just gray scale.
if ndims(im) == 3
    im = rgb2gray(im);
end

% Now make the pattern circular
[X,Y] = meshgrid((1:imageSize) - centerPoint(1),(1:imageSize) - centerPoint(2));
imRadius = sqrt(X.^2 + Y.^2);
% ieNewGraphWin; imagesc(imRadius); colormap(gray); colorbar; axis image
idx = (imRadius > radius);
im(idx) = 0;

%another circular
if(rand()>0.5)
    temp1 = round(imageSize/12 * 5);
    temp2 = round(imageSize/6);
    centerPoint = [temp1 + 1+randi(temp2), temp1 + 1+randi(temp2)];
    [X,Y] = meshgrid((1:imageSize) - centerPoint(1),(1:imageSize) - centerPoint(2));
    imRadius = sqrt(X.^2 + Y.^2);
    idx = (imRadius > radius);
    im(idx) = 0;
end
im = imrotate(im,randi(30));
im = imresize(im,[imageSize,imageSize]);
%%
if nargout == 2
    % Fill in params
    params.dotMean = dotMean;
    params.dotSD = dotSD;
    params.dotOpacity = dotOpacity;
    params.dotRadius = dotRadius;
    params.lineMean = lineMean;
    params.lineSD = lineSD;
    params.lineOpacity = lineOpacity;
    params.lineWidth = lineWidth;
    params.segmentLength = segmentLength;
    params.nsides = nSides;
end

end

%% Utility function
function xy = RandomPointsInUnitCircle(num_points)
% Random point generation within the unit circle

r = rand(1, num_points, 'single');   % Between 0 and 1

theta = rand(1, num_points, 'single') * 2 * pi; % Between 0 and 2pi

xy = [r .* cos(theta); r .* sin(theta)];  % Convert r,theta to xy

end
