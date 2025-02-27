function scene = sceneHDRLights(varargin)
% Create a HDR chart of circles, squares and lines (lights) on a black
% background
%
%   scene = sceneHDRLights()
%
% Inputs
%
% Returns
%   scene:         HDR chart as a scene
%
% See also
%   sceneHDRChart, sceneReflectanceChart, macbethChartCreate

%Examples:
%{
  scene = sceneHDRLights();
  sceneWindow(scene);
%}
%{
  scene = sceneHDRLights('n circles',4,'radius',repmat(0.01,1,4),'circle colors',{'white'});
  sceneWindow(scene);
%}
%{
 params = hdrlightsP;
 scene = sceneCreate('hdr lights',params);
%}

%%
varargin = ieParamFormat(varargin);
p = inputParser;

p.addParameter('imagesize',384,@isnumeric);

p.addParameter('ncircles',4,@isnumeric);
p.addParameter('radius',[0.01,0.035,0.07,0.1],@isvector);
p.addParameter('circlecolors',{'white','green','blue','yellow','magenta','white'},@iscell);

p.addParameter('nlines',4,@isnumeric);
p.addParameter('linelength',0.02,@isnumeric);
p.addParameter('linecolors',{'white','green','blue','yellow','magenta','white'},@iscell);

p.parse(varargin{:});

%% Spatial pattern
imSize   = p.Results.imagesize;
img = zeros(imSize,imSize);

%% Put circles with different sizes across the top third
nCircles = p.Results.ncircles;
radius   = p.Results.radius;
cColors  = p.Results.circlecolors;
% y = round(imSize * 0.25);
% xvals = round(linspace(0.2,0.8,nCircles)*imSize);
y = round(imSize * 0.5);
xvals = round(linspace(0.5,0.5,nCircles)*imSize);
radius = radius*imSize;
for ii = 1:numel(xvals)
    %cc = mod(ii,numel(cColors)) + 1;    
    %img = insertShape(img,'filled-circle',[xvals(ii),y,radius(ii)],'Color',cColors{cc});
    %red = 0.05+0.95*rand();
    %green = 0.05+0.95*rand();
    %blue = 0.05+0.95*rand();
    %more case that will really happen in real dataset
    light_type = randi([1, 5]);  % choose those appear more in real world

    switch light_type
        case 1  % 白光
            red = 0.7 + 0.3 * rand();
            green = 0.7 + 0.3 * rand();
            blue = 0.7 + 0.3 * rand();
        case 2  % 蓝光
            red = 30/255 + 0.1 * rand();
            green = 144/255 + 0.1 * rand();
            blue = 0.95 + 0.05 * rand();
        case 3  % 蓝光
            red = 135/255 + 0.1 * rand();
            green = 206/255 + 0.1 * rand();
            blue = 0.95 + 0.05 * rand();
        % case 3  % 红光
        %     red = 0.5 + 0.5 * rand();
        %     green = 0.05 + 0.45 * rand();
        %     blue = 0.05 + 0.45 * rand();
        % case 4  % 紫光
        %     red = 0.5 + 0.5 * rand();
        %     green = 0.05 + 0.1 * rand();
        %     blue = 0.5 + 0.5 * rand();
        case 4  % 太阳光
            red = 0.8 + 0.2 * rand();
            green = 0.7 + 0.3 * rand();
            blue = 0.5 + 0.5 * rand();
        case 5  % 太阳光
            red = 0.8 + 0.2 * rand();
            green = 0.7 + 0.3 * rand();
            blue = 0.5 + 0.5 * rand();
    end
    % 生成椭圆，随机调整长短轴
    %aspect_ratio = 1 + (rand - 0.5) * 0.5;  % 随机生成一个比例，控制长短轴比例，1表示圆形，<1表示横向拉长，>1表示纵向拉长
    
    % 根据随机生成的比例调整长短轴
    diff = 0.1*radius(ii)+0.9*radius(ii)*rand();
    long_axis = radius(ii)+diff;  % 长轴
    short_axis = radius(ii)-diff;  % 短轴
    
    % 随机生成一个旋转角度，范围是 0 到 360 度
    angle = rand() * 360;
    
    % 插入椭圆，确保传入五列矩阵，添加旋转角度
    img = insertShape(img, 'filled-ellipse', [xvals(ii), y, long_axis, short_axis, angle], 'Color', [red, green, blue]);


end

%% Put lines of different thickness and orientation around the middle
nLines       = p.Results.nlines;
lineLength   = p.Results.linelength;
lColors      = p.Results.linecolors;

y = round(imSize * 0.5);
xvals = round(linspace(0.1,0.8,nLines)*imSize);
lineLength = round(lineLength*imSize);
hw = round([1,7*lineLength; 1,3*lineLength; 3*lineLength,1; 8*lineLength,1]);
for ii = 1:numel(xvals)
    cc = mod(ii,numel(lColors)) + 1;
    img = insertShape(img,'filled-rectangle',[xvals(ii),y,hw(ii,1),hw(ii,2)],'Color',lColors{cc});
end

%% Squares
%%I have deleted it, can't control
% y = round(imSize * 0.75);
% xvals = round(linspace(0.1,0.7,3)*imSize);
% squareEdge = imSize/64;
% hw = [2 2; 5 5; 9 9]*squareEdge;
% for ii = 1:numel(xvals)
%     img = insertShape(img,'filled-rectangle',[xvals(ii),y,hw(ii,1),hw(ii,2)],'Color','white');
% end
% ieNewGraphWin; imagesc(img); axis image

% Add a uniform low level to set the dynamic range
wave = 400:10:700;
scene = sceneFromFile(img,'rgb',1e5,displayCreate,wave);
sceneU = sceneCreate('uniform',imSize,wave);
sceneU = sceneSet(sceneU,'mean luminance',1e-2);
scene = sceneAdd(scene,sceneU);
end




