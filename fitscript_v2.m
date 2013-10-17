
clear all
close all

vid = mmreader('extract4a.avi');
F = get(vid, 'NumberOfFrames');
Q = []; % rows = sequences of least squares coefficients in x
V = []; % rows = sequences of least squares coefficients in y
blur = fspecial('gaussian',7,7);

image = read(vid,1);
imshow(image)
disp('Click the origin.');
[x0 y0]=ginput(1);
disp('Click the left boundary.');
[leftx lefty]=ginput(1);
disp('Click the right boundary');
[rightx righty]=ginput(1);
disp('Click the upper boundary.');
[upx upy]=ginput(1);
disp('CLick the lower boundary.');
[lowx lowy]=ginput(1);

for ff = 1:F    
    image = read(vid,ff);
    I = rgb2gray(image);
    for ii=1:9
        I = imfilter(I, blur);
    end
    bw = edge(I, 'canny');
    bw(:,1:leftx) = 0; % remove pixels to the left of eye
    bw(:,rightx:end) = 0; % remove pixels to the right of eye
    bw(lowy:end,:) = 0; % remove pixels from below the eye
    bw(1:upy,:) = 0; % remove pixels from above the eye
    
    [x y] = find(bw==1);
    points = [x y];
    hull = convhull(x,y);
    points = [points(hull,2), points(hull,1)];
    bwHull = zeros(size(bw));
    for ii=1:length(points)
        bwHull(points(ii,2),points(ii,1))=1;
    end
    bw = bwHull;
% imshow(bw)
% hull(ff)=getframe;

    bwLine=zeros(size(bw));
    for ii=1:length(points)-1
        x0 = points(ii,2);
        y0 = points(ii,1);
        x1 = points(ii+1,2);
        y1 = points(ii+1,1);
        bwLine = drawLine(bwLine, x0, y0, x1, y1, 1);
    end
    bw = bwLine;
    % Re-order points with respect to angle
    [x y] = find(bw==1);
    points = [x y];
%     x0 = mean(x);
%     y0 = mean(y);
    dx = zeros(length(x),1); dy = zeros(length(y),1);
    for ii=1:length(points)
        dx(ii) = x(ii) - x0; dy(ii) = y(ii) - y0;
    end
    [theta rho] = cart2pol(dx,dy);
    sortedPol = sortrows([theta rho]);
    [x y] = pol2cart(sortedPol(:,1),sortedPol(:,2));

% Trigonometric least squares fitting in space
    theta = (0:(length(x)-1)).*(2*pi / length(x));
    C=cos(theta(:)*(0:10));
    S=sin(theta(:)*(1:10));
    A=[C S];
    v = A\x; % y coefficients: a,b
    q = A\y; % x coefficients: c,d
    Q(:,end+1) = q; % rows are sequences of coefficients
    V(:,end+1) = v; % rows are sequences of coefficients
        plot(A*q,-A*v)
        axis('equal')
%         fit(ff) = getframe;   
end

% Trigonometric least sqaures fitting in time
% numSequences = 2*numModes + 1;
numSequencs = 21;
phi = (0:(F-1)).*(2*pi / F);
newXCoeffs = [];
newYCoeffs = [];
cModes = cos(phi(:)*(0:5));
sModes = sin(phi(:)*(1:5));
modeMatrix = [cModes sModes];
xCoeffMat=[];
yCoeffMat=[];
for tt = 1:length(phi)
    for ii=1:21
        smoothingCoeffs = modeMatrix \ Q(ii,:)';
        coeff = smoothingCoeffs(1);
        for jj=2:6
            coeff = coeff+smoothingCoeffs(jj)*cos((jj-1)*phi(tt));
        end
        for jj=7:11
            coeff = coeff+smoothingCoeffs(jj)*sin((jj-6)*phi(tt));
        end
        newXCoeffs(ii)=coeff;
    end
    xCoeffMat(:,tt) = newXCoeffs;
end
for tt = 1:length(phi)
    for ii=1:21
        smoothingCoeffs = modeMatrix \ V(ii,:)';
        coeff = smoothingCoeffs(1);
        for jj=2:6
            coeff = coeff+smoothingCoeffs(jj)*cos((jj-1)*phi(tt));
        end
        for jj=7:11
            coeff = coeff+smoothingCoeffs(jj)*sin((jj-6)*phi(tt));
        end
        newYCoeffs(ii)=coeff;
    end    
    yCoeffMat(:,tt) = newYCoeffs;
end

% smoothed video generation
for ff=1:F %F
    qq = xCoeffMat(:,ff);
    vv = yCoeffMat(:,ff);
    xx = A*qq;
    yy = A*vv;
%     plot(xx, -yy)
%     axis('equal')  
% Optional: overlay on original video
% Convert back to row/col coords    
    xx = floor(xx+y0);
    yy = floor(yy+x0);
    frame = read(vid,ff);
    frame=frame(:,:,1);
    imshow(frame);
    hold on
    plot(xx,yy)
    axis('equal');
    hold off
    mov(ff)=getframe(gcf);
end