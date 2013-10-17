function image = preprocess(image)

    % If image is colored, take average of color channels
    image = (image(:,:,1)+image(:,:,2)+image(:,:,3))/3;
    
    % Apply Gaussian filter
    blur = fspecial('gaussian',7,7);
    iterations = 3;
    for ii = 1:iterations
        imfilter(image, blur);
    end
    
    % Optional Weighted Sobel Filtering
%     hSobel = fspecial('sobel');
%     vSobel = hSobel';
%     hSobelImage = imfilter(image, hSobel);
%     vSobelImage = imfilter(image, vSobel);
    
end