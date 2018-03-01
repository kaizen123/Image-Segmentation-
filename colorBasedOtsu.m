% MACHINE VISION - Project


clear all;
close all;
clc;


%------------------------------------------------------------------------------------------

                                    % Image Segmentation                                                                                                                                                                                                                                                                                                                                                                                                           

%------------------------------------------------------------------------------------------
                                        % (c)

% Color-Based Segmentation with HSV Color Model Using Otsu's Method of Thresholding.
%------------------------------------------------------------------------------------------

tic

% Read image.
input = imread('MVI_4117_frame_0677.bmp');
if (size(input, 3) ~= 3)
    error('Input image must be RGB.')
else
    % Display the original image.
    figure; imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image', 'FontSize', 17);

    % Perform morphological opening on the image with structuring element strel('disk',10)
    background = imopen(input,strel('disk',10)); % disk with radius 10
    figure; imshow(background);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('background','FontSize',17)

    % Cut the background so we can get better threshold.
    I=input-background;
    figure; imshow(I);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image without the background','FontSize',17)

    % Convert RGB colormap to HSV colormap.
    hsvI = rgb2hsv(I);
    H=hsvI(:,:,1);
    S=hsvI(:,:,2);
    V=hsvI(:,:,3);

    % Display all the components of HSV color model.
    figure;
    subplot(2,2,1), imshow(hsvI);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image in HSV colormap','FontSize',17)
    subplot(2,2,2), imshow(H);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('hue','FontSize',17)
    subplot(2,2,3), imshow(S); 
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('saturation','FontSize',17)
    subplot(2,2,4), imshow(V);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('value','FontSize',17)
    suptitle('Components of HSV color model.')

    % Find the histogram counts and the bin locations.
    counts = imhist(V);
    nbinsV = length(counts);
    p = counts / sum(counts);
    sigma_b = zeros;

    % Calculating best threshold using Otsu's method.
    for t = 1 : nbinsV
        q_L = sum(p(1 : t));
        q_H = sum(p(t + 1 : end));
        m_L = sum(p(1 : t) .* (1 : t)') / q_L;
        m_H = sum(p(t + 1 : end) .* (t + 1 : nbinsV)') / q_H;
        sigma_b(t) = q_L * q_H * (m_L - m_H)^2;
    end
    [~,threshold_otsuV] = max(sigma_b(:));

    % Use the traditional method of thresholding rather than im2bw.
    % For this method we use Otsu thresholding from previous step.
    bwV=V;
    [m,n]=size(V);
    threshold = threshold_otsuV/nbinsV;
    for i=1:m
        for j=1:n
            if( bwV(i,j)>threshold )
                % 255(1)->white 
                bwV(i,j)=255;
            else
                % 0->black
                bwV(i,j)=0;
            end;
        end;
    end;

    %Display image after Otsu thresholding.
    figure; imshow(bwV);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image after Otsu thresholding','FontSize',12)

    % Remove objects with less than 50 connected pixels.
    bwV = bwareaopen(bwV, 50);
    figure; imshow(bwV);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image after removing unwanted objects','FontSize',12)

    % Return the connected components found in binary image. Connectivity
    % is 4, meaning two-dimensional four-connected neighborhood of
    % components.
    cc = bwconncomp(bwV, 4);
    leaf = false(size(bwV));
    leaf(cc.PixelIdxList{30}) = true;
    leaf(cc.PixelIdxList{33}) = true;
    leaf(cc.PixelIdxList{36}) = true;
    figure; imshow(leaf);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('binary image - leaves', 'FontSize', 17);
    
    % Color the white leaves in the binary image.
    lastimage=zeros(size(input)); 
    lastimage(:,:,1) = leaf.*double(input(:,:,1));
    lastimage(:,:,2) = leaf.*double(input(:,:,2));
    lastimage(:,:,3) = leaf.*double(input(:,:,3));

    %  Display the segmented image.
    figure, imagesc(uint8(lastimage)), colormap gray
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image', 'FontSize', 17);
    
    % The result is colored leaves in the foreground with black background
    % Display the result side by side with the original image.
    figure
    subplot(1,2,1), imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image','FontSize',17)
    subplot(1,2,2), imshow(uint8(lastimage));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image','FontSize',17)
    suptitle('Image segmentation using Otsu thresholding')
end

% Calculate elapsed time using tic toc.
otsutime=toc;
display(['Color-Based Segmentation with HSV Color Model Using Otsu Method of Thresholding: Elapsed Time = ',num2str(otsutime),' seconds'])
