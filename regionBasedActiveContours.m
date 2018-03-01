% MACHINE VISION - Project


clear all;
close all;
clc;


%------------------------------------------------------------------------------------------

                                    % Image Segmentation                                                                                                                                                                                                                                                                                                                                                                                                           

%------------------------------------------------------------------------------------------
                                        % (e)

% Color image segmentation using active contours.
%------------------------------------------------------------------------------------------

tic

% Read image.
input = imread('MVI_4117_frame_0677.bmp');     

if (size(input, 3) ~= 3)
    error('Input image must be RGB.')
else
    % Display the original image.
    figure, imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image', 'FontSize', 17);
    
    % Perform morphological opening on the image with structuring element
    % strel('disk',10).
    background = imopen(input,strel('disk',10));
    figure; imshow(background);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('background','FontSize',17)
    
    % Image pre-processing to bring the leaves to the foreground while 
    % cutting the background.
    forgr=input-background;
    figure; imshow(forgr);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image without the background','FontSize',17)
    
    % Convert RGB image to grayscale.
    imgG = rgb2gray(forgr);
    figure;imshow(imgG);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image in grayscale','FontSize',17)
    
    % Adjust image intensity values.
    imgAdj = imadjust(imgG);
    figure;imshow(imgAdj);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('intensity values adjustment','FontSize',17)

    % Region Selection
    % Create mask so that we focus on a specific region (3 leaves region).
    mask = zeros(size(imgAdj));           % zeros = black
    mask(352:end-60,245:end-325) = 1;     % ones = white
    figure;imshow(mask);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('initial contour location','FontSize',17);
    
    % Region boundaries on image.
    figure; imshow(imgAdj)
    hold on
    visboundaries(mask) %Introduced in Matlab R2015a
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('region boundaries on image','FontSize',17)
    hold off

    % Segment image into foreground and background using active contours.
    bw = activecontour(imgAdj,mask,400);
    figure; imshow(bw);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image using active contour','FontSize',17)
    hold off
    
    % Remove objects with less than 50 connected pixels.
    bw = bwareaopen(bw, 50);
    figure; imshow(bw);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image after removing unwanted objects','FontSize',17)
    
    % Color the white leaves in the binary image.
    segimg=zeros(size(input)); 
    segimg(:,:,1) = bw.*double(input(:,:,1));
    segimg(:,:,2) = bw.*double(input(:,:,2));
    segimg(:,:,3) = bw.*double(input(:,:,3));

    %  Display the segmented image.
    figure, imagesc(uint8(segimg)), colormap gray
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image', 'FontSize', 17);
    
     % Original image.
    figure
    subplot(1,2,1), imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image','FontSize',17)
    % Segmented image.
    subplot(1,2,2), imshow(uint8(segimg));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image','FontSize',17)
    suptitle('Color image segmentation using active contours')
end

% Calculate elapsed time using tic toc.
activecontourtime=toc;
display(['Color image segmentation using active contour: Elapsed Time = ',num2str(activecontourtime),' seconds'])
