% MACHINE VISION - Project


clear all;
close all;
clc;


%------------------------------------------------------------------------------------------

                                    % Image Segmentation                                                                                                                                                                                                                                                                                                                                                                                                           

%------------------------------------------------------------------------------------------
                                           % (b)                                                                                                                                                                                                                                                                                                                                                                                                                                    

% Color-Based Segmentation Using the L*a*b* Color Space.
%------------------------------------------------------------------------------------------

tic

% Read image
input = imread('MVI_4117_frame_0677.bmp');     

if (size(input, 3) ~= 3)
    error('Input image must be RGB.')
else
    % Display the original image.
    figure, imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image', 'FontSize', 17);
    
    load regioncoordinates;
    % Let's assume we can see 6 major colors in the image: dark/medium/light green, dark/light yellow and white.  
    % We can visually distinguish these colors from one another using the L*a*b* colorspace (also known as CIELAB or CIE L*a*b*).
    % Our approach is to choose a small sample region for each color and to calculate each sample region's average color in 'a*b*' space. 
    % We will use these color markers to classify each pixel.
    nColors = 6;
    sample_regions = false([size(input,1) size(input,2) nColors]);
    for count = 1:nColors
        sample_regions(:,:,count) = roipoly(input,region_coordinates(:,1,count), region_coordinates(:,2,count));
    end
    
    % Convert the RGB image to L*a*b* color space using makecform & applycform.
    colorTransform = makecform('srgb2lab');
    lab_input = applycform(input, colorTransform);
    % Calculate the mean 'a*' and 'b*' value for each area that we extracted with roipoly. 
    % These values serve as our color markers in 'a*b*' space.
    a = lab_input(:,:,2);
    b = lab_input(:,:,3);
    color_markers = zeros([nColors, 2]);
    for count = 1:nColors
        color_markers(count,1) = mean2(a(sample_regions(:,:,count)));
        color_markers(count,2) = mean2(b(sample_regions(:,:,count)));
    end
   
    % Each color marker now has an 'a*' and a 'b*' value. 
    % We can classify each pixel in the lab_input image by calculating the Euclidean distance between that pixel and each color marker. 
    % The smallest distance will tell us that the pixel most closely matches that color marker.
    color_labels = 0:nColors-1;
    % Initialize matrices to be used in the nearest neighbor classification.
    a = double(a);
    b = double(b);
    distance = zeros([size(a), nColors]);
    % Perform classification.
    for count = 1:nColors
        distance(:,:,count) = ( (a - color_markers(count,1)).^2 + (b - color_markers(count,2)).^2 ).^0.5;
    end
    [~, label] = min(distance,[],3);
    label = color_labels(label);
    clear distance;
    
    % The label matrix contains a color label for each pixel in the input image. 
    % Use the label matrix to separate objects in the original image by color.
    rgb_label = repmat(label,[1 1 3]);
    segmented_images = zeros([size(input), nColors],'uint8');
    for count = 1:nColors
     leaves = input;
     leaves(rgb_label ~= color_labels(count)) = 0;
     segmented_images(:,:,:,count) = leaves;
    end
    % Region 5 is the one that includes the leaves.
    figure, imshow(segmented_images(:,:,:,5));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('yellow objects', 'FontSize', 17);
    
    % Convert image to binary image.
    binary_picture=im2bw(segmented_images(:,:,:,5));
    figure, imshow(binary_picture);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('binary image', 'FontSize', 17);
    % Removes all connected components (objects) that have fewer than 100 pixels from the binary image.
    binary_picture_edit = bwareaopen(binary_picture, 100);
    figure, imshow(binary_picture_edit);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('binary image - remove unwanted regions', 'FontSize', 17);
    
    % Find connected components in binary image.
    CC = bwconncomp(binary_picture_edit, 4);
    leaf = false(size(binary_picture_edit));
    leaf(CC.PixelIdxList{17}) = true;
    leaf(CC.PixelIdxList{11}) = true;
    leaf(CC.PixelIdxList{3}) = true;
    figure; imshow(leaf);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('binary image - leaves', 'FontSize', 17);
    
    [l, m, n]=size(input);
    % Recoloring the binary image using the original.
    leaves=input;
    for i=1:l 
        for j=1:m
            for k=1:n
                 if leaf(i,j)==0
                      leaves(i,j,k)=0;
                 end
            end
        end
    end

    % Display the segmented image.
    figure, imshow(leaves);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image', 'FontSize', 17);
    
    % Original image.
    figure
    subplot(1,2,1);
    imshow(input)
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image', 'FontSize', 17); 
    % Segmented image.
    subplot(1,2,2);
    imshow(leaves)
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image', 'FontSize', 17);
    suptitle('Color-Based Segmentation Using the L*a*b* Color Space');
end

% Calculate elapsed time using tic toc.
labtime=toc;
display(['Color-Based Segmentation Using the L*a*b* Color Space: Elapsed Time = ',num2str(labtime),' seconds'])
