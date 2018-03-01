% MACHINE VISION - Project


clear all;
close all;
clc;


%------------------------------------------------------------------------------------------

                                    % Image Segmentation                                                                                                                                                                                                                                                                                                                                                                                                           

%------------------------------------------------------------------------------------------
                                           % (d)                                                                                                                                                                                                                

% Color-Based Segmentation Using k-Means Clustering.
%------------------------------------------------------------------------------------------

tic

% Read image.
input = imread('MVI_4117_frame_0715.bmp');     

if (size(input, 3) ~= 3)
    error('Input image must be RGB.')
else
    % Display the original image.
    figure, imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image', 'FontSize', 17);
    
    % Convert Image from RGB Color Space to L*a*b* Color Space.
    % Convert the image to L*a*b* color space using makecform and applycform.
    cform = makecform('srgb2lab');
    lab_he = applycform(input,cform);
    
    % Since the color information exists in the 'a*b*' space, objects are pixels with 'a*' and 'b*' values.
    % Use kmeans to cluster the objects into six clusters using the Euclidean distance metric.
    ab = double(lab_he(:,:,2:3));
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab,nrows*ncols,2);

    nColors = 6; % determined experimentally
    % Repeat the clustering 3 times to avoid local minima.
    [cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', 'Replicates',3);
    
    % For every object in input, kmeans returns an index corresponding to a cluster. 
    % Label every pixel in the image with its cluster_index.
    pixel_labels = reshape(cluster_idx,nrows,ncols);
    figure, imshow(pixel_labels,[]);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image labeled by cluster index', 'FontSize', 17);

    segmented_images = cell(1,6);
    rgb_label = repmat(pixel_labels,[1 1 3]);

    for k = 1:nColors
        color = input;
        color(rgb_label ~= k) = 0;
        segmented_images{k} = color;
    end

%     figure, imshow(segmented_images{1}), title('objects in cluster 1');
%     figure, imshow(segmented_images{2}), title('objects in cluster 2');
%     figure, imshow(segmented_images{3}), title('objects in cluster 3');
%     figure, imshow(segmented_images{4}), title('objects in cluster 4');
%     figure, imshow(segmented_images{5}), title('objects in cluster 5');
%     figure, imshow(segmented_images{6}), title('objects in cluster 6');
    
    % Recall that the 'L*' layer contains the brightness values of each color. 
    % Find the cluster that contains the yellow objects. 
    % Extract the brightness values of the pixels in this cluster and threshold them.
    % We programmatically determine the index of the cluster containing the blue objects because kmeans will not return the same cluster_idx value every time. 
    % We do this using the cluster_center value, which contains the mean 'a*' and 'b*' value for each cluster. 
    mean_cluster_value = mean(cluster_center,2); 
    [tmp, idx] = sort(mean_cluster_value);
    yellow_cluster_num = idx(4); % Determined experimentally.

    L = lab_he(:,:,1);
    yellow_idx = find(pixel_labels == yellow_cluster_num);
    L_yellow = L(yellow_idx);
    isleaf= im2bw(L_yellow,1); % Determine the level experimentally. In this case, graythresh not helpful.
    % Use the mask is_leaf to label which pixels belong to the leaves. Then display the leaves in a separate image.
    
    yellowobjects = repmat(uint8(0),[nrows ncols]);
    yellowobjects(yellow_idx(isleaf==false)) = 1;
    yellowobjects = repmat(yellowobjects,[1 1 3]);
    yellow = input;
    yellow(yellowobjects ~= 1) = 0;
    figure, imshow(yellow);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('yellow objects', 'FontSize', 17);
    
    % Convert image to binary image.
    binary_yellow=im2bw(yellow);
    figure, imshow(binary_yellow);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('binary image', 'FontSize', 17);
    % Removes all connected components (objects) that have fewer than 100 pixels from the binary image.
    binary_yellow_edit = bwareaopen(binary_yellow, 100);
    figure, imshow(binary_yellow_edit);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('binary image - remove unwanted regions', 'FontSize', 17);
    
    % Find connected components in binary image.
    CC = bwconncomp(binary_yellow_edit, 4);
    leaf = false(size(binary_yellow_edit));
    leaf(CC.PixelIdxList{17}) = true;
    leaf(CC.PixelIdxList{13}) = true;
    leaf(CC.PixelIdxList{4}) = true;
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
    suptitle('Color-Based Segmentation Using k-Means Clustering');
end

% Calculate elapsed time using tic toc.
kmeanstime=toc;
display(['Color-Based Segmentation Using K-Means Clustering: Elapsed Time = ',num2str(kmeanstime),' seconds'])
