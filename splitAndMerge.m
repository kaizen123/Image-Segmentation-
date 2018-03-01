% MACHINE VISION - Project


close all;
clear all;
clc;


%------------------------------------------------------------------------------------------

                                    % Image Segmentation                                                                                                                                                                                                                                                                                                                                                                                                           

%------------------------------------------------------------------------------------------
                                        % (f)

% Image Segmentation using Split and Merge Method.
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
    forgr=input-background;
    figure; imshow(forgr);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image without the background','FontSize',17)

    % Convert RGB image to grayscale.
    imgG = rgb2gray(forgr);
    figure; imshow(imgG);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image in grayscale','FontSize',17)
    
    % Adjust image intensity values.
    imAdj = imadjust(imgG);
    figure;imshow(imAdj);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('intensity values adjustment','FontSize',17)

    % Image padding in a size that corresponds to the next power of two
    % to prepare the image for qtdecomp
    Q = 2^nextpow2(max(size(imAdj)));
    [M, N] = size(imAdj);
    pad = padarray(imAdj, [Q - M, Q - N], 'post');
    figure;imshow(pad)
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image after padding','FontSize',17)

    % Divide a square image into four equal-sized square blocks, and
    % then test each block to see if meets some criterion of homogeneity.
    S = qtdecomp(pad,.8);

    % Now merge by looking at each quadregion and setting all its 
    % elements to 1 if the block satisfies the predicate.
    
    % Get the size of the largest block. Use full because S is sparse.
    Lmax = full(max(S(:)));
    % Set the output image initially to all ones. The MARKER array is
    % used later to establish connectivity.
    imconn = ones(size(pad));
    MARKER = ones(size(pad));
    % Begin the merging stage.
    for K = 1:Lmax
        % Returns in VALS an array containing the DIM-by-DIM
        % blocks in the quadtree decomposition of I.
        [vals, r, c] = qtgetblk(pad, S, K);
        if ~isempty(vals)
            % Check the predicate for each of the regions
            % of size K-by-K with coordinates given by vectors
            % r and c.
            for I = 1:length(r)
                xlow = r(I); ylow = c(I);
                xhigh = xlow + K - 1; yhigh = ylow + K - 1;
                region = pad(xlow:xhigh, ylow:yhigh);
                % Finding predicate (flag).
                sd = std2(region);
                m = mean2(region);
                flag = (sd > 10) & (m > 0) & (m < 125);

                if flag 
                    imconn(xlow:xhigh, ylow:yhigh) = 0;
                    MARKER(xlow, ylow) = 0;
                end
            end
        end
    end
    % Finally, obtain each connected region and label it with a
    % different integer value using function bwlabel.
    imconn = bwlabel(imreconstruct(MARKER, imconn));
    figure; imshow(imconn)
    axis equal; axis tight; axis off;  
    title('image with connected regions (padded)','FontSize',17)
    
    % Crop the padding.
    imcrop = imconn(1:M, 1:N);
    figure; imshow(imcrop)
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image with connected regions (cropped)','FontSize',17)

    % Remove objects with less than 50 connected pixels.
    bw = bwareaopen(imcrop, 50);
    figure; imshow(bw);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image after removing unwanted objects','FontSize',12)
    
    % Return the connected components found in binary image. Connectivity
    % is 4, meaning two-dimensional four-connected neighborhood of
    % components.
    cc = bwconncomp(bw, 4);
    leaf = false(size(bw));
    leaf(cc.PixelIdxList{40}) = true;
    leaf(cc.PixelIdxList{42}) = true;
    leaf(cc.PixelIdxList{44}) = true;
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
    
    % Original image.
    figure
    subplot(1,2,1), imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('original image','FontSize',17)
    % Segmented image.
    subplot(1,2,2), imshow(uint8(lastimage));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image','FontSize',17)
    suptitle('Image Segmentation using Split and Merge Method')
end

% Calculate elapsed time using tic toc.
splitmergetime=toc;
display(['Image Segmentation using Split and Merge Method: Elapsed Time = ',num2str(splitmergetime),' seconds'])
