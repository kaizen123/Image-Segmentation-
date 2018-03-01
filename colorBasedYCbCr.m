% MACHINE VISION - Project


clear all;
close all;
clc;


%------------------------------------------------------------------------------------------

                                    % Image Segmentation                                                                                                                                                                                                                                                                                                                                                                                                           

%------------------------------------------------------------------------------------------
                                           % (a)                                                                                                                                                                                                                                                                                                                                                                                                                                    

% Color Image Segmentation Considering the YCbCr Color Space Based on (Color) Thresholding.
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
    
    figure;
    subplot(2,2,1), imshow(input);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image in RGB colormap','FontSize',17)
    subplot(2,2,2), imshow(input(:,:,1));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('Red','FontSize',17)
    subplot(2,2,3), imshow(input(:,:,2));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('Green','FontSize',17)
    subplot(2,2,4), imshow(input(:,:,3));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('Blue','FontSize',17)
    suptitle('Components of RGB color model')
   
    % Convert the RGB values in map to the YCbCr color space.
    imgy = rgb2ycbcr(input);
    
    figure;
    subplot(2,2,1), imshow(imgy);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('image in YCbCr colormap','FontSize',17)
    subplot(2,2,2), imshow(imgy(:,:,1));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('Y (Luminance)','FontSize',17)
    subplot(2,2,3), imshow(imgy(:,:,2));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('Cb (Blue-difference)','FontSize',17)
    subplot(2,2,4), imshow(imgy(:,:,3));
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('Cr (Red-difference)','FontSize',17)
    suptitle('Components of YCbCr color model')
    
    % Convert image to double precision.
    dinput=double(input);   
   
    % Converting RGB space into YCbCr space using following formula.
    % Cb: Chrominance-Blue and Cr: Chrominance-Red.
    Cb = 0.001 * dinput(:,:,1) - 0.2 * dinput(:,:,2) + 0.4 * dinput(:,:,3) + 128;
    Cr =  0.45 * dinput(:,:,1) - 0.35 * dinput(:,:,2) - 0.09 * dinput(:,:,3) + 128;
    [rows,columns]=size(dinput(:,:,1));
    
    figure, hist(Cb);
    figure, hist(Cr);
    
    % Preallocate memory.
    leaves=zeros(rows,columns);

    for i=1:rows
        for j=1:columns 
            % Segmenting the yellow colored areas.     
            if  Cb(i,j)>=140 && Cr(i,j)>=140      % determined experimentally
                leaves(i,j)=1;            
            else       
                leaves(i,j)=0;    
            end    
        end
    end
    
    figure, imshow(leaves);
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('binary image', 'FontSize', 17);
    
    % Use the colored object mask to mask out the colored-only portions of the rgb image.
    output(:,:,1) = dinput(:,:,1).*leaves;   
    output(:,:,2) = dinput(:,:,2).*leaves; 
    output(:,:,3) = dinput(:,:,3).*leaves; 

    %  Display the segmented image.
    figure, imagesc(uint8(output)), colormap gray
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
    imagesc(uint8(output)), colormap gray
    axis equal; axis tight; axis off;  
    set(gcf, 'Color', 'White'); title('segmented image', 'FontSize', 17);
    suptitle('Image segmentation considering the YCbCr color space based on (color) thresholding');
end

% Calculate elapsed time using tic toc.
ycbcrtime=toc;
display(['Color Image Segmentation Considering the YCbCr Color Space Based on (Color) Thresholding: Elapsed Time = ',num2str(ycbcrtime),' seconds'])
