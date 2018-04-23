%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The implementation of the BM3D paper (Image Denoising by Sparse 3-D 
% Transform-Domain Collaborative Filtering)
% For more information, you can see the following two papers:
% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4271520
% http://www.ipol.im/pub/art/2012/l-bm3d/article.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
sigma=0.1225; window_size= 8; search_width= 19; 
l2= 0; selection_number = 8; l3= 2.7;

org_img = (imread('lena.jpg'));
org_img = org_img(500:900,470:770, :);
org_img = imresize(org_img,[256-search_width*2,256-search_width*2]);
imwrite(org_img,'results/input_image.jpg');

for channel = 1:3
    img = padarray(org_img(:,:,channel),[search_width search_width], ...
        'symmetric','both');
    noisy_image = imnoise(img,'gaussian', 0, sigma*sigma);
    noisy_img(:,:,channel) = noisy_image;
    basic_result(:,:,channel) = first_step(noisy_image, sigma, ...
        window_size, search_width, l2, l3, selection_number);
    basic_padded = padarray(basic_result(:,:,channel), ...
        [search_width search_width],'symmetric','both');
    final_result(:,:,channel) = second_step(noisy_image,basic_padded, ...
        sigma, window_size, search_width, l2, selection_number);
end
noisy_img = noisy_img(search_width+1:end-search_width, ...
    search_width+1:end-search_width,:);
imwrite(noisy_img,'results/noisy_image.jpg');
imwrite(uint8(basic_result),'results/res_phase1.jpg');
imwrite(uint8(final_result),'results/res_phase2.jpg');
img = org_img;
mypsnr(1) = psnr(uint8(double(noisy_img)*mean(img(:)) ... 
    /mean(noisy_img(:))),img);
mypsnr(2) = psnr(uint8(final_result*mean(img(:)) ...
    /mean(final_result(:))),(img));
disp(mypsnr);