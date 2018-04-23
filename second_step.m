function result = second_step(noisy_image, basic_res, sigma, ws,sw, l2, sn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first step of the BM3D Image Denoising Method
% Inputs : noisy_image, basic_res(The result of the previous step), sigma, 
% ws(window_size), sw(search_window_width), l2(lambda_2D), l3(lambda_3d), 
% sn(number of selections)
% Output : basic_result (The result image of the first step
% For more information, you can see the following two papers:
% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4271520
% http://www.ipol.im/pub/art/2012/l-bm3d/article.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    image_size = size(noisy_image);
    numerator = zeros(size(noisy_image));
    denomerator = zeros(size(noisy_image));
    bpr = (sw*2+1) - ws + 1; % number of blocks per row/col
    center_block = bpr^2 / 2 + bpr/2 + 1;
    for i = sw+1:image_size(1)-sw
        for j = sw+1:image_size(2)-sw
            window  = noisy_image (i-sw:i+sw, j-sw:j+sw);
            window2 = basic_res(i-sw:i+sw, j-sw:j+sw);
            blocks  = double(im2col(window, [ws ws], 'sliding'));
            blocks2 = double(im2col(window2, [ws ws], 'sliding'));
            dist = zeros(size(blocks,2),1);
            for k = 1:size(blocks,2)
                tmp = wthresh(blocks2(:,center_block),'h',sigma*l2)- ...
                    wthresh(blocks2(:,k),'h',sigma*l2);
                tmp = reshape(tmp, [ws ws]);
                tmp = norm(tmp,2)^2;
                dist(k) = tmp/(ws^2);
            end
            [~, I] = sort(dist);
            I = I(1:sn);
            blocks = blocks(:, I);
            blocks2 = blocks2(:, I);
            p = zeros([ws ws sn]);
            basic_p = zeros([ws ws sn]);
            for k = 1:sn
                p(:,:,k) = reshape(blocks(:,k), [ws ws]);
                basic_p(:,:,k) = reshape(blocks2(:,k), [ws ws]);
            end
            basic_p = dct3(basic_p);
            wp = zeros(sn,1);
            for k = 1:sn
                tmp = basic_p(:,:,k);
                tmp = norm(tmp,1).^2;
                wp(k) = tmp/(tmp+(sigma^2));
            end
            p = dct3(p);
            for k = 1:sn
                p(:,:,k) = p(:,:,k)*wp(k);
            end
            p = dct3(p,'inverse');
            wp = 1/sum(wp(:).^2);
            for k = 1:sn
                q = p(:,:,k);
                x = max(1,i-sw) + floor((center_block-1)/bpr);
                y = max(1,j-sw) + (mod(center_block-1,bpr));
                numerator(x:x+ws-1 , y:y+ws-1) = ...
                    numerator(x:x+ws-1 , y:y+ws-1) + (wp * q);
                denomerator(x:x+ws-1 , y:y+ws-1) = ...
                    wp + denomerator(x:x+ws-1 , y:y+ws-1);
            end
        end
    end 
    result = numerator./denomerator;
    result = result(sw+1:end-sw,sw+1:end-sw);
end