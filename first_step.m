function basic_result = first_step(noisy_image, sigma, ws, sw, l2, l3, sn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first step of the BM3D Image Denoising Method
% Inputs : noisy_image, sigma, ws(window_size), sw(search_window_width)
% l2(lambda_2D), l3(lambda_3d), sn(number of selections)
% Outputs : basic_result (The result image of the first step
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
            window = noisy_image(i-sw:i+sw , j-sw:j+sw);
            blocks = double(im2col(window, [ws ws], 'sliding'));
            dist = zeros(size(blocks,2),1);
            for k = 1:size(blocks,2)
                tmp = wthresh(blocks(:,center_block),'h',sigma*l2)- ...
                    wthresh(blocks(:,k),'h',sigma*l2);
                tmp = reshape(tmp, [ws ws]);
                tmp = norm(tmp,2)^2;
                dist(k) = tmp/(ws^2);
            end
            [~, inds] = sort(dist);
            inds = inds(1:sn);
            blocks = blocks(:, inds);
            p = zeros([ws ws sn]);
            for k = 1:sn
                p(:,:,k) = reshape(blocks(:,k), [ws ws]);
            end
            p = dct3(p);
            p = wthresh(p,'h',sigma*l3);
            wp = 1/sum(p(:)>0);
            p = dct3(p,'inverse');
            for k = 1:sn
                x = max(1,i-sw) + floor((center_block-1)/bpr);
                y = max(1,j-sw) + (mod(center_block-1,bpr));
                numerator(x:x+ws-1 , y:y+ws-1) = ...
                    numerator(x:x+ws-1 , y:y+ws-1) + (wp * p(:,:,k));
                denomerator(x:x+ws-1 , y:y+ws-1) = ...
                    denomerator(x:x+ws-1 , y:y+ws-1) + wp;
            end
        end
    end 
    basic_result = numerator./denomerator;
    basic_result = basic_result(sw+1:end-sw,sw+1:end-sw);
end