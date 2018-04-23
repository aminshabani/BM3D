function res = dct3(data, mode)
if(nargin > 2 && mode == 'inverse')
    res = idct(data,'dim',3);
    for i = 1:size(data,3)
        res(:,:,i)=idct2(res(:,:,i));
    end
else
    res=zeros(size(data));
    for i=1:size(data,2)
        res(:,:,i)=dct2(data(:,:,i));
    end
    dct(res,'dim',3);
end
end