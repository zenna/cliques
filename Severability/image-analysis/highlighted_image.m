function [ out ] = highlighted_image( I, pixels )
%HIGHLIGHTED_IMAGE Generates a new image with particular pixels highlighted
%   I = original image
%   pixels = pixels to highlight in node format



[u,v,w]=size(I);

dim2 = ceil(pixels./u);
dim1 = mod(pixels,u);
dim1(dim1==0)=u; % Make sure that we have no 0 indices.

mask=false(u,v);
for i=1:length(pixels)
    mask(dim1(i),dim2(i))=1;
end

if w==1
    out_red = I;
    out_green = I;
    out_blue = I;
elseif w==3
    out_red = I(:,:,1);
    out_green=I(:,:,2);
    out_blue= I(:,:,3);
else
    error('Input error in image variable');
end

diceroll=randi(5);
c_choice=[255 255 0; 255 0 255; 0 255 255; 0 255 0; 0 0 255];

%diceroll=4;
out_red(mask) = (c_choice(diceroll,1)/2+out_red(mask)/2);
out_green(mask) = (c_choice(diceroll,2)/2+out_green(mask)/2);
out_blue(mask) = (c_choice(diceroll,3)/2+out_blue(mask)/2);
%out_red(mask)=c_choice(diceroll,1);
%out_green(mask)=c_choice(diceroll,2);
%out_blue(mask)=c_choice(diceroll,3);

out = cat(3, out_red, out_green, out_blue);
%imshow(out);

out = cast(out,'uint8');

end

