function [ key ] = keygen(c_vector)
% Converts a vector containing community identification into a string
% key for use in a struct "hash table"
% c_vector *must* be one-dimensional, and the value x of each entry
% corresponds to the node x being a member of the community. There is no
% node "0".
    v = zeros(1,max(c_vector));
    for i=1:length(c_vector)
        if c_vector(i)~=0
            v(c_vector(i))=1;
        end
    end
    if  mod(length(v),4)~=0
        v(4*ceil(length(v)/4))=0;
    end
    
    key=char(ones(1,length(v)/4)*48);
    for i=0:(length(v)/4 - 1)
        key(end-i)=dec2hex(v(4*i+1)+v(4*i+2)*2+v(4*i+3)*4+v(4*i+4)*8);
    end
    
    
end