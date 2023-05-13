function [vec_out] = convert2anglephase(vec_in,opt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

vec_out = vec_in;
for j = 1:size(vec_in,3)
    for i = 1:opt.Np
        vec_out(i,1,j) = norm([vec_in(i,1,j) vec_in(i,2,j)]);
        vec_out(i,2,j) = angle(vec_in(i,1,j)+1i*vec_in(i,2,j));
        vec_out(i,3,j) = vec_in(i,3,j);
    end
end

