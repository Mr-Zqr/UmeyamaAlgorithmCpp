% generate 3 random 3D point coordinates
clear all
x = [rand(1,3)', rand(1,3)', rand(1,3)', rand(1,3)']
x(4,:) = [1 1 1 1];
T = [rotz(157)*rotx(60)*roty(65),rand(1,3)';0 0 0 1]
y = T*x
