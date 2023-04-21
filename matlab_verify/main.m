x_test = x(1:3,:);
y_test = y(1:3,:);

mean_x = [mean(x_test(1,:)); mean(x_test(2,:)); mean(x_test(3,:))];
mean_y = [mean(y_test(1,:)); mean(y_test(2,:)); mean(y_test(3,:))];

H = zeros(3);
for i = 1:4
    H = H + (y_test(:,i)-mean_y)*(x_test(:,i)-mean_x)';
end
H = H/4;
[U,D,V] = svd(H);

S = eye(3);
if det(H) < 0
    S(3,3) = -1;
end

test_R = U*S*V'
delta_x = 0;
for i = 1:4
    delta_x = delta_x + (x_test(:,i) - mean_x)'*(x_test(:,i) - mean_x);
end
delta_x = delta_x/4;
t = mean_y - test_R*mean_x

c = trace(D*S)/delta_x



