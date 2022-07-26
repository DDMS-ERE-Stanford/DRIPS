clear all;

rng(1);
rand_p = rand(3,100);
rand_p(1,:) = 0.8+0.2*rand_p(1,:);
rand_p(2,:) = 0.9+0.2*rand_p(2,:);
rand_p(3,:) = 500+200*rand_p(3,:);
e_dmd = zeros(1,100);
for i = 1:100
    e_dmd(i) = PROM(rand_p(1,i),rand_p(2,i),rand_p(3,i));i
end

save('test100_error_dmd.mat',"e_dmd","rand_p");



