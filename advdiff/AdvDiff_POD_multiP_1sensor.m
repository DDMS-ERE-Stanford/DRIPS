close all;
clear all;

dt = 0.01;
t = 0:dt:1;
N_t = length(t);
u = [1000,2000,3000]; kappa = [100,200,300]; ybar = [0.36,0.39,0.42];
N_par = 3;
p = zeros(length(u),length(kappa),length(ybar),N_par);
[p(:,:,:,1),p(:,:,:,2),p(:,:,:,3)] = ndgrid(u,kappa,ybar);
size_p = size(p);
size_p1 = size(p(:,:,:,1));

% test data
rng(4);
rand_p = rand(3,1);
u_test = 1000+2000*rand_p(1);kappa_test = 100+200*rand_p(2);ybar_test = 0.36+0.06*rand_p(3);
[w_test,bc_test,K_test,B_test] = adv_diff([u_test;kappa_test;ybar_test],t);

% loc = [20,38];
% timewindow = 1:101;
% obs_test =  sensor_assign(vec_to_imag(w_test,bc_test),loc,timewindow);

[N_w,~] = size(w_test);
N_q = 10;
w_train = zeros([N_w,N_t,size_p1]);
C = zeros(1,N_w);C(37*floor(sqrt(N_w))+19) = 1;

wr_train = zeros([N_w,N_t,size_p1]);
obs_r_train =  zeros([1,N_t,size_p1]);
K = zeros([N_w,N_w,size_p1]);
B = zeros([N_w,1,size_p1]);

ROB_train = zeros([N_w,N_q,size_p1]);
Kr_train = zeros([N_q,N_q,size_p1]);
Br_train = zeros([N_q,size_p1]);
Cr_train = zeros([1,N_q,size_p1]);

tic
for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            [w_train(:,:,i1,i2,i3),bc,K(:,:,i1,i2,i3),B(:,:,i1,i2,i3)] = adv_diff([u(i1);kappa(i2);ybar(i3)],t);
%             visualize(vec_to_imag(w_train(:,:,i1,i2,i3),bc));
        end
    end
end
toc

tic
for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            [U,~,~] = svd(w_train(:,:,i1,i2,i3),'econ');
            ROB_train(:,:,i1,i2,i3) = U(:,1:N_q);
            Kr_train(:,:,i1,i2,i3) = ROB_train(:,:,i1,i2,i3)'*K(:,:,i1,i2,i3)*ROB_train(:,:,i1,i2,i3);
            Br_train(:,i1,i2,i3) = ROB_train(:,:,i1,i2,i3)'*B(:,:,i1,i2,i3);
            Cr_train(:,:,i1,i2,i3) = C*ROB_train(:,:,i1,i2,i3);
            obs_r_train(1,:,i1,i2,i3) = ROM_construction(w_train(:,:,i1,i2,i3),ROB_train(:,:,i1,i2,i3),...
   Kr_train(:,:,i1,i2,i3) ,Br_train(:,i1,i2,i3),Cr_train(:,:,i1,i2,i3),N_t);
        end
    end
end
%% interpolate ROBs
i0 = [1,1,1]; % reference point
Gamma_train = compute_Gamma_train(ROB_train,i0,p);
toc

tic
ROB_test_interp = interpolate_ROB(Gamma_train,ROB_train,i0,p,u_test,kappa_test,ybar_test,'linear');
[ROB_test_interp_general,Kr_test_interp,Br_test_interp,Cr_test_interp] = interpolate_PROM...
    (ROB_train,ROB_test_interp,Kr_train,Br_train,Cr_train,i0,p,u_test,kappa_test,ybar_test,'spline');
obs_r_test_interpPROM = ROM_construction(w_test,ROB_test_interp_general,...
    Kr_test_interp,Br_test_interp,Cr_test_interp,N_t);
toc

figure
width = 5;     % Width in inches
height = 4;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
plot(t,C*w_test,'k-','LineWidth',1);
plot(t,obs_r_test_interpPROM,'r-.','LineWidth',1);

e_pod = norm(C*w_test-obs_r_test_interpPROM)./norm(C*w_test);
% end

function [T_vec,BC,K,bias] = adv_diff(p,t)
dt= t(2)-t(1);
U1 = p(1);
kappa = p(2);
ybar = p(3);
N = 75;
x = linspace(0,1,N);y = x;
h = x(2)-x(1);

internalPoints = N-2;
%% Diffusion matrix
d   = -kappa*ones(internalPoints,1);
spd = spdiags([d -2*d d], -1:1,internalPoints,internalPoints);
Iz  = speye(internalPoints);
D  = kron(Iz,spd)+kron(spd,Iz);
%% Advection matrix
a   = U1*h*ones(internalPoints,1);
spa = spdiags([-a a], -1:0,internalPoints,internalPoints);
A = kron(Iz,spa);

M = A+D;
%% Modification for T_{:,1}, T_{:,N-2}
M(1:internalPoints,1:internalPoints) = ...
    M(1:internalPoints,1:internalPoints) - kappa*Iz;
M(end+1-internalPoints:end,end+1-internalPoints:end) = ...
    M(end+1-internalPoints:end,end+1-internalPoints:end)-kappa*Iz;
%% No Modification for T_{1,:}
%% Modification for T_{N-2,:}
M(internalPoints:internalPoints:end,internalPoints:internalPoints:end) = ...
    M(internalPoints:internalPoints:end,internalPoints:internalPoints:end)...
    -kappa*Iz;

%% Dirichlet BC
BC = 300*ones(internalPoints,1);
BC(25:49) = 300+325*(sin(3*pi*abs([25:49]*h-ybar))+1);

b = sparse(1:internalPoints:internalPoints^2+1-internalPoints,1,...
    (U1*h+kappa)*BC,internalPoints^2,1);

T_vec = zeros((N-2)*(N-2),length(t));
T_vec(:,1) = 300* ones((N-2)*(N-2),1);
K = inv(eye((N-2)*(N-2))+dt*M);
bias = (eye((N-2)*(N-2))+dt*M)\(dt*b);
for iter = 1:length(t)-1 
    T_vec(:,iter+1) = K*T_vec(:,iter)+bias;
end    
end

% function Tr = ROM_construction(T,ROB,Kr,Br,N_t)
% q = ROB'*T;
% for n = 2:N_t
%    q(:,n) =  Kr*q(:,n-1)+Br;
% end
% Tr = ROB*q;
% end

function obs_r = ROM_construction(T,ROB,Kr,Br,Cr,N_t)
q = ROB'*T;
for n = 2:N_t
   q(:,n) =  Kr*q(:,n-1)+Br;
end
obs_r = Cr*q;
end


function Gamma_train = compute_Gamma_train(ROB_train,i0,p)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,1));
    Gamma_train = zeros([N_w,N_q,size_p1]);
    matrix_pre1 = (eye(N_w)-ROB_train(:,:,i0(1),i0(2),i0(3))*ROB_train(:,:,i0(1),i0(2),i0(3))');
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                [U,Sigma,Z] = svd(matrix_pre1*...
            ROB_train(:,:,i1,i2,i3)*inv(ROB_train(:,:,i0(1),i0(2),i0(3))'*ROB_train(:,:,i1,i2,i3)),'econ');
                Gamma_train(:,:,i1,i2,i3) = U*diag(atan(diag(Sigma)))*Z';
            end
        end
    end
end

function ROB_test = interpolate_ROB(Gamma_train,ROB_train,i0,p,u_test,kappa_test,ybar_test,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,1));
%     Gamma_train = zeros([N_w,N_q,size_p1]);
%     
%     
%     matrix_pre1 = (eye(N_w)-ROB_train(:,:,i0(1),i0(2),i0(3))*ROB_train(:,:,i0(1),i0(2),i0(3))');
%     for i1 = 1:size_p1(1)
%         for i2 = 1:size_p1(2)
%             for i3 = 1:size_p1(3)
%                 [U,Sigma,Z] = svd(matrix_pre1*...
%             ROB_train(:,:,i1,i2,i3)*inv(ROB_train(:,:,i0(1),i0(2),i0(3))'*ROB_train(:,:,i1,i2,i3)),'econ');
%                 Gamma_train(:,:,i1,i2,i3) = U*diag(atan(diag(Sigma)))*Z';
%             end
%         end
%     end

    Gamma_test = zeros(N_w,N_q);
    parfor k = 1:N_w
        for s = 1:N_q
            Gamma_test(k,s) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_train(k,s,:,:,:),size_p1),u_test,kappa_test,ybar_test,interpolate_method);
        end
    end
%     Gamma_test = zeros(N_w*N_q,1);
%     parfor k = 1:N_w*N_q
%         Gamma_test(k) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
%                 reshape(Gamma_train(k,:,:,:),size_p1),u_test,kappa_test,ybar_test,interpolate_method);
%     end
%     Gamma_test = reshape(Gamma_test,[N_w,N_q]);
    [U,Sigma,Z] = svd(Gamma_test,'econ');
    ROB_test = ROB_train(:,:,i0(1),i0(2),i0(3))*Z*diag(cos(diag(Sigma)))+U*diag(sin(diag(Sigma)));
end

function [ROB_test,Kr_test,Br_test,Cr_test] = interpolate_PROM(ROB_train,ROB_test,Kr_train,Br_train,Cr_train,i0,p,u_test,kappa_test,ybar_test,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,1));
    %% Step A
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                P = ROB_train(:,:,i1,i2,i3)'* ROB_train(:,:,i0(1),i0(2),i0(3));
                [U,~,Z] = svd(P,'econ');
                Q = U*Z';
                Kr_train(:,:,i1,i2,i3) = Q'*Kr_train(:,:,i1,i2,i3)*Q;
                Br_train(:,i1,i2,i3) = Q'*Br_train(:,i1,i2,i3);
                Cr_train(:,:,i1,i2,i3) = Cr_train(:,:,i1,i2,i3)*Q;
            end
        end
    end
    P = ROB_test'* ROB_train(:,:,i0(1),i0(2),i0(3));
    [U,~,Z] = svd(P,'econ');
    Q = U*Z';
    ROB_test = ROB_test*Q;

    %% Step B
    Gamma_K_train = manifold_log(Kr_train(:,:,i0(1),i0(2),i0(3)),Kr_train,'real');
    Gamma_K_test = zeros(N_q,N_q);
    for k = 1:N_q
        for s = 1:N_q
            Gamma_K_test(k,s) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_K_train(k,s,:,:,:),size_p1),u_test,kappa_test,ybar_test,interpolate_method);
        end
    end
    Kr_test = manifold_exp(Kr_train(:,:,i0(1),i0(2),i0(3)),Gamma_K_test,'real');

    Gamma_B_train = manifold_log(Br_train(:,i0(1),i0(2),i0(3)),Br_train,'real');
    Gamma_bias_test = zeros(N_q,1);
    for k = 1:N_q
        Gamma_bias_test(k) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_B_train(k,:,:,:),size_p1),u_test,kappa_test,ybar_test,interpolate_method);
    end
    Br_test = manifold_exp(Br_train(:,i0(1),i0(2),i0(3)), Gamma_bias_test, 'real');
    
    Gamma_C_train = manifold_log(Cr_train(:,:,i0(1),i0(2),i0(3)),Cr_train,'real');
    Gamma_C_test = zeros(1,N_q);
        for s = 1:N_q
            Gamma_C_test(1,s) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_C_train(1,s,:,:,:),size_p1),u_test,kappa_test,ybar_test,interpolate_method);
        end
    Cr_test = manifold_exp(Cr_train(:,:,i0(1),i0(2),i0(3)),Gamma_C_test,'real');
end


function Gamma = manifold_log(X,Y,matrix_type)
if strcmp(matrix_type,'real')
    Gamma = Y-X;
end
if strcmp(matrix_type,'nonsingular')
    Gamma = log(Y*inv(X));
end
if strcmp(matrix_type,'SPD')
    R = chol(X);
    Gamma = log(inv(R)*Y*inv(R));
end
end

function Y = manifold_exp(X,Gamma,matrix_type)
if strcmp(matrix_type,'real')
    Y = X+Gamma;
end
if strcmp(matrix_type,'nonsingular')
    Y = exp(Gamma)*X;
end
if strcmp(matrix_type,'SPD')
    R = chol(X);
    Y = R*exp(Gamma)*R;
end
end

function [ROB, Kr,Br] = DMD(X,r)
    X1 = X(:,1:end-1);
    X1_tilde = [X1;ones(1,size(X1,2))];
    X2 = X(:,2:end);
    [V,Sigma,Z] = svd(X1_tilde,'econ');
    V = V(:,1:r+1); Sigma = Sigma(1:r+1,1:r+1); Z = Z(:,1:r+1);
    A_tilde = X2*Z/Sigma*V';
    K = A_tilde(:,1:end-1);
    B = A_tilde(:,end);
    [ROB,~,~] = svd(X1,'econ');
    ROB = ROB(:,1:r);
    Kr = ROB'*K*ROB;
    Br = ROB'*B;
end



function T_imag = vec_to_imag(T_vec,BC)
    [internalPoints_sq,Nt] = size(T_vec);
    N = sqrt(internalPoints_sq)+2;
    T_imag = zeros(N,N,Nt);
    T_in = reshape(T_vec,[N-2,N-2,Nt]);
    T_imag(2:N-1,2:N-1,:) = T_in;
    T_imag(1,2:end-1,:) = reshape(repmat(BC',Nt,1)',[1,N-2,Nt]);
    T_imag(end,:,:) = T_imag(end-1,:,:);
    T_imag(:,1,:) = T_imag(:,2,:);
    T_imag(:,end,:) = T_imag(:,end-1,:);
end


function obs = sensor_assign(T_imag,loc,timewindow)
obs = T_imag(loc(1),loc(2),timewindow);
obs = reshape(obs,[size(timewindow),1]);
end






