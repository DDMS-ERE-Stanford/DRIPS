% clear all;
% function [e_train_rate_dmd,e_rate_dmd,t_offline,t_online] = linear_stanford_rate(order)
order = 8;
N_t = 500;
t = linspace(0,5000,N_t);
 
c = [1,1.5,2];sL = [1,2,3];sR = [1,2,3];
N_par = 3;
p = zeros(length(c),length(sL),length(sR),N_par);
[p(:,:,:,1),p(:,:,:,2),p(:,:,:,3)] = ndgrid(c,sL,sR);
size_p = size(p);
size_p1 = size(p(:,:,:,1));
 
%% test data
rng(4);
rand_p = rand(3,1);
c_test = 1+rand_p(1);sL_test = 1+2*rand_p(2);sR_test = 1+2*rand_p(3);
w_test = stanford(c_test,sL_test,sR_test,t);
y_test = observable(w_test,order);


[N_w,~] = size(w_test);
w_train = zeros([N_w,N_t,size_p1]);
w_train_dmd = zeros([N_w,N_t,size_p1]);
e_train_rate_dmd = zeros(size_p1);

% fprintf('generate data time:');
% tic
% for i1 = 1:size_p1(1)
%     for i2 = 1:size_p1(2)
%         for i3 = 1:size_p1(3)
%             w_train(:,:,i1,i2,i3) = stanford(c(i1),sL(i2),sR(i3),t);
%         end
%     end
% end
% toc
load('training_data.mat');

y_train = zeros([order,N_t,size_p1]);
y_train_dmd = zeros([order,N_t,size_p1]);
N_q = order;
ROB_train_dmd = zeros([order,N_q,size_p1]);
Kr_train_dmd = zeros([N_q,N_q,size_p1]);
Br_train_dmd = zeros([N_q,size_p1]);

fprintf('offline time:');
tic
for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            y_train(:,:,i1,i2,i3) = observable(w_train(:,:,i1,i2,i3),order);
            [ROB_train_dmd(:,:,i1,i2,i3), Kr_train_dmd(:,:,i1,i2,i3),Br_train_dmd(:,i1,i2,i3)] = DMD(y_train(:,:,i1,i2,i3),N_q);
        end
    end
end
toc
t_offline = toc;

for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            y_train_dmd(:,:,i1,i2,i3) = ROM_construction(y_train(:,:,i1,i2,i3),ROB_train_dmd(:,:,i1,i2,i3),...
                Kr_train_dmd(:,:,i1,i2,i3),Br_train_dmd(:,i1,i2,i3),N_t);
            w_train_dmd(:,:,i1,i2,i3) = y_train_dmd(1,:,i1,i2,i3)/2;
            e_train_rate_dmd(i1,i2,i3) = norm(w_train(:,:,i1,i2,i3)-w_train_dmd(:,:,i1,i2,i3))./norm(w_train(:,:,i1,i2,i3));
        end
    end
end
    

%% interpolate dmd PROMs
i0 = [1,1,1]; % reference point
Gamma_train_dmd = compute_Gamma_train(ROB_train_dmd,i0,p);
toc
t_offline = toc;

fprintf('online time:');
tic
ROB_test_interp_dmd = interpolate_ROB(Gamma_train_dmd,ROB_train_dmd,i0,p,c_test,sL_test,sR_test,'linear');
[ROB_test_interp_general_dmd,Kr_test_interp_dmd,Br_test_interp_dmd] = interpolate_PROM...
    (ROB_train_dmd,ROB_test_interp_dmd,Kr_train_dmd,Br_train_dmd,i0,p,c_test,sL_test,sR_test,'spline');
yr_test_interpPROM_dmd = ROM_construction(y_test,ROB_test_interp_general_dmd,...
    Kr_test_interp_dmd,Br_test_interp_dmd,N_t);
wr_test_interpPROM_dmd = yr_test_interpPROM_dmd(1,:)/2;
toc
t_online = toc;
e_rate_dmd = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);
% end



 
function w = stanford(s_h,sL,sR,tlist)
    x = [167,167,349,349,462,462,260,167,167,258,561,652,652,455,455,341,341,562,652,652,547,269];
    y = [146,295,295,226,226,327,327,417,691,782,782,690,542,542,603,603,508,508,418,146,42,42];
    y = y-10;
 
    thermalmodel = createpde('thermal','transient');
    R1 =  [3;4;0;800;800;0;0;0;800;800];
    n = length(x);
    P1 = [2;n;x';y'];
    R1 = [R1;zeros(length(P1)-length(R1),1)];
    gd = [R1,P1]; 
    sf = 'R1-P1';
    ns = char('R1','P1');
    ns = ns';
    g = decsg(gd,sf,ns);
    geometryFromEdges(thermalmodel,g);
    thermalProperties(thermalmodel,'ThermalConductivity',1,...
                                   'MassDensity',1,...
                                   'SpecificHeat',s_h);
    internalHeatSource(thermalmodel,0);
 
    %% BC
    thermalBC(thermalmodel,'Edge',1,'Temperature',sL);
    thermalBC(thermalmodel,'Edge',23,'Temperature',sR);
    thermalBC(thermalmodel,'Edge',3:22,'Temperature',3);
    thermalBC(thermalmodel,'Edge',25:26,'Temperature',3);
 
    %% IC
    thermalIC(thermalmodel,1);
 
    msh = generateMesh(thermalmodel,'Hmax',20);
 
    pos1 = [0.08 0.2 0.19 0.7];
    pos2 = [0.41 0.2 0.19 0.7];
    pos3 = [0.74 0.2 0.19 0.7];
    c1_pos = [0.29 0.2 0.01 0.7];
    c2_pos = [0.62 0.2 0.01 0.7];
    c3_pos = [0.95 0.2 0.01 0.7];
 
 
%     thermalmodel.SolverOptions.ReportStatistics ='on';
    result = solve(thermalmodel,tlist);
    Q = zeros(length(tlist),26);
    for i = 1:26
        Q(:,i) = evaluateHeatRate(result,'Edge',i);
    end
    w = Q(:,3)';
end
 
function Tr = ROM_construction(T,ROB,Kr,Br,N_t)
q = ROB'*T;
for n = 2:N_t
   q(:,n) =  Kr*q(:,n-1)+Br;
end
Tr = ROB*q;
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
function [ROB_test,Kr_test,Br_test] = interpolate_PROM(ROB_train,ROB_test,Kr_train,Br_train,i0,p,u_test,kappa_test,ybar_test,interpolate_method)
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
    [ROB,~,~] = svd(X1,'econ');
    ROB = ROB(:,1:r);
    [V,Sigma,Z] = svd(X1_tilde,'econ');
    V = V(:,1:r+1); Sigma = Sigma(1:r+1,1:r+1); Z = Z(:,1:r+1);
    V1 = V(1:size(X1,1),:);V2 = V(end,:);
 
    Kr = ROB'*X2*Z/Sigma*V1'*ROB;
    Br = ROB'*X2*Z/Sigma*V2';
    
%     A_tilde = X2*Z/Sigma*V';
%     K = A_tilde(:,1:end-1);
%     B = A_tilde(:,end);
%     Kr = ROB'*K*ROB;
%     Br = ROB'*B;
end
 
function y = observable(x,order)
y = hermiteH(1,x);
for n = 2:order
    y = [y;hermiteH(n,x)];
end
end
