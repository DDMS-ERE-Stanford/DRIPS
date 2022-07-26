% clear all;
% close all;
% function e_dmd = PROM(velxi_test,rho_test,mu_test)
velxi = [0.8,0.9,1];rho = [0.9,1,1.1];mu = [500,600,700];mu = 1./mu;
N_par = 3;
p = zeros(length(velxi),length(rho),length(mu),N_par);
[p(:,:,:,1),p(:,:,:,2),p(:,:,:,3)] = ndgrid(velxi,rho,mu);
size_p = size(p);
size_p1 = size(p(:,:,:,1));

%% test data
load('../test_data.mat');
mu_test = 1./mu_test;
% w_test = Main(velxi_test,rho_test,mu_test);


[N_w,N_t] = size(w_test);
N_q = 10;

load('../training_data.mat');
% w_train_dmd = zeros([N_w,N_t,size_p1]);
% e_train_dmd = zeros(size_p1);
ROB_train_dmd = zeros([N_w,N_q,size_p1]);
Kr_train_dmd = zeros([N_q,N_q,size_p1]);
Br_train_dmd = zeros([N_q,size_p1]);

tic
for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            [ROB_train_dmd(:,:,i1,i2,i3), Kr_train_dmd(:,:,i1,i2,i3),Br_train_dmd(:,i1,i2,i3)] = DMD(w_train(:,:,i1,i2,i3),N_q);
%             w_train_dmd(:,:,i1,i2,i3) = ROM_construction(w_train(:,:,i1,i2,i3),...
%                 ROB_train_dmd(:,:,i1,i2,i3), Kr_train_dmd(:,:,i1,i2,i3),Br_train_dmd(:,i1,i2,i3),N_t);
%             e_train_dmd(i1,i2,i3) = norm(w_train(:,:,i1,i2,i3)-w_train_dmd(:,:,i1,i2,i3))./norm(w_train(:,:,i1,i2,i3));
        end
    end
end
toc

%% interpolate dmd PROMs
tic
i0 = [1,1,1]; % reference point
ROB_test_interp_dmd = interpolate_ROB(ROB_train_dmd,i0,p,velxi_test,rho_test,mu_test,'linear');
[ROB_test_interp_general_dmd,Kr_test_interp_dmd,Br_test_interp_dmd] = interpolate_PROM...
    (ROB_train_dmd,ROB_test_interp_dmd,Kr_train_dmd,Br_train_dmd,i0,p,velxi_test,rho_test,mu_test,'spline');
wr_test_interpPROM_dmd = ROM_construction(w_test,ROB_test_interp_general_dmd,...
    Kr_test_interp_dmd,Br_test_interp_dmd,N_t);
toc
e_dmd = norm(w_test-wr_test_interpPROM_dmd)./norm(w_test);

% end
% visualize(w_test);
visualize(wr_test_interpPROM_dmd);
% visualize(log(abs(w_test-wr_test_interpPROM_dmd)));




function Tr = ROM_construction(T,ROB,Kr,Br,N_t)
q = ROB'*T;
for n = 2:N_t
   q(:,n) =  Kr*q(:,n-1)+Br;
end
Tr = ROB*q;
end
function ROB_test = interpolate_ROB(ROB_train,i0,p,u_test,kappa_test,ybar_test,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,1));
    Gamma_train = zeros([N_w,N_q,size_p1]);
    
    matrix_pre1 = (eye(N_w)-ROB_train(:,:,i0(1),i0(2),i0(3))*ROB_train(:,:,i0(1),i0(2),i0(3))');
    matrix_pre2 = ROB_train(:,:,i0(1),i0(2),i0(3))';
    
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                [U,Sigma,Z] = svd( matrix_pre1*...
            ROB_train(:,:,i1,i2,i3)*inv(matrix_pre2*ROB_train(:,:,i1,i2,i3)),'econ');
                Gamma_train(:,:,i1,i2,i3) = U*diag(atan(diag(Sigma)))*Z';
            end
        end
    end
    Gamma_test = zeros(N_w,N_q);
    parfor k = 1:N_w
        for s = 1:N_q
            Gamma_test(k,s) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_train(k,s,:,:,:),size_p1),u_test,kappa_test,ybar_test,interpolate_method);
        end
    end
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

function visualize(T)
    T = reshape(T,[100,200,size(T,2)]);
    pos1 = [0.08 0.2 0.19 0.7];
    pos2 = [0.41 0.2 0.19 0.7];
    pos3 = [0.74 0.2 0.19 0.7];
    c1_pos = [0.29 0.2 0.01 0.7];
    c2_pos = [0.62 0.2 0.01 0.7];
    c3_pos = [0.95 0.2 0.01 0.7];

    figure
    width = 6.6;     % Width in inches
    height = 2;    % Height in inches
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);

    subplot('Position',pos1);
    imagesc(T(:,:,1))
    colormap jet;
    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
%     title('Reference $t = 0$','FontUnits','points','interpreter','latex',...
%         'FontSize',10);
%     lim = caxis;
    hcb = colorbar;
    title(hcb,{'$Q_{DMD}$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c1_pos,'FontSize',7);

    subplot('Position',pos2)
    imagesc(T(:,:,floor(size(T,3)/2)))
    colormap jet;
    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
%     title('Reference $t$ half','FontUnits','points','interpreter','latex',...
%         'FontSize',10);
%     caxis(lim);
    hcb = colorbar;
    title(hcb,{'$Q_{DMD}$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c2_pos,'FontSize',7);

    subplot('Position',pos3)
    imagesc(T(:,:,end));
    colormap jet;
    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
%     title('Reference $t = 10000$','FontUnits','points','interpreter','latex',...
%         'FontSize',10);
%     caxis(lim);
    hcb = colorbar;
    title(hcb,{'$Q_{DMD}$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c3_pos,'FontSize',7);
    print(gcf,'test_dmd.png','-dpng','-r1500');  
end

