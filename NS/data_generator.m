% % clear all;
% % close all;
% 
% velxi = [0.8,0.9,1];rho = [0.9,1,1.1];mu = [500,600,700];
% N_par = 3;
% p = zeros(length(velxi),length(rho),length(mu),N_par);
% [p(:,:,:,1),p(:,:,:,2),p(:,:,:,3)] = ndgrid(velxi,rho,mu);
% size_p = size(p);
% size_p1 = size(p(:,:,:,1));
% 
% %% test data
% rng(4);
% rand_p = rand(3,1);
% velxi_test = 0.8+0.2*rand_p(1);rho_test = 0.9+0.2*rand_p(2);mu_test = 500+200*rand_p(3);
% w_test = Main(velxi_test,rho_test,1/mu_test);
% % velxi in [0.8,1], rho in [0.9,1.1] ,mu in [500,700] 
% 
% 
% 
% 
% [N_w,N_t] = size(w_test);
% w_train = zeros([N_w,N_t,size_p1]);
% w_train_dmd = zeros([N_w,N_t,size_p1]);
% e_train_rate_dmd = zeros(size_p1);
% 
% tic
% for i1 = 1:size_p1(1)
%     for i2 = 1:size_p1(2)
%         for i3 = 1:size_p1(3)
%             [velxi(i1),rho(i2),1/mu(i3)]
%             w_train(:,:,i1,i2,i3) = Main(velxi(i1),rho(i2),1/mu(i3));
%             visualize(w_train(:,:,i1,i2,i3));
%         end
%     end
% end
% toc

% visualize(w_train(:,:,1,1,1))
visualize(w_train(:,:,3,3,3))




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

%     lim = caxis;
    hcb = colorbar;
    title(hcb,{'$Q$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c1_pos,'FontSize',7);

    subplot('Position',pos2)
    imagesc(T(:,:,floor(size(T,3)/2)))
    colormap jet;
    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);

%     caxis(lim);
    hcb = colorbar;
    title(hcb,{'$Q$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c2_pos,'FontSize',7);

    subplot('Position',pos3)
    imagesc(T(:,:,end));
    colormap jet;
    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);

%     caxis(lim);
    hcb = colorbar;
    title(hcb,{'$Q$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c3_pos,'FontSize',7);
    
     print(gcf,'train2.png','-dpng','-r1500');  
end
