% calculate dG 20200220 by chtung
clear
close all

TC = [25 250 300 350 400];
for t = 1:length(TC)
    T = TC(t) + 273; % Temperature
    element = 'VNbMoTaW';
    fraction_compound = 0:0.01:1.01;
    fraction_dissolve = 0:0.01:1.01;
    for i = 1:length(fraction_compound)
        for j = 1:length(fraction_compound)
            %% input
            
            % ------------------------------------------------------------
            % composition
            concentration = zeros(6,1);
            concentration(1) =  0.24; % V
            concentration(2) =  0.24; % Nb
            concentration(3) =  0.24; % Mo   input concentration here!
            concentration(4) =  0.24; % Ta
            concentration(5) =  0.24; % W
            concentration(6) = 98.8; % Cu
            % ------------------------------------------------------------
            
            X_c = fraction_compound(i); % compound fraction (0~1)
            X_ss = fraction_dissolve(j); % fraction of element not forming compound dissolve in Cu
            
            disp(['X_c = ',num2str(X_c)])
            disp(['X_ss = ',num2str(X_ss)])
            
            % ------------------------------------------------------------
            c_element(1) = 3;
            c_element(2) = 4;
            c_Stoichio(1) = 1; % e.g. Nb1Mo1, c_element(1) = 2 c_element(2) = 3
            c_Stoichio(2) = 1; %              c_Stoichio(1) = 1 c_Stoichio(2) = 1;
            % ------------------------------------------------------------
            
            concentration_c = zeros(6,1);
            concentration_c(c_element) = ...
                min(concentration(c_element)'./c_Stoichio).*c_Stoichio...
                *X_c; % compound composition
            concentration_ss = [(concentration(1:5) - concentration_c(1:5))*X_ss;concentration(6)]; % solid solution composition
            concentration_p = concentration-concentration_ss-concentration_c;
                        
            % thermodynamics
            
            delta_H = zeros(6);
            delta_H(:,1) = [0;-1;0;-1;-1;5];
            delta_H(:,2) = [-1;0;-6;0;-8;3];
            delta_H(:,3) = [0;-6;0;-5;0;19];
            delta_H(:,4) = [-1;0;-5;0;-7;2];
            delta_H(:,5) = [-1;-8;0;-7;0;22];
            delta_H(:,6) = [5;3;19;2;22;0];
            
            delta_H_c = [0.00, 0.00, 0.00, 0.00, 0.00;...
                0.00, 0.00, -9.40, 0.00, 0.00;...
                0.00, -9.40, 0.00,	-11.00, 0.00;...
                0.00, 0.00, -11.00, 0.00, -6.70;...
                0.00, 0.00,	0.00, -6.70, 0.00];
            
            delta_S_c = [0.00, 0.00, 0.00, 0.00, 0.00;...
                0.00, 0.00,	5.00, 0.00,	0.00;...
                0.00, 5.00,	0.00, 4.90,	0.00;...
                0.00, 0.00,	4.90, 0.00,	4.60;...
                0.00, 0.00,	0.00, 4.60,	0.00];
            
            %% enthalpy
            
            % precipitate
            concentration_p_n = concentration_p/sum(concentration_p,1);
            concentration_p_n(isnan(concentration_p_n)) = 0;
            concentration_p_n(isinf(concentration_p_n)) = 0;
            concentration_p_sum = sum(concentration_p,1);
            dH_p = sum(concentration_p_n'.*concentration_p_n.*delta_H,'all')/2*4*concentration_p_sum/100;
            
            % compound
            concentration_c_n = concentration_c/sum(concentration_c,1);
            concentration_c_n(isnan(concentration_c_n)) = 0;
            concentration_c_sum = sum(concentration_c,1);
            dH_c = min(concentration_c(c_element))/100*delta_H_c(c_element(1),c_element(2));
            
            % dissolve in cu
            concentration_ss_n = concentration_ss/sum(concentration_ss,1);
            concentration_ss_n(isnan(concentration_ss_n)) = 0;
            concentration_ss_sum = sum(concentration_ss,1);
            dH_ss = sum(concentration_ss_n'.*concentration_ss_n.*delta_H,'all')/2*4*concentration_ss_sum/100;
            % dH_ss = 0;
            
            dH(i,j) = dH_p + dH_c + dH_ss; % (KJ/mol)
            disp(['dH = ',num2str(dH(i,j))])
            %% entropy
            
            % precipitate
            concentration_p_n_log = concentration_p_n;
            concentration_p_n_log(concentration_p_n==0) = 1;
            dS_p = -8.314*sum(concentration_p_n_log.*log(concentration_p_n_log),'all')*concentration_p_sum/100;
            % dS_p = 0;
            
            % compound
            dS_c = min(concentration_c(c_element))/100*delta_S_c(c_element(1),c_element(2));
            
            % solid solution
            concentration_ss_n_log = concentration_ss_n;
            concentration_ss_n_log(concentration_ss_n==0) = 1;
            dS_ss = -8.314*sum(concentration_ss_n_log.*log(concentration_ss_n_log),'all')*concentration_ss_sum/100;
            
            dS(i,j) = dS_p + dS_c + dS_ss; % (J/mol/K)
            disp(['dS = ',num2str(dS(i,j))])
            
            dG(i,j) = dH(i,j)*1000 - T*dS(i,j);
            
            c_sum(i,j) = concentration_c_sum/sum(concentration(1:5));
            ss_sum(i,j) = (concentration_ss_sum-concentration(6))/sum(concentration(1:5));
            p_sum(i,j) = concentration_p_sum/sum(concentration(1:5));
            all_sum(i,j) = c_sum(i,j)+ss_sum(i,j)+p_sum(i,j);
        end
    end
    
    %% plot
    [XX0,YY0] = meshgrid(fraction_compound, fraction_compound);
    XX = XX0+YY0/2;
    YY = YY0*sqrt(3)/2;
    
    Xc_sum = c_sum+ss_sum/2;
    Yss_sum = ss_sum*sqrt(3)/2;
    
    figure(t);
    dG = real(dG);
    dG(dG==0) = NaN;
    
    hold on
    pcolor(Xc_sum,Yss_sum,dG);
    shading interp
    daspect([1 1 1])
    xlim([-0.1 1.1]);
    ylim([-0.01 sqrt(3)/2]);
    
    TX0 = [Xc_sum(end-1,1) Xc_sum(end-1,end-1) 1];
    TY0 = [Yss_sum(end-1,1) Yss_sum(end-1,end-1) 0];
    fill(TX0,TY0,'k','LineStyle','none')
    
    % grid
    xg1 = [0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0];
    yg1 = [0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9]*sqrt(3)/2/10;
    cg1 = [xg1;yg1];
    plot(cg1(1,:),cg1(2,:),'-','Color','k')
    
    cg2 = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)]*[1 0;0 -1]*cg1;
    plot(cg2(1,:),cg2(2,:),'-','Color','k')
    
    cg3 = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)]*(cg2)+[cos(pi/3); -sin(pi/3)]/2;
    plot(cg3(1,:),cg3(2,:),'-','Color','k')
    
    % hide    
    TXR = [0.5 1.5 1.5];
    TYR = [sqrt(3)/2 sqrt(3)/2 -sqrt(3)/2];
    fill(TXR,TYR,'w','LineStyle','none')
    
    TXL = [0.5 -0.5 -0.5];
    TYL = [sqrt(3)/2 sqrt(3)/2 -sqrt(3)/2];
    fill(TXL,TYL,'w','LineStyle','none')
    
    TXD = [0 1 1 0];
    TYD = [0.001 0.001 -1 -1];
    fill(TXD,TYD,'w','LineStyle','none')
    
    axis off
    set (gcf,'Color','#FFFFFF')
    set (gcf,'Position',[0,0,900,750])

    c = colorbar;
    c.LineWidth = 2;
    c.FontSize = 22;
    set(get(c,'title'),'string','\DeltaG (J/mol)','fontsize',28);
    
    TX2 = [0.5 0 1];
    TY2 = [sqrt(3)/2 0 0];
    f_triangle = fill(TX2,TY2,'w','LineWidth',2);
    set(f_triangle,'facealpha',0)
    
    %min dG
    dGr = dG;
    for d = 1:length(fraction_compound)
    dGr(d,length(fraction_compound)-d+1)=NaN;
    end
    dGmin(t) = min(dGr(1:end-1,1:end-1),[],'all');
    xgmin0 = c_sum(dGr==dGmin(t));
    ygmin0 = ss_sum(dGr==dGmin(t));
    xgmin = xgmin0+ygmin0/2;
    ygmin = ygmin0*sqrt(3)/2;
    
    plot(xgmin,ygmin,'kx','MarkerSize',22,'LineWidth',9)
    plot(xgmin,ygmin,'wx','MarkerSize',16,'LineWidth',4)
    c = ['(',num2str(xgmin0,'%0.2f'),', ',num2str(ygmin0,'%0.2f'),', ',num2str((1-xgmin0-ygmin0)*ceil(1-xgmin0-ygmin0),'%0.2f'),')'];
    
    if (xgmin+ygmin*sqrt(3)/2)<0.72
        text(xgmin+0.05,ygmin+0.025,c,'Color','w','FontSize',22,'fontname','Arial')
    elseif ygmin<0.6
        text(xgmin-0.45,ygmin+0.015,c,'Color','w','FontSize',22,'fontname','Arial')
    else
        t60 = text(xgmin-0.03,ygmin-0.08,c,'Color','w','FontSize',22,'fontname','Arial');
        set(t60,'Rotation',-60);
    end
    
    % tick label 1
    t1x = 0:0.1:1;
    t1y = zeros(length(t1x),1);
    t1c = arrayfun(@num2str, 0:0.1:1, 'UniformOutput', false);
    text(t1x,t1y-0.03,t1c,'Color','k','FontSize',22, 'HorizontalAlignment', 'center','fontname','Arial')
    
    % tick label 2
    t2x = (0:0.1:1)/2;
    t2y = (0:0.1:1)*sqrt(3)/2;
    t2c = arrayfun(@num2str, 1:-0.1:0, 'UniformOutput', false);
    text(t2x-0.014,t2y+0.025,t2c,'Color','k','FontSize',22, 'HorizontalAlignment', 'right','fontname','Arial')
    
    % tick label 3
    t3x = 1-(0:0.1:1)/2;
    t3y = (0:0.1:1)*sqrt(3)/2;
    t3c = arrayfun(@num2str, 0:0.1:1, 'UniformOutput', false);
    text(t3x+0.014,t3y+0.025,t3c,'Color','k','FontSize',22, 'HorizontalAlignment', 'left','fontname','Arial')
    
    % axis label 1
    a1x = 0.5;
    a1y = 0;
    a1c = 'Cu+5B''';
    text(a1x,a1y-0.15,a1c,'Color','k','FontSize',28,'HorizontalAlignment','center','VerticalAlignment', 'bottom','fontname','Arial')
    
    % axis label 2
    a2x = 0.25;
    a2y = sqrt(3)/4;
    a2c = 'Cu+5B';
    text(a2x-0.2,a2y+0.05,a2c,'Color','k','FontSize',28,'HorizontalAlignment','center','VerticalAlignment', 'bottom','fontname','Arial')
    
    % axis label 3
    a3x = 0.75;
    a3y = sqrt(3)/4;
    a3c = 'Cu-5B';
    text(a3x+0.2,a3y+0.05,a3c,'Color','k','FontSize',28,'HorizontalAlignment','center','VerticalAlignment', 'bottom','fontname','Arial')
    
    % temperature & element
    text(-0.05,0.925,[num2str(T),' K'],'FontSize',28,'HorizontalAlignment','left','fontname','Arial')
    text(-0.05,0.84,element,'FontSize',24,'HorizontalAlignment','left','fontname','Arial')
    
    colormap((parula+0.0)/1.0)
    
    dG_dissolve_no_compound(t) = dG(1,length(fraction_compound)-1);
    dG_precipitate_no_compound(t) = dG(1,1);
    dG_dissolve_compound(t) = dG(length(fraction_compound)-1,length(fraction_compound)-1);
    dG_precipitate_compound(t) = dG(length(fraction_compound)-1,1);
    
    %% relative yield
    concentration_c_th = zeros(6,1);
    concentration_c_th(c_element) = ...
        min(concentration(c_element)'./c_Stoichio).*c_Stoichio...
        *1; % theoretical yield concentration
    X_ty = sum(concentration_c_th)/sum(concentration(1:5)); % theoretical yield fraction
    ry = c_sum(1:101,1)/X_ty;
    dG_ry = dG(1:101,1);
    dG_ry_relative = dG_ry - dG_ry(1);
    
    % loop
    dG(isnan(dG)) = 0;
    dG_ry_1 = dG(1:101,1);
    dG_ry_2 = dG(101,1:101)';
    dG_ry_3 = dG(101:-1:1,101);
    dG_ry_4 = dG(1,101:-1:1)';
    dG_ry_loop(:,t) = [dG_ry_1(1:100);dG_ry_2(1:100);dG_ry_3(1:100);dG_ry_4];
    dG_ry_loop_relative(:,t) = dG_ry_loop(:,t) - dG_ry_loop(1);
end

% %% dG(T)
% figure(6);
% hold on
% box on
% plot(TC,[dG_dissolve_no_compound;dG_precipitate_no_compound;dG_dissolve_compound;dG_precipitate_compound]','-s','MarkerSize',10)
% plot(TC,dGmin','-ok')
% hold off
% xlabel('T (\circC)','fontsize',18);
% ylabel('\DeltaG (J/mol)','fontsize',18);
% set (gcf,'Position',[0,0,800,600])
% set (gca,'fontsize',15)
% legend('dissolve in Cu', 'HEA+Cu', 'dissolve in Cu+compound', 'HEA+Cu+compound', '\DeltaG_m_i_n','fontsize',12)
% xlim([0 425])

% plot_dg

%%
disp(' ')
disp(element)
disp(['dissolve in Cu,          dG(T) = ',num2str(dG_dissolve_no_compound,'%0.4f ')])
disp(['HEA+Cu,                  dG(T) = ',num2str(dG_precipitate_no_compound,'%0.4f ')])
disp(['dissolve in Cu+compound, dG(T) = ',num2str(dG_dissolve_compound,'%0.4f ')])
disp(['HEA+Cu+compound,         dG(T) = ',num2str(dG_precipitate_compound,'%0.4f ')])
disp(['                     dG_min(T) = ',num2str(dGmin,'%0.4f ')])