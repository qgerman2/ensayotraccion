close all, clear all, clc,

%% ===========Lectura de los datos========================

% Lectura de los seis documentos en el archivo xslx
% muestras{i} contiene los datos de la muestra i-esima
muestra = cell(1,6);
for i = 1:6
    % Lectura datos iniciales:
    % Ancho(mm) | Espesor(mm) | Largo calibrado(mm) |
    muestra{i}.initdata = xlsread('MS1_Laboratorio_DP1000.xlsx',i,'B24:B26','basic');
    % Lectura informacion ensayos
    muestra{i}.trans = xlsread('MS1_Laboratorio_DP1000.xlsx',i,'G40:G298','basic'); %Deformacion transversal (mm)
    muestra{i}.long = xlsread('MS1_Laboratorio_DP1000.xlsx',i,'H40:H298','basic');  % Deformacion longitudinal (mm)
    muestra{i}.force = xlsread('MS1_Laboratorio_DP1000.xlsx',i,'I40:I298','basic'); % Fuerza (N)
end

%% ============Calculos Preliminares===================

% Calculos area inicial(mm^2)
for i = 1:6
    muestra{i}.area = prod(muestra{i}.initdata(1:2)); % Ancho*Espesor(mm^2)
end
% Calculo esfuerzo nominal(MPa)
for i=1:6
    muestra{i}.stress = muestra{i}.area\muestra{i}.force(:,end); % Fuerza(N)/Area inicial(mm^2)
end

%% =========== Calculo constantes elasticas =====================
%Limite proporcional y modulo de young
epsilon = 20000;
for i = 1:6
    %dominio x para evaluar polinomio
    paso = 0.00001;
    x_dom = muestra{i}.long(1):paso:muestra{i}.long(end);
    %polinomio de curva esfuerzo - deformacion longitudinal
    [C, iu] = unique(muestra{i}.long);
    p = polyfit(muestra{i}.long(iu, :), muestra{i}.stress(iu, :), 6);
    %derivar polinomio
    p_d = polyder(p);
    p_d2 = polyder(p_d);
    %encontrar y = 0 en segunda derivada
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    xder2_cero = x_dom(zci(polyval(p_d2, x_dom)));
    %encontrar cuando la primera derivada se desvía mucho (limite proporcional)
    yder = polyval(p_d, xder2_cero(1));
    xder_desvio = find(polyval(p_d, x_dom) > yder - epsilon, 1, 'last');
    muestra{i}.proplim = polyval(p, x_dom(xder_desvio));
    muestra{i}.proplim_x = x_dom(xder_desvio);
    %encontrar pendiente en rango hasta limite proporcional (modulo de young)
    x_dom_proporcional = x_dom(1:xder_desvio);
    c = polyfit(x_dom_proporcional, polyval(p, x_dom_proporcional), 1);
    muestra{i}.young = c(1);
    %encontrar punto mas cercano entre recta con pendiente de young y datos (limite de fluencia)
    offset = 0.002;
    dist_corta = -1;
    for x = x_dom
        r = [x, muestra{i}.young * (x - offset)];
        pv = [x, polyval(p, x)];
        dist = norm(r - pv);
        if (dist_corta == -1 || dist < dist_corta)
            dist_corta = dist;
            muestra{i}.yieldlim = pv(2);
            muestra{i}.yieldlim_x = pv(1);
        end
    end
    %polinomio de curva deformacion transversal - deformacion longitudinal
    p = polyfit(muestra{i}.long(iu, :), muestra{i}.trans(iu, :), 6);
    %encontrar pendiente en rango hasta limite proporcional (modulo de poisson)
    c = polyfit(x_dom_proporcional, polyval(p, x_dom_proporcional), 1);
    muestra{i}.poisson = -c(1);
    %modulo de corte
    muestra{i}.corte = muestra{i}.young/(2*(1+muestra{i}.poisson));
end

%% ===========Graficos=============================

% Deformacion vs Esfuerzo en el mismo gráfico
color = [
    0.1, 0.1, 1;
    1, 0.1, 0.1;
    0.5, 0.5, 1;
    1, 0.5, 0.5;
    0.7, 0.7, 1;
    1, 0.7, 0.7;
];
figure('Name','Deformación vs Esfuerzo','Position',[0 0 1366 768])
for i = 1:6
    p = plot(muestra{i}.long, muestra{i}.stress,'LineWidth', 2, 'Color', color(i,:));
    h(i) = p(1);
    xlabel('$\epsilon_l $','Interpreter','latex','FontSize',16),ylabel('$\sigma [MPa]$','Interpreter','latex','FontSize',16),
    axis([0 ,1.1e-2,0, 800])
    grid on;
    hold on;
end
legend(h, {'Muestra 1 a $0^{\circ}$','Muestra 2 a $0^{\circ}$','Muestra 1 a $45^{\circ}$',...
    'Muestra 2 a $45^{\circ}$','Muestra 1 a $90^{\circ}$','Muestra 2 a $90^{\circ}$'},...
    'Interpreter','latex','Location','southeast','FontSize',14);
hold off;

% Deformacion longitudinal vs Esfuerzo
figure('Name','Deformación longitudial vs Esfuerzo','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    % Grafica principal 
    plot(muestra{i}.long,muestra{i}.stress);
    % Graficar puntos de importancia
    hold on
    p1 = plot(muestra{i}.yieldlim_x, muestra{i}.yieldlim, 'o');
    p2 = plot(muestra{i}.proplim_x, muestra{i}.proplim, 'o');
    p3 = line([0.002, muestra{i}.yieldlim_x], [0, muestra{i}.yieldlim],'LineStyle','--');
    title({'Deformacion longitudinal vs Esfuerzo';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\epsilon_l $','Interpreter','latex','FontSize',14),ylabel('$\sigma [MPa]$','Interpreter','latex','FontSize',14),
    legend([p1,p2,p3],{'Limite de fluencia','Limite proporcional','Criterio del Offset'},'Location','southeast','FontSize',14),
    axis([0 ,1.1e-2,0, 800])
    hold off
    grid on
end

% Deformacion longitudinal vs Deformacion transversal
figure('Name','Deformación longitudial vs Deformacion transversal','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    plot(muestra{i}.long,muestra{i}.trans);
    title({'Deformacion longitudinal vs Deformacion transversal';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\epsilon_l $','Interpreter','latex','FontSize',14),ylabel('$\epsilon_t $','Interpreter','latex','FontSize',14),
    axis([0 ,1.1e-2, -5e-3 ,0])
    grid on
end

% Desplegar constates elasticas

Constantes = {'Modulo de Young';'Modulo de Poisson';'Modulo de corte';'Limite proporcional';'Limite de fluencia'};
Muestra1_0grad = [muestra{1}.young ; muestra{1}.poisson ; muestra{1}.corte ; muestra{1}.proplim ; muestra{1}.yieldlim];
Muestra2_0grad = [muestra{2}.young ; muestra{2}.poisson ; muestra{2}.corte ; muestra{2}.proplim ; muestra{2}.yieldlim];
Muestra1_45grad = [muestra{3}.young ; muestra{3}.poisson ; muestra{3}.corte ; muestra{3}.proplim ; muestra{3}.yieldlim];
Muestra2_45grad = [muestra{4}.young ; muestra{4}.poisson ; muestra{4}.corte ; muestra{4}.proplim ; muestra{4}.yieldlim];
Muestra1_90grad = [muestra{5}.young ; muestra{5}.poisson ; muestra{5}.corte ; muestra{5}.proplim ; muestra{5}.yieldlim];
Muestra2_90grad = [muestra{6}.young ; muestra{6}.poisson ; muestra{6}.corte ; muestra{6}.proplim ; muestra{6}.yieldlim];
T = table(Constantes,Muestra1_0grad,Muestra1_0grad,Muestra1_45grad,Muestra2_45grad,Muestra1_90grad,Muestra2_90grad);
disp(T)

%% ========= Graficos Esfuerzo vs Ctes elasticas ==============================

% Calculo constantes
for i = 1:6
    ind = find(muestra{i}.stress > muestra{i}.yieldlim,1); 
    muestra{i}.youngsec = muestra{i}.stress(2:ind)./muestra{i}.long(2:ind);             % Modulo de Young Secante
    muestra{i}.youngtg = (muestra{i}.stress(3:ind) - muestra{i}.stress(2:ind-1))./...   % Modulo de Young Tangente
                          (muestra{i}.long(3:ind) - muestra{i}.long(2:ind-1));
    muestra{i}.poissonsec = -muestra{i}.trans(2:ind)./muestra{i}.long(2:ind);            % Modulo de Poisson secante
end
% Graficos
% Esfuerzo vs Modulo de Young secante
figure('Name','Esfuerzo vs Modulo de Young secante','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    plot(muestra{i}.stress(1:length(muestra{i}.youngsec)),muestra{i}.youngsec,'o','MarkerSize',2)
    title({'Esfuerzo vs Modulo de Young Secante';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\sigma[MPa] $','Interpreter','latex','FontSize',14),ylabel('$E[MPa]$','Interpreter','latex','FontSize',14),
    axis([0,muestra{i}.yieldlim,1e5,3e5])
    grid on
end

% Esfuerzo vs Modulo de Young tangente
figure('Name','Esfuerzo vs Modulo de Young Tangente','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    plot(muestra{i}.stress(1:length(muestra{i}.youngtg)),muestra{i}.youngtg,'o','MarkerSize',2)
    title({'Esfuerzo vs Modulo de Young Tangente';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\sigma[MPa] $','Interpreter','latex','FontSize',14),ylabel('$E[MPa]$','Interpreter','latex','FontSize',14),
    axis([0,muestra{i}.yieldlim,0,3e5])
    grid on
end

% Esfuerzo vs Modulo de Young discretos
figure('Name','Esfuerzo vs Modulo de Young discretos','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    hold on
    p1 = plot(muestra{i}.stress(1:length(muestra{i}.youngtg)),muestra{i}.youngtg,'ob','MarkerSize',2);
    p2 = plot(muestra{i}.stress(1:length(muestra{i}.youngsec)),muestra{i}.youngsec,'or','MarkerSize',2);
    title({'Esfuerzo vs Modulo de Young Discretos';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\sigma[MPa] $','Interpreter','latex','FontSize',14),ylabel('$E[MPa]$','Interpreter','latex','FontSize',14),
    legend([p1,p2],{'Modulo de Young tangente','Modulo de Young secante'},'FontSize',14)
    axis([100,muestra{i}.yieldlim,0,3e5])
    grid on
end

%% Esfuerzo vs Modulo de Poisson
figure('Name','Esfuerzo vs Modulo de Poisson','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    plot(muestra{i}.stress(1:length(muestra{i}.poissonsec)),muestra{i}.poissonsec,'o','MarkerSize',2);
    title({'Esfuerzo vs Modulo de Poisson';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\sigma[MPa] $','Interpreter','latex','FontSize',14),ylabel('$\nu[MPa]$','Interpreter','latex','FontSize',14),
    axis([0,muestra{i}.yieldlim,0,.5])
    grid on
end