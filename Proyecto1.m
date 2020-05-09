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

%% ===========Graficos=============================

% Deformacion longitudinal vs Esfuerzo
figure('Name','Deformacion longitudial vs Esfuerzo','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    plot(muestra{i}.long,muestra{i}.stress);
    title({'Deformacion longitudinal vs Esfuerzo';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\epsilon_l $','Interpreter','latex','FontSize',14),ylabel('$\sigma [MPa]$','Interpreter','latex','FontSize',14),
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
    grid on
end

%% =========== Calculo constantes elasticas =====================

% Aqui se considera el criterio del offset, para calcular las constantes por minimos cuadrados
for i = 1:6
    muestra{i}.proplim = find(muestra{i}.long > (0.002),1) - 1;                                     % Indice para el proplime de fluencia
    c = polyfit(muestra{i}.long(1:muestra{i}.proplim),muestra{i}.stress(1:muestra{i}.proplim),1);   % Ajuste polinomial lineal
    muestra{i}.young = c(1);                                                                        % Módulo de Young
    c = polyfit(muestra{i}.long(1:muestra{i}.proplim),muestra{i}.trans(1:muestra{i}.proplim),1);    % Ajuste polinomial lineal
    muestra{i}.poisson = -c(1);                                                                     % Módulo de Poisson
    muestra{i}.corte = muestra{i}.young/(2*(1+muestra{i}.poisson));                                 % Módulo de corte
    muestra{i}.proplim = muestra{i}.stress(muestra{i}.proplim(1));                                  % Limite proporcional
end

% Desplegar constates elasticas

Constantes = {'Modulo de Young';'Modulo de Poisson';'Modulo de corte';'proplime proporcional'};
Muestra1_0grad = [muestra{1}.young ; muestra{1}.poisson ; muestra{1}.corte ; muestra{1}.proplim];
Muestra2_0grad = [muestra{2}.young ; muestra{2}.poisson ; muestra{2}.corte ; muestra{2}.proplim];
Muestra1_45grad = [muestra{3}.young ; muestra{3}.poisson ; muestra{3}.corte ; muestra{3}.proplim];
Muestra2_45grad = [muestra{4}.young ; muestra{4}.poisson ; muestra{4}.corte ; muestra{4}.proplim];
Muestra1_90grad = [muestra{5}.young ; muestra{5}.poisson ; muestra{5}.corte ; muestra{5}.proplim];
Muestra2_90grad = [muestra{6}.young ; muestra{6}.poisson ; muestra{6}.corte ; muestra{6}.proplim];
T = table(Constantes,Muestra1_0grad,Muestra1_0grad,Muestra1_45grad,Muestra2_45grad,Muestra1_90grad,Muestra2_90grad);
disp(T)

%% ========= Graficos Esfuerzo vs Ctes elasticas ==============================

yieldlim = [400,400,400,400,400,400]; % Limite de fluencia
for i = 1:6
    muestra{i}.yieldlim = yieldlim(i);
end
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
    grid on
end

% Esfuerzo vs Modulo de Poisson
figure('Name','Esfuerzo vs Modulo de Poisson','Position',[0 0 1366 768])
grad = [0,0,45,45,90,90];
m = [1,2,1,2,1,2];
j = 1;
for i = [1,3,5,2,4,6]
    subplot(2,3,j)
    j = j+1;
    plot(muestra{i}.stress(1:length(muestra{i}.poissonsec)),muestra{i}.poissonsec,'o','MarkerSize',2);
    title({'Esfuerzo vs Modulo de Posson';strcat('DP-1000 a $',num2str(grad(i)),'^{\circ}$ Muestra $',num2str(m(i)),'$')},...
        'Interpreter','latex','FontSize',14),
    xlabel('$\sigma[MPa] $','Interpreter','latex','FontSize',14),ylabel('$\nu[MPa]$','Interpreter','latex','FontSize',14),
    grid on
end