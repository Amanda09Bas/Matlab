%% Simulated Annealing
clear all; clc; close all;
tic()
%% INICIALIZAÇÃO
RO=0.3;
x_o=[0;0];
NS = 100; %número de simulações
alfa = 0.99; %[0.8 a 0.99]
%% FUNÇÃO
%Função Unimodal(f=-0,25; x1=-1,5 e x2=2,5)
%f = @(x1,x2)(2*x1^2 + x2^2 + 2*x1*x2 + x1 - 2*x2 + 3);
% Função Multimodal
%Beale's function (f=0;x1=3 e x2=0,5)
 %f = @(x1,x2)(1.5-x1+x1*x2)^2+(2.25-x1+x1*x2^2)^2+(2.625-x1+x1*x2^3)^2;%(sin(1/2*x1^2-1/4*x2^2+3)*cos(2*x1+1-exp(x2)));
%Lévi function (f=0; x1=1 e x2=1)
 % f=@(x1,x2)(sin(3*pi*x1))^2+(x1-1)^2*(1+(sin(3*pi*x2))^2)+(x2-1)^2*(1+(sin(2*pi*x2))^2);
%Easom function (f=-1, x1=pi e x2=pi)
f=@(x1,x2)(-cos(x1)*cos(x2)*exp(-((x1-pi)^2+(x2-pi)^2)));
for m = 1:NS
    x = x_o;
    N = 3;
    for i = 1:N
        gama = randn(length(x),1);
        gama_arq(:,i) = gama;
        x_o = x_o + (RO*gama);
        xo_arq(:,i) = x_o;
        x1 = x_o(1); x2 = x_o(2);
        Fo(i) = f(x1,x2);
    end
    
    j = find(Fo==(min(Fo)));
    x_min = xo_arq(:,j);
    f_x_min = f(x_min(1),x_min(2));
    y = x_min;
    
    for i = 1:(length(Fo)-1)
        Fo_abs(i) = abs(Fo(i)-Fo(i+1));
    end
    H = max(Fo_abs);
    q = 1;
    p = 1;
    T = 1;
    while T > 0.0001
        T = 3*(alfa^q)*H;
        x = y + (alfa^q)*RO*randn(length(x),1);
        f_x = f(x(1),x(2));
        f_y = f(y(1),y(2));
        if f_x <= f_y
            y = x;
            if f_x < f_x_min
                x_min = x;
            end
        else
            lambda = exp((f_y-f_x)/T);
            fi = rand(1,1);
            if fi <= lambda
                y = x;
            end
        end
        p = p+1;
        if p > 100
            p = 1;
            q = q+1;
        end
        T = alfa*T;
    end
    x_aux(:,m) = x;
    f_aux(:,m) = f(x(1),x(2));
    %MAT_DADOS=[linha s x1 x2 f(x1,x2)]
    mat_dados(m,:) = [m x(1) x(2) f_aux(m)];
end

%% RESULTADO
toc()
f_med=min(f_aux);
lin_ot = find(f_aux==(min(f_aux)));
res_x1 =mat_dados(lin_ot,2);
res_x2 =mat_dados(lin_ot,3);
x_otimo = [res_x1; res_x2];
disp('----------------------------------------------------------------');
disp('*****************************************************************');
disp(' ');
disp(' --------------- Simulated Annealing --------------------');
disp(['O valor mínimo da função-objetivo é : f(x1,x2)= ',num2str(f_med)]);
disp('');
disp('Onde:');
disp(['x1 = ', num2str(x_otimo(1))]);
disp(['x2 = ', num2str(x_otimo(2))]);
disp(' ');
disp('*****************************************************************');
disp('----------------------------------------------------------------');
%% Gráfico
% % Domínio da imagem
% x_1=linspace(-4.5,4.5,NS);
% x_2=linspace(-4.5,4.5,NS);
% z=zeros(length(x_2),length(x_2));
% for i=1:length(x_1)
%     for j=1:length(x_2)
%          aux_z=f(x_1(i), x_2(j));
%         z(i,j)=aux_z;
%     end
% end
% surfc(x_1,x_2,z);
