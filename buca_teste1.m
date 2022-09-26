%% TRABALHO DE OTIMIZAÇÃO
% Direções Aleatórias
clc; clear all; close all;
%% INICIALIZAÇÃO DOS DADOS

num_var=2;
x0 = [0 0]';
x1 = x0(1);
x2 = x0(2);
sig=0.0001;
fx=(2*(x1)^2+x2^2+2*x1*x2+x1-2*x2+3);
fx_min = 2.9;

% Variáveis de interesse
criterio_parada = sig*(fx - fx_min);
delta_fx = 1;
aux = eye(length(x0));
delta = 0.0001;
criterio_par = 1;
k = 1;
var = zeros;
fx_novo=0;
%% ALGORITMO
while delta_fx >criterio_parada
    dk = randn(num_var,1);
    % Sistema
    a=0;
    b=1;
    % A
    alfa_a = b - 0.618*(b - a);
    aux_fa = x0 + alfa_a*dk;
    fa = (2*(aux_fa(1)^2)) + (aux_fa(2)^2) + (2*aux_fa(1)*aux_fa(2)) + aux_fa(1) - (2*aux_fa(2)) + 3;
    % B
    alfa_b = a + 0.618*(b - a);
    aux_fb = x0 + alfa_b*dk;
    fb = (2*(aux_fb(1)^2)) + (aux_fb(2)^2) + (2*aux_fb(1)*aux_fb(2)) + aux_fb(1) - (2*aux_fb(2)) + 3;

    epsilon = 0.001;
    cont = 1;
    %Análise
    while b-a > epsilon
        if fa < fb
            b = alfa_b;
            alfa_b = alfa_a;
            alfa_a = b - 0.618*(b - a);
            fb = fa;
            aux_fa = x0 + alfa_a*dk;
            fa = (2*(aux_fa(1)^2)) + (aux_fa(2)^2) + (2*aux_fa(1)*aux_fa(2)) + aux_fa(1) - (2*aux_fa(2)) + 3;
        else
            a = alfa_a;
            alfa_a = alfa_b;
            alfa_b = a + 0.618*(b - a);
            fa = fb;
            aux_fb = x0 + alfa_b*dk;
            fb = (2*(aux_fb(1)^2)) + (aux_fb(2)^2) + (2*aux_fb(1)*aux_fb(2)) + aux_fb(1) - (2*aux_fb(2)) + 3;
        end
        cont = cont+ 1;
    end
    alfa = (a+b)/2;
    aux_x = x0;
    x0 = x0 + alfa*dk;
    %% Estabilização da Função Objetivo
    fx_novo(k) = (2*(x0(1)^2)) + (x0(2)^2) + (2*x0(1)*x0(2)) + x0(1) -(2*x0(2)) + 3;
    if length(fx_novo) > 5
        fx = fx_novo;
        var = [fx_novo(k) fx_novo(k-1) fx_novo(k-2) fx_novo(k-3) fx_novo(k-4) fx_novo(k-5)];
        delta_fx = abs(max(var) - min(var));
    end
    k = k+1;
end
fx=(2*(x0(1)^2))+ (x0(2)^2) + (2*x0(1)*x0(2)) + x0(1) - (2*x0(2)) + 3;
%% erro
x1_calc=-1.5;
x2_calc=2.5;
erro_x1=((x1_calc-x0(1))/x1_calc)*100;
erro_x2=((x2_calc-x0(2))/x2_calc)*100;
%% RESULTADOS

disp(' **************************************************************************** ');
disp('                          Método de Direções Aleatórias  ');
disp(' --------------------------------------------------------------------------- ');
disp(' '); 
disp(['O mínimo da função: ' num2str(fx)]);
disp('Os valores de x1 e x2 que minimizam a função objetivo são:');
disp(['x1 = ' num2str(x0(1))]);
disp(['x2 = ' num2str(x0(2))]);
disp(' '); 
disp([' Foram necessárias ',  num2str(k) ,    'interações']);
disp(' ');
disp(' --------------------------------------------------------------------------- ');
disp(' **************************************************************************** ');


