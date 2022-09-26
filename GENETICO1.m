%% Algoritmo Genético - SSGA
clear all; clc; close all;
tic()
%% INICIALIZAÇÃO
x_o = [-(pi-1) (pi+1); -(pi-1) (pi+1)];
N = 1000;
RO_1 = 0.02;
k=0;
RO_2 = 0.99;% RO_2 -> [0.8 a 0.99]
y = 1.8;
%% FUNÇÃO
%Função Unimodal(f=-0,25; x_o=-1,5 e x2=2,5)
%f = @(x_o,x2)(2*x_o^2 + x2^2 + 2*x_o*x2 + x_o - 2*x2 + 3);
% Função Multimodal
%Beale's function (f=0;x_o=3 e x2=0,5)
 %f = @(x_o,x2)(1.5-x_o+x_o*x2)^2+(2.25-x_o+x_o*x2^2)^2+(2.625-x_o+x_o*x2^3)^2;%(sin(1/2*x_o^2-1/4*x2^2+3)*cos(2*x_o+1-exp(x2)));
%Lévi function (f=0; x_o=1 e x2=1)
 % f=@(x_o,x2)(sin(3*pi*x_o))^2+(x_o-1)^2*(1+(sin(3*pi*x2))^2)+(x2-1)^2*(1+(sin(2*pi*x2))^2);
%Easom function (f=-1, x_o=pi e x2=pi)
f=@(x_o,x2)(-cos(x_o)*cos(x2)*exp(-((x_o-pi)^2+(x2-pi)^2)));
%% Algoritmo
while k < 150
    k=k+1;
    %% Passo1: Criação da População
   %principal( modo aleatório)
    if k==1
        for n = 1:N
            for j = 1:length(x_o)
                for i = 1:length(x_o)
                    pop = x_o(i,j) - x_o(i,j)*(rand(1)*(x_o(i,2)-x_o(i,1))+x_o(i,1));
                end
                populacao(j,n) = pop;
            end
        end
        % Avalia Função-Objetivo
        for i = 1:N
            Fk(i) = f(populacao(1,i),populacao(2,i));
        end
    else
        % Cria a População
        populacao2 = zeros(m-1,length(pais_elite));
        for n = num_elite:length(pais_elite)
            for j = 1:length(x_o)
                for i = 1:length(x_o)
                    pop2 = x_o(i,j) - x_o(i,j)*(rand(1)*(x_o(i,2)-x_o(i,1))+x_o(i,1));
                end
                populacao2(j,n) = pop2;
            end
        end

        % Função-Objetivo
        for i = num_elite:length(pais_elite)
            Fk_2(i) = f(populacao2(1,i),populacao2(2,i));
        end
       polulacao_2 = [populacao2; Fk_2];
        polulacao_1 = pais_elite +polulacao_2;
    end
    if k==1
        polulacao_1 = [populacao; Fk];
        Fk_aux = Fk;
        populacao_geral = populacao;
        num_pais = 50;
    else
        
        Fk_aux = polulacao_1(m,:);
        populacao_geral = polulacao_1(1:m-1,:);
    end
     
    %% Passo 2: Seleciona Pais

    for i = 1:num_pais
        menor_fk(i) = find(Fk_aux==min(Fk_aux));
        melhor_fk(i) = min(Fk_aux);
        menor_x(:,i) = populacao_geral(:,menor_fk(i));
        melhor_x(:,i) = menor_x(:,i);
        Fk_aux(:,menor_fk(i)) = [];
        populacao_geral(:, menor_fk(i)) = [];
    end
    seleciona_pais = [melhor_x; melhor_fk];
    % Arquivando melhores pais
    melhor_pais = seleciona_pais(:,1);
%    
    %%  cruzamento,mutação,avaliação,função de ajuste,seleção.
    %% Passo 3: Cruzamento
   [seleciona_pais_2_aux]=fcruz(melhor_x,num_pais,seleciona_pais,RO_2);
    % Avaliando a Função-Objetivo
    for i = 1:length(seleciona_pais_2_aux)
        Fk2(i) = f(seleciona_pais_2_aux(1,i),seleciona_pais_2_aux(2,i));
    end
    Fk_xf = [seleciona_pais_2_aux; Fk2];
    pos_melhor_Fk_xf = find(Fk2==min(Fk2));
    melhor_pais_cruzamento = Fk_xf(:,pos_melhor_Fk_xf);
    % Comparando pais 
    if melhor_pais_cruzamento(length(melhor_pais)) < melhor_pais(length(melhor_pais))
        melhor_pais = melhor_pais_cruzamento;
    end
    %% Passo 4: Mutação
    [maximos_x,Fk_2max]=fmut(Fk2,Fk_xf,seleciona_pais_2_aux,RO_1);
    %% Passo 5: Avaliação
    % Mutação dos Piores Representantes
    aux_maximos_x = maximos_x;
    aux_maximos_x(3,:) = [];
    xr = x_o(:,2) - x_o(:,1);
    deltai = 0.05*randn*xr;
    aux_mut1 = aux_maximos_x(:,1)+deltai;
    aux_mut2 = aux_maximos_x(:,2)+deltai;
    f_mut1 = f(aux_mut1(1),aux_mut1(2));
    f_mut2 = f(aux_mut2(1),aux_mut2(2));
    mut1 = [aux_mut1; f_mut1];
    mut2 = [aux_mut2; f_mut2];
    mut_total = [mut1 mut2];
    % Rearranjando a população, piores são subtituidos pelas mutações
    for i=1:length(Fk_2max)
        if i==1
            Fk_xf(:,Fk_2max(i)) = mut_total(:,i);
        else
            Fk_xf(:,Fk_2max(i)+1) = mut_total(:,i);
        end
    end
    for i=1:(length(mut_total)-1)
        if mut_total(length(mut_total),i) < melhor_pais(length(melhor_pais),1)
            melhor_pais = mut_total(:,i);
        end
    end
    %% Passo 6: Função de Ajuste (slide1)
    [FT,m,n ]=fajuste(Fk_xf,y);
    %% Passo 7: Seleção
    %Busca dos melhores
[seleciona_pais_sel]=fseleciona(Fk_xf,m);
    %% Passo 8: Elitização
    %REarranjando a população
    num_elite = 26;
    aux_elite=zeros(m,length(seleciona_pais_sel));
    for i=num_elite:length(seleciona_pais_sel)
        seleciona_pais_sel(:,i) = 0;
    end
    pais_elite = seleciona_pais_sel;
    melhores_populacoes(:,k) = pais_elite(:,1);
    
end
toc()
%% RESULTADOS
f_melhores = melhores_populacoes(m,:);%menor pop_best
melhor_P = find(f_melhores==min(f_melhores));
pop_best = melhor_P(randi(length(melhor_P),1));
min_global = melhores_populacoes(1:m-1,pop_best);
funcao_final = melhores_populacoes(m,pop_best);
disp('----------------------------------------------------------------');
disp('*****************************************************************');
disp(' ');
disp(' --------------- Algoritmo Genético --------------------');
disp(['O pop_best mínimo da função-objetivo é : f(x1,x2)= ',num2str(funcao_final)]);
disp('');
disp('Onde:');
disp(['x1 = ', num2str(min_global(1))]);
disp(['x2 = ', num2str(min_global(2))]);
disp(' ');
disp('*****************************************************************');
disp('----------------------------------------------------------------');

