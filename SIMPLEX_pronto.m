%% MÉTODO SIMPLEX
%Amanda Botelho Amaral
%%
clc;clear all; close all
format short;
%% ENTRADA DE DADOS
%Exemplo Apostila
% F_Z=[3  1  3];
% S_A=[2  1  1;
%     1  2  3;
%     2  2  1];
% b=[2 5 6];
%Segundo Exemplo
F_Z=[3  5];
S_A=[1  0;
     0  1;3 2];
b=[4 6 18];

%GErando a matriz Z
MAX=(-1);
MIN=(1);
f_z=F_Z*MAX;

%% Gerando a matriz ampliada
[equa,variaveis]=size(S_A);
folga=eye(equa,equa); %matriz de variáveis de folga
for i=variaveis+1:(variaveis+equa+1)
    f_z(i)=0;
end
mat_sol=[S_A,folga,b';f_z]; % Matriz do sistema a ser resolvido
mat_dados=mat_sol;
[A,B]=size(mat_sol);

%% Método Simplex
cont=0;
d=0;
while cont<500
    cont=cont+1;
    %% Condição de Parada
    d=0;
    mat_sol=mat_dados;
    for i=1:B
        x=mat_sol(A,i);
        if x>=0
            d=d+1;
        end
    end
    % Se não houver mais nenhum coeficiente negativo
    if d==B
        disp('Solução Ótima Encontrada');
        break
    end
    %% PROGRAMA
    %  Identificar a variável que entra
    %Escolher o menor valor negativo de z,
    [p1,var_entra]=find(mat_sol(A,:)==(min(min(mat_sol(A,:)))));
    %CAso tenha dois valores iguais
    if length(p1)>1
        p1=p1(1);
        var_entra=var_entra(1);
    end
    coluna_pivo=(mat_sol(:,var_entra));
    %Identificar a linha que sai( o Pivô)
   bbb=length(mat_dados);
   b=mat_dados(1:equa,bbb);
    caso1=b./ coluna_pivo(1:equa);
    for i=1:length(caso1)
        if caso1(i)==inf
            caso_aux(i)=100;
        else 
            caso_aux(i)=  caso1(i);
        end
    end
    
    [linha_sai,p1_1]=find(caso1==(min(min(caso_aux))));
    pivo =mat_dados(linha_sai,var_entra);
    %Coluna pivo: var_entra
    %Linha pivo: linha_sai
    %Atualizar a Matriz ampliada
    LP=mat_dados(linha_sai,:)/pivo;
    for i=1:A
        multiplicador=LP*(-mat_dados(i,var_entra));
        if i==linha_sai
            mat_dados(i,:)=LP;
        else
            mat_dados(i,:)=mat_dados(i,:)+multiplicador;
        end
    end
 

end

%% RESULTADOS
%Montando o vetor com as variáveis do sistema
for i=1:variaveis+equa
    dados(i)={[' A_',num2str(i)]};
end
dados(variaveis+equa+1)={ ' b'};
disp('----------------------------------------------------------------------');
disp('Sistema Final');
disp(dados);
disp(mat_dados)
var1=0;
var2=0;
for i=1:variaveis+equa
    if sum(mat_dados(:,i))==1
        var1=var1+1;
        %sol(var1)=i;
        %Montando o vetor com os valores ótimos da solução
        %Variáveis básicas
        resp_mat(:,var1)=mat_dados(1:A-1,i);
        vetor_sol(:,var1)=resp_mat(:,var1)'*mat_dados(1:A-1,B);
        resp_var(:,var1)=dados(i);
        x(i)=vetor_sol(:,var1);
    else
        %Variáveis não básicas
        var2=var2+1;
        resp_var0(:,var2)=dados(i);
        p_var0(:,var2)=0;
        x(i)=0;
    end

end
disp('----------------------------------------------------------------------');
disp('Solução do ótima do sistema')
disp('');
disp(['x=[', num2str(x),']']);
disp('');
disp('----------------------------------------------------------------------');
disp('Variáveis básicas');
disp(resp_var);
disp(vetor_sol);
disp('Variáveis não básicas');
disp(resp_var0);
disp(p_var0)
z_otimo=mat_dados(A,B)*(-1);
disp(['O valor ótimo de Z : ',num2str(z_otimo)]);
disp('----------------------------------------------------------------------');