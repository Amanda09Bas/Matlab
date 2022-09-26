function [maximos_x,Fk_2max]=fmut(Fk2,Fk_xf,seleciona_pais_2_aux,RO_1)

    Fk2_aux = Fk2;
    Fk_xf_aux = Fk_xf;
    c = length(seleciona_pais_2_aux)*RO_1;
    % Verificação das piores representantes
    for i=1:ceil(c)
        Fk_2max(i) = find(Fk2_aux==max(Fk2_aux));
        maximos_fk(i) = max(Fk2_aux);
        if i==1
            x_2max(:,i) = Fk_xf_aux(:,Fk_2max(i));
        else
            x_2max(:,i) = Fk_xf_aux(:,Fk_2max(i)+1);
        end
        maximos_x(:,i) = x_2max(:,i);
        Fk2_aux(:,Fk_2max(i))= [];
    end

end