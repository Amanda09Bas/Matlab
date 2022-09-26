function [seleciona_pais_2_aux]=fcruz(melhor_x,num_pais,seleciona_pais,RO_2)

   for i = 1:length(melhor_x)
        if i < (num_pais/2)+1
            pop_1(:,i) = melhor_x(:,i);
        else
            for j = 1:(num_pais/2)
                pop_2(:,j) = melhor_x(:,(num_pais/2)+j);
            end
        end
    end
    for i = 1:(length(seleciona_pais)/2)
        xg(:,i) = RO_2*pop_1(:,i) + (1-RO_2)*pop_2(:,i);
    end
    seleciona_pais_2_aux = [melhor_x xg];
end