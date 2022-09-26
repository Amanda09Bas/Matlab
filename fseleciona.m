function [seleciona_pais_sel]=fseleciona(Fk_xf,m)

    Fk_xf_12 = Fk_xf(1:(m-1),:);
    Fk_xf_3 = Fk_xf(m,:);
    for i = 1:length(Fk_xf)
        Fk_5min_sel(i) = find(Fk_xf_3==min(Fk_xf_3));
        minimos_sel(i) = min(Fk_xf_3);
        x_5min_sel(:,i) = Fk_xf_12(:,Fk_5min_sel(i));
        minimos_x_sel(:,i) = x_5min_sel(:,i);
        Fk_xf_3(:,Fk_5min_sel(i)) = [];
        Fk_xf_12(:, Fk_5min_sel(i)) = [];
    end
    seleciona_pais_sel = [minimos_x_sel; minimos_sel];
end