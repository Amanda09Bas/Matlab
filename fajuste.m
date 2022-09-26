function [FT,m,n ]=fajuste(Fk_xf,y)

     [m,n] = size(Fk_xf);
    aux_J = Fk_xf(m,:);
    Jmedia = (sum(Fk_xf(m,:)))/n;
    Jmax = max(aux_J);
    Jmin = min(aux_J);

    v = ((y*Jmin)-Jmax)/(y-1);
    if Jmin >= v
        alfa_J = Jmedia*((y-1)/(Jmax-Jmedia));
        beta_J = Jmedia*((Jmax-(y*Jmedia))/(Jmax-Jmedia));
    else
        alfa_J = Jmedia/(Jmax-Jmin);
        beta_J = -(Jmedia*Jmin)/(Jmax-Jmin);
    end
    FT = (alfa_J*aux_J)+beta_J;
end