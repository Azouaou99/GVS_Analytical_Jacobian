function lambda_dom = Valeur_propre_dominante(vx_norm, lambda)
    % Seuils pour filtrer les pulsations trop faibles ou trop élevées
    puls_min = 1e-2;   % seuil inférieur à ajuster selon le système
    puls_max = 200;     % seuil supérieur à ajuster selon l2échelle typique

    % Fréquence naturelle (pulsation) en rad/s
    omega = lambda;

    % Masque pour ne garder que les valeurs raisonnables
    masque = (omega >= puls_min) & (omega <= puls_max);

    if ~any(masque)
        warning('Aucune valeur propre dans la plage dynamique pertinente. Retour de NaN.');
        lambda_dom = NaN;
        return;
    end

    % Appliquer le filtre
    vx_filt = vx_norm(masque);
    lambda_filt = lambda(masque);
    omega_filt = omega(masque);

    % Critère pondéré : norme^2 / pulsation (faible pulsation = peu dynamique)
    dominance_critere = vx_filt;%(vx_filt.^2) ./ omega_filt;

    % Trouver la dominante
    [~, idx_dom] = max(dominance_critere);
    lambda_dom = lambda_filt(idx_dom);
end

