function [f, amplitude] = fft_signal(signal, Fs, show_plot)
% fft_signal - Calcule et affiche la FFT d’un signal
%
% Entrées :
%   signal     : vecteur du signal temporel
%   Fs         : fréquence d’échantillonnage (en Hz)
%   show_plot  : booléen (true/false) pour afficher le spectre
%
% Sorties :
%   f          : vecteur des fréquences (Hz)
%   amplitude  : amplitude du spectre (unilatéral)

    % Vérification et mise en forme
    signal = signal(:);  % s’assurer qu’il est en colonne
    N = length(signal);  % nombre d1échantillons
    t = (0:N-1)/Fs;      % vecteur temps (utile pour affichage)

    % Appliquer la FFT
    Y = fft(signal);

    % Calcul du spectre unilatéral (positif)
    P2 = abs(Y / N);          % spectre double
    P1 = P2(1:floor(N/2)+1);  % spectre simple (réel)
    P1(2:end-1) = 2 * P1(2:end-1);  % correction d'amplitude

    % Fréquences associées
    f = Fs * (0:floor(N/2)) / N;

    % Affichage optionnel
    if nargin < 3 || show_plot
        figure;
        plot(f, P1, 'b', 'LineWidth', 1.5);
        grid on;
        title('Spectre en fréquence du signal');
        xlabel('Fréquence (Hz)');
        ylabel('Amplitude');
    end

    % Sorties
    amplitude = P1;
end
