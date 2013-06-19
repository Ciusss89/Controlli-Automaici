function [ ] = Graph1( Hzp,L,Tp,Sp )
%GRAPH1 Di seguito verranno stampanti a video tutti i grafici che possono
% ritornare utili. I grafici utili al progetto a tempo continuo usano il
% tratto VERDE.

% Bode Chart di L  
    figure,bode(Hzp,'g')     
    title('Diagramma di Bode di L')
    grid on
    
    % Asintotic Bode Chart di L
    [num,den]=tfdata(Hzp,'v');
    %figure,asbode(num,den)
    title('Diagramma di Bode Asintotico')
    
    % Polar Chart
    [re,im]=nyquist(Hzp);   
    figure,plot(squeeze(re), squeeze(im),'g')
    title('Diagramma Polare di L')
    grid on
    
    % Nyquist Chart nell'intorno di zero
    % WARNING : Se dovesse risultare girato al contrario, e i diagrammi di
    % nichols sono sfasati di molto, CAMBIARE IL SEGNO a kc_v. (vedi 5.0)
    [numH,denH]=tfdata(Hzp,'v');
    figure,nyquist1(numH,denH)
    title('Diagramma di Nyquist attorno allo 0')
    grid on
    
    % Nyquist Chart Moderno
    figure,nyquist(Hzp,'g'); 
    title('Diagramma di Nyquist moderno')
    grid on
    
    % Bode Chart di T, S, L
    T=L/(1+L);
    S=1/(1+L);
    
    figure,bode(T,'b')
    title('Diagramma di Bode di L - T - S')
    hold on
    bode(S,'r')
    hold on
    bode(L,'k')
    grid on
    
    % Nichols Chart Sistema Analogico
    figure,ngridcustom(Tp,Sp) % si traccia il diagramma di Tp0 e Sp0
    hold on 
    nichols(Hzp,'g') %si traccia il diagramma di Nichols del sistema analogico
    title('Diagramma di Nichols del Sistema Analogico')
    grid on
    set(gca,'ylim',[-50,50])
    %set(gca,'xlim',[-360,180]) % allarga l'asse x
    hold on
end

