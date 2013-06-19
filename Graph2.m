function [ ris_num_Gcz,ris_den_Gcz ] = Graph2( Lz,Gcz,Tp,Sp )
%GRAPH2 Di seguito verranno stampanti a video tutti i grafici che possono
% ritornare utili. I grafici utili al progetto a tempo discreto usano il
% tratto ROSSO.
    % Bode Chart Sistema Digitale
    figure('Name',' Diagramma di Bode di Lz per LoopShaping ');
    bode(Lz,'r')
    grid on;
    
    % Nichols Chart Sistema Digitale  
    figure('Name',' Diagramma di Nichols del Sistema Digitale ');
    [ris_num_Gcz,ris_den_Gcz]= tfdata(Gcz,'v');
    [num_Lz,den_Lz]= tfdata(Lz,'v');
    ngridcustom(Tp,Sp) % disegna le curve a modulo costante  
    hold on
    nichols(num_Lz,den_Lz,'r') % disegna la funzione ad anello in nichols
    grid on
    set(gca,'ylim',[-50,50])
end

