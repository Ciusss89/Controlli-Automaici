%% CONTROLLI AUTOMATICI,  revisione 4, 30/1/2013. 
%  AUTORI: Andrea Rizzo - Giuseppe Tipaldi
%  IMPOSTAZIONE GRAFICA E ULTIME MODIFICHE A OPERA DI GIUSEPPE.
%  PER USARE IL FILTRO ANTI ALIASING BISOGNA RECUPERARE ALTRA FASE (ALTRI 25 ALMENO)
%  OLTRE QUELLA NECESSARIA A STABILIZZARE IL SISTEMA.
%  PRENDI VISIONE DEL FILE ngridcustom.m PER RISOLVERE EVENTUALI ANOMALIE
%  DI NICHOLS. 

    clc; clear all; close all; s=tf('s');
    scrsz = get(0,'ScreenSize');
        
%% Funzione di tra sferimento dell'impianto Gp:
% Inserire la funzione di trasferimento dell'impianto in esame.
    Gp= FUNZIONE_IMPIANTO ;
    fprintf('Funzione di trasferimento, espressa in forma ZERO-POLI-GUADAGNO\n\t');
    Gp_zpk=zpk(Gp)
    
%% Specifiche che descrivono il sistema:
% Inserire i parametri che descrivono l'impianto in esame.
    Kd= ; Gs=  ; Ga=0e3; Gr=0; Gf=1/(Kd*Gs);
    
    % Specificare la costante di guadagno del blocco Gda, se diversa da una
    % costante, porre a 0 Gda, e introdurre la funzione Gda_s.
    % Nel caso in cui non risulti specificata porre Gda=1.
    Gda=1; 
    Gda_s=s*0;
    [num_Gda,den_Gda]=tfdata(Gda_s,'v');

    % Specificare la costante di guadagno del blocco Gdp, se diversa da una
    % costante, porre a 0 Gdp, e introdurre la funzione Gdp_s.
    % Nel caso in cui non risulti specificata porre Gdp=1.
    Gdp=0;
    Gdp_s=0;
    [num_Gdp,den_Gdp]=tfdata(Gdp_s,'v');

    p=   0 ;  % n.ro di poli nell'origine dell'impianto. 
    h_r= 0 ; % indicare il tipo di sistema scelto (slide part4-35).
    
%% Specifiche del Gradino: 
% Inserire i parametri delle specifche da soddisfare per la risposta al
% gradino [Overshoot, RiseTime, SettlingTime, assestamento % ]. Se non
% presenti assegnare 0 alla variabile.
    over =  0.0; tr =0.00 ; ts = 0.0; aph=0.00;
    
    %Inserire la specifica sul riferimento che deve usare il comando.
    % WARING, nella specifica del riferimento viene specificato il tipo di
    % riferimeto, riferimenti disponibili:
    % Gradino unitario [ R_u ] || Rampa r(t) [R_0]
    R_0=1;   % Rampa r(t)=t, esempio se r(t)=t/4  -> R_O=0.25
    R_u=1;   % Gradino unitario, valore 1  
    R=R_0; 
    
%% Errori a regime :
% Specificare gli errori tollerati a regime:
    e_r = 0e-4;  % Errore quando agisce il riferimento.
    e_da = 0e-2;   % Errore quando agisce il disturbo sull'attuatore.
    e_dp = 0e-4;   % Errore quando agisce il disturbo sull'impianto.
    e_ds = 0e-4; % Errore quando agisce il disturbo sul sensore.
        
%% [ 1.0 ] Calcolo Kp :
% Il kp è il k stazionario.  Lim [s^p*Gp(s)] per s->0
% p poli nell'orgine di Gp.
    fprintf('Guadagno a Bassa Frequenza\n');
    kp = dcgain((s^p)*Gp); fprintf('\tKp = %.3f\n',kp);
    
%% [ 2.0 ] RIFERIMENTO. 
% In relazione al tipo di riferimento (Gradino - Rampa - Parabola) del
% sistema e alla specifica riguardante l'errore massimo tollerato ( Nullo -
% finito non nullo - inf ), determino il valore di h per il riferimento.
% Un errore sul riferimento NULLO, non aggiunge vincoli al kc. 
% Un errore sul riferimento RIFINITO NON NULLO, richide di determinare il 
% kc che verifica la disequazione.     
    kc_r=0; % NO EDIT
    u_r=-1; % NO EDIT
  
    if e_r~=0     
        if h_r == 0
        	kc_r = abs( (R_0*Kd^2-Kd*e_r)/(e_r*kp*Ga) );
        elseif h_r == 1 
            kc_r = abs( (R_0*Kd^2) / (e_r*kp*Ga) );
        else
            kc_r = abs( (R_0*Kd^2)/(e_r*kp*Ga) );
        end
        u_r = h_r - p; 
        if u_r < 0
            kc_r = 0;
        else
            fprintf('Le specifiche sul riferimento impongono i seguenti limiti:');
            fprintf('\n\t|Kc| >= %.3f \n\tu >= %d',kc_r,u_r);
        end
    end
    % Il risultato qui determinato deve trovare conferma con la risoluzione
    % manuale del limite. Il limite deve essere svolto su carta con i
    % passaggi algebrici ben chiari per imposizione da parte del prof.

   

%% [ 2.1 ] Disturbi polinomiali TIPS:
% Ricorda che Gp=kp/s^p e Gc/s^u. 
% Ramo diretto, E_da
% e_da<=>limit[e_da(t),x->inf]<=>limit[s*Gp(s)/(1+Gp*Gc*Ga*Gf*Gs)*Da/s^h+1,s->0]
% Ramo feedback, E_dp
% e_dp<=>limit[e_dp(t),x->inf]<=>limit[s*Gd(s)/(1+Gp*Gc*Ga*Gf*Gs)*Dp/s^h+1,s->0]

    
%% [ 2.1 ] Disturbo sull'attuatore da : 
% Il disturbo può essere polinomiale o sinusoidale, bisogna specificare
% quale disturbo prendere in considerazione, si specifichi la configurazione 
% del disturbo presente e si annulli settando a 0 quella non presente.

%%    -    Disturbo da polinomiale.
% Specificare l'ampiezza, nel caso sia una costante (Da_c), o una rampa, (Da_r)
% Il disturbo polimoniale, richiede il calcolo del limite, la deteminazaione 
% del kc_da, e il grado u per cui tale condizione risulti verificata.
    Da_c=0e-3; Da_r=0;               
    
    % Kc_x è il valore di kc corrispondente a questa specica, e il u_x è il
    % u per cui questa verifica risulta specifica.
    kc_a=abs(0);
    u_a=-inf;       
         
%%    -    Disturbo da sinusoidale 
% Specificare l'ampiezza (a_da), la sua pulsazione, (pulse_da). ATTENZIONE: 
% La Wl e la conseguente Wc ≥ 2Wl DEVONO essere calcolate a MANO!
    a_da = 0;         % Ampiezza distubo.
    Pulse_da = 0;     % Pulsazione Disturbo.
    w_l_letta = 0;    % La Wl letta tramite maschera su carta semiLog.
    
    if a_da ~=0
        MS_lf=20*log10(e_da/a_da); % Modulo della funzione S a Low_Frequency 
        w_c=w_l_letta*2;  % Wc calcolata in funzione del disturbo presente
        fprintf('Disturbo sul ramo diretto attuatore:\n\tMS_LF = %.3fdB, Wc ≥ %.3f rad/s \n',MS_lf,w_c);    
    end

    
%% [ 2.2 ] Disturbo sull'impianto dp : 
% Il disturbo può essere polinomiale o sinusoidale, bisogna specificare
% quale disturbo prendere in considerazione, si specifichi la configurazione 
% del disturbo presente e si annulli settando a 0 quella non presente.   

%%    -    Disturbo dp polinomiale.
% Specificare l'ampiezza, nel caso sia una costante (Dp_c), o una rampa (Dp_r)
% Il disturbo polimoniale, richiede il calcolo del limite, la deteminazaione 
% del kc_dp, e il grado u per cui tale condizione risulti verificata.
    Dp_c=0e-4; Dp_r=0;
    
    % Kc_x è il valore di kc corrispondente a questa specica, e il u_x è il
    % u per cui questa verifica risulta specifica.
    kc_p=abs(0); 
    u_p=-inf;       

%%    -    Disturbo dp sinusoidale 
% Specificare l'ampiezza (a_dp), la sua pulsazione, (pulse_dp). ATTENZIONE: 
% La Wl e la conseguente Wc ≥ 2Wl DEVONO essere calcolate a MANO!
    a_dp=0e-2;         % Ampiezza distubo
    Pulse_dp =0e-0;      % Pulsazione Disturbo
    w_l_letta=0e-1; % La Wl letta tramite maschera su carta semiLog
    
    if a_dp ~=0
        MS_lf=20*log10(e_dp/a_dp); % Modulo della funzione S a Low_Frequency 
        w_c=w_l_letta*2; % Wc calcolata in funzione del disturbo presente
        fprintf('Disturbo sul ramo diretto:\n\tMS_LF = %.3fdB, Wc ≥ %.3f rad/s \n',MS_lf,w_c);
    end

        
%% [ 2.3 ] Disturbo sul sensore ds  : 
% Il disturbo può essere polinomiale o sinusoidale, bisogna specificare
% quale disturbo prendere in considerazione, si specifichi la configurazione 
% del disturbo presente e si annulli settando a 0 quella non presente.   

%%    -    Disturbo ds polinomiale :
% Specificare l'ampiezza, nel caso sia una costante (Ds_c), o una rampa (Ds_r)    
% Il disturbo polimoniale, richiede il calcolo del limite, la deteminazaione 
% del kc_ds, e il grado u per cui tale condizione risulti verificata.
    Ds_c=0; Ds_r=0;
       
    % Kc_x è il valore di kc corrispondente a questa specica, e il u_x è il
    % u per cui questa verifica risulta specifica.
    kc_s=abs(0); % |Kc| calcolato mediante il limite.
    u_s=-inf;    % Valore che verifica il limite  
    
%%    -    Disturbo ds sinusoidale :
% Specificare l'ampiezza (a_ds), la sua pulsazione, (pulse_ds). ATTENZIONE: 
% La Wh e la conseguente Wc ≤ (Wh/2) DEVONO essere calcolate a MANO!
    a_ds=0e-2;       % Ampiezza disturbo
    Pulse_ds=0;      % Pulsazione Disturbo
    w_h_letta=0e2;      % Wh calcolata in funzione del disturbo presente
    
    if a_ds ~=0
        MT_hf=20*log10(e_ds*Gs/a_ds); % Modulo della funzione T ad High_Frequency 
    	w_c_s=w_h_letta/2;
        fprintf('Disturbo sul ramo di feedback:\n\tMT_HF = %.3f dB, Wc ≤ %.3f rad/s \n',MT_hf,w_c_s);
    end


%% [ 3.0 ] Sovraelongazione, Smorzamento, Picchi risonanza ( Sp, Tp), margini (mG_T, mPh_T, mG_S, mPh_S).
% Lo smorzamento può essere calcolato analitcamente (soluzione proposta)
% oppure dedotto graficamente in funzione della sovraelongazione.
    
    eps = (abs(log(over)))/(sqrt((pi^2)+(log(over))^2)); % Smorzamento;
    Tp=1/(2*eps*sqrt(1-eps^2));  Tp_dB=20*log10(Tp);     % Picco di Risonanza di T
    Sp=(2*eps*sqrt(2+4*eps^2+2*sqrt(1+8*eps^2)))/(sqrt(1+8*eps^2)+4*eps^2-1); Sp_dB=20*log10(Sp);    % Picco di Risonanza di S
    mG_T=(1/Tp)+1; mGdB_T=20*log10(mG_T);
    mG_S=(1/(Sp-1))+1; mGdB_S=20*log10(mG_S);
    mPh_T=2*asin(1/(2*Tp)); mPhT_g=(mPh_T*180)/pi;
    mPh_S=2*asin(1/(2*Sp)); mPhS_g=(mPh_S*180)/pi;
    fprintf('Smorzamento - Tp - Sp e relativi margini di Guadagno e Fase\n');
    fprintf('\tTp = %.3f [ %.2fdB ]\n\tSp = %.3f. [ %.2fdB ]\n',Tp,Tp_dB,Sp,Sp_dB);
    fprintf('Margine di guadagno di T\n\tmG_T = %.3f [ %.3fdB ]\nMargine di guadagno di S\n\tmG_S = %.3f [ %.3fdB ]\nMargine di fase di T\n\tmPh_T = %.0f\nMargine di fase di S\n\tmPh_S = %.0f\n',mG_T,mGdB_T,mG_S,mGdB_S,mPhT_g,mPhS_g);
    mg_des = max(mGdB_T,mGdB_S);   % trovo margine di guadagno desiderato
    mPh_des = max( mPhT_g,mPhS_g); % trovo margine di fase desiderato
    mPh_dig = mPh_des + 10;        % guadagno 10° in più
    grad_stable = -180 + mPh_dig; % limite minimo di stabilità
    fprintf('I margini desiderati sono i seguenti:\n\tmargine di guadagno = %.3fdB\n\tmargine di fase analogico = %.0f\n\tmargine di fase digitale = %.0f\n',mg_des,mPh_des,mPh_dig);
    fprintf('Minima fase da rispettare per il digitale: %.0f\n',grad_stable);
    fprintf('Smorzamento: \n\t%.3f\n',eps);
    
    
%% [ 3.1 ] Rise time e Settling time:
% In funzione dello smorzamento, determinare dalle curve i valori di Tr e Ts
% I valori vanno letttti dalle curve.
    tr_smo = 0;   % Rise Time corrispondente allo smorzamento
    ts_smo = 0;   % Settling Time corrispondente allo smorzamento
    
    xtr_c=tr_smo/tr;
    xts_c=ts_smo/ts;
    fprintf('Con un tempo di salita:\n\ttr = %.3fs si ha una Wc ≥ %.2f rad/s\n',tr,xtr_c);
    fprintf('Con un tempo di assestamento:\n\tts = %.3fs si ha una Wc ≥ %.2f rad/s\n',ts,xts_c);

%% [ 4.0 ] Basi del progetto:
% Kc: Il guadagno del controllore (Kc) deve essere il massimo fra quelli 
% fin qui calcolati, MA SOLO A PARITA' DI u PIU' ALTO.
% u : Il grado del controllore deve essere il massimo di quelli fin qui
% determinati. 

% Wc: Il disturbo sinusoidale sul sensore fornisce la massima Wc mentre, tutte le restanti 
%     pulsazioni calcolate forniscono il limite inferiore. Scegliere
%     aribitrariamente un valore. 
%     NOTA il sitema risulta tanto veloce, quanto maggiore è l'Wc scelta.
    
    % WARNING, trova solo il massimo dei kc, senza capire se è il kc
    % corretto o meno. prestare attenzione, ed eventulmente specificarlo di
    % seguito.
    gain_vect=[kc_r,kc_a,kc_p,kc_s]; kc=(max(gain_vect)+0.4*max(gain_vect));
    if kc==0
        kc=1;
    end
    %kc=; % Kc manuale!

    u_c_vect=[u_r,u_a,u_p,u_s];    u=max(u_c_vect);
	
    if u < 0
        u=0; % u manuale!
    end
    
    h=u+p;
    wprog=0;
    
    fprintf('Dalla traduzione delle specifiche risulta:');
    fprintf('\n\tkc = %.3f.\n\th=%1.0f.\n\tu=%1.0f.\n\tWc = %.2f. rad/s.\n',kc,h,u,wprog);
    
    %Funzione di trasferimento del filtro antialising. 
    Glp=1; %NO EDIT
    
%% [ 5.0 ] Elaborazione & Progetto del Sistema Analogico 
% Si procede alla stabilizzazione e al progetto delle reti. Conviene
% stabilire inizialmente quanta fase si deve recuperare, e valutare in un
% secondo momomento il guadagno in eccesso.
% La stabilità e il rispetto delle specifche si raggiunge inserendo delle
% reti dinamiche e bilanciandone il loro contributo nella funzione del
% controllore.
    
    %                    ------ WARING ------- 
    % Quando cambiare segno al Kc:
    % Nichols è girato al contratio o non ruota corettamente attorno al punto critico
    % La fase di L risulta posiva 
    Kc_v = kc;
    %Kc_v=-kc; 
    
    %Il sistema risulta stabile -N=Pol, questo non sarà mai vero 
    % (perchè se lo fosse il sistema sarebbe già stabile), quindi il
    % sistema risulta stabilizzabile mediante l'aggiunta di reti dinamiche.
    % N rotazioni attorno il punto critico.
    % Pol poli a parte reale positiva della funzione L.
    
%% [ 5.1 ] Rete derivativa (rete anticipatrice a rucupero fase) Rd : 
% In md si deve specificare il livello di normalizzazione della rete usata,

    Rd = 1;
   	md = 1;
    md2= 1;
    norm_rd=0;
    norm_rd2=0;
    zd = wprog/norm_rd;
    zd2 = wprog/norm_rd2;
    %Rd=(1+s/zd)/(1+s/(zd*md));
    %Rd=(1+(s/zd))/(1+(s/(zd*md))) * (1+(s/zd2))/(1+(s/(zd2*md2)))  ;
        
%% [ 5.2 ] Rete integrativa (rete attenuatrice a perdita di guadagno) Ri:
% In mi si deve specificare il livello di normalizzazione della rete usata,
% La rete integrativa, in generale tende a far crescere l'effetto coda,
% quanto più è grande la normalizzazione.
    Ri = 1;
    mi = 1;
    mi2= 1;
    norm_ri=0;
    norm_ri2=0;
    pii = wprog/norm_ri;
    pii2 = wprog/norm_ri2;
    %Ri = (1+(s/(mi*pii)))/(1+(s/pii));
    %Ri = (1+(s/(mi*pii)))/(1+(s/pii)) * (1+(s/(mi2*pii2)))/(1+(s/pii2));

%% [ 5.3 ] Rete Pi :
% Questo tipo di rete fornisce un controllore migliore.
% Uso: posso inserire uno zero, per ogni polo nell'origine del controllore.
% Quindi con un u=1 posso inserire un solo zero, mediante la rete pi, con u=2
% è possibilie inserire 2 reti pi e cosi via. Questo perchè ogni zero
% introdotto dalla rete pi inserita risulta 'compensato' dal polo
% nell'orgine del controllore.
    Rpi = 1;
    if u>0
        normi_pi=0;
        z_pi=wprog/normi_pi;
        Rpi=(1+(s/z_pi));
        
        normi_pi2=0;
        z2_pi=wprog/normi_pi2;
        %Rpi=(1+(s/z_pi))*(1+(s/z2_pi));      
    end

%% [ 6.0 ] Sistema Analogico:   
% Verranno di seguito calcolate le funzioni del blocco controllore, la 
% funzione ad anello, segue la visualizzazione grafica. 
    
    Gc=(Kc_v/(s^u))*Rd*Ri*Rpi; % controllore
    md_n= md *md2 ;   % prodotto dei vari md usati (se più di una)  
    mi_n= mi *mi2;   % prodotto dei vari mi usati (se più di una) 
    if u == 0
        Umax = R*Gr*Ga*Kc_v*(md_n/mi_n);
        fprintf('Comando massimo: \n\tUmax = %.10f\n',Umax);
    else 
        fprintf('Per avere una valutazone analitica dell''azione del comando, u=0.');
    end
    
    L=Gc*Ga*Gp*Gs*Gf; % Funzione ad anello (tempo continuo)
    [num_Gp,den_Gp]=tfdata(Gp,'v');
    [num_Gc,den_Gc]=tfdata(Gc,'v');    
    Hzp = zpk(L);   

%% [ 6.1 ] GRAFICI:  
% In verde i tratti relativi al tempo continuo.
    Graph1(Hzp,L,Tp,Sp)
    
%% [ 7.0 ] Sistema Digitale :
% Nella parte seguente bisogna specificare i parametri per il controllore
% digitale. Il tempo di sampling (Ts) conviene lasciarlo a 0.1, per avere una
% perdita di fase minima. 
% La Wceff dedotta mediante lettura grafica sul grafico di nichols, 
% corrisponde alla Wc per cui si ha un guadagno unitorio (intersezione della Lz con
% l'asse a 0db). NOTA Wceff va letta prima di lanciare la simulazione e
% dopo aver stabilizzato il sistema.
    Wceff=wprog; 
    %Wceff=0; 
    Ts=0.1/(Wceff);    
    [Gcz,Lz]=TDmode(Ts,Glp,Gc,Ga,Gp,Gs,Gf);
    
%% [ 7.1 ] GRAFICI:
% In rosso i tratti relativi al tempo discreto.
    [num_Gcz,den_Gcz]=Graph2(Lz,Gcz,Tp,Sp);
    
%% [ 8.0 ] Paramentri di simulazione per simulink :
% Generalmente non si ha necessità di modificare questo paragrafo, tranne
% per anticipare o ritardare la durata della simulazione.
% Costante di tempo [ tau ], Passo di campionamento [sample_t], durata [stop_t]
    tau=max(abs(1/(real(pole(Gp))))); sample_t=tau/100; stop_t=chop(50*tau,1); 
    fprintf('Parametri di Simulazione:\n');
    fprintf('\ttau = %.10f\n\tstop_t = %d\n\tsample_t = %.10f\n',tau,stop_t,sample_t);
    
     
%% [ 9.0 ] Gc, funzione di trasferimento del controllore, e reti usate:
% Reti usate per la stabilizzazione (prendere nota solo di quelle usate) e
% funzione fdt di Gc
    fprintf('Funzione di trasferimento di L:'); zpk(L)
    fprintf('Funzione di trasferimento del controllore:'); zpk(Gc)
    fprintf('Rete integrativa:'); zpk(Ri)
    fprintf('Rete derivativa:'); zpk(Rd)
    fprintf('Rete pi:'); zpk (Rpi)
    
%% [ 10.0 ] Filtro antialiasing :
% Progetto filtro antialising.La frequenza di taglio del filtro deve essere
% posta un'ottava prima della pulsazione effettiva.
% Usando un filtro del secondo ordine, sul disturbo si giunge con
% un'attenuazione di 12 db su ottava. Il filtro usato è un filtro BESSEL.
	lp_bssl_d=1; lp_bssl_n=1; AA_order=2; % NO EDIT
    
    Pulse_sempling=(2*pi/Ts);
    fprintf('Verifica del teorema del shannon, Ws>2Wds')
    fprintf('\n\t W_sempling %.3f rad/s.\n\t W_ds %.3f rad/s.\n',Pulse_sempling,Pulse_ds);
    
    if Pulse_sempling > (2*Pulse_ds)
        fprintf('Non si rileva utile un filtro antialising, il sistema è ben campionato.\n')
    else
        pulse_d=(Pulse_sempling-Pulse_ds);
        if ( pulse_d < 0)
            fprintf('Caso strano, ommetto il progetto.\n')
        else
            Pulse_AAF=ceil((pulse_d/2));
            [lp_bssl_n,lp_bssl_d]=besself(AA_order,Pulse_AAF);
            Glp = TF(lp_bssl_n,lp_bssl_d);
            figure('Name',' Comportamento del filtro:'); freqs(lp_bssl_n,lp_bssl_d);
            fprintf('Filtro antialiasing: \n\t Pulsazione disturbo w=%.3f \n\t Pulsazione FILTRO AA ω=%.3f \n',pulse_d,Pulse_AAF);
            close all
            % L'uso del filtro, impone il ricalcolo della funzione L, e dei
            % grafici. CI SI ASPETTA UNA PERDITA DI FASE.
            L=Glp*Gc*Ga*Gp*Gs*Gf;
            Hzp = zpk(L);
            [Gcz,Lz]=TDmode(Ts,Glp,Gc,Ga,Gp,Gs,Gf);
            Graph2(Lz,Gcz,Tp,Sp);
        end
    end

   
%% [ 11.0 ] Generazione segnali :  
% Blocco dedicato alla generazione di riferimenti particolari, inserire
% nelle variabili sig_gen_amp e sin_gen_freq, ampiezza e frequenza
% desiderate.
    sig_gen_amp= 1;
    sig_gen_freq= 0.5;
    
%close all;

%% [ 12.0 ] Output della simulazione :  
% Output generato dalla simulazione.
    figure('Name',' Simout','Position',[1 scrsz(4) scrsz(3) scrsz(4)]); 
    plot(simout.time, simout.signals.values);
