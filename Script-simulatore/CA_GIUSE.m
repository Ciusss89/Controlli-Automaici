%% CONTROLLI AUTOMATICI,  revisione 5, 7/6/2013. 
%  AUTORI: Andrea Rizzo - Giuseppe Tipaldi
%  IMPOSTAZIONE GRAFICA E ULTIME MODIFICHE A OPERA DI GIUSEPPE.
%  PRENDI VISIONE DEL FILE ngridcustom.m PER RISOLVERE EVENTUALI ANOMALIE
%  DI NICHOLS.

    bdclose all;slCharacterEncoding('Windows-1252'); 
    clc; clear all; close all; s=tf('s'); scrsz = get(0,'ScreenSize');
    fprintf('r5.1 Giuseppe Tipaldi %s \r',date);
    fprintf('Setup_ \n');
            
%% DATI PRINCIPALI DEL PROBLEMA:
    % Blocchi costituienti il tipo di controllore.
    Gp =  0;
    fprintf(' Gp(s), espressa in forma ZERO-POLI-GUADAGNO;' ); zpk(Gp)
    fprintf('******************************************************\n*');
    p  = 0;  % n.ro di poli nell'origine dell'impianto (Funzione Gp). 
    Kd = 0; 
    Gs = 0e-0; 
    Ga = 0.e0; 
    Gr = 1e0; 
    Gf = 1/(Kd*Gs);
    fprintf(' Poli nell''origine %d\n* Guadagno Kd %.3g\n* Funzione Sensore %.3g\n* Funzione Attenuatore %.3g\n* Funzione Gr %.3g\n* Funzione Gf %.3g\n\n*',p,Kd,Gs,Ga,Gr,Gf);
    
    % Specificare la costante di guadagno del blocco Gda. Se diversa da una
    % costante, porre a 0 Gda, e introdurre la funzione Gda_s.
    % Nel caso in cui non risulti specificata porre Gda=1.
    Gda = 1; 
    Gda_s = 0*s;
    [num_Gda,den_Gda]=tfdata(Gda_s,'v');

    % Specificare la costante di guadagno del blocco Gdp. Se diversa da una
    % costante, porre a 0 Gdp, e introdurre la funzione Gdp_s.
    % Nel caso in cui non risulti specificata porre Gdp=1.
    Gdp = 1;
    Gdp_s = s*0;
    [num_Gdp,den_Gdp]=tfdata(Gdp_s,'v');
    
    % Specificare i paramentri della risposta al gradino.
    over = 0.0;   % Overshoot.
    tr = 0;    % Rise Time.
    ts = 0;      % Settling Time.
    aph = 0.05;    % Volore % di settling.
    % Calcolo lo Smorzamento
    eps = (abs(log(over)))/(sqrt((pi^2)+(log(over))^2));
    tr_smo = 0;   % Rise Time corrispondente allo smorzamento
    ts_smo = 0;   % Settling Time corrispondente allo smorzamento 
    fprintf(' Overshoot %d%%\n* Rise Time %.3g Sec\n* Settling Time %.3g Sec\n* Assestamento %% %.3g\n* Smorzamento %.3g\n* Lettura Rise Time smorzato %.3g\n* Lettura Settling Time smorzato %.3g\n\n*',over*100,tr,ts,aph,eps,tr_smo,ts_smo);
    
    % Specificare gli errori tollerati a regime:
    e_r  = 0e-0;   % Errore quando agisce il riferimento.
    e_da = 0e-0;   % Errore quando agisce il disturbo sull'attuatore.
    e_dp = 0e-0;   % Errore quando agisce il disturbo sull'impianto.
    e_ds = 0e-0;   % Errore quando agisce il disturbo sul sensore.
    fprintf(' Errore sul riferimento %.3g\n* Errore sull''attenuatore %.3g\n* Errore sull''impianto %.3g\n* Errore sul sensore %.3g\n\n*',e_r,e_da,e_dp,e_ds);
        
    % Gradino unitario [ R_u ] || Rampa r(t) [R_0]
    R_0=1/1;   % Rampa r(t)=t, esempio se r(t)=t/4  -> R_O=0.25
    R_u=1;   % Gradino unitario, valore 1.  
    R=R_u;   % Specificare quale dei precedenti riferimenti usare.
    
%% Disturbi che agiscono sul sistema:
% Ricorda che per dsturbi polinomiali Gp=kp/s^p e Gc/s^u. 
% Il ramo diretto, è riferito all'attuatore.
% e_da<=>limit[e_da(t),x->inf]<=>limit[s*Gp(s)*Gda(s)/(1+Gp*Gc*Ga*Gf*Gs)*Da/s^h+1,s->0]
% Il ramo ramo di feedback è riferito al sensore.
% e_dp<=>limit[e_dp(t),x->inf]<=>limit[s*Gdp(s)/(1+Gp*Gc*Ga*Gf*Gs)*Dp/s^h+1,s->0]  

     % Attuatore [ da ]
     Da_c = 0e-2;     % Inserire qui il valore se risulta una COSTANTE.
     Da_r = 0e-2;      % Inserire qui il valore se risulta una RAMPA.          
     a_da = 0;          % Inserire qui l'ampiezza distubo se SINUSOIDALE.
     Pulse_da = 0;      % Inserire qui la w del disturbo se SINUSOIDALE.
     % Se determinati in fase traduzione delle secifiche inserire:
        w_l_letta_a = 0;     % La Wl letta tramite maschera su carta semiLog.
        kc_a = abs(0);       % Se calcolato specificare il valore. 
        u_a = -inf;          % Se determinato specificare il valore.
    
        
     % Impianto [ dp ]
     Dp_c = 0e-2;      % Inserire qui il valore se risulta una COSTANTE.
     Dp_r = 0e-2;      % Inserire qui il valore se risulta una RAMPA.
     a_dp = 0e-2;         % Inserire qui l'ampiezza distubo se SINUSOIDALE.
     Pulse_dp =0e-3;    % Inserire qui la w del disturbo se SINUSOIDALE.
     % Se determinati in fase traduzione delle secifiche inserire:
        w_l_letta_p = 0; % La Wl letta tramite maschera su carta semiLog
        kc_p = abs(0);      % Se calcolato specificare il valore.
        u_p = -inf;         % Se determinato specificare il valore.
    
        
    % Sensore [ ds ]
    Ds_c = 00e-2;      % Inserire qui il valore se risulta una COSTANTE. 
    Ds_r = 00e-2;      % Inserire qui il valore se risulta una RAMPA.
    a_ds = 00e-2;          % Inserire qui l'ampiezza distubo se SINUSOIDALE.
    Pulse_ds = 0;      % Inserire qui la w del disturbo se SINUSOIDALE.
    % Se determinati in fase traduzione delle secifiche inserire:
        w_h_letta_s = 0e0; % La Wl letta tramite maschera su carta semiLog
        kc_s=abs(0);        % Se calcolato specificare il valore.
        u_s=-inf;           % Se determinato specificare il valore.  
  
        
%% [ 1.0 ] Calcolo Kp    
% Il kp è il k stazionario.  Lim [s^p*Gp(s)] per s->0, 
    kp = dcgain((s^p)*Gp); fprintf(' Guadagno a Bassa Frequenza: \n\tKp = %.3f\n*',kp);
    
    
%% [ 2.0 ] DETERMINAZIONE DEI Kc. 
    % RIFERIMENTO:
    % A meno di stravolgimenti della struttura del controllore, questo è il
    % solo dei 4 ipotetici kc che non cambia espressione di calcolo.
    
    h_r= 0 ;   % Indicare il tipo di sistema scelto per la specifica sul riferimento.
   
    kc_r=0;  ignore_kc_r=0;% NO EDIT
    u_r = h_r - p; 
    if u_r < 0
    	kc_r = 0;
        ignore_kc_r=1;
    end
    if ( e_r~=0 && ignore_kc_r==0)     
        if h_r == 0
        	kc_r = abs( (R_0*Kd^2)/(e_r*kp*Ga+Kd) );
        elseif h_r == 1 
            kc_r = abs( (R_0*Kd^2)/(e_r*kp*Ga) );
        else
            kc_r = abs( (R_0*Kd^2)/(e_r*kp*Ga) );
        end
    end
    % I kc relativi all'attuatore, impianto, sensore, vanno calcolati a
    % mano.
    gain_vect=[kc_r,kc_a,kc_p,kc_s]; tmp=[0,0,0,0];
    omega_min=[0,0,0,0]; w_c_s=0;

%%    -    Disturbi sinusoidali
    % Calcolo disturbi sinusoidali.

    % Impianto, disturbo sul ramo diretto.
    if a_da ~=0
        MS_lf=20*log10(e_da/a_da); % Modulo della funzione S a Low_Frequency 
        w_c_a=w_l_letta_a*2;  % Wc calcolata in funzione del disturbo presente
        omega_min(3)=w_c_a;
        fprintf(' Disturbo sul ramo diretto attuatore:\n\tMS_LF = %.3fdB, Wc > %.3f rad/s \n*',MS_lf,w_c_a);    
    end
 
    % Impianto, disturbo sul ramo diretto.
    if a_dp ~=0
        MS_lf=20*log10(e_dp/a_dp); % Modulo della funzione S a Low_Frequency 
        w_c_p=w_l_letta_p*2; % Wc calcolata in funzione del disturbo presente
        omega_min(4)=w_c_p;
        fprintf(' Disturbo sul ramo diretto:\n\tMS_LF = %.3fdB, Wc > %.3f rad/s \n*',MS_lf,w_c_p);
    end
    
    % Sensore, disturbo sul ramo di feedback
    if a_ds ~=0
        MT_hf=20*log10(e_ds*Gs/a_ds); % Modulo della funzione T ad High_Frequency 
    	w_c_s=w_h_letta_s/2;
        fprintf(' Disturbo sul ramo di feedback:\n\tMT_HF = %.3fdB, Wc < %.3f rad/s \n*',MT_hf,w_c_s);
    end
    
%% [ 3 ] Rise time e Settling time:
% In funzione dello smorzamento, determinare dalle curve i valori di Tr e Ts
% I valori vanno letttti dalle curve.
    xtr_c=0;xts_s=0; % NO EDIT
    if tr ~=0
        xtr_c=tr_smo/tr; omega_min(1)=xtr_c;
        fprintf(' Con un tempo di salita:\n\ttr = %.6fs si ha una Wc > %.2f rad/s \n*',tr,xtr_c);
    end
    if ts ~=0
        xts_c=ts_smo/ts; omega_min(2)=xts_c;
        fprintf(' Con un tempo di assestamento:\n\tts = %.6fs si ha una Wc > %.2f rad/s\n\n*',ts,xts_c);    
    end
 %% [ 3.1 ] TIPO DI SISTEMA - Kc - Pulsazione di progetto.
    u_c_vect=[u_r,u_a,u_p,u_s];  u=max(u_c_vect);
 
    if u < 0
        u=0; % Warning u negativi no! Specificare un u manuale.  
    end
    for j=1:4,
        if u_c_vect(j) == u
            tmp(j)=gain_vect(j);
        end
    end
    Wc_min=max(omega_min);
    kc=max(tmp);
    h=u+p;
    fprintf(' Specifica riferimento kc=%.3g, u=%d.\n* Specifica attenuatore kc=%.3g, u=%d.\n* Specifica impianto kc=%.3g, u=%d.\n* Specifica sensore kc=%.3g, u=%d. \n\n*',gain_vect(1),u_c_vect(1),gain_vect(2),u_c_vect(2),gain_vect(3),u_c_vect(3),gain_vect(4),u_c_vect(4));
    if kc == 0
       kc=1; % Kc Libero.
       fprintf(' VINCOLI DETERMINATI: \n\tGuadagno controllore Kc LIBERO. \n\tGrado controllore u=%d. \n\tGrado Tipo sistema h=%d\n',u,h);
    else
       fprintf(' VINCOLI DETERMINATI: \n\tGuadagno controllore Kc>%.3g. \n\tGrado controllore u=%d. \n\tGrado Tipo sistema h=%d\n',kc,u,h); 
    end
    

    if w_c_s == 0
       fprintf(' \tW massima : LIBER0\n'); 
       fprintf(' \tW minima  : %.3g rad/s\n*',Wc_min);
    else 
       if w_c_s < Wc_min
            %fprintf(' \tW massima : %.3g rad/s\n',w_c_s);
            fprintf(' \tW minima  : %.3g rad/s\n*',Wc_min); 
       else
            fprintf(' \tRange  W  : %.3g<W<%.3g  rad/s\n*',Wc_min,w_c_s); 
       end
    end 
   

%% [ 3.2 ] Sovraelongazione, Smorzamento, Picchi risonanza ( Sp, Tp), margini (mG_T, mPh_T, mG_S, mPh_S).
% Lo smorzamento può essere calcolato analitcamente (soluzione proposta)
% oppure dedotto graficamente in funzione della sovraelongazione.
    
    Tp=1/(2*eps*sqrt(1-eps^2));  Tp_dB=20*log10(Tp);     % Picco di Risonanza di T
    Sp=(2*eps*sqrt(2+4*eps^2+2*sqrt(1+8*eps^2)))/(sqrt(1+8*eps^2)+4*eps^2-1); Sp_dB=20*log10(Sp);    % Picco di Risonanza di S
    mG_T=(1/Tp)+1; mGdB_T=20*log10(mG_T);
    mG_S=(1/(Sp-1))+1; mGdB_S=20*log10(mG_S);
    mPh_T=2*asin(1/(2*Tp)); mPhT_g=(mPh_T*180)/pi;
    mPh_S=2*asin(1/(2*Sp)); mPhS_g=(mPh_S*180)/pi;
    fprintf(' Smorzamento - Tp - Sp e relativi margini di Guadagno e Fase\n');
    fprintf('\tTp = %.3f [ %.2fdB ]\n\tSp = %.3f. [ %.2fdB ]\n*',Tp,Tp_dB,Sp,Sp_dB);
    fprintf(' Margine di guadagno di T\n\tmG_T = %.3f [ %.3fdB ]\n* Margine di guadagno di S\n\tmG_S = %.3f [ %.3fdB ]\n* Margine di fase di T\n\tmPh_T = %.0f\n* Margine di fase di S\n\tmPh_S = %.0f\n*',mG_T,mGdB_T,mG_S,mGdB_S,mPhT_g,mPhS_g);
    mg_des = max(mGdB_T,mGdB_S);   % trovo margine di guadagno desiderato
    mPh_des = max( mPhT_g,mPhS_g); % trovo margine di fase desiderato
    mPh_dig = mPh_des + 10;        % guadagno 10° in più
    grad_stable = -180 + mPh_dig;  % limite minimo di stabilità
    fprintf(' I margini desiderati sono i seguenti:\n\tmargine di guadagno = %.3fdB\n\tmargine di fase analogico = %.0f\n\tmargine di fase digitale = %.0f\n*',mg_des,mPh_des,mPh_dig);
    fprintf(' Minima fase da rispettare per il digitale: %.0f\n',grad_stable); 
    fprintf('******************************************************\n*');
    

%% [ 4.0 ] Basi del progetto:
% Kc: Il guadagno del controllore (Kc) deve essere il massimo fra quelli 
% fin qui calcolati, MA SOLO A PARITA' DI u PIU' ALTO.
% u : Il grado del controllore deve essere il massimo di quelli fin qui
% determinati. 

% Wc Kc  
%     NOTA il sitema risulta tanto veloce, quanto maggiori Wc Kc  scelta.
% Wceff 
%     Letta sul grafico di nichols, pulsazione per cui |Lz|=0db

    Rd=1; Ri=1; Rpi=1; Glp=1;  % NO EDIT
    FAA_ON = 0;         % ABILITA IL PROGETTO DEL FILTRO AA SE POSTO A 1
    SAVE_SIMU_DATA = 0; % ABILITA IMPORTAZIONE GRAFICI SE SIMULATO. 
    kc=kc+kc*0.0;   
    wprog=0;
     
    
    
    fprintf('Progetto:');
    fprintf('\n\tkc = %.3f.\n\th=%1.0f.\n\tu=%1.0f.\n\tWc = %.2f. rad/s.\n',kc,h,u,wprog);
  
    % Reti derivative, introduce un anticipo di fase. L'entità della fase 
    % recuperata dipende dal volere di md. NON E' OPPORTUNO UTILIZZARE RETI
    % CON Md ECCESSIVAMENTE ELEVATI, in questa ipotesi conviene aggingere 
    % una seconda rete per raggiungere il recupero di fase richiesto.
    % L'attività del comando cresce al crescere di Md.
    
    % Definendo un md diverso da 1, viene automaticamente determinata la
    % rete che garantisce il miglior recupero di fase. 
    md=1; md2=1; 
                     % Valore md, e fase massima recuperata :
    % Md: 2 Fase: 19.5 | Md: 4 Fase: 36.9 | Md: 6 Fase: 45.6 | Md: 8 Fase: 51.1 
    % Md: 10 Fase: 54.9 | Md: 12 Fase: 57.8 | Md: 14 Fase: 60.1 | Md: 16 Fase: 61.9
    
    % Reti integrative. 
    % Definendo un mi diverso da 1, viene abilitata la rete. Per non
    % rallentare il sistema è opportuno scengliere una normalizzazione non
    % troppo alta.
    mi = 1; norm_ri=0; 
    mi2= 1; norm_ri2=0; 
       
    % Rete Pi :
    % Questo tipo di rete fornisce un controllore più economico. 
    % Uso: posso inserire uno zero, per ogni polo nell'origine del
    % controllore, solo se u>0.
    % Inserndo il valore di nomalizzazione diverso da 0, verrà abilitata la rete.
    norm_Pi=0;
       
%% [ 5.0 ] Elaborazione & Progetto del Sistema Analogico 

%                    ------ WARING ------- 
    % Quando cambiare segno al Kc:
    % Nichols è girato al contratio o non ruota corettamente attorno al punto critico
    % La fase di L risulta posiva 
    Kc_v = kc;
    %Kc_v=-kc; 
    
    % Criterio di Nyquest :
    %     Ipotesi per avere N ben definito : il digramma di nyquist non
    %     taglia in -1,0.
    %     *N   Numero di rotazioni compiute in SENSO ORARIO ( se kc positivo )
    %          dalla funzione L intorno al punto critico. Per kc negativi le
    %          rotazioni si contano come negative.
    %     *PiL Numero di poli instabili ( parte reale positiva) della
    %          funzione L
    %     Riusciamo a stabilizzare un sistema mediante l'uso di reti
    %     dinamiche correttive se :
    %                       N=-PiL
%% [ 5.1 ] Rete correttive 
    
    % Rete anticipatrice a rucupero fase: 
    if  md ~=1
        zd=wprog/(chop(sqrt(md),2,1));
        % norm_rd=1.2; zd = wprog/norm_rd; Decommenta per progetto manuale
        Rd=(1+s/zd)/(1+s/(zd*md)); fase_Rd=(asin((md-1)/(md+1)))*((180)/(pi));
        gain_rd=-20*log10(1/md);
        fprintf('Rete derivativa 1:');
        fprintf('\n\tMd: %.3g.\n\tFase recuperata: %.3g.\n\tAumento massimo del modulo per md scelto (worst case) : %.3gdB.\n\tZd: %.5g.\n',md,fase_Rd,gain_rd,zd);
    end
    if  md2 ~=1
        zd2=wprog/(chop(sqrt(md2),2,1)); 
        % norm_rd2=1.2; zd2 = wprog/norm_rd2; Decommenta per progetto manuale
        Rd=Rd*(1+(s/zd2))/(1+(s/(zd2*md2))); fase_Rd2=(asin((md2-1)/(md2+1)))*((180)/(pi));
        gain_rd2=-20*log10(1/md2);
        fprintf('Rete derivativa 2:');
        fprintf('\n\tMd: %.3g.\n\tFase recuperata: %.3g.\n\tAumento massimo del modulo per md scelto (worst case) : %.3gdB.\n\tZd: %.5g.\n',md2,fase_Rd2,gain_rd2,zd2);
    end
    
    % Rete attenuatrice a perdita di guadagno Ri:
    if mi ~=1
        %pii=1/(max((100/wprog),((10*mi)/(wprog))));  Autodeterminazione pii
        pii = wprog/norm_ri; 
        Ri = (1+(s/(mi*pii)))/(1+(s/pii));
        loss_ri=20*log10(1/mi);
        fprintf('Rete integrativa :');
        fprintf('\n\tMi: %.3g.\n\tPi: %.5g.\n\tGuadagno perso per mi scelto (best case) : %.3gdB.\n',mi,pii,loss_ri);
    end 
    if mi2 ~=1
        %pii2=1/(max((100/wprog),((10*mi2)/(wprog)))); Autodeterminazione pii2.
        pii2 = wprog/norm_ri2;
        Ri = Ri*(1+(s/(mi2*pii2)))/(1+(s/pii2));
        loss_ri2=20*log10(1/mi2);
        fprintf('Rete integrativa 2:');
        fprintf('\n\tMi: %.3g.\n\tPi: %.5g.\n\tGuadagno perso per mi scelto (best case) : %.3gdB.\n',mi2,pii2,loss_ri2);
    end
    
    % Rete pi 
    if ( u>0 && norm_Pi ~=0 )        
        z_pi=wprog/norm_Pi; Rpi=(1+(s/z_pi));
        %normi_pi2=0; z2_pi=wprog/normi_pi2; Rpi=Rpi*(1+(s/z2_pi));  
        fprintf('Rete proporzionale integrativa:');
        fprintf('\n\tZero: %.5g.\n',z_pi);
    end
    fprintf('******************************************************\n');

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
    
    % GRAFICI: In verde i tratti relativi al tempo continuo.
    Graph1(Hzp,L,Tp,Sp)
    
%% [ 7.0 ] Sistema Digitale :
% Nella parte seguente bisogna specificare i parametri per il controllore
% digitale. Il tempo di sampling (Ts) conviene lasciarlo a 0.1, per avere una
% perdita di fase minima. 
    [Gm,Pm,Wcg,Wcp]=margin(L); Wceff=Wcp;
    Ts=0.1/(Wceff); [Gcz,Lz]=TDmode(Ts,Glp,Gc,Ga,Gp,Gs,Gf);
    %GRAFICI: In rosso i tratti relativi al tempo discreto.
    [num_Gcz,den_Gcz]=Graph2(Lz,Gcz,Tp,Sp);
         
    
%% [ 8.0 ] Paramentri di simulazione per simulink :
% Generalmente non si ha necessità di modificare questo paragrafo, tranne
% per anticipare o ritardare la durata della simulazione.
% Costante di tempo [ tau ], Passo di campionamento [sample_t], durata [stop_t]
    tau=max(abs(1/(real(pole(Gp))))); sample_t=tau/100; stop_t=chop(50*tau,1); 
    fprintf('Parametri di Simulazione:\n');
    fprintf('\ttau = %.10f\n\tstop_t = %d\n\tsample_t = %.10f\n',tau,stop_t,sample_t);
    
    
%% [ 9.0 ] Filtro antialiasing :
% Frequenza di taglio del filtro posta poco prima della pulsazione effettiva fc.
% Usando un filtro del secondo ordine, sul disturbo si giunge con
% un'attenuazione di 12 db su ottava. Il filtro usato è un filtro BESSEL.
	lp_bssl_d=1; lp_bssl_n=1; AA_order=2; % NO EDIT
    
    Pulse_sempling=(2*pi/Ts);
    fprintf('Verifica del teorema del shannon, Ws>2Wds')
    fprintf('\n\t W_sempling %.3f rad/s.\n\t W_ds %.3f rad/s.\n',Pulse_sempling,Pulse_ds);
    
    if Pulse_sempling > (2*Pulse_ds)
        fprintf('Non si rileva utile un filtro antialising, sistema ben campionato.\n')
    else
        if FAA_ON == 1
            pulse_d=(Pulse_sempling-Pulse_ds);
            if ( pulse_d < 0)
                fprintf('Caso strano, ommetto il progetto.\n')
            else
                Pulse_AAF=ceil((pulse_d*0.8)); [lp_bssl_n,lp_bssl_d]=besself(AA_order,Pulse_AAF);
                Glp = tf(lp_bssl_n,lp_bssl_d);
                close all
                figure('Name',' Comportamento del filtro:'); freqs(lp_bssl_n,lp_bssl_d);
                fprintf('Filtro antialiasing: \n\t Pulsazione disturbo w=%.3f \n\t Pulsazione FILTRO AA �=%.3f \n',pulse_d,Pulse_AAF);
                % L'uso del filtro, impone il ricalcolo della funzione L, e dei
                % grafici. CI SI ASPETTA UNA PERDITA DI FASE.
                L=Glp*Gc*Ga*Gp*Gs*Gf; Hzp = zpk(L); [Gcz,Lz]=TDmode(Ts,Glp,Gc,Ga,Gp,Gs,Gf);
                Graph2(Lz,Gcz,Tp,Sp); [Gm,Pm,Wcg,Wcp]=margin(Lz);
            end
        else
            fprintf('Filtro necessario, calcolo non abilitato. FAA_ON=0.\n')
        end
    end

%% [ 10.0 ] Indcatori stabilità :   

   fprintf('Indicatori di stabilita: \n\t Margine di GUADAGNO Gm=%.3g dB, alla pulsazione Wgm=%.3g rad/s \n\t Margine di FASE Pm=%.3g , alla pulsazione Wpm=%.3g rad/s \n',Gm,Wcg,Pm,Wcp);

   
%% [ 11.0 ] Generazione segnali :  
% Blocco dedicato alla generazione di riferimenti particolari, inserire
% nelle variabili sig_gen_amp e sin_gen_freq, ampiezza e frequenza
% desiderate.
    sig_gen_amp= 0.5;
    sig_gen_freq= 0.25;
    
%close all;

%% [ 12.0 ] Output della simulazione :  
% Importa nello spazione di lavoro di matlab output generato da simulink.
if SAVE_SIMU_DATA == 1
    figure('Name',' Output','Position',[1 scrsz(4) scrsz(3) scrsz(4)]); 
    plot(simout.time, simout.signals.values);
    figure('Name',' Output error','Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    plot(simout1.time, simout1.signals.values);
end