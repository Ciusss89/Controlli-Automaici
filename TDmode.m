function [ risGcz,risLz] = TDmode( Ts,Glp,Gc,Ga,Gp,Gs,Gf)
%TDMODE Discretizzazione controllore.
    s=tf('s');
    risGcz=c2d(Gc,Ts,'matched'); % discretizzazione Gc
    Gczoh=1/(1+s*Ts/2);
    risLz=Gczoh*Glp*Gc*Ga*Gp*Gs*Gf;

end

