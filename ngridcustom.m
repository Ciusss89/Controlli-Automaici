function[]=ngridcustom(Tp0,Sp0)
status=ishold;
limit=axis;
if nargin==1 | nargin==2
    hold off
else
    axis(axis)
    hold on
    if limit(3)>=0;
        disp('Questo non � un Nichols Chart!!!');
        hold off
        return
    end;
end

px=linspace(-359.99,-0.01,500);
i=sqrt(-1);
[p,m]=meshgrid(px,Tp0);
[p,s]=meshgrid(px,Sp0);
z=m.*exp(i*p/180*pi);
z1=s.*exp(i*p/180*pi);
g=z./(1-z);
g1=(1-z1)./z1;
gain=20*log10(abs(g));
gain1=20*log10(abs(g1));
 

%% Warning.
% Se si incappa in problemi grafici come quelli dati dal tema di settembre
% 2012, abilitare la riga 31,32 disabilitando la 34-35
% phase=rem(angle(g)/pi*180+360,360);
% phase1=rem(angle(g1)/pi*180+360,360);   

 phase=rem(angle(g)/pi*180+360,360)-360;
 phase1=rem(angle(g1)/pi*180+360,360)-360; 

plot(phase,gain,'-b')
hold on
grid on
plot(phase1,gain1,'-b')

if nargin==1 | nargin==2
    set(gca,'xlim',[-360,0]);
    set(gca,'ylim',[-50,50]);
end
plot([-180-180],[-100 100],'b--')
