 close all
t=botezat.X.Data;
u=double(botezat.Y(1,3).Data');
w=double(botezat.Y(1,2).Data')%% viteza 
y=double(botezat.Y(1,1).Data')%% pozitie
subplot(3,1,1)
plot(t, u*200), grid, shg
title('Factorul de umplere u')
xlabel('Timp (s)')
ylabel('u(Factorul de umplere [%])')
subplot(3,1,2)
plot(t,y), grid, shg
title('poziția unghiulară')
xlabel('Timp (s)')
ylabel('y(impulsuri)')
subplot(3,1,3)
plot(t,w), grid, shg
title('viteza unghiulară')
xlabel('Timp (s)')
ylabel('w([rad/sec])')
figure
plot(t, [u*200, y, w]), grid, shg
title('Datele proiectului')
xlabel('Timp (s)')
ylabel('Amplitudine')
i1=750
i2=879
i3=1263
i4=1373

yst=mean(w(i3:i4))
ust=mean(u(i3:i4))
y0=mean(w(i1:i2))
u0=mean(u(i1:i2))
k=(yst)/(ust)
i5=1003
i6=1029
y63=0.63*(yst-y0)+y0
hold on
plot(t,y63*ones(1,length(t)),'black','LineStyle','--')
T=t(i6)-t(i5)
A=[-1/T];
B=[k/T];
C=[1];
D=0;
i7=963
i8=1002
Tm=t(i8)-t(i7)
N=round(Tm/(t(2)-t(1))) % timp mort de calculat 

uN=[u(1)*ones(N,1);u(1:end-N)]
ysim = lsim(A,B,C,D,uN,t,w(1));%% 
Emp=sqrt(sum((w-ysim).^2)/length(w));
Q=norm(w-ysim)/norm(w-mean(w))*100;

[num,den]=ss2tf(A,B,C,D)
H_wu=tf(num,den,'iodelay',Tm)
h=tf(k/C,[T 1],'iodelay',Tm)
%% viteza ->pozitie

iv1 = 4566;
iv2 = 4650;
delta_y = y(iv2)-y(iv1);
delta_t = t(iv2)-t(iv1);
w_med = mean(w(iv1:iv2));
Ky = delta_y/(delta_t*w_med) %k integrator = 4.6261

% H = 4.6261/s;

A_int = 0;
B_int = Ky;
C_int = 1;
D_int = 0;
sys = ss(A_int,B_int,C_int,D_int);
ysim2 = lsim(sys,w,t,y(1));
plot(t,[u.*200,y,ysim2]);

Empn1 = norm(y-ysim2)/norm(y-mean(y));
[num,den]=ss2tf(A_int,B_int,C_int,D_int)
Hintegrator=tf(num,den)




%%
i9=1582;
i10=2955;
i11=4666;
i12=5869;
Te=t(2)-t(1);
data_id=iddata(w(i9:i10),u(i9:i10),Te);
data_vd=iddata(w(i11:i12),u(i11:i12),Te);

data_id_p=iddata(y(i9:i10),w(i9:i10),Te);
data_vd_p=iddata(y(i11:i12),w(i11:i12),Te);

%% ARX-Identificare viteza
m_arx=arx(data_id,[1,1,1])

figure
resid(data_vd,m_arx,10)

figure
compare(data_vd,m_arx)
Hyu = tf(m_arx.B, m_arx.A, Te, 'variable', 'z^-1')
Hyu_c = d2c(Hyu, 'zoh')

%% metoda OE  ptr intrare-viteză

m_oe=oe(data_id,[1,1,1])

figure
resid(data_vd,m_oe)

figure
compare(data_vd,m_oe)
Hyu = tf(m_oe.B, m_oe.F, Te, 'variable', 'z^-1')
Hyu_c = d2c(Hyu, 'zoh')
%% ARX Pozitie
p_arx=arx(data_id_p,[1,1,0])

figure 
resid(data_vd_p,p_arx)

figure
compare(data_vd_p,p_arx)



%% OE  de la viteza la pozitie  
p_oe=oe(data_id_p,[1,1,0])

figure
resid(data_vd_p,p_oe)

figure
compare(data_vd_p,p_oe)
Hyu = tf(p_oe.B, p_oe.F, Te, 'variable', 'z^-1')
Kwy=p_oe.B/0.00084
Hyu_c=tf(Kwy,[1 0])
%Hyu_c = d2c(Hyu, 'zoh')
%% decimarea datelor 

u1=u(i9:9:i10)
y1=w(i9:9:i10)
v1=y(i9:9:i10)

u2=u(i11:9:i12)
y2=w(i11:9:i12)
v2=y(i11:9:i12)

% date indentificare intrare->viteza 
data_id_d=iddata(y1,u1,9*Te)
data_vd_d=iddata(y2,u2,9*Te)

% date decimare viteza->pozitie
data_id_pd=iddata(v1,y1,9*Te)
data_vd_pd=iddata(v2,y2,9*Te)

%% ARX cu decimare ptr intrare-viteza
d_arx=arx(data_id_d,[1,1,1])

figure
resid(data_vd_d,d_arx)

figure
compare(data_vd_d,d_arx)

Hyu = tf(d_arx.B, d_arx.A, 9*Te, 'variable', 'z^-1')
Hyu_c = d2c(Hyu, 'zoh')

%% metoda OE cu decimare intrare-viteza
d_oe=oe(data_id_d,[1,1,1])

figure
resid(data_vd_d,d_oe)

figure
compare(data_vd_d,d_oe)

Huw = tf(d_oe.B, d_oe.F, 9*Te, 'variable', 'z^-1')
Huw_c = d2c(Hyu, 'zoh')

%% arx viteza - pozitie cu decimare
close all
dp_arx=arx(data_id_pd,[1,1,0])
figure
resid(data_vd_pd,dp_arx)
figure
compare(data_vd_pd,dp_arx)
Hwu = tf(dp_arx.B, dp_arx.A, 9*Te, 'variable', 'z^-1')
Ki = dp_arx.B/0.00756; % Perioada de achizitie
Hwu_c = tf(Ki, [1 0])
%Hwu_c = d2c(Hyu, 'zoh')
%% OE cu decimare viteza->pozitie

pd_oe=oe(data_id_pd,[1,1,0])%% ultimu e 0 deoarece nu este tact de intarziere intre viteza si pozitie 

figure
resid(data_vd_pd,pd_oe)


figure
compare(data_vd_pd,pd_oe)

Hwy = tf(pd_oe.B, pd_oe.F, 9*Te, 'variable', 'z^-1')
Kwy=pd_oe.B/0.00756
Hwy_c=tf(Kwy,[1 0])
%Hyu_c = d2c(Hyu, 'zoh')
%%
d_oe=oe(data_id_d,[1,1,1])%% ultimu e 0 deoarece nu este tact de intarziere intre viteza si pozitie 

figure
resid(data_vd_d,d_oe)


figure
compare(data_vd_d,d_oe)

Hyu = tf(d_oe.B, d_oe.F, Te, 'variable', 'z^-1')
Hyu_c = d2c(Hyu, 'zoh')
%% armax pozitie cu decimare

dp_armax=armax(data_id_pd,[1,1,1,0])
figure
resid(data_vd_pd,dp_armax)
figure
compare(data_vd_pd,dp_armax)
Hwu = tf(dp_armax.B, dp_armax.A, 9*Te, 'variable', 'z^-1')
Ki = dp_armax.B/0.00756; % Perioada de achizitie
Hwu_c = tf(Ki, [1 0])