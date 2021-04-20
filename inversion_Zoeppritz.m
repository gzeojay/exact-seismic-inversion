% c：反演结果
% c1：平滑之后的初始模型
% ct：原始测井曲线
%%
clear all;  close all;  clc;
tic
% return
%% 参数选择
dt=0.002;   %采样间隔
nwt=65;     %子波采样点数
nt0=32;     %子波中心数

ntheta=9;  dtheta=5;  stheta=5;   %入射角： 个数，间隔，初始角度

% 未动参数

% beta=2*1e2;   %% PP-PS modified cauchy

alpha1=1.0;

alpha2=20;       %%  信噪比为2时，alpha1=1.0;alpha2=20; 不含噪声时，alpha1=1.0;alpha2=20;

beta=1e-2;   %% PP-PS cauchy 3e-1 无噪 5630； 有噪 2-0 8272;   PP 6E-2 7665

% zoeppritz正演参数
ipol=0;
%          ipol = 1, the reflection or transmission
%          coefficient "coef" (see below) is returned as output
%          in polar form, i.e., the real part of "coef" is the
%          amplitude (i.e., magnitude) of the complex coefficient
%          and the imaginary part of "coef" is the phase angle
%          of the coefficient in radians. If ipol .ne. 1, then
%          the coefficient is returned in rectangular form, i.e.,
%          the real and imaginary parts of "coef" are the real
%          and imaginary parts of the coefficient.

incwav=1;
%          incwav = 1 for an incident P wave from above.
%          incwav = 2 for an incident S wave from above.
%          incwav = 3 for an incident P wave from below.
%          incwav = 4 for an incident S wave from below.
% irfwav=1;
%          irfwav = 1 for a reflected/transmitted P wave in upper medium.
%          irfwav = 2 for a reflected/transmitted S wave in upper medium.
%          irfwav = 3 for a transmitted/reflected P wave in lower medium.
%          irfwav = 4 for a transmitted/reflected S wave in lower medium.
%%  测井byy2_output
% 输入初始模型

% AAA=load('byy2_output_quwenzi.txt');      %导入测井数据
% [m,n]=size(AAA);
% 
% %25594-22794=2800，选取2850m~3200m深度(350m)的地层,共2800行数据（包含目的层）
% %均匀选取A中的2800行抽取为140行，20行抽取为1行，便于显示
% AA=AAA(22794:25594,:);
% nt=139;                        %纵向采样点
% 
% c=zeros(nt+1,5);
% for i=1:nt+1               %注意多取一行，便于计算最底层界面的反射系数
%     vpt=1000./AA(:,2);     %c中第一列为纵波速度，转换单位μs/m→km/s    Vp+0.7(获得更好看的模型数据)
%     c(i,1)=mean(vpt((i-1)*20+1:i*20));       %20行取平均
%         vps=304.8./AA(:,7);    %c中第七列为横波速度，转换单位μs/ft→km/s，1000000×0.3048=304800   Vs+0.5(获得更好看的模型数据)
%         c(i,2)=mean(vps((i-1)*20+1:i*20));       %20行取平均
%         rhot=AA(:,6);          %c中第六列为密度，g/cm3
%         c(i,3)=mean(rhot((i-1)*20+1:i*20));      %20行取平均
% end
% % c(:,2)=-0.055*c(:,1).^2+1.017*c(:,1)-1.031; % 横波速度
% % c(:,3)=1.74*c(:,1).^0.25; %密度
% %
% for i=1:3      %平滑2次
%     c(:,1)=smooth(c(:,1),5);
%     c(:,2)=smooth(c(:,2),5);
%     c(:,3)=smooth(c(:,3),5);
% end
% %%
% % 1. 经验公式 -----------------------
% Por=1.16-0.41*c(:,3);  % 1.16-0.41*c(:,3)-0.05
% 
% % 2. 
% 
% %%
% c(:,4)=Por;
% 
% % F
% gamma_dry=mean(c(:,1))^2/mean(c(:,2))^2-0.5; % mean(c(:,1))^2/mean(c(:,2))^2-0.3
% M=c(:,1).^2.*c(:,3);
% U=c(:,2).^2.*c(:,3);
% f=M-gamma_dry*U;
% F=f./Por;
% c(:,5)=F;
% % dry vp
% c(:,1)=((M-f)./c(:,3)).^(1/2);

%%  测井 bz212.las

fid_well=fopen('bz212.las','r');
for i=1:81
   tmp=fgets(fid_well);
end
nsample=510-81;
tmpd=fscanf(fid_well,'%f %f %f %f %f %f',[6,nsample]);
tmpd=tmpd';
% return
nsample=151;
nstart=250;
nt=nsample-1;
c=zeros(nt+1,5);
c(:,1)=tmpd(nstart:nstart+nsample-1,5)/1000;
c(:,2)=tmpd(nstart:nstart+nsample-1,3)/1000;
c(:,3)=tmpd(nstart:nstart+nsample-1,6)/1000;
c(:,4)=tmpd(nstart:nstart+nsample-1,4);
%%
c(:,2)=-0.055*c(:,1).^2+1.017*c(:,1)-1.031; % 横波速度
c(:,3)=1.74*c(:,1).^0.25; %密度
smoothstep=5;
for i=1:5      %平滑2次
    c(:,1)=smooth(c(:,1),smoothstep);
    c(:,2)=smooth(c(:,2),smoothstep);
    c(:,3)=smooth(c(:,3),smoothstep);
    c(:,4)=smooth(c(:,4),smoothstep);
end

%%
% % 1. 经验公式 -----------------------
Por=1.16-0.41*c(:,3);  % 1.16-0.41*c(:,3)-0.05
% 
% % 2. 
f=ones(nt+1,1);
for i=1:nt+1
    [vv,cc,rho1]=bkus(f(i),c(i,3),c(i,1),c(i,2));
    c11=cc(1);c33=cc(2);c55=cc(3);c66=cc(4);c13=cc(5);
    c12=c11-2*c66;
    epsilon(i)=(c11-c33)/(2*c33);
end
smooth_p=5;
epsilon=smooth(epsilon,smooth_p);
% Por=1e14*epsilon-Por;
Por=Por-1e14*epsilon;
% Por=5e14*Por-(c(:,2)+1*c(:,3))/2;
% % return
% % 归一化
% xmax=max(Por);xmin=min(Por);
% rx=(0.2-0.1)/(xmax-xmin);
% Por=(Por-xmin)*rx+0.005;
% %%
c(:,4)=Por;
%%
% F
% gamma_dry=mean(c(:,1))^2/mean(c(:,2))^2-0.05; % 躁
gamma_dry=mean(c(:,1))^2/mean(c(:,2))^2-0.35; % 躁
% gamma_dry=mean(c(:,1))^2/mean(c(:,2))^2;
M=c(:,1).^2.*c(:,3);
U=c(:,2).^2.*c(:,3);
f=M-gamma_dry*U;
F=f./c(:,4);
c(:,5)=F;
% dry vp
c(:,1)=((M-f)./c(:,3)).^(1/2);
% return
%%
figure();
subplot(1,5,1);plot(c(:,1),'k');view(90,90);ylabel('Vp(km.s^-^1)');xlabel('Time (ms)');xlim([0 150]);%ylim([3.3 5]);
subplot(1,5,2);plot(c(:,2),'k');view(90,90);ylabel('Vs(km.s^-^1)');xlim([0 150]);%ylim([2 2.6]);
subplot(1,5,3);plot(c(:,3),'k');view(90,90);ylabel('rho(g.cm^-^3)');xlim([0 150]);%ylim([2.4 2.6]);
subplot(1,5,4);plot(c(:,4),'k');view(90,90);ylabel('Por');xlim([0 150]);%ylim([0.03 0.07]);
subplot(1,5,5);plot(c(:,5),'k');view(90,90);ylabel('F');xlim([0 150]);%ylim([0 inf]);

% return

%% % % % % % % % wavelet%%%%%%%%
tlen=(nwt-1)*dt;    fm=30;
wl=ricker1(dt,fm,tlen);
figure;plot(wl,'linewidth',1.5);grid on;title(['主频fm=' num2str(fm)]);
w=zeros(nwt,nt);
for ii=1:nt
    w(:,ii)=wl;
end

ww=zeros(nt);
for it1=1:nt
    for it2=1:nt
        tt0=it2-it1+nt0+1;
        if tt0<=nwt&&tt0>=1
            ww(it1,it2)=w(tt0,it2);
        end
    end
end

%% % % % % reflectivity%%%%%%%%%%% zoeppritz正演,合成地震记录
rpp=zeros(nt,ntheta);
rps=zeros(nt,ntheta);
for itheta=1:ntheta;
    theta=pi*itheta*dtheta/180;
    for i=1:nt;
        vs1=c(i,2);vs2=c(i+1,2);
        vp1=c(i,1);vp2=c(i+1,1);
        rho1=c(i,3);rho2=c(i+1,3);
        Por1=c(i,4);Por2=c(i+1,4);
        F1=c(i,5);F2=c(i+1,5);
        
        
        A=zeros(4,4);
        B=zeros(4,1);
        %%
        H12=(vs1^2)/(vp1^2+F1*Por1/rho1);
        H13=(vp2^2+F2*Por2/rho2)^(1/2)/(vp1^2+F1*Por1/rho1)^(1/2);
        H14=(vs2^2)/(vp1^2+F1*Por1/rho1);
        %没问题
        H22=vs1/(vp1^2+F1*Por1/rho1)^(1/2);
        H23=(vp2^2+F2*Por2/rho2)/(vp1^2+F1*Por1/rho1);
        H24=vs2/(vp1^2+F1*Por1/rho1)^(1/2);
        %没问题
        H32a=(vp1^2+F1*Por1/rho1)^(1/2)/vs1;
        H32b=H22;
        H33a=vs2^2/vs1^2;
        H33b=H23;
        H34a=vs2*(vp1^2+F1*Por1/rho1)^(1/2)/vs1^2;
        H34b=H14;
        %没问题
        H41=H12;
        H42=H12;
        H43a=H13;
        H43b=H14;
        H44=H14;
        % 没问题
        G41=H12;
        %% A(1,:)
        A(1,1)=-sin(theta);
        A(1,2)=-sqrt(1-H12*sin(theta)^2);
        A(1,3)=H13*sin(theta);
        A(1,4)=-sqrt(1-H14*sin(theta)^2);
        %没问题
        %% A(2,:)
        A(2,1)=cos(theta);
        A(2,2)=-sin(theta)*H22;
        A(2,3)=sqrt(1-H23*sin(theta)^2);
        A(2,4)=sin(theta)*H24;
        %没问题
        %% A(3,:)
        A(3,1)=sin(2*theta);
        A(3,2)=H32a-2*sin(theta)^2*H32b;
        A(3,3)=2*(rho2/rho1)*H33a*sqrt(sin(theta)^2-H33b*sin(theta)^4);
        A(3,4)=-(rho2/rho1)*H34a*(1-2*H34b*sin(theta)^2);
        %没问题
        %% A(4,:)
        A(4,1)=2*H41*sin(theta)^2-1;
        A(4,2)=2*H42*sqrt(sin(theta)^2-H42*sin(theta)^4);
        A(4,3)=(rho2/rho1)*H43a*(1-2*H43b*sin(theta)^2);
        A(4,4)=2*(rho2/rho1)*H44*sqrt(sin(theta)^2-H44*sin(theta)^4);
        %没问题
        %% B(:,1)
        B(1,1)=sin(theta);
        B(2,1)=cos(theta);
        B(3,1)=sin(2*theta);
        B(4,1)=1-2*G41*sin(theta)^2;
        %没问题
        RR=(inv(A))*B;
        rpp(i,itheta)=RR(1);
        rps(i,itheta)=RR(2);
    end
end
% % % % % synseis%%%%%%%%%%%%%
d=zeros(1,2*nt*ntheta);% 地震数据存为一行
dpp=zeros(nt,ntheta);dps=zeros(nt,ntheta);
for itheta=1:ntheta;
    dpp(:,itheta)=ww*rpp(:,itheta);
    dps(:,itheta)=ww*rps(:,itheta);
end
for itheta=1:ntheta;
    for it=1:nt;
        d(it+(itheta-1)*nt)=dpp(it,itheta);
        d(it+(itheta-1+ntheta)*nt)=dps(it,itheta);
    end;
end;
% return
%d是纵横波联合数据

% figure;
% subplot(141);wigb(dpp);
% title('PP');
% xlabel('Degree');
% ylabel('Time(ms)');
% set(gca,'XTick',0:5:10);set(gca,'XTickLabel',{'0','25','50'});
% set(gca,'YTick',0:40:120);set(gca,'YTickLabel',{'0','40','80','120'});
% subplot(142);wigb(dps);
% title('PS');
% xlabel('Degree');
% set(gca,'XTick',0:5:10);set(gca,'XTickLabel',{'0','25','50'});
% set(gca,'YTick',0:40:120);set(gca,'YTickLabel',{'0','40','80','120'});

%%%%%%display 含噪数据显示
% load noisesnr3.mat
% d=d+bbb(1:2700);
% 
% 
% for itheta=1:ntheta;
%     for it=1:nt;
%         dpp1(it,itheta)=d(it+(itheta-1)*nt);
%         dps1(it,itheta)=d(it+(itheta-1+ntheta)*nt);
%     end;
% end;


% subplot(143);wigb(dpp1);
% title('PP with noise ');
% xlabel('Degree');
% set(gca,'XTick',0:5:10);set(gca,'XTickLabel',{'0','25','50'});
% set(gca,'YTick',0:40:120);set(gca,'YTickLabel',{'0','40','80','120'});
% % ylabel('Time(s)','FontSize',10,'FontWeight','bold');
% subplot(144);
% wigb(dps1);
% title('PS with noise ');
% xlabel('Degree');
% set(gca,'XTick',0:5:10);set(gca,'XTickLabel',{'0','25','50'});
% set(gca,'YTick',0:40:120);set(gca,'YTickLabel',{'0','40','80','120'});
% return


%% % % % % % % % % % % % % % % % % % % % 反演
% invesion W is grendient matrix, Wt is transpose of grendient matrix 正演算子关于模型参数的一阶偏导数
% H is second order derivative matrix, wis wavelet matrix varies with time
%%%%initial model
ct=c;  % 原始测井曲线
c1=c;  % 平滑之后的初始模型

nn=40;  % 平滑参数
% for i=1:nn
%    c1(:,1)=smooth(c(:,1),20,'lowess');
%    c1(:,2)=smooth(c(:,2),20,'lowess');
%    c1(:,3)=smooth(c(:,3),20,'lowess');
% %    c=c1;
% end
% ct=c1;

% for i=1:nn
%    c1(:,1)=smooth(c(:,1),20);
%    c1(:,2)=smooth(c(:,2),20);
%    c1(:,3)=smooth(c(:,3),20);
% %    c=c1;
% end

%  window=101;
% % window=51;
% % window=31;
% c1(:,1)=Filter(c(:,1),window);
% c1(:,2)=Filter(c(:,2),window);
% c1(:,3)=Filter(c(:,3),window);
% c=c1;

for i=1:20      %平滑2次
    c1(:,1)=smooth(c(:,1),30,'lowess');
    c1(:,2)=smooth(c(:,2),30,'lowess');
    c1(:,3)=smooth(c(:,3),30,'lowess');
    c1(:,4)=smooth(c(:,4),30,'lowess');
    c1(:,5)=smooth(c(:,5),30,'lowess');
    c=c1;                                     %这样使得下面第一次反演从初始模型开始计算, 平滑之后的初始模型
end
u=[c1(1:nt,1);c1(1:nt,2);c1(1:nt,3);c1(1:nt,4);c1(1:nt,5)];

mean_c1(:,1)=mean(c1(:,1))*ones(nt,1);   %求Vp的均值，存成1列
mean_c1(:,2)=mean(c1(:,2))*ones(nt,1);   %求Vs的均值，存成1列
mean_c1(:,3)=mean(c1(:,3))*ones(nt,1);   %求ρ的均值，存成1列
mean_c1(:,4)=mean(c1(:,4))*ones(nt,1);   %求invQp的均值，存成1列
mean_c1(:,5)=mean(c1(:,5))*ones(nt,1);   %求invQs的均值，存成1列
mean_c1=[mean_c1; mean_c1(end,:)];       %将最后一行复制一遍
%  Cov_c=cov(c1);

% Construct Wt % % % % % % % % % % % % % % % % % % % % % % % % % % % %
Wt=zeros(5*nt,ntheta*nt*2);  %  Wt is transpose of grendient matrix 正演算子关于模型参数的一阶偏导数
for it1=1:nt
    for itheta=1:ntheta
        angle=itheta*dtheta;
        vp1=c(it1,1);vp2=c(it1+1,1);
        vs1=c(it1,2);vs2=c(it1+1,2);
        rho1=c(it1,3);rho2=c(it1+1,3);
        Por1=c(it1,4);Por2=c(it1+1,4);
        F1=c(it1,5);F2=c(it1+1,5);
        if it1>1
            rho0=c(it1-1,3);vp0=c(it1-1,1);vs0=c(it1-1,2);Por0=c(it1-1,4);F0=c(it1-1,5);
            rp2=Rp2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);       % p波反射系数 下行波 正演算子关于下层纵波速度的一阶偏导数
            rs2=Rs2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
            rr2=Rr2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
            rpor2=Rpor2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
            rf2=Rf2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
        end
        
        rp1=Rp1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);         % p波反射系数 上行波 正演算子关于上层纵波速度的一阶偏导数
        rs1=Rs1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
        rr1=Rr1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
        rpor1=Rpor1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
        rf1=Rf1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
        
        for it2=1:nt
            tt0=it2-it1+nt0+1;
            if tt0<=nwt&&tt0>=1
                %             in c program tt0=it2-it1+nt0;
                %             in c program if tt0<nwt&&tt0>=;
                Wt(it1,it2+itheta*nt-nt)=w(tt0,it1)*rp1(1);
                Wt(it1+nt,it2+itheta*nt-nt)=w(tt0,it1)*rs1(1);
                Wt(it1+2*nt,it2+itheta*nt-nt)=w(tt0,it1)*rr1(1);
                Wt(it1+3*nt,it2+itheta*nt-nt)=w(tt0,it1)*rpor1(1);
                Wt(it1+4*nt,it2+itheta*nt-nt)=w(tt0,it1)*rf1(1);
                
                Wt(it1,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rp1(2);
                Wt(it1+nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rs1(2);
                Wt(it1+2*nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rr1(2);
                Wt(it1+3*nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rpor1(2);
                Wt(it1+4*nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rf1(2);
            end
            tt1=it2-it1+1+nt0+1;
            if tt1<=nwt&&tt1>=1&&it1>1
                %  in c program tt0=it2-it1+1+nt0;
                %  in c program if tt1<nwt&&tt1>=0&&it1>0;
                
                Wt(it1,it2+itheta*nt-nt)=Wt(it1,it2+itheta*nt-nt)+w(tt1,it1-1)*rp2(1);
                Wt(it1+nt,it2+itheta*nt-nt)=Wt(it1+nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rs2(1);
                Wt(it1+2*nt,it2+itheta*nt-nt)=Wt(it1+2*nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rr2(1);
                Wt(it1+3*nt,it2+itheta*nt-nt)=Wt(it1+3*nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rpor2(1);
                Wt(it1+4*nt,it2+itheta*nt-nt)=Wt(it1+4*nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rf2(1);
                
                Wt(it1,it2+(itheta-1+ntheta)*nt)=Wt(it1,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rp2(2);
                Wt(it1+nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rs2(2);
                Wt(it1+2*nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+2*nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rr2(2);
                Wt(it1+3*nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+3*nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rpor2(2);
                Wt(it1+4*nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+4*nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rf2(2);
            end
        end
    end
end



%% PP inversion
wp=zeros(5*nt,ntheta*nt);
for k=1:5*nt
    for l=1:ntheta*nt;
        wp(k,l)=Wt(k,l);
    end
end
WtW1=wp*wp';


%%  PS inversion
ws=zeros(5*nt,ntheta*nt);
for k=1:5*nt
    for l=1:ntheta*nt;
        ws(k,l)=Wt(k,l+ntheta*nt);
    end
end
WtW2=ws*ws';

%% PP inversion
%   GTG=WtW1*WtW1';
%%  joint inversion
  % Wt is transpose of grendient matrix 目标函数关于模型参数的二阶偏导数的第一项

  % GTG=Wt*Wt';
GTG=alpha1*WtW1+alpha2*WtW2; 

%sigma=1e-12;

I=eye(5*nt,5*nt)*1e0;

Cov_c=cov(ct);      % 模型参数的方差矩阵
unitm1=eye(nt);
invCm1=inv(Cov_c);
Cm1=kron(invCm1,unitm1);

% inversion
dcal=zeros(1,2*nt*ntheta);   %一行


%% 开始迭代求解
nter=5;
for iter=1:nter
    fprintf('第 %d 次迭代:\n',iter);
    % Construct Wtd % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for itheta=1:ntheta;
        theta=pi*itheta*dtheta/180;
        for i=1:nt;
            vs1=c(i,2);  vs2=c(i+1,2);
            vp1=c(i,1);vp2=c(i+1,1);
            rho1=c(i,3);rho2=c(i+1,3);
            Por1=c(i,4);Por2=c(i+1,4);
            F1=c(i,5);F2=c(i+1,5);
            
            
                    A=zeros(4,4);
        B=zeros(4,1);
        %%
        H12=(vs1^2)/(vp1^2+F1*Por1/rho1);
        H13=(vp2^2+F2*Por2/rho2)^(1/2)/(vp1^2+F1*Por1/rho1)^(1/2);
        H14=(vs2^2)/(vp1^2+F1*Por1/rho1);
        %没问题
        H22=vs1/(vp1^2+F1*Por1/rho1)^(1/2);
        H23=(vp2^2+F2*Por2/rho2)/(vp1^2+F1*Por1/rho1);
        H24=vs2/(vp1^2+F1*Por1/rho1)^(1/2);
        %没问题
        H32a=(vp1^2+F1*Por1/rho1)^(1/2)/vs1;
        H32b=H22;
        H33a=vs2^2/vs1^2;
        H33b=H23;
        H34a=vs2*(vp1^2+F1*Por1/rho1)^(1/2)/vs1^2;
        H34b=H14;
        %没问题
        H41=H12;
        H42=H12;
        H43a=H13;
        H43b=H14;
        H44=H14;
        % 没问题
        G41=H12;
        %% A(1,:)
        A(1,1)=-sin(theta);
        A(1,2)=-sqrt(1-H12*sin(theta)^2);
        A(1,3)=H13*sin(theta);
        A(1,4)=-sqrt(1-H14*sin(theta)^2);
        %没问题
        %% A(2,:)
        A(2,1)=cos(theta);
        A(2,2)=-sin(theta)*H22;
        A(2,3)=sqrt(1-H23*sin(theta)^2);
        A(2,4)=sin(theta)*H24;
        %没问题
        %% A(3,:)
        A(3,1)=sin(2*theta);
        A(3,2)=H32a-2*sin(theta)^2*H32b;
        A(3,3)=2*(rho2/rho1)*H33a*sqrt(sin(theta)^2-H33b*sin(theta)^4);
        A(3,4)=-(rho2/rho1)*H34a*(1-2*H34b*sin(theta)^2);
        %没问题
        %% A(4,:)
        A(4,1)=2*H41*sin(theta)^2-1;
        A(4,2)=2*H42*sqrt(sin(theta)^2-H42*sin(theta)^4);
        A(4,3)=(rho2/rho1)*H43a*(1-2*H43b*sin(theta)^2);
        A(4,4)=2*(rho2/rho1)*H44*sqrt(sin(theta)^2-H44*sin(theta)^4);
        %没问题
        %% B(:,1)
        B(1,1)=sin(theta);
        B(2,1)=cos(theta);
        B(3,1)=sin(2*theta);
        B(4,1)=1-2*G41*sin(theta)^2;
        %没问题
            
            RR=(inv(A))*B;
            rpp(i,itheta)=RR(1);
            rps(i,itheta)=RR(2);
        end
    end
    for itheta=1:ntheta;
        dpp(:,itheta)=ww*rpp(:,itheta);
        dps(:,itheta)=ww*rps(:,itheta);
    end
    for itheta=1:ntheta;
        for it=1:nt;
            dcal(it+(itheta-1)*nt)=dpp(it,itheta);
            dcal(it+(itheta-1+ntheta)*nt)=dps(it,itheta);
        end;
    end;
    
    dmiff=d-dcal;  % 差异数据
    %dmiff=smooth(dmiff);
    %
    
    %% pp inversion
    dmiffpp=zeros(nt*ntheta,1);
    for l=1:nt*ntheta
        dmiffpp(l)=dmiff(l);
    end
    Wtd1=wp*dmiffpp;
    misfit=sum(dmiffpp.*dmiffpp);
    fprintf('misfit=\t%f\n',misfit);
    
    %% ps inversion
    dmiffps=zeros(nt*ntheta,1);
    for l=1:nt*ntheta
        dmiffps(l)=dmiff(l+nt*ntheta);
    end
    Wtd2=ws*dmiffps;
    misfit=sum(dmiffps.*dmiffps);
    fprintf('misfit=\t%f\n',misfit);
    
        %% PP inversion
%   Wtd=alpha1*Wtd1;
     %% joint inversion
 Wtd=alpha1*Wtd1+alpha2*Wtd2;
 
%     Wtd=Wt*dmiff';   % 目标函数对模型参数一阶偏导的第一项
    
    misfit=sum(dmiff.*dmiff);
    fprintf('misfit=\t%f\n',misfit);
    
    
    aaa=[c(1:nt,1);c(1:nt,2);c(1:nt,3);c(1:nt,4);c(1:nt,5)];
       
       
    Q1=modificauchy1testDCm1(c,mean_c1,nt,Cov_c);       % 正则化项
    %     Q1=modificauchy1testCm1(c,nt,Cov_c);
    %     Q1=modifycauchy1test2(c,nt,Cov_c);
    %     Q1=cauchy1testCm1(c,nt,Cov_c);
    %      Q1=cauchy1testn(c,nt,Cov_c);
    %Q4=cauchy1testnew2(c,nt,Cov_c);
    %     Q1=cauchy1testnew(c,nt,Cov_c);
    %Wtd1=Wtd+beta*Q1;
    %Q=cauchy2testCm(c,nt,Cov_c);
    %        for it111=1:3;
    %Wtd1=Wtd-beta*Q4;
    Hessian=GTG+beta*Q1;            % 目标函数对模型参数的二阶偏导数
    condH=cond(Hessian);
    fprintf('海森矩阵条件数=\t%f\n',condH);
    

% beta=1e-3;   %% PP-PS cauchy 3e-1 无噪 5630； 有噪 2-0 8272;   PP 6E-2 7665
%         r=(GTG+beta*I)\(Wtd);
%     r=(GTG+beta*Q1)\Wtd;          % B\A等价于inv(B)*A    m(n+1)=m(n)-p*inv(H)*Wtd;  模型参数的步长 beta=5e-2;
     r=(GTG+beta*Cm1)\(Wtd);    % beta=1e-3

%     r=(GTG+beta*Q1)\(Wtd-beta*Q1*([c(1:nt,1);c(1:nt,2);c(1:nt,3);c(1:nt,4);c(1:nt,5)]-u)); % beta=1e-2; 

%       r=(GTG+beta*Cm1+beta*Q1)\(Wtd); %beta=1e-3; 
%       r=(GTG+beta*Cm1+beta*Q1)\(Wtd-beta*(Cm1+Q1)*([c(1:nt,1);c(1:nt,2);c(1:nt,3);c(1:nt,4);c(1:nt,5)]-u)); % beta=5e-4;


%     r=(GTG+beta*(Cm1+Q1))\(Wtd-beta*(Cm1+Q1)*(aaa-u)); % beta=5e-4; 

    %     c2=c;
    %
    %     c2(1:nt,1)=c2(1:nt,1)+r(1:nt);%%最后一个点不更新
    %     c2(1:nt,2)=c2(1:nt,2)+r(nt+1:nt*2);
    %     c2(1:nt,3)=c2(1:nt,3)+r(nt*2+1:nt*3);
    %
    %     Q1=cauchy1testCm1(c2,nt,Cov_c);
    %      end
    
    %% withnoise && PS
    %     nnn=3;
    % for i=1:nnn;
    %    r(1:nt)=smooth(r(1:nt),5,'lowess');
    %    r(nt+1:2*nt)=smooth(r(nt+1:2*nt),5,'lowess');
    %    r(2*nt+1:3*nt)=smooth(r(2*nt+1:3*nt),5,'lowess');
    % end;
    
    
    %%
    
    % for i2=1:nnn+1
    %      r(2*nt+1:3*nt)=smooth(r(2*nt+1:3*nt),7,'lowess');
    % end
    
    c(1:nt,1)=c(1:nt,1)+r(1:nt);%%最后一个点不更新        % 平滑之后的初始模型
    c(1:nt,2)=c(1:nt,2)+r(nt+1:nt*2);
    c(1:nt,3)=c(1:nt,3)+r(nt*2+1:nt*3);
    c(1:nt,4)=c(1:nt,4)+r(nt*3+1:nt*4);
    c(1:nt,5)=c(1:nt,5)+r(nt*4+1:nt*5);
    
    %     c(1:nt,3)=c(1:nt,3)+r(nt*2+1:nt*3)*1.4;   %%%%%  withnoise
    % %
    
    %%  withnoise
    %    nnn=4;
    % for i=1:nnn;
    %    c(1:nt,1)=smooth(c(1:nt,1),5,'lowess');
    %    c(1:nt,2)=smooth(c(1:nt,2),5,'lowess');
    %    %c(1:nt,3)=smooth(c(1:nt,3),10,'lowess');
    % end;
    % for i2=1:nnn-1
    %     c(1:nt,3)=smooth(c(1:nt,3),7,'lowess');
    % end
    %
    %
    % %% ----------往下计算模型参数更新后的目标函数关于模型参数的二阶偏导数的第一项
    Wt=zeros(5*nt,ntheta*nt*2);   % 正演算子关于模型参数的一阶偏导数
    for it1=1:nt
        for itheta=1:ntheta
            angle=itheta*dtheta;
            vp1=c(it1,1);vp2=c(it1+1,1);
            vs1=c(it1,2);vs2=c(it1+1,2);
            rho1=c(it1,3);rho2=c(it1+1,3);
            Por1=c(it1,4);Por2=c(it1+1,4);
            F1=c(it1,5);F2=c(it1+1,5);
            if it1>1
                rho0=c(it1-1,3);vp0=c(it1-1,1);vs0=c(it1-1,2);Por0=c(it1-1,4);F0=c(it1-1,5);
                rp2=Rp2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);       % p波反射系数 下行波 正演算子关于下层纵波速度的一阶偏导数
                rs2=Rs2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
                rr2=Rr2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
                rpor2=Rpor2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
                rf2=Rf2(vp0,vs0,rho0,Por0,F0,vp1,vs1,rho1,Por1,F1,angle);
            end
            
            rp1=Rp1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);         % p波反射系数 上行波 正演算子关于上层纵波速度的一阶偏导数
            rs1=Rs1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
            rr1=Rr1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
            rpor1=Rpor1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
            rf1=Rf1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle);
            
            for it2=1:nt
                tt0=it2-it1+nt0+1;
                if tt0<=nwt&&tt0>=1
                    %             in c program tt0=it2-it1+nt0;
                    %             in c program if tt0<nwt&&tt0>=;
                    Wt(it1,it2+itheta*nt-nt)=w(tt0,it1)*rp1(1);
                    Wt(it1+nt,it2+itheta*nt-nt)=w(tt0,it1)*rs1(1);
                    Wt(it1+2*nt,it2+itheta*nt-nt)=w(tt0,it1)*rr1(1);
                    Wt(it1+3*nt,it2+itheta*nt-nt)=w(tt0,it1)*rpor1(1);
                    Wt(it1+4*nt,it2+itheta*nt-nt)=w(tt0,it1)*rf1(1);
                    
                    Wt(it1,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rp1(2);
                    Wt(it1+nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rs1(2);
                    Wt(it1+2*nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rr1(2);
                    Wt(it1+3*nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rpor1(2);
                    Wt(it1+4*nt,it2+(itheta-1+ntheta)*nt)=w(tt0,it1)*rf1(2);
                end
                tt1=it2-it1+1+nt0+1;
                if tt1<=nwt&&tt1>=1&&it1>1
                    %  in c program tt0=it2-it1+1+nt0;
                    %  in c program if tt1<nwt&&tt1>=0&&it1>0;
                    
                    Wt(it1,it2+itheta*nt-nt)=Wt(it1,it2+itheta*nt-nt)+w(tt1,it1-1)*rp2(1);
                    Wt(it1+nt,it2+itheta*nt-nt)=Wt(it1+nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rs2(1);
                    Wt(it1+2*nt,it2+itheta*nt-nt)=Wt(it1+2*nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rr2(1);
                    Wt(it1+3*nt,it2+itheta*nt-nt)=Wt(it1+3*nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rpor2(1);
                    Wt(it1+4*nt,it2+itheta*nt-nt)=Wt(it1+4*nt,it2+itheta*nt-nt)+w(tt1,it1-1)*rf2(1);
                    
                    Wt(it1,it2+(itheta-1+ntheta)*nt)=Wt(it1,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rp2(2);
                    Wt(it1+nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rs2(2);
                    Wt(it1+2*nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+2*nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rr2(2);
                    Wt(it1+3*nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+3*nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rpor2(2);
                    Wt(it1+4*nt,it2+(itheta-1+ntheta)*nt)=Wt(it1+4*nt,it2+(itheta-1+ntheta)*nt)+w(tt1,it1-1)*rf2(2);
                end
            end
        end
    end
    
        %% PP inversion
    wp=zeros(5*nt,ntheta*nt);
for k=1:5*nt
    for l=1:ntheta*nt;
        wp(k,l)=Wt(k,l);
    end
end
WtW1=wp*wp';
%% PS inversion
ws=zeros(5*nt,ntheta*nt);
for k=1:5*nt
    for l=1:ntheta*nt;
        ws(k,l)=Wt(k,l+ntheta*nt);
    end
end
WtW2=ws*ws';

%% PP inversion

% WtW=alpha1*WtW1;
%% PP PS joint inversion

WtW=alpha1*WtW1+alpha2*WtW2; 
%     WtW=Wt*Wt';      % 目标函数关于模型参数的二阶偏导数的第一项
    
%%
    GTG=WtW;     % 目标函数关于模型参数的二阶偏导数的第一项
end
%     c(1:nt,1)=c(1:nt,1)+r(1:nt);%%最后一个点不更新
%     c(1:90,2)=c(1:90,2)+r(nt+1:nt+90)*8;
% %     c(1:nt,3)=c(1:nt,3)+r(nt*2+1:nt*3);
%     c(1:50,3)=c(1:50,3)+r(nt*2+1:nt*2+50)*5;
%

fprintf('迭代结束\n');
rmse=sum((c-ct).*(c-ct))/sum(ct.*ct);
fprintf('结果误差=\t%f\n',rmse);

%    nnn=6;
% for i=1:nnn;
%    c(1:nt,1)=smooth(c(1:nt,1),5,'lowess');
%    c(1:nt,2)=smooth(c(1:nt,2),5,'lowess');
%    c(1:nt,3)=smooth(c(1:nt,3),7,'lowess');
% end;

%%
% x=(1:nt)*dt;
% figure(3);
% subplot(1,3,1);h=plot(ct(1:nt,1)/1000,x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c(1:nt,1)/1000,x,'b');set(h,'linewidth',1.5);
% h=plot(c2(1:nt,1)/1000,x,'k');set(h,'linewidth',1.5);
% %xlim([3 3.7]);
% xlim([2.5 3.1]);%ylim([0 0.6]);
% ylabel('Time(s)','fontsize',12);set(gca,'fontsize',12);
% xlabel('Vp(km/s)','fontsize',12);set(gca,'yDir','rev');
%
% subplot(1,3,2);h=plot(ct(1:nt,2)/1000,x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c(1:nt,2)/1000,x,'b');set(h,'linewidth',1.5);
% h=plot(c2(1:nt,2)/1000,x,'k');set(h,'linewidth',1.5);
% xlabel('Vs(km/s)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([1.05 1.41]);%ylim([0 0.6]);
% %xlim([1.5 2.2]);
%
% subplot(1,3,3);h=plot(ct(1:nt,3)/1000,x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c(1:nt,3)/1000,x,'b');set(h,'linewidth',1.5);
% h=plot(c2(1:nt,3)/1000,x,'k');set(h,'linewidth',1.5);
% xlabel('Rho(g/cc)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([2.2 2.38]);

% c_cau=zeros(400,1);
% c_cau(:,1)=c(101:nt,1)-ct(101:nt,1);
% c_cau(:,2)=c(101:nt,2)-ct(101:nt,2);
% c_cau(:,3)=c(101:nt,3)-ct(101:nt,3);

Cor=zeros(3,1);

Cor1=corrcoef(c(101:nt,1),ct(101:nt,1));
Cor(1,1)=Cor1(2,1);

Cor2=corrcoef(c(101:nt,2),ct(101:nt,2));
Cor(2,1)=Cor2(2,1);

Cor3=corrcoef(c(101:nt,3),ct(101:nt,3));
Cor(3,1)=Cor3(2,1);




c_mcau(:,1)=c(:,1)-ct(:,1);
c_mcau(:,2)=c(:,2)-ct(:,2);
c_mcau(:,3)=c(:,3)-ct(:,3);

% c_cau(:,1)=c(:,1)-ct(:,1);
% c_cau(:,2)=c(:,2)-ct(:,2);
% c_cau(:,3)=c(:,3)-ct(:,3);

dc=ct(1:1:nt,:)-c(1:1:nt,:);
%

figure;


subplot(151);h=plot(ct(1:1:nt,1),'k');set(h,'linewidth',1);hold on;
% h=plot(c(1:1:nt,1),x,'bh');set(h,'Markersize',1.5);
h=plot(c(1:1:nt,1),'k-.');set(h,'linewidth',2);
h=plot(c1(1:1:nt,1),'k:');set(h,'linewidth',1.0);
% ylim([2.5 3.5]);
xlabel('Time(ms)','fontsize',12);
% set(gca,'XTick',1:53:107);set(gca,'XTickLabel',{'2052','2158','2264'});
ylabel('Vp(km/s)','fontsize',12);%set(gca,'yDir','rev');
view(90,90)

% x=(1:nt)*dt;
subplot(152);h=plot(ct(1:1:nt,2),'k');set(h,'linewidth',1);hold on;
% h=plot(c(1:1:nt,2),x,'bh');set(h,'Markersize',1.5);
h=plot(c(1:1:nt,2),'k-.');set(h,'linewidth',2);
h=plot(c1(1:1:nt,2),'k:');set(h,'linewidth',1.0);
% set(gca,'XTick',1:53:107);set(gca,'XTickLabel',{'2052','2158','2264'});
ylabel('Vs(km/s)');%set(gca,'yDir','rev');
% ylim([1.5 2.5]);%ylim([0 0.6]);
view(90,90)

subplot(153);h=plot(ct(1:1:nt,3),'k');set(h,'linewidth',1);hold on;
% h=plot(c(1:1:nt,3),x,'bh');set(h,'Markersize',1.5);
h=plot(c(1:1:nt,3),'k-.');set(h,'linewidth',2);
h=plot(c1(1:1:nt,3),'k:');set(h,'linewidth',1.0);
% set(gca,'XTick',1:53:107);set(gca,'XTickLabel',{'2052','2158','2264'});
ylabel('Rho(g/cm^3)');%set(gca,'yDir','rev');
% ylim([2.4 3]);%set(gca,'XTick',[2.0 2.2 2.4 2.6 2.8]);
% ylim([2.3 2.5])
view(90,90)

subplot(154);h=plot(ct(1:1:nt,4),'k');set(h,'linewidth',1);hold on;
% h=plot(c(1:1:nt,3),x,'bh');set(h,'Markersize',1.5);
h=plot(c(1:1:nt,4),'k-.');set(h,'linewidth',2);
h=plot(c1(1:1:nt,4),'k:');set(h,'linewidth',1.0);
% set(gca,'XTick',1:53:107);set(gca,'XTickLabel',{'2052','2158','2264'});
ylabel('Por');%set(gca,'yDir','rev');
% ylim([0 0.1]);%set(gca,'XTick',[2.0 2.2 2.4 2.6 2.8]);
% ylim([0.1 0.2])
view(90,90)

subplot(155);h=plot(ct(1:1:nt,5),'k');set(h,'linewidth',1);hold on;
% h=plot(c(1:1:nt,3),x,'bh');set(h,'Markersize',1.5);
h=plot(c(1:1:nt,5),'k-.');set(h,'linewidth',2);
h=plot(c1(1:1:nt,5),'k:');set(h,'linewidth',1);
% set(gca,'XTick',1:53:107);set(gca,'XTickLabel',{'2052','2158','2264'});
ylabel('F');%set(gca,'yDir','rev');
% ylim([50 150]);%set(gca,'XTick',[2.0 2.2 2.4 2.6 2.8]);
% ylim([0 0.3])
view(90,90)
%%
% c4=zeros(nt+1,3);
%
%
% for i=1:nt+1;
%     c4(i,2)=(c(i,1).^2-2*c(i,2).^2)/(2*(c(i,1).^2-c(i,2).^2));
%     c4(i,1)=((3*c(i,1).^2-4*c(i,2).^2)/(c(i,1).^2-c(i,2).^2)).*c(i,3).*c(i,2).^2;
% end
% c4(:,1)=c4(:,1)/(1e9);
% c4(:,3)=c(:,3);
%
%
% c11=zeros(nt+1,3);
%
%
% for i=1:nt+1;
%     c11(i,2)=(c1(i,1).^2-2*c1(i,2).^2)/(2*(c1(i,1).^2-c1(i,2).^2));
%     c11(i,1)=((3*c1(i,1).^2-4*c1(i,2).^2)/(c1(i,1).^2-c1(i,2).^2)).*c1(i,3).*c1(i,2).^2;
% end
%
% c11(:,1)=c11(:,1)/(1e9);
%
% c11(:,3)=c1(:,3);
%
%
% ct1=zeros(nt+1,3);
%
%
% for i=1:nt+1;
%     ct1(i,2)=(ct(i,1).^2-2*ct(i,2).^2)/(2*(ct(i,1).^2-ct(i,2).^2));
%     ct1(i,1)=((3*ct(i,1).^2-4*ct(i,2).^2)/(ct(i,1).^2-ct(i,2).^2)).*ct(i,3).*ct(i,2).^2;
% end
%
% ct1(:,1)=ct1(:,1)/(1e9);
%
% ct1(:,3)=ct(:,3);
%
%
% x=(1:nt)*dt;
% figure(10);
% %gtext('a)')
% text(1834.37061403509,-0.0386382623224727,'\fontsize{15}\rm c)');
% subplot(1,3,1);h=plot(ct1(1:nt,1)/10,x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c4(1:nt,1)/10,x,'b');set(h,'linewidth',1.5);
% h=plot(c11(1:nt,1)/10,x,'k');set(h,'linewidth',1.5);
% xlim([0.76 1.3]);
% ylim([0.2 1])
%
%  text(1.145,0.860701754385965+0.2,'\times 10^{10}');
%
% title('Young'' modulus/(N/m^{2})','fontsize',9);
% ylabel('Time(s)','fontsize',12);set(gca,'fontsize',12);
% % xlabel('Young modulus/(N/m^{2})','fontsize',10);
% set(gca,'yDir','rev');
%
%
% subplot(1,3,2);h=plot(ct1(1:nt,2),x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c4(1:nt,2),x,'b');set(h,'linewidth',1.5);
% h=plot(c11(1:nt,2),x,'k');set(h,'linewidth',1.5);
% % xlabel('Poisson ratio（-）','fontsize',10);
% set(gca,'yDir','rev');set(gca,'fontsize',12);
% title('Poisson'' ratio（-）','fontsize',9);
% % set(gca,'XAxisLocation','top');
% % xlim([1.05 1.41]);%ylim([0 0.6]);
% xlim([0.355 0.39]);
% ylim([0.2 1])
%
% subplot(1,3,3);h=plot(ct1(1:nt,3),x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c4(1:nt,3),x,'b');set(h,'linewidth',1.5);
% h=plot(c11(1:nt,3),x,'k');set(h,'linewidth',1.5);
% % xlabel('Density(g/cc)','fontsize',10);
% title('Density(kg/m^3)','fontsize',9);
% set(gca,'yDir','rev');set(gca,'fontsize',12);
% % set(gca,'XAxisLocation','top');
% xlim([2250 2380]);
% set(gca,'XTick',[2280 2350]);
% ylim([0.2 1])
%%


% x=(1:1:nt)*dt;
% figure(4);
%
% subplot(131);h=plot(ct(1:1:nt,1)/1000,x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c(1:1:nt,1)/1000,x,'b-');set(h,'linewidth',1.5);
% h=plot(c1(1:1:nt,1)/1000,x,'k');set(h,'linewidth',1.5);
% %xlim([3 3.7]);
% xlim([2.5 3.1]);%ylim([0 0.6]);
% ylabel('Time(s)','fontsize',12);set(gca,'fontsize',12);
% xlabel('Vp(km/s)','fontsize',12);set(gca,'yDir','rev');
% ylim([0.2 1])
%
% % x=(1:nt)*dt;
% subplot(132);h=plot(ct(1:1:nt,2)/1000,x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c(1:1:nt,2)/1000,x,'b-');set(h,'linewidth',1.5);
% h=plot(c1(1:1:nt,2)/1000,x,'k');set(h,'linewidth',1.5);
% xlabel('Vs(km/s)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([1.1 1.41]);%ylim([0 0.6]);
% %xlim([1.5 2.2]);
% ylim([0.2 1])
%
% subplot(133);h=plot(ct(1:1:nt,3)/1000,x,'r');set(h,'linewidth',1.5);hold on;
% h=plot(c(1:1:nt,3)/1000,x,'b-');set(h,'linewidth',1.5);
% h=plot(c1(1:1:nt,3)/1000,x,'k');set(h,'linewidth',1.5);
% xlabel('Rho(g/cm^3)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([2.25 2.38]);set(gca,'XTick',[2.3 2.38]);
% ylim([0.2 1])


% x=(1:nt)*dt;
% figure(5);
% subplot(131);h=plot(ct(1:nt,1)/1000,x,'k');set(h,'linewidth',1.5);hold on;
% h=plot(c1(1:nt,1)/1000,x,'k--');set(h,'linewidth',1.5);
% %xlim([3 3.7]);
% xlim([2.5 3.1]);%ylim([0 0.6]);
% ylabel('Time(s)','fontsize',12);set(gca,'fontsize',12);
% xlabel('Vp(km/s)','fontsize',12);set(gca,'yDir','rev');
% ylim([0.2 1])
% xlim([2.5 3.1]);%ylim([0 0.6]);
%
%
% subplot(132);h=plot(ct(1:nt,2)/1000,x,'k');set(h,'linewidth',1.5);hold on;
% h=plot(c1(1:nt,2)/1000,x,'k--');set(h,'linewidth',1.5);
% xlabel('Vs(km/s)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([1.1 1.41]);%ylim([0 0.6]);
% %xlim([1.5 2.2]);
% ylim([0.2 1])
%
% subplot(133);h=plot(ct(1:nt,3)/1000,x,'k');set(h,'linewidth',1.5);hold on;
% h=plot(c1(1:nt,3)/1000,x,'k--');set(h,'linewidth',1.5);
% xlabel('Rho(g/cm^3)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([2.25 2.38]);set(gca,'XTick',[2.3 2.38]);
% ylim([0.2 1])

% %xlim([1.8 2.6]);
% %ylim([0 0.6]);
%
% % c_cau=c_mcau;
% % clear c_mcau;
% load c_cau23.mat
% % load c_mcau22.mat
% c_cau=abs(c_cau);
% c_mcau=abs(c_mcau);
%
% x=(1:nt)*dt;
% figure();
% subplot(1,3,1);h=plot(c_cau(1:nt,1)/1000,x,'b');set(h,'linewidth',0.5);hold on;
% h=plot(c_mcau(1:nt,1)/1000,x,'r');set(h,'linewidth',0.5);
% %xlim([3 3.7]);
% xlim([0 0.09]);
% ylabel('Time(s)','fontsize',12);set(gca,'fontsize',12);
% xlabel('Vp(km/s)','fontsize',12);set(gca,'yDir','rev');
% ylim([0.2 1])
%
%
% subplot(1,3,2);h=plot(c_cau(1:nt,2)/1000,x,'b');set(h,'linewidth',0.5);hold on;
% h=plot(c_mcau(1:nt,2)/1000,x,'r');set(h,'linewidth',0.5);
% xlabel('Vs(km/s)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([0 0.06]);
% %xlim([1.5 2.2]);
% ylim([0.2 1])
% % Markersize
%
% subplot(1,3,3);h=plot(c_cau(1:nt,3)/1000,x,'b');set(h,'linewidth',0.5);hold on;
% h=plot(c_mcau(1:nt,3)/1000,x,'r');set(h,'linewidth',0.5);
% xlabel('Rho(g/cm^3)','fontsize',12);set(gca,'yDir','rev');set(gca,'fontsize',12);
% xlim([0 0.045]);
% ylim([0.2 1])



toc
