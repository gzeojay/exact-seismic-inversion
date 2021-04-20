function r=Rf1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle)
theta=pi*angle/180;
Af1=zeros(4,4);
Bf1=zeros(4,1);
A=zeros(4,4);
B=zeros(4,1);
%%
H12=(vs1^2)/(vp1^2+F1*Por1/rho1);
H13=(vp2^2+F2*Por2/rho2)^(1/2)/(vp1^2+F1*Por1/rho1)^(1/2);
H14=(vs2^2)/(vp1^2+F1*Por1/rho1);
H22=vs1/(vp1^2+F1*Por1/rho1)^(1/2);
H23=(vp2^2+F2*Por2/rho2)/(vp1^2+F1*Por1/rho1);
H24=vs2/(vp1^2+F1*Por1/rho1)^(1/2);
H32a=(vp1^2+F1*Por1/rho1)^(1/2)/vs1;
H32b=H22;
H33a=vs2^2/vs1^2;
H33b=H23;
H34a=vs2*(vp1^2+F1*Por1/rho1)^(1/2)/vs1^2;
H34b=H14;
H41=H12;
H42=H12;
H43a=H13;
H43b=H14;
H44=H14;
%
G41=H12;
%% A(1,:)
A(1,1)=-sin(theta);
A(1,2)=-sqrt(1-H12*sin(theta)^2);
A(1,3)=H13*sin(theta);
A(1,4)=-sqrt(1-H14*sin(theta)^2);
%% A(2,:)
A(2,1)=cos(theta);
A(2,2)=-sin(theta)*H22;
A(2,3)=sqrt(1-H23*sin(theta)^2);
A(2,4)=sin(theta)*H24;
%% A(3,:)
A(3,1)=sin(2*theta);
A(3,2)=H32a-2*sin(theta)^2*H32b;
A(3,3)=2*(rho2/rho1)*H33a*sqrt(sin(theta)^2-H33b*sin(theta)^4);
A(3,4)=-(rho2/rho1)*H34a*(1-2*H34b*sin(theta)^2);
%% A(4,:)
A(4,1)=2*H41*sin(theta)^2-1;
A(4,2)=2*H42*sqrt(sin(theta)^2-H42*sin(theta)^4);
A(4,3)=(rho2/rho1)*H43a*(1-2*H43b*sin(theta)^2);
A(4,4)=2*(rho2/rho1)*H44*sqrt(sin(theta)^2-H44*sin(theta)^4);
%% B(:,1)
B(1,1)=sin(theta);
B(2,1)=cos(theta);
B(3,1)=sin(2*theta);
B(4,1)=1-2*G41*sin(theta)^2;



% -------- A对上层流体F1求偏导
H12_f1=-(Por1/rho1)*(vs1^2)/(vp1^2+F1*Por1/rho1)^2;
H13_f1=-(1/2)*(vp1^2+F1*Por1/rho1)^(-3/2)*(vp2^2+F2*Por2/rho2)^(1/2)*(Por1/rho1);
H14_f1=-(Por1/rho1)*(vs2^2)/(vp1^2+F1*Por1/rho1)^2;

H22_f1=-(1/2)*(vp1^2+F1*Por1/rho1)^(-3/2)*vs1*(Por1/rho1);
H23_f1=-(Por1/rho1)*(vp2^2+F2*Por2/rho2)/(vp1^2+F1*Por1/rho1)^2;
H24_f1=-(1/2)*(vp1^2+F1*Por1/rho1)^(-3/2)*vs2*(Por1/rho1);

H32a_f1=(1/2)*(Por1/rho1)*(vp1^2+F1*Por1/rho1)^(-1/2)/vs1;
H32b_f1=H22_f1;
H33b_f1=H23_f1;
H34a_f1=(1/2)*vs2*(Por1/rho1)*(vp1^2+F1*Por1/rho1)^(-1/2)/vs1^2;
H34b_f1=H14_f1;

H41_f1=H12_f1;
H42_f1=H12_f1;
H43a_f1=H13_f1;
H43b_f1=H14_f1;
H44_f1=H14_f1;
%
G41_f1=H12_f1;

%% Af1(1,:)
Af1(1,1)=0;
Af1(1,2)=(1/2)*(1-H12*sin(theta)^2)^(-1/2)*sin(theta)^2*H12_f1;
Af1(1,3)=H13_f1*sin(theta);
Af1(1,4)=(1/2)*(1-H14*sin(theta)^2)^(-1/2)*sin(theta)^2*H14_f1;

%% Af1(2,:)
Af1(2,1)=0;
Af1(2,2)=-H22_f1*sin(theta);
Af1(2,3)=-(1/2)*(1-H23*sin(theta)^2)^(-1/2)*sin(theta)^2*H23_f1;
Af1(2,4)=H24_f1*sin(theta);

%% Af1(3,:)
Af1(3,1)=0;
Af1(3,2)=H32a_f1-2*sin(theta)^2*H32b_f1;
Af1(3,3)=-(rho2/rho1)*H33a*(sin(theta)^2-H33b*sin(theta)^4)^(-1/2)*sin(theta)^4*H33b_f1;
Af1(3,4)=-(rho2/rho1)*H34a_f1*(1-2*H34b*sin(theta)^2)+2*(rho2/rho1)*sin(theta)^2*H34b_f1*H34a;

%% Af1(4,:)
Af1(4,1)=2*sin(theta)^2*H41_f1;
Af1(4,2)=2*H42_f1*(sin(theta)^2-H42*sin(theta)^4)^(1/2)-H42*(sin(theta)^2-H42*sin(theta)^4)^(-1/2)*sin(theta)^4*H42_f1;
Af1(4,3)=(rho2/rho1)*H43a_f1*(1-2*H43b*sin(theta)^2)-2*(rho2/rho1)*sin(theta)^2*H43b_f1*H43a;
Af1(4,4)=2*(rho2/rho1)*H44_f1*(sin(theta)^2-H44*sin(theta)^4)^(1/2)-(rho2/rho1)*H44*(sin(theta)^2-H44*sin(theta)^4)^(-1/2)*sin(theta)^4*H44_f1;

% -------- B对上层流体F1求偏导
%%
Bf1(1,1)=0;
Bf1(2,1)=0;
Bf1(3,1)=0;
Bf1(4,1)=-2*sin(theta)^2*G41_f1;

%% 

r=(inv(A))*(-Af1*(inv(A))*B+Bf1);    % 关于纵波速度vs1反射系数的一阶偏导数
end