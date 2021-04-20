function r=Rp1(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle)
theta=pi*angle/180;
Ap1=zeros(4,4);
Bp1=zeros(4,1);
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


% -------- A对上层纵波速度vp1求偏导
H12_p1=(-2*vp1*vs1^2)/(vp1^2+F1*Por1/rho1)^2;
H13_p1=-(vp1^2+F1*Por1/rho1)^(-3/2)*vp1*(vp2^2+F2*Por2/rho2)^(1/2);
H14_p1=(-2*vp1*vs2^2)/(vp1^2+F1*Por1/rho1)^2;

H22_p1=-(vp1^2+F1*Por1/rho1)^(-1/2)*vp1*vs1/(vp1^2+F1*Por1/rho1);
H23_p1=-2*(vp2^2+F2*Por2/rho2)*vp1/(vp1^2+F1*Por1/rho1)^2;
H24_p1=-(vp1^2+F1*Por1/rho1)^(-1/2)*vp1*vs2/(vp1^2+F1*Por1/rho1);

H32a_p1=(vp1^2+F1*Por1/rho1)^(-1/2)*vp1/vs1;
H32b_p1=H22_p1;
H33b_p1=H23_p1;
H34a_p1=(vp1^2+F1*Por1/rho1)^(-1/2)*vp1*vs2/vs1^2;
H34b_p1=H14_p1;

H41_p1=H12_p1;
H42_p1=H12_p1;
H43a_p1=H13_p1;
H43b_p1=H14_p1;
H44_p1=H14_p1;
%
G41_p1=H12_p1;
%% Ap1(1,:)
Ap1(1,1)=0;
Ap1(1,2)=(1/2)*(1-H12*sin(theta)^2)^(-1/2)*sin(theta)^2*H12_p1;
Ap1(1,3)=H13_p1*sin(theta);
Ap1(1,4)=(1/2)*(1-H14*sin(theta)^2)^(-1/2)*sin(theta)^2*H14_p1;
%% Ap1(2,:)
Ap1(2,1)=0;
Ap1(2,2)=-H22_p1*sin(theta);
Ap1(2,3)=-(1/2)*(1-H23*sin(theta)^2)^(-1/2)*sin(theta)^2*H23_p1;
Ap1(2,4)=H24_p1*sin(theta);
%% Ap1(3,:)
Ap1(3,1)=0;
Ap1(3,2)=H32a_p1-2*sin(theta)^2*H32b_p1;
Ap1(3,3)=-(rho2/rho1)*H33a*(sin(theta)^2-H33b*sin(theta)^4)^(-1/2)*sin(theta)^4*H33b_p1;
Ap1(3,4)=-(rho2/rho1)*H34a_p1*(1-2*H34b*sin(theta)^2)+2*(rho2/rho1)*sin(theta)^2*H34b_p1*H34a;
%% Ap1(4,:)
Ap1(4,1)=2*sin(theta)^2*H41_p1;
Ap1(4,2)=2*H42_p1*(sin(theta)^2-H42*sin(theta)^4)^(1/2)-H42*(sin(theta)^2-H42*sin(theta)^4)^(-1/2)*sin(theta)^4*H42_p1;
Ap1(4,3)=(rho2/rho1)*H43a_p1*(1-2*H43b*sin(theta)^2)-2*(rho2/rho1)*sin(theta)^2*H43b_p1*H43a;
Ap1(4,4)=2*(rho2/rho1)*H44_p1*(sin(theta)^2-H44*sin(theta)^4)^(1/2)-(rho2/rho1)*H44*(sin(theta)^2-H44*sin(theta)^4)^(-1/2)*sin(theta)^4*H44_p1;

% -------- B对上层纵波速度vp1求偏导
%%
Bp1(1,1)=0;
Bp1(2,1)=0;
Bp1(3,1)=0;
Bp1(4,1)=-2*sin(theta)^2*G41_p1;

%%

r=(inv(A))*(-Ap1*(inv(A))*B+Bp1);    % 关于纵波速度vp1反射系数的一阶偏导数
end
