function r=Rp2(vp1,vs1,rho1,Por1,F1,vp2,vs2,rho2,Por2,F2,angle)
theta=pi*angle/180;
Ap2=zeros(4,4);
Bp2=zeros(4,1);
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



% -------- A对下层纵波速度vp2求偏导

H13_p2=(vp2^2+F2*Por2/rho2)^(-1/2)*vp2/(vp1^2+F1*Por1/rho1)^(1/2);
H23_p2=2*vp2/(vp1^2+F1*Por1/rho1);
H33b_p2=H23_p2;
H43a_p2=H13_p2;

%% Ap2(1,:)
Ap2(1,1)=0;
Ap2(1,2)=0;
Ap2(1,3)=H13_p2*sin(theta);
Ap2(1,4)=0;
%% Ap2(2,:)
Ap2(2,1)=0;
Ap2(2,2)=0;
Ap2(2,3)=-(1/2)*(1-H23*sin(theta)^2)^(-1/2)*sin(theta)^2*H23_p2;
Ap2(2,4)=0;
%% Ap2(3,:)
Ap2(3,1)=0;
Ap2(3,2)=0;
Ap2(3,3)=-(rho2/rho1)*H33a*(sin(theta)^2-H33b*sin(theta)^4)^(-1/2)*sin(theta)^4*H33b_p2;
Ap2(3,4)=0;
%% Ap2(4,:)
Ap2(4,1)=0;
Ap2(4,2)=0;
Ap2(4,3)=(rho2/rho1)*H43a_p2*(1-2*H43b*sin(theta)^2);
Ap2(4,4)=0;

% -------- B对下层纵波速度vp2求偏导
%%
Bp2(1,1)=0;
Bp2(2,1)=0;
Bp2(3,1)=0;
Bp2(4,1)=0;

%%

r=(inv(A))*(-Ap2*(inv(A))*B+Bp2);    % 关于纵波速度vp2反射系数的一阶偏导数
end
