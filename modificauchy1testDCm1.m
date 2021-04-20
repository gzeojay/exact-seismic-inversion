function Q1=modificauchy1testDCm1(c,mean_c,N,pusai)  %%柯西约束的一阶导数
%lamda=0.01;
%s=norm(pusai);
%inpusai=inv(pusai+lamda*s);
%inpusai=inv(pusai);
unitm=eye(N);      % 单位矩阵
invCm=inv(pusai);         % 模型参数扰动的协方差矩阵的逆
Cm=kron(invCm,unitm);      % 模型参数协方差      % D'*invPusai*D
%invCm=inv(Cm);
inv_model=zeros(5*N,1);
for k=1:N
        inv_model(k,1)=c(k,1)-mean_c(k,1);
        inv_model(N+k,1)=c(k,2)-mean_c(k,2);
        inv_model(2*N+k,1)=c(k,3)-mean_c(k,3);  
        inv_model(3*N+k,1)=c(k,4)-mean_c(k,4);
        inv_model(4*N+k,1)=c(k,5)-mean_c(k,5);
 end

    
    
    ffai=1+inv_model'*Cm*inv_model;     % 1+(m-u)'*Cm*(m-u)
    
    

% Q1=zeros(3*N,1);

% Q1=2*Cm/ffai^2;         %%改进型柯西约束的一阶导数

Q1=4*Cm/ffai;         %%柯西约束的一阶导数

