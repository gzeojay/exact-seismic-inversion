function Q1=modificauchy1testDCm1(c,mean_c,N,pusai)  %%����Լ����һ�׵���
%lamda=0.01;
%s=norm(pusai);
%inpusai=inv(pusai+lamda*s);
%inpusai=inv(pusai);
unitm=eye(N);      % ��λ����
invCm=inv(pusai);         % ģ�Ͳ����Ŷ���Э����������
Cm=kron(invCm,unitm);      % ģ�Ͳ���Э����      % D'*invPusai*D
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

% Q1=2*Cm/ffai^2;         %%�Ľ��Ϳ���Լ����һ�׵���

Q1=4*Cm/ffai;         %%����Լ����һ�׵���

