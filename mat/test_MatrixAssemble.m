function err=test_MatrixAssemble(N)
h=1/N;
tmpM1=diag(-2*ones(N-1,1))+diag(ones(N-2,1),1)+diag(ones(N-2,1),-1);
tmpM2=eye(N);
Au=kron(tmpM2,tmpM1);
tmpM1=diag(-2*ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
%Boundary Condition
% tmpM1(1,1:3)=[-5 2 -1/5];
% tmpM1(end,end:-1:end-2)=[-5 2 -1/5];
tmpM1(1,1)=-3;
tmpM1(end,end)=-3;
tmpM2=eye(N-1);
Au=Au+kron(tmpM1,tmpM2);
Au=-1/h/h*Au;
tmpM1=zeros(N-1,N); tmpM1(1:N:end)=-1;tmpM1(N:N:end)=1;
tmpM2=eye(N);
Bu=kron(tmpM2,tmpM1);
Bu=1/h*Bu;

tmpM1=diag(-2*ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
%Boundary Condition
% tmpM1(1,1:3)=[-5 2 -1/5];
% tmpM1(end,end:-1:end-2)=[-5 2 -1/5];
tmpM1(1,1)=-3;
tmpM1(end,end)=-3;
tmpM2=eye(N-1);
Av=kron(tmpM2,tmpM1);
tmpM1=diag(-2*ones(N-1,1))+diag(ones(N-2,1),1)+diag(ones(N-2,1),-1);
tmpM2=eye(N);
Av=Av+kron(tmpM1,tmpM2);
Av=-1/h/h*Av;
tmpM1=zeros(N-1,N); tmpM1(1:N:end)=-1;tmpM1(N:N:end)=1;
tmpM2=eye(N);
Bv=kron(tmpM1,tmpM2);
Bv=1/h*Bv;

A=blkdiag(Au,Av);
B=[Bu;Bv];
M=[A,B;-B',zeros(N*N)];

% test with true solution
func_f=@(x,y)-4*pi^2*(2*cos(2*pi*x)-1)*sin(2*pi*y)+x^2;
func_g=@(x,y)4*pi^2*(2*cos(2*pi*y)-1)*sin(2*pi*x);
func_u=@(x,y)(1-cos(2*pi*x))*sin(2*pi*y);
func_v=@(x,y)-(1-cos(2*pi*y))*sin(2*pi*x);
func_p=@(x,y)x^3/3-1/12;

u=zeros(N*N-N,1);
f=zeros(N*N-N,1);
for i=1:N-1
    for j=1:N
        u(i+(j-1)*(N-1))=func_u(i*h,(j-0.5)*h);
        f(i+(j-1)*(N-1))=func_f(i*h,(j-0.5)*h);
    end
end
v=zeros(N*N-N,1);
g=zeros(N*N-N,1);
for i=1:N
    for j=1:N-1
        v(i+(j-1)*N)=func_v((i-0.5)*h,j*h);
        g(i+(j-1)*N)=func_g((i-0.5)*h,j*h);
    end
end
p=zeros(N*N,1);
for i=1:N
    for j=1:N
        p(i+(j-1)*N)=func_p((i-0.5)*h,(j-0.5)*h);
    end
end

b=[f;g;zeros(N*N,1)];
x=[u;v;p];
err=max(abs(M*x-b));

% upU=zeros(N*N-N,N*N);
% for j=1:N
%     org=[(N-1),N]*(j-1);
%     if j==1 || j==N
%         upU(1+org(1),1+org(2))=1/2; upU(N-1+org(1),N+org(2))=-1/2;
%     else
%         upU(1+org(1),1+org(2))=1/3; upU(N-1+org(1),N+org(2))=-1/3;
%     end
%     for i=2:N-1
%        if j==1 || j==N
%             upU(i+org(1),i+org(2))=1/3; upU(i-1+org(1),i+org(2))=-1/3;
%        else
%             upU(i+org(1),i+org(2))=1/4; upU(i-1+org(1),i+org(2))=-1/4;
%        end
%     end
% end
% upV=zeros(N*N-N,N*N);
% for j=1:N
%     org=[N,N]*(j-1);
%     if j==1
%         upV(org(1)+1,org(2)+1)=1/2; upV(org(1)+N,org(2)+N)=1/2;
%         for i=2:N-1
%             upV(org(1)+i,org(2)+i)=1/3;
%         end
%     elseif j==N
%         upV(org(1)+1-N,org(2)+1)=-1/2; upV(org(1),org(2)+N)=-1/2;
%         for i=2:N-1
%             upV(org(1)+i-N,org(2)+i)=-1/3;
%         end
%     else
%         upV(org(1)+1,org(2)+1)=1/3; upV(org(1)+N,org(2)+N)=1/3;
%         upV(org(1)+1-N,org(2)+1)=-1/3; upV(org(1),org(2)+N)=-1/3;
%         for i=2:N-1
%             upV(org(1)+i,org(2)+i)=1/4; upV(org(1)+i-N,org(2)+i)=-1/4;
%         end
%     end
% end
% u=rand(N*N-N,1);
% v=rand(N*N-N,1);
% r=B'*[u;v];
% norm(r)
% u=u+h*upU*r;
% v=v+h*upV*r;
% r=-B'*[u;v];
% norm(r);



end