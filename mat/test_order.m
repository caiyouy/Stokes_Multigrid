N=[10,20,40];
Err=ones(3,1);
for i=1:3
    Err(i)=test_MatrixAssemble(N(i));
end
loglog(N,Err);
hold on;
loglog(N,1./N.^2);
legend("Err Inf norm","Standard Order 2")