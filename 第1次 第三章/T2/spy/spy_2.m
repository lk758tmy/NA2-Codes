a=input("a:");
n=input("n:");
A=zeros(n,n);
A(1,1)=1;
for i=2:n
    A(i,i)=1;
    A(1,i)=a;
    A(i,1)=a;
end
[L,U,P]=lu(A);
spy(L+U)