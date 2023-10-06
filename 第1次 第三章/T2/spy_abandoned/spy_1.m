a=input("a:");
n=input("n:");
A=zeros(n,n);
A(n,n)=1;
for i=1:n
    A(i,i)=1;
    A(n,i)=a;
    A(i,n)=a;
end
[L,U,P]=lu(A);
spy(L+U)