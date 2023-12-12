format long;
for n=100:101
    for q=2:3
        A=zeros(n,n);
        for i=1:n-1
            A(i,i)=2-q+0.001;
            A(i,i+1)=-1;
            A(i+1,i)=-1;
        end;
        A(n,n)=2-q+0.001;
%         v1=rand(n);
        v1=zeros(n);
        for i=1:n
            v1(i)=1;
        end;
        v2=zeros(n);

        cnt=0;
        m1=1;
        m2=0;
        while abs(1/m1-1/m2)>1e-8
            m2=m1;
            v2=A\v1;
            m1=norm(v2,"inf");
            v1=v2/m1;
            cnt=cnt+1;
        end;
        cnt
        1/m1+q-0.001
    end;
end;