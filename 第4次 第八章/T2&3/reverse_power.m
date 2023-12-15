format long;
for n=100:101
    for q=2:3
        u=zeros(n);
        if(q==2)
            if(n==100)
                l=2-2*cos(51*pi/101);
                for j=1:n
                    u(j)=sin(j*51*pi/101);
                end
            else
                l=2;
                for j=1:n
                    u(j)=sin(j*pi/2);
                end
            end
        else
            if(n==100)
                l=2-2*cos(67*pi/101);
                for j=1:n
                    u(j)=sin(j*67*pi/101);
                end
            else
                l=3;
                for j=1:n
                    u(j)=sin(j*pi*2/3);
                end
            end
        end

        A=zeros(n,n);
        for i=1:n-1
            A(i,i)=2-q+0.001;
            A(i,i+1)=-1;
            A(i+1,i)=-1;
        end
        A(n,n)=2-q+0.001;
%         v1=rand(n);
        v1=zeros(n);
        for i=1:n
            v1(i)=1;
        end
        v2=zeros(n);

        cnt=0;
        m1=1;
        m2=0;
        while abs(1/m1-1/m2)>1e-8
            e=dot(u,v1);
            e=e(1)/(norm(v1)*norm(u));
            e=sqrt(1-e*e);
            ll=1/m1+q-0.001;
            b=[cnt;abs(ll-l);e]

            m2=m1;
            v2=A\v1;
            m1=norm(v2,"inf");
            v1=v2/m1;
            cnt=cnt+1;
        end
        e=dot(u,v1);
        e=e(1)/norm(v1)/norm(u);
        e=sqrt(1-e*e);
        ll=1/m1+q-0.001;
        b=[cnt;abs(ll-l);e]
    end
end