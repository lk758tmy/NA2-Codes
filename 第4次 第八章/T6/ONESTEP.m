format long;
B1=[1.018011831;1.072672941;1.128231160;1.184632771;1.241823189;1.299747102; 
1.358348467;1.417570569;1.477356158;1.537647359;1.598385848;1.659512885;
1.720969327;1.782695703;1.844632305;1.906719215;1.968896382];
B2=[1.053812884;1.108523287;1.164079316;1.220428251;1.277516671;1.335290395;
1.393694647;1.452674024;1.512172557;1.572133832;1.632500969;1.693216689;
1.754223414;1.815463282;1.876878195;1.938409887; 1.999999993];
for n=100:101
    for j=1:17
        if(n==100)
            q=B1(j);
        else
            q=B2(j);
        end;
        A=zeros(n,n);
        for i=1:n-1
            A(i,i)=2-q;
            A(i,i+1)=-1;
            A(i+1,i)=-1;
        end;
        A(n,n)=2-q;
        v1=rand(n);
        %v1=zeros(n);
        %for i=1:n
        %    v1(i)=1;
        %end;
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
        1/m1+q
    end;
end;