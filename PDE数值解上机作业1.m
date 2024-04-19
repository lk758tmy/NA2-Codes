format shortG
for J=[18,36,72,144,288]
%for J=286
    d_x=2*pi/J;    d_t=0.4*d_x*d_x;
    x_array=transpose((-pi+d_x):d_x:pi);
    N=ceil(1/d_t);  T=N*d_t;
    ERR=0.0000001/sqrt(d_x);
    %Advancing Matrixes
    B0=diag(ones(1,J)*0.2)+diag(ones(1,J-1)*0.4,1)+diag(ones(1,J-1)*0.4,-1);
    B0(1,J)=0.4; B0(J,1)=0.4;
    B1=diag(ones(1,J)*1.8)-diag(ones(1,J-1)*0.4,1)-diag(ones(1,J-1)*0.4,-1);
    B1(1,J)=-0.4; B1(J,1)=-0.4;
    B0_2=diag(ones(1,J)*0.6)+diag(ones(1,J-1)*0.2,1)+diag(ones(1,J-1)*0.2,-1);
    B0_2(1,J)=0.2; B0_2(J,1)=0.2;
    B1_2=diag(ones(1,J)*1.4)-diag(ones(1,J-1)*0.2,1)-diag(ones(1,J-1)*0.2,-1);
    B1_2(1,J)=-0.2; B1_2(J,1)=-0.2;
    for m=1:2
        if m==1     %Initial Value 1 + Real Solution
            u0=zeros(J,1);
            for j=1:J
                if(j>=J/4 && j<=3*J/4)
                    u0(j)=1;
                end            
            end
            u1=0.5*ones(J,1);
            u_tmp=ones(J,1);
            i=0;
            while norm(u_tmp)>ERR
                %让级数的相邻误差在10^-7量级，比pde数值解的误差小至少2个量级
                i=i+1; nn=2*i-1;
                u_tmp=2*(2*mod(i,2)-1)*cos(nn*x_array)*exp(-nn*nn*T)/nn/pi;
                u1=u1+u_tmp;
            end
        else        %Initial Value 2 + Real Solution
            u0=pi*ones(J,1)-abs(x_array);
            u1=pi*ones(J,1)/2;
            u_tmp=ones(J,1);
            nn=-1;
            while norm(u_tmp)>ERR
                nn=nn+2;
                u_tmp=4*cos(nn*x_array)*exp(-nn*nn*T)/pi/nn/nn;
                u1=u1+u_tmp;
            end
        end
        for k=1:3
            u=u0;
            %m=3;  u=cos(x_array);  u1=exp(-T)*cos(x_array); %光滑解
            if k==1       %Fully Explicit Scheme
                for n=1:N
%                     if n==N    %改变最后一层的时间步长，注意真解也需要相应变化
%                         miu=(1-(N-1)*d_t)/d_x/d_x;
%                         B0=diag(ones(1,J)*(1-2*miu))+diag(ones(1,J-1)*miu,1)+diag(ones(1,J-1)*miu,-1);
%                         B0(1,J)=miu; B0(J,1)=miu;
%                     end
                    u=B0*u;
                end
            elseif k==2   %Fully Implicit Scheme
                for n=1:N
                    u=B1\u;
                end
            else          %Crank-Nicolson Scheme
                for n=1:N
                    u=B1_2\(B0_2*u);
                end
            end
            [J,m,k,norm(u-u1)*sqrt(d_x),norm(u-u1,"inf")]
            %m=1 J=18时的误差反而小，因为4不整除18，避过了初值的间断点
        end
    end
end
