X=zeros(21,200000);
Y=zeros(200000,1);
w=zeros(21,1);
u0=-K*x;
for i=1:200000
    if abs(u0(1,i))>3
        u0(1,i)=3*sign(u0(1,i));
    end
    if abs(u0(2,i))>20
        u0(2,i)=20*sign(u0(2,i));
    end
end
for k=1:20
for i=1:200000
    xi=x(:,i);
    dsL=[2*xi(1) 0 0;
        0 2*xi(2) 0;
        0 0 2*xi(3);
        xi(2) xi(1) 0;
        xi(3) 0 xi(1);
        0 xi(3) xi(2);
        4*xi(1)^3 0 0;
        0 4*xi(2)^3 0;
        0 0 4*xi(3)^3;
        2*xi(1)*xi(2)^2 2*xi(1)^2*xi(2) 0;
        2*xi(1)*xi(3)^2 0 2*xi(3)*xi(1)^2;
        0 2*xi(2)*xi(3)^2 2*xi(3)*xi(2)^2;
        2*xi(1)*xi(2)*xi(3) xi(1)^2*xi(3) xi(1)^2*xi(2);
        xi(2)^2*xi(3) 2*xi(1)*xi(2)*xi(3) xi(2)^2*xi(1);
        xi(3)^2*xi(2) xi(3)^2*xi(1) 2*xi(1)*xi(2)*xi(3);
        3*xi(1)^2*xi(2) xi(1)^3 0;
        3*xi(1)^2*xi(3) 0 xi(1)^3;
        xi(2)^3 3*xi(2)^2*xi(1) 0;
        xi(3)^3 0 3*xi(3)^2*xi(1);
        0 xi(3)^3 3*xi(3)^2*xi(2);
        0 3*xi(2)^2*xi(3) xi(2)^3];
    X(:,i)=dsL*[2*xi(1)+xi(2)+xi(3);xi(1)-xi(2)+u0(2,i);xi(3)+u0(1,i)];
    if abs(u0(1,i))==3
        W1=18*log(6)-9*log(9);
    else 
        W1=6*u0(1,i)*atanh(u0(1,i)/3)+9*log(1-u0(1,i)^2/9);
    end
    if abs(u0(2,i))==20
        W2=800*log(40)-400*log(400);
    else 
        W2=40*u0(2,i)*atanh(u0(2,i)/20)+400*log(1-u0(2,i)^2/400);
    end
    Y(i)=xi'*xi+W1+W2;
end
w=-X'\Y;
unew=zeros(2,200000);
for i=1:200000
    xi=x(:,i);
    dsL=[2*xi(1) 0 0;
        0 2*xi(2) 0;
        0 0 2*xi(3);
        xi(2) xi(1) 0;
        xi(3) 0 xi(1);
        0 xi(3) xi(2);
        4*xi(1)^3 0 0;
        0 4*xi(2)^3 0;
        0 0 4*xi(3)^3;
        2*xi(1)*xi(2)^2 2*xi(1)^2*xi(2) 0;
        2*xi(1)*xi(3)^2 0 2*xi(3)*xi(1)^2;
        0 2*xi(2)*xi(3)^2 2*xi(3)*xi(2)^2;
        2*xi(1)*xi(2)*xi(3) xi(1)^2*xi(3) xi(1)^2*xi(2);
        xi(2)^2*xi(3) 2*xi(1)*xi(2)*xi(3) xi(2)^2*xi(1);
        xi(3)^2*xi(2) xi(3)^2*xi(1) 2*xi(1)*xi(2)*xi(3);
        3*xi(1)^2*xi(2) xi(1)^3 0;
        3*xi(1)^2*xi(3) 0 xi(1)^3;
        xi(2)^3 3*xi(2)^2*xi(1) 0;
        xi(3)^3 0 3*xi(3)^2*xi(1);
        0 xi(3)^3 3*xi(3)^2*xi(2);
        0 3*xi(2)^2*xi(3) xi(2)^3];
    unew(1,i)=-3*tanh(1/2*[0 0 1]*dsL'*w/3);
    unew(2,i)=-20*tanh(1/2*[0 1 0]*dsL'*w/20);
end
u0=unew;
end