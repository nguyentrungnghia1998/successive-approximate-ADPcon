X=zeros(24,200000);
Y=zeros(200000,1);
w=zeros(24,1);
u0=satlin(-5*x2(1,:)-3*x2(2,:));
for k=1:20
for i=1:200000
    xi=x2(:,i);
    dsL=[2*xi(1) 0;
        0 2*xi(2);
        xi(2) xi(1);
        4*xi(1)^3 0;
        0 4*xi(2)^3;
        3*xi(1)^2*xi(2) xi(1)^3;
        2*xi(1)*xi(2)^2 2*xi(1)^2*xi(2);
        xi(2)^3 3*xi(1)*xi(2)^2;
        6*xi(1)^5 0;
        0 6*xi(2)^5;
        5*xi(1)^4*xi(2) xi(1)^5;
        4*xi(1)^3*xi(2)^2 2*xi(1)^4*xi(2);
        3*xi(1)^2*xi(2)^3 3*xi(1)^3*xi(2)^2;
        2*xi(1)*xi(2)^4 4*xi(1)^2*xi(2)^3;
        xi(2)^5 5*xi(1)*xi(2)^4;
        8*xi(1)^7 0;
        7*xi(1)^6*xi(2) xi(1)^7;
        6*xi(1)^5*xi(2)^2 2*xi(1)^6*xi(2);
        5*xi(1)^4*xi(2)^3 3*xi(1)^5*xi(2)^2;
        4*xi(1)^3*xi(2)^4 4*xi(1)^4*xi(2)^3;
        3*xi(1)^2*xi(2)^5 5*xi(1)^3*xi(2)^4;
        2*xi(1)*xi(2)^6 6*xi(1)^2*xi(2)^5;
        xi(2)^7 7*xi(1)*xi(2)^6;
        0 8*xi(2)^7];
    X(:,i)=dsL*[xi(1)+xi(2)-xi(1)*(xi(1)^2+xi(2)^2);-xi(1)+xi(2)-xi(2)*(xi(1)^2+xi(2)^2)+u0(i)];
    if abs(u0(i))==1
        W1=2*log(2);
    else 
        W1=2*u0(i)*atanh(u0(i))+log(1-u0(i)^2);
    end
    Y(i)=xi'*xi+W1;
end
w=-X'\Y;
unew=zeros(1,200000);
for i=1:200000
    xi=x2(:,i);
    dsL=[2*xi(1) 0;
        0 2*xi(2);
        xi(2) xi(1);
        4*xi(1)^3 0;
        0 4*xi(2)^3;
        3*xi(1)^2*xi(2) xi(1)^3;
        2*xi(1)*xi(2)^2 2*xi(1)^2*xi(2);
        xi(2)^3 3*xi(1)*xi(2)^2;
        6*xi(1)^5 0;
        0 6*xi(2)^5;
        5*xi(1)^4*xi(2) xi(1)^5;
        4*xi(1)^3*xi(2)^2 2*xi(1)^4*xi(2);
        3*xi(1)^2*xi(2)^3 3*xi(1)^3*xi(2)^2;
        2*xi(1)*xi(2)^4 4*xi(1)^2*xi(2)^3;
        xi(2)^5 5*xi(1)*xi(2)^4;
        8*xi(1)^7 0;
        7*xi(1)^6*xi(2) xi(1)^7;
        6*xi(1)^5*xi(2)^2 2*xi(1)^6*xi(2);
        5*xi(1)^4*xi(2)^3 3*xi(1)^5*xi(2)^2;
        4*xi(1)^3*xi(2)^4 4*xi(1)^4*xi(2)^3;
        3*xi(1)^2*xi(2)^5 5*xi(1)^3*xi(2)^4;
        2*xi(1)*xi(2)^6 6*xi(1)^2*xi(2)^5;
        xi(2)^7 7*xi(1)*xi(2)^6;
        0 8*xi(2)^7];
    unew(i)=-tanh(1/2*[0 1]*dsL'*w);
end
u0=unew;
end