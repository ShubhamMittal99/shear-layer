a1=[0.01,0.01,0.02,0.02,0.04];
a2=[0,0.02,0.02,0.04,0.02];
phi=[0,0,pi/4,pi/2,-pi/4];
c=2*pi;
t=10;
y=0;
x= linspace(0,25,4000);
z=zeros(3,4000);

for i=1:5
    for j=1:4000
        z(i,j)= (4*a1(i)*sech(y)*sech(y)*sech(y)*sin(x(j)-c*t)) + (4*a2(i)*sech(y)*((sech(y)*sech(y)) - 3/8)*sin((x(j)-c*t + phi(i))/2)) + sech(y)*sech(y) ;
        
    end
    figure()
    plot(x,z(i,:))
end 
