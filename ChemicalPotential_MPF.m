clf;
u=4;
del=0.05;
x=[-del:0.01:1+del];
y=[-del:0.01:1+del];

f=zeros(length(x),length(y));

for i=1:length(x)
    for j=1:length(y)
        x1=x(i); y1=y(j); z1=1-x(i)-y(j);
        f(i,j)=x1^2*y1^2+x1^2*z1^2+y1^2*z1^2;
        f(i,j)=f(i,j)*u;
        if ( x(i)+y(j) )>(1+del)
            f(i,j)=0;
        end;
    end;
end;
figure(1); 
[xa,ya]=meshgrid(x,y);
surf(x,y,f);
set(gca,'ztick',[]);
title('multi-well');
set(gcf, 'color', 'w');

%--------------------------------------------------------------------------parabolic
u=1;
for i=1:length(x)
    for j=1:length(y)
        a=abs(x(i)); b=abs(y(j)); c=abs(1-x(i)-y(j));
        f(i,j)=a*b+a*c+b*c;
        f(i,j)=f(i,j)*u;
        if ( x(i)+y(j) )>(1+del)
            f(i,j)=0;
        end;
    end;
end;
figure(3);
[xa,ya]=meshgrid(x,y);
surf(x,y,f);
set(gca,'ztick',[]);
set(gcf, 'color', 'w');
title('multi-parabolic');

for i=1:length(x) % penalty of 3-phase coexistence
    for j=1:length(y)
        f(i,j)=10*u*x(i)*y(j)*(1-x(i)-y(j));
        if ( x(i)+y(j) )>1
            f(i,j)=0;  % 0.12
        end;
    end;
end;

figure(4);
surf(x,y,f);
set(gca,'ztick',[]);
set(gcf, 'color', 'w');
title('penalty term');