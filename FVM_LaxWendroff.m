clear
close all
clc

L=15;
N=500; %% points in space
dx=2*L/N; %% step size in space
x=-L:dx:L; %% domain
x=x';

c=7; %% Wave speed

dt=1/600; %% Step size in time
T=1; %% Final time
t=0:dt:T;
Nt=numel(t);

alpha=c*dt/dx; %% CFL


%% Initial condition

u0=zeros(1,numel(x));
for(jj=1:numel(x)-1)
    if(x(jj)<-L/5)
        u0(jj)=2;
    elseif(x(jj)<=L/5 && x(jj)>=-L/5)
        u0(jj)=(3-5*x(jj)/L)/2;
    elseif(x(jj)>L/5)
        u0(jj)=1;
    end
end

%% Analytical solution

ue=u0';
for(jj=1:numel(t))
    xt= x + c*ue*dt;
    xt= [x(1)-dx; xt ; x(end)+dx];
    ut=[0;ue;0];
    ue= interp1(xt, ut, x);
end

%% Lax wendroff method


u=u0;
for(i=2:N+1)
    K(i,i)=-1/2;
    K(i-1,i)=1/2;

    KP(i,i)=1/2;
    KP(i,i-1)=-1/2;
end
K(1,1)=-0.5;
KP(1,1)=0.5;

for(i=1:Nt)
    if(c<0)
        unp=u-((alpha*(K*(u.*u)'))+(alpha^2*(K*(u.*(K*(u.*u)')')')))';
    else
        unp=u-((alpha*(KP*(u.*u)'))+(alpha^2*(KP*(u.*(KP*(u.*u)')')')))';
    end

    u=unp;

end

 error_lax=sum((u-ue')/ue')^2;
 
plot(x,unp,x,u0,x,ue,"LineWidth",2)
legend("U_{n+1}","U_{0}","Exact solution")
xlabel("Distance in m")
ylabel("Amplitude")
fontsize(gca,15,'points')