clearvars, close all

% Computation of the numerical solution with different FD
% schemes

f=@(x)1+0.*x; 
uex=@(x,mu,alfa)-(x + exp((alfa.*x)./mu) - x.*exp(alfa./mu) - 1)./(alfa.*(exp(alfa./mu) - 1));

L=1; % b-a
u0=0; uL=0; % Dirichlet BCs.

hh=[0.1 0.01];
alfaa=[1 -1];
muu=[0.1 0.01];

fprintf('\nFD: Backward Scheme\n');
image=0;
figure()

for i=1:numel(hh)
    for j=1:numel(alfaa)
        for k=1:numel(muu)
            
            h=hh(i);
            alfa=alfaa(j);
            mu=muu(k);
            
            image=image+1;

            n=L/h+1; % total number of nodes. L/h-1 internals + 2 boundaries
            xnodes=linspace(0,L,n); % x-coordinates of grid points

            Uh=FDbackward(alfa,h,mu,n,u0,uL);

            subplot(4,2,image)
            plot(xnodes,Uh,'b-',xnodes,uex(xnodes,mu,alfa),'k--')
            xlabel('x')

            str=sprintf('alfa=%d\nh=%3.2f\nmu=%3.2f',alfa,h,mu);
            err=norm(Uh'-uex(xnodes,mu,alfa),'inf');

            text(xnodes(ceil(numel(xnodes)/2)),Uh(ceil(numel(xnodes)/2)),str)
            legend('Uh(x)','uex(x)')
            fprintf('\nalfa=%d, h=%4.2f, mu=%4.2f\nerr=%12.10e\n',alfa,h,mu,err);
            
        end
    end
end
sgtitle('Finite Difference method: backward scheme')


fprintf('\n\nFD: Centred Scheme\n');
image=0;
figure()

for i=1:numel(hh)
    for j=1:numel(alfaa)
        for k=1:numel(muu)
            
            h=hh(i);
            alfa=alfaa(j);
            mu=muu(k);
            
            image=image+1;

            n=L/h+1; % total number of nodes. L/h-1 internals + 2 boundaries
            xnodes=linspace(0,L,n); % x-coordinates of grid points

            Uh=FDcentr(alfa,h,mu,n,u0,uL);

            subplot(4,2,image)
            plot(xnodes,Uh,'b-',xnodes,uex(xnodes,mu,alfa),'k--')
            xlabel('x')

            str=sprintf('alfa=%d\nh=%3.2f\nmu=%3.2f',alfa,h,mu);
            err=norm(Uh'-uex(xnodes,mu,alfa),'inf');

            text(xnodes(ceil(numel(xnodes)/2)),Uh(ceil(numel(xnodes)/2)),str)
            legend('Uh(x)','uex(x)')
            fprintf('\nalfa=%d, h=%4.2f, mu=%4.2f\nerr=%12.10e\n',alfa,h,mu,err);
            
        end
    end
end
sgtitle('Finite Difference method: centered scheme')

fprintf('\n\nFD: Forward Scheme\n');
image=0;
figure()

for i=1:numel(hh)
    for j=1:numel(alfaa)
        for k=1:numel(muu)
            
            h=hh(i);
            alfa=alfaa(j);
            mu=muu(k);
            
            image=image+1;

            n=L/h+1; % total number of nodes. L/h-1 internals + 2 boundaries
            xnodes=linspace(0,L,n); % x-coordinates of grid points

            Uh=FDforward(alfa,h,mu,n,u0,uL);

            subplot(4,2,image)
            plot(xnodes,Uh,'b-',xnodes,uex(xnodes,mu,alfa),'k--')
            xlabel('x')

            str=sprintf('alfa=%d\nh=%3.2f\nmu=%3.2f',alfa,h,mu);
            err=norm(Uh'-uex(xnodes,mu,alfa),'inf');

            text(xnodes(ceil(numel(xnodes)/2)),Uh(ceil(numel(xnodes)/2)),str)
            legend('Uh(x)','uex(x)')
            fprintf('\nalfa=%d, h=%4.2f, mu=%4.2f\nerr=%12.10e\n',alfa,h,mu,err);
            
        end
    end
end
sgtitle('Finite Difference method: forward scheme')
