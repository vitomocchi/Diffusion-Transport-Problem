clearvars, close all

L=1; % b-a
ro0=0; roL=0; % Dirichlet BCs.

hh=[0.1 0.01];
alfaa=[1 -1];
muu=[0.1 0.01];
image=0;

f=@(x)1+0.*x; 
a=@(x,mu,alfa)mu.*exp(alfa.*x./mu);

roex=@(x,mu,alfa)-(exp(-(alfa.*x)./mu).*(x + exp((alfa.*x)./mu) - x.*exp(alfa./mu)...
    - 1))./(alfa.*(exp(alfa./mu) - 1));

uex=@(x,mu,alfa)-(exp(-(alfa.*x)./mu).*(x + exp((alfa.*x)./mu) - x.*exp(alfa./mu)...
    - 1))./(alfa.*(exp(alfa./mu) - 1)).*exp(alfa*x/mu); % uex = roex * exp(alfa*x/mu)

for i=1:numel(hh)
    for j=1:numel(alfaa)
        for k=1:numel(muu)
            
            h=hh(i);
            alfa=alfaa(j);
            mu=muu(k);
            
            image=image+1;

            n=L/h-1; % internal nodes
            xnodes=linspace(0,L,n+2); 
            
            %%%%%%% ALGORITHM FOR THE DISCRETIZAITON %%%%%%%

            diag = a((xnodes(2:end-1)-1/2*h),mu,alfa)+a((xnodes(2:end-1)+1/2*h),mu,alfa);
            offdiaginf = [-a((xnodes(3:end-1)-1/2*h),mu,alfa), 0]; % spdiags ignores last value.
            offdiagsup = [0, -a((xnodes(2:end-2)+1/2*h),mu,alfa)]; % spdiags ignores first value.

            % Create the sparse matrix using spdiags().
            A = spdiags([offdiaginf', diag', offdiagsup'],-1:1, n, n);

            A=A/h^2;

            b=f(xnodes(2:end-1)); b=b'; % known term (column).
            b(1)=b(1)+ro0*a((xnodes(2)-1/2*h),mu,alfa)/h^2;  % modifies with Dirichlet BC in 0.
            b(end)=b(end)+roL*a((xnodes(end-1)+1/2*h),mu,alfa)/h^2; % modifies with Dirichlet BC in 1.

            roh=A\b; % dim Uh: (nx1)
            roh=[ro0;roh;roL];
            
            %%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%
    
            figure()
            str=sprintf('alfa=%d h=%3.2f mu=%3.2f',alfa,h,mu);
            
            % plotting of roh e roex
            
            subplot(1,2,1)
            plot(xnodes,roh,'b*',xnodes,roex(xnodes,mu,alfa), 'k-')
            xlabel('x')

            legend('roh(x)','roex(x)')
            
            % plotting of Uh e uex using the change of variables

            Uh=roh.*exp((alfa/mu).*xnodes)';
            
            subplot(1,2,2)
            plot(xnodes,Uh,'r*',xnodes,uex(xnodes,mu,alfa), 'g-')
            xlabel('x')

            err=norm(Uh'-uex(xnodes,mu,alfa),'inf');
            
            legend('Uh(x)','uex(x)')
            
            fprintf('\nalfa=%d, h=%4.2f, mu=%4.2f\n',alfa,h,mu);
            fprintf('err(u)=%6.4e\n',err);
            
            sgtitle(str)
            
            clear roh
            clear Uh
        end
    end
end
