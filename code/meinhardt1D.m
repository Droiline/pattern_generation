function meinhardt1D

% implements simple numerics for Meinhardt model
% http://www.eb.tuebingen.mpg.de/fileadmin/uploads/pdf/Emeriti/Hans_Meinhardt/Old_paper_pdf/87-shells.pdf

% numerical parameters
    nCols = 260;
    nRows = 6500;  

% model parameters that give a travelling wave
%     dx = 1;
%     dx2 = dx*dx;
%     dt = 0.5;
%     
%     Da = 0.01;
%     Di = 0.0;
%     params.eqn = 2;
%     params.r = 0.02;
%     params.k = 0.004;
%     params.mu = 0.05;
%     params.nu = 0.03;
%     params.r_a = 0.05;
%     params.r_i = 0.0075;
%     params.rho = zeros(1,nCols);
  
% model parameters that give a branching travelling wave
    dx = 1;
    dx2 = dx*dx;
    dt = 0.5;
    
    Da = 0.015;
    Di = 0.0;
    params.eqn = 3;
    params.r = 0.02;
    params.k = 0.004;
    params.mu = 0.05;
    params.nu = 0.03;
    params.r_a = 0.05;
    params.r_i = 0.0075;
    params.rho = zeros(1,nCols);    
    
% stability test
    stabilityA = Da * dt/dx2;
    stabilityI = Di * dt/dx2;

    disp(sprintf('Activator stability is %f, should be < 0.5',stabilityA));
    disp(sprintf('Inhibitor stability is %f, should be < 0.5',stabilityI));
    
% initialise model
    activator = zeros(nRows,nCols);
    activator(1,125:135) = 1;
    % activator(1,80:120) = 0.0;
    % activator(1,:) = randn(1,nCols)*0.25 + 1;
    for col = 1:nCols
        if rand(1,1) > 0.25
            activator(1,col) = 0.5;
        else
            activator(1,col) = 0;
        end
        % activator(1,col) = (sin(2*pi*col/(nCols/4))+1.0)/2;
    end
    
    inhibitor = zeros(nRows,nCols) + 0.5;
    % inhibitor(1,:) = randn(1,nCols)*0.25 + 1;
    
    R = 0.5;
    
% main loops
    for t = 1:nRows-1
        
        % params.rho = 0.05 + randn(1,nCols).*0.15;   
        
        oldArow = activator(t,:);
        oldIrow = inhibitor(t,:);
        
        R = R + dt * 0.1 * (sum(oldArow)/nCols - R);
        Rt(t) = R;
        
        for x = 1:nCols
            
            %newRho = 0.05 + randn(1,1).*0.15;
            newRho = params.r_a + randn(1,1).*params.r_a/2;

            reactionA = activatorReaction(oldArow(x), oldIrow(x), params, newRho, R);            
            diffusionA = Da * diffusion( oldArow, x, nCols, dx2 );
            activator(t+1,x) = activator(t,x) + dt * (reactionA + diffusionA);
            
            if reactionA > 10000
                disp('exploding ....');
            end
            
            newRho = params.r_i + randn(1,1).*params.r_i/2;

            reactionI = inhibitorReaction(oldArow(x), oldIrow(x), params, newRho, R);            
            diffusionI = Di * diffusion( oldIrow, x, nCols, dx2 );
            inhibitor(t+1,x) = inhibitor(t,x) + dt * (reactionI + diffusionI);

        end
    end
    
    figure
    pcolor(activator);shading flat;
    colorbar;
    colormap hot
    caxis([0 5])
    
    figure
    pcolor(inhibitor);shading flat
    colorbar
    colormap jet
    
end

function dA = activatorReaction(oldA, oldI, params, rho, R)

    if params.eqn == 1
    % implement equation 1 from meinhardt paper
        aStar2 = ((oldA*oldA)/(1 + params.k * oldA*oldA)) + params.r_a;
        dA = (rho * oldI * aStar2) - (params.mu * oldA);
    
    elseif params.eqn == 2
    % implement equation 2 from meinhardt paper
        aStar2 = (oldA*oldA)/(1 + params.k * oldA*oldA);
        dA = ((rho * (aStar2 + params.r))/oldI) - (params.mu * oldA);
    else
    % implement equation 2 from meinhardt paper
        aStar2 = (oldA*oldA)/(1 + params.k * oldA*oldA);
        dA = (rho * (aStar2 + params.r))/(0.1 + oldI) - params.mu * oldA;        
    end
    
    
end

function dI = inhibitorReaction(oldA, oldI, params, rho, R)

    if params.eqn == 1
    % implement equation 1 from meinhardt paper
        aStar2 = ((oldA*oldA)/(1 + params.k * oldA*oldA)) + params.r_a;
        dI = params.s - rho * oldI * aStar2 - params.nu * oldI;

    elseif params.eqn == 2
    % implement equation 2 from meinhardt paper
        aStar2 = (oldA*oldA)/(1 + params.k * oldA*oldA);
        dI = (rho * aStar2) - (params.nu * oldI) + params.r_i;
        
    else
        % implement equation 2 from meinhardt paper
        aStar2 = (oldA*oldA)/(1 + params.k * oldA*oldA);
        dI = rho * aStar2 - params.nu/R * oldI + params.r_i;
        
    end
    
end

function diff = diffusion( row, x, nX, dx2 )
    if (x == 1)
        diff = (row(x) + row(x+1) - 2 * row(x) )/dx2;
    elseif(x < nX)
        diff = (row(x-1) + row(x+1) - 2 * row(x) )/dx2;
    else
        diff = (row(x-1) + row(x) - 2 * row(x) )/dx2;
    end
end



