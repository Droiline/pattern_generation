function meinhardt1D

% implements simple numerics for Meinhardt model
% http://www.eb.tuebingen.mpg.de/fileadmin/uploads/pdf/Emeriti/Hans_Meinhardt/Old_paper_pdf/87-shells.pdf

% numerical parameters
    nCols = 200;
    nRows = 5000;

% model parameters for Fig.3 in paper
%     dx = 0.05;
%     dx2 = dx*dx;
%     dt = 0.001;
%     
%     Da = 0.002;
%     Di = 0.4;
%     
%     params.eqn = 1;
%     params.r = 0.01;
%     params.k = 0.0;
%     params.mu = 0.01;
%     params.nu = 0.0;
%     params.s = 0.015;
%     params.r_a = 0.001;
%     params.r_i = 0.0075;
%     params.rho = params.r + rand(1,nCols).*0.025; 
    
    
% model parameters for Fig.5 in paper
    dx = 0.5;
    dx2 = dx*dx;
    dt = 0.1;
    
    Da = 0.1;
    Di = 0.0;
    params.eqn = 2;
    params.r = 0.05;
    params.k = 0.0004;
    params.mu = 0.05; 
    params.nu = 0.03;
    params.r_a = 0.02;
    params.r_i = 0.0075;
    params.rho = params.r + rand(1,nCols).*0.000000005;
    
% model parameters for Fig.4 in paper
%     Da = 0.01;
%     Di = 0.4;
%     params.eqn = 2;
%     params.r = 0.2;
%     params.k = 0.15;
%     params.mu = 0.02;
%     params.nu = 0.02;
%     params.r_a = 0.01;
%     params.r_i = 0.0075;
%     params.rho = params.r + rand(1,nCols).*0.01;   
    
% stability test
    stabilityA = Da * dt/dx2;
    stabilityI = Di * dt/dx2;

    disp(sprintf('Activator stability is %f, should be < 0.5',stabilityA));
    disp(sprintf('Inhibitor stability is %f, should be < 0.5',stabilityI));
    
% initialise model
    activator = zeros(nRows,nCols) + 0.5;
    activator(1,80:120) = 0.0;
    %activator(1,:) = randn(1,nCols)*0.25 + 2;
    
    inhibitor = zeros(nRows,nCols) + 0.5;
    %inhibitor(1,:) = randn(1,nCols)*0.25 + 1;
    
% main loops
    for t = 1:nRows-1
        oldArow = activator(t,:);
        oldIrow = inhibitor(t,:);
        
        for x = 1:nCols
            reactionA = activatorReaction(oldArow(x), oldIrow(x), params, params.rho(x));            
            diffusionA = Da * diffusion( oldArow, x, nCols, dx2 );
            activator(t+1,x) = activator(t,x) + dt * (reactionA + diffusionA);
            
            if reactionA > 10000
                disp('exploding ....');
            end

            reactionI = inhibitorReaction(oldArow(x), oldIrow(x), params, params.rho(x));            
            diffusionI = Di * diffusion( oldIrow, x, nCols, dx2 );
            inhibitor(t+1,x) = inhibitor(t,x) + dt * (reactionI + diffusionI);

        end
    end
    
    figure
    pcolor(activator);shading flat;
    colorbar;
    colormap hot
    
    figure
    pcolor(inhibitor);shading flat
    colorbar
    colormap jet
    
end

function dA = activatorReaction(oldA, oldI, params, rho)

    if params.eqn == 1
    % implement equation 1 from meinhardt paper
        aStar2 = ((oldA*oldA)/(1 + params.k * oldA*oldA)) + params.r_a;
        dA = (rho * oldI * aStar2) - params.mu * oldA;
    
    else
    % implement equation 2 from meinhardt paper
        aStar2 = (oldA*oldA)/(1 + params.k * oldA*oldA);
        dA = (rho * (aStar2 + params.r))/oldI - params.mu * oldA;
    end
    
end

function dI = inhibitorReaction(oldA, oldI, params, rho)

    if params.eqn == 1
    % implement equation 1 from meinhardt paper
        aStar2 = ((oldA*oldA)/(1 + params.k * oldA*oldA)) + params.r_a;
        dI = params.s - rho * oldI * aStar2 - params.nu * oldI;

    else
    % implement equation 2 from meinhardt paper
        aStar2 = (oldA*oldA)/(1 + params.k * oldA*oldA);
        dI = rho * aStar2 - params.nu * oldI + params.r_i;
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



