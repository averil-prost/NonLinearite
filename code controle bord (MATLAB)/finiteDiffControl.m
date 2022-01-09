%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Contrôle par différences finies, deuxième version
%   Mémoire MFA - 2021-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
clf
global xx tt nu EX

fprintf("//////////////////////////////////////////////////////\n");
fprintf("//////////////////////////////////////////////////////\n");

%{
    Convention de stockage : 
    ------> x, i=1,..,N
    |
    |
    v t, j=1,..,M
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paramètres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maillages
l = pi/4; % Omega = [0,l]
N = 20; % points en espace
T = 1.0; % temps final
M = 60; % points en temps

% pondération du laplacien
nu = 1;

% nombre d'itérations maximal de la minimisation
Nit_max = 50;

% choix de l'emplacement des contrôles (logical)
Gin1 = 0; % sur la composante {x=0}x[0,T]
GinN = 1; % sur la composante {x=l}x[0,T]

% pour les fonctions, voir la classe associée.
auto = 1;
stencil = 1;

% EX = "test";
% EX = "LIN";
EX = "TRO";
% EX = "DISC";

affichage = false;
comparaison = true;

MNlist = round(20*1.4.^[1:6]);
% MNlist = [60];
erreurs = zeros(length(MNlist),4);

for imn=1:length(MNlist)
    N=MNlist(imn);
    M=N;
    erreurs(imn,1) = 1/N;
    fprintf("//////////////////////////////////////////////////////\nN=M=%d\n",M);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = T/(M-1); % bords inclus
    h = l/(N-1); % bords inclus
    xx = linspace(0,l,N);
    tt = linspace(0,T,M);

    % variables
    %{ 
    Choix d'indices : 
        - pour l'état, u_{1,1},u_{2,1},...,u{N,1},u_{1,2},... en espace d'abord
               u_0 est dans UnetGn
               PUIS les contrôles : si uniquement l'un des deux bords,
               g_2,...,g_m. Sinon, g^1_2,g^2_2,g^1_3,g^2_3,... (pas de contrôle initial)
               avec g^1 composante x=0, g^2 composante x=l.
        - pour l'adjoint, même ordre que u 
    %}
    UnetGn = 0.5*ones(M*N+(Gin1+GinN)*(M-1),Nit_max+1); % état et contrôles
    Wn = 0.5*ones(M*N,Nit_max+1); % état adjoint

    [lb, ub] = fdcfunc.boxes (Gin1,GinN,M*N);

    %%% pour test : solution exacte de TRO
    if (EX == "TRO")
        UetG_sol = zeros(M*N+(Gin1+GinN)*(M-1),1);
        W_sol = zeros(M*N,1);
        for j=1:M
            i=N; gline = fdcfunc.gligne(i,j,N,M,Gin1,GinN);
            UetG_sol (gline) = min(1,max(0,(exp(tt(j))-exp(1/3))/(exp(2/3)-exp(1/3))));
            for i=1:N
                line = fdcfunc.ligne (i,j,N);
                UetG_sol (line) = exp(-tt(j)) * cos(xx(i));
                W_sol (line) = -exp(tt(j)) * cos(xx(i));
            end
        end
    end

    % solution par matlab d'un système parabolique 
    % POUR UNE FONCTION COUT DE LA FORME 1/2 \int_\Sigma (g - bobj)^2
    if (EX == "LIN")
        UetG_sol = zeros(M*N+(Gin1+GinN)*(M-1),1);
        sol = pdepe (0, @fdcfunc.pdefun, @fdcfunc.icfun, @fdcfunc.bcfun, xx, tt);
        UetG_sol(1:M*N) = reshape(sol(:,:,1)',[],1);
        if (Gin1 * GinN == 0) % un seul contrôle
            UetG_sol((M*N+1):end) = fdcfunc.bobj(l*GinN,tt(2:end)');
        else
            UetG_sol((M*N+1):end) = reshape([fdcfunc.bobj(0,tt(2:end));fdcfunc.bobj(l,tt(2:end))],[],1);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % itérations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oldJvalue = 10e7; newJvalue = 10e6; stop = 0; n=1;
    while (n <= Nit_max && ~stop)

        % actualisation des contraintes / coût
        [Cmat, Cvect] = fdcfunc.system (M,N,h,dt,UnetGn(:,n),Gin1,GinN,auto,stencil,"const");
%         [Cmat, Cvect] = fdcfunc.constraints (M,N,h,dt,UnetGn(:,n),Gin1,GinN);
%         [Cmat, Cvect] = fdcfunc.constraints (M,N,h,dt,UnetGn(:,n),Gin1,GinN,auto,stencil);
        [Jmat, Jvect] = fdcfunc.cost (M,N,h,dt,UnetGn(:,n),Wn(:,n),Gin1,GinN);

        % minimisation & actualisation de UnetGn
        UnetGn (:,n+1) = quadprog (Jmat,Jvect,[],[],Cmat,Cvect,lb,ub);

        % calcul de la fonctionnelle 
        oldJvalue = newJvalue;
        newJvalue = fdcfunc.Jvalue (UnetGn (:,n+1),dt,h,N,M,Gin1,GinN);
        ecart = abs(newJvalue - oldJvalue)/(oldJvalue+1e-6);
        fprintf("Valeur de la fonctionnelle à l'étape %d : %e, écart : %e\n", n, newJvalue, ecart);
        stop = (newJvalue < 1e-8 || ecart < 1e-8); % Nan < qqch = false. SALE

        % calcul de W_n
%         [admat, advect] = fdcfunc.system (M,N,h,dt,UnetGn(:,n+1),Gin1,GinN,auto,stencil,"adj");
        [admat, advect] = fdcfunc.adjointSys (M,N,h,dt,UnetGn(:,n+1),Gin1,GinN,0,stencil);
        Wn(:,n+1) = admat \ advect;

    %     % calcul de W_n
    %     [admat, advect] = fdcfunc.adjointSys2 (M,N,h,dt,UnetGn(:,n),UnetGn(:,n+1),Wn(:,n),Gin1,GinN);
    %     Wn(:,n+1) = admat \ advect;

        if (EX == "TRO")
            erreur_u = sqrt(h*dt)*norm(UnetGn(1:M*N,n+1) - UetG_sol(1:M*N));
            erreur_g = sqrt(h*dt)*norm(UnetGn(M*N+1:end,n+1) - UetG_sol(M*N+1:end));
            erreur_w = sqrt(h*dt)*norm(Wn(:,n+1) - W_sol);
            fprintf("erreur : u %f, g %f, w %f\n", erreur_u, erreur_g, erreur_w);
        end
        if (EX == "LIN")
            erreur_u = sqrt(h*dt)*norm(UnetGn(1:M*N,n+1) - UetG_sol(1:M*N));
            erreur_g = sqrt(h*dt)*norm(UnetGn(M*N+1:end,n+1) - UetG_sol(M*N+1:end));
            erreur_w = 1.0;
            fprintf("erreur : u %e, g %e\n", erreur_u, erreur_g);
        end
        
        n = n+1;
    end % itérations de minimisation
    
    if (comparaison)
        % enregistrement des erreurs finales
        erreurs (imn,2) = erreur_u;
        erreurs (imn,3) = erreur_g;
        erreurs (imn,4) = erreur_w;
    end % if comparaison
    
end % for different mesh size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% affichage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (affichage)
    hold on
%     subplot(121);
    figure
    utils.ploot (UnetGn(1:M*N,n), xx, tt);
    title("État");
%     if (EX=="LIN") utils.ploot (UetG_sol(1:M*N), xx, tt); end
%     subplot(122);
    figure
    fdcfunc.plotControl(Gin1,GinN,UnetGn(M*N+1:end,n));
    title("Contrôle");
%     if (EX=="LIN") fdcfunc.plotControl(Gin1,GinN,UetG_sol(M*N+1:end)); end
end % if affichage

if (comparaison && (length(erreurs(1,1:end))>1))
    figure
    for i=1:length(MNlist)
       fprintf("h=%.6f\terr_u=%e\n", erreurs(i,1),erreurs(i,2));
    end
    utils.plotErreur(erreurs(:,1),erreurs(:,2:4)',["u","g","w"]);
end % if comparaison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}










