%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fonctions associées à finiteDiffControl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef fdcfunc
    methods (Static) % uniquement pour le rangement

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assemblage du système à l'étape n
        %   Entrées : 
        %       * M,N int, dimensions (temps/espace)
        %       * h, dt real, pas du schéma
        %       * nu la pondération du laplacien 
        %       * UetG vecteur de taille (M*N+M) approximation courante
        %       * Gin1, GinN logicals : a-t-on un contrôle en i=1, i=N ?
        %       * auto : utiliser les schémas automatisés
        %       * sl : int, (stencil) rayon de la boule centrée en x pour une approx 
        %       * which : string, soit contraintes ("const") soit adjoint ("adj")
        %   Sorties : 
        %       * Cmat matrice M*N * (M*N+<nb de côtés avec contrôle>*(M-1)) 
        %       * Cvect vecteur M*N + <nb de côtés avec contrôle>*(M-1)
        %   les contraintes sont Cmat*X = Cvect
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [mat, vect] = system (M,N,h,dt,UnetGn,Gin1,GinN,auto,sl,sys)
            global xx tt nu
            mat = sparse (M*N,M*N +(Gin1+GinN)*(M-1)*(sys=="const")); 
            vect = sparse (M*N,1);
            if (sys=="const")
                % condition initiale pour les constraintes
                mat(1:N,1:N) = eye (N);
                vect (1:N) = fdcfunc.u_0(xx)';
            else
                % condition finale pour le système adjoint
                for i=1:N
                    line = fdcfunc.ligne (i,M,N);
                    mat(line,line) = 1.0;
                    vect (line) = vect (line) + fdcfunc.phi_u(xx(i),UnetGn(line));
                end
            end
            for j=((sys=="adj")+2*(sys==("const"))):(M-(sys=="adj")) % temps > 0
                for i=[1,N] % bords
                    line = fdcfunc.ligne (i,j,N); 
                    bu = fdcfunc.b_u(xx(i),tt(j),UnetGn(line));
                    if (auto) 
                        localScheme = fdcfunc.getLocalScheme(i,N,h,"partial_nu",sl,[],[]);
                        localScheme(sl+2,sl+2) = localScheme(sl+2,sl+2) + bu;
                    else
                        localScheme = fdcfunc.partial_nu (i,h) + diag([0,0,bu,0,0]);
                    end
                    mat(line,1:M*N) = utils.flatten(localScheme,j,i,M,N);
                    if (sys=="const")
                        vect (line) = vect(line) + bu * UnetGn(line) - fdcfunc.b(xx(i),tt(j),UnetGn(line));
                        if ((Gin1 && i==1) || (GinN && i==N))
                            gline = fdcfunc.gligne (i,j,N,M,Gin1,GinN);
                            mat(line,gline) = mat(line,gline) - 1.0; 
                        end 
                    else
                        gline = fdcfunc.gligne (i,j,N,M,Gin1,GinN);
                        vect (line) = vect (line) + fdcfunc.psi_u(xx(i),tt(j),UnetGn(line),UnetGn(gline));
                    end
                end
                for i=2:N-1 % intérieur
                    line = fdcfunc.ligne (i,j,N); 
                    du = fdcfunc.d_u(xx(i),tt(j),UnetGn(line));
                    if (auto)
                        if (sys=="const")
                            [locDt, scdmdt] = fdcfunc.getLocalScheme(j,M,dt,"partial_t",sl,[vect(fdcfunc.ligne(i,1,N))],[]);
                        else
                            [locDt, scdmdt] = fdcfunc.getLocalScheme(j,M,dt,"-partial_t",sl,[],[vect(fdcfunc.ligne(i,M,N))]);                            
                        end
                        loclap = fdcfunc.getLocalScheme(i,N,h,"lap",sl,[],[]);
                        localScheme = utils.mypad(locDt,round(0.5*(length(loclap(1,:))-length(locDt(1,:)))),0) ...
                            - nu * loclap + utils.mypad(du,round(0.5*(length(loclap(1,:)))-1),0);
                    else
                        scdmdt = 0; 
                        localScheme = fdcfunc.Euler(j,M,dt,1) - nu * fdcfunc.Lap (h,i,N) + diag([0,0,du,0,0]);
                    end
                    mat(line,1:M*N) = utils.flatten(localScheme,j,i,M,N); 
                    
                    if (sys=="const")
                        vect (line) = vect (line) + du * UnetGn(line) - fdcfunc.d(xx(i),tt(j),UnetGn(line)) - scdmdt;
                    else
                        vect (line) = vect (line) + fdcfunc.varphi_u(xx(i),tt(j),UnetGn(line));
                    end
                end % for i
            end % for j
        end % constraints
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assemblage du système à l'étape n
        %   Entrées : 
        %       * M,N int, dimensions (temps/espace)
        %       * h, dt real, pas du schéma
        %       * nu la pondération du laplacien 
        %       * UetG vecteur de taille (M*N+M) approximation courante
        %       * Gin1, GinN logicals : a-t-on un contrôle en i=1, i=N ?
        %   Sorties : 
        %       * Cmat matrice M*N * (M*N+<nb de côtés avec contrôle>*(M-1)) 
        %       * Cvect vecteur M*N + <nb de côtés avec contrôle>*(M-1)
        %   les contraintes sont Cmat*X = Cvect
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [mat, vect] = constraints (M,N,h,dt,UnetGn,Gin1,GinN)
            global xx tt nu
            mat = sparse (M*N,M*N +(Gin1+GinN)*(M-1)); 
            vect = sparse (M*N,1);
            % condition initiale pour les constraintes
            mat(1:N,1:N) = eye (N);
            vect (1:N) = fdcfunc.u_0(xx)';
            for j=2:M % temps > 0
                for i=[1,N] % bords
                    line = fdcfunc.ligne (i,j,N); 
                    bu = fdcfunc.b_u(xx(i),tt(j),UnetGn(line));
                    localScheme = fdcfunc.partial_nu (i,h) + diag([0,0,bu,0,0]);
                    mat(line,1:M*N) = utils.flatten(localScheme,j,i,M,N);
                    vect (line) = vect(line) + bu * UnetGn(line) - fdcfunc.b(xx(i),tt(j),UnetGn(line));
                    if ((Gin1 && i==1) || (GinN && i==N))
                        gline = fdcfunc.gligne (i,j,N,M,Gin1,GinN);
                        mat(line,gline) = mat(line,gline) - 1.0; 
                    end 
                end
                for i=2:N-1 % intérieur
                    line = fdcfunc.ligne (i,j,N); 
                    du = fdcfunc.d_u(xx(i),tt(j),UnetGn(line));
                    localScheme = fdcfunc.Euler(j,M,dt,1) - nu * fdcfunc.Lap (h,i,N) + diag([0,0,du,0,0]);
                    mat(line,1:M*N) = utils.flatten(localScheme,j,i,M,N); 
                    vect (line) = vect (line) + du * UnetGn(line) - fdcfunc.d(xx(i),tt(j),UnetGn(line));
                end % for i
            end % for j
        end % system2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assemblage du système adjoint à l'étape n
        %   Entrées : 
        %       * M,N int, dimensions (temps/espace)
        %       * h, dt real, pas du schéma
        %       * UnetGn vecteur de taille (M*N+M)
        %       * auto : utiliser les schémas automatisés
        %     < * sl : int, rayon de la boule centrée en x pour une approx >
        %   Sorties : 
        %       * admat matrice (M*N)^2
        %       * advect vecteur M*N 
        %   le système est admat * w = advect
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [admat, advect] = adjointSys (M,N,h,dt,UnetGn,Gin1,GinN,auto,sl)
            global xx tt nu
            admat = sparse(M*N,M*N);
            advect = sparse (M*N,1);
            for j=1:M-1 % temps < T
                for i=[1,N] % bords
                    line = fdcfunc.ligne (i,j,N); 
                    gline = fdcfunc.gligne (i,j,N,M,Gin1,GinN);
                    if (auto) 
                        localScheme = fdcfunc.selectOp(h,i,N,"partial_nu",sl)...
                            + diag([zeros(1,sl),fdcfunc.b_u(xx(i),tt(j),UnetGn(line)),zeros(1,sl)]);
                    else
                        localScheme = fdcfunc.partial_nu (i,h) + diag([0,0,fdcfunc.b_u(xx(i),tt(j),UnetGn(line)),0,0]);
                    end
                    admat (line,1:M*N) = utils.flatten(localScheme,j,i,M,N);
                    
                    advect (line) = advect (line) + fdcfunc.psi_u(xx(i),tt(j),UnetGn(line),UnetGn(gline));
                end 
                for i=2:N-1 % intérieur
                    line = fdcfunc.ligne (i,j,N);
                    if (auto)
                        localScheme = - fdcfunc.selectOp(dt,j,M,"partial_t",sl) ...
                            - nu * fdcfunc.selectOp(h,i,N,"lap",sl) ...
                            + diag([zeros(1,sl),fdcfunc.d_u(xx(i),tt(j),UnetGn(line)),zeros(1,sl)]);
                    else
                        localScheme = fdcfunc.Euler(j,M,dt,0) - nu * fdcfunc.Lap (h,i,N) ...
                            + diag([0,0,fdcfunc.d_u(xx(i),tt(j),UnetGn(line)),0,0]);
                    end
                    admat(line,1:M*N) = utils.flatten(localScheme,j,i,M,N);
                    
                    advect (line) = advect (line) + fdcfunc.varphi_u(xx(i),tt(j),UnetGn(line));
                end % for i
            end % for j
            % condition finale
            j=M;
            for i=1:N
                line = fdcfunc.ligne (i,j,N);
                admat(line,line) = 1.0;
                advect (line) = advect (line) + fdcfunc.phi_u(xx(i),UnetGn(line));
            end
        end % adjointSys
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assemblage du système adjoint à l'étape n, version 2
        %   Entrées : 
        %       * M,N int, dimensions (temps/espace)
        %       * h, dt real, pas du schéma
        %       * UnetGn vecteur de taille (M*N+M)
        %   Sorties : 
        %       * admat matrice (M*N)^2
        %       * advect vecteur M*N 
        %   le système est admat * w = advect
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [admat, advect] = adjointSys2 (M,N,h,dt,UnetGn,UetG,Wn,Gin1,GinN)
            global xx tt nu
            admat = sparse(M*N,M*N);
            advect = sparse (M*N,1);
            for j=1:M-1 % temps < T
                for i=[1,N] % bords
                    line = fdcfunc.ligne (i,j,N); 
                    gline = fdcfunc.gligne (i,j,N,M,Gin1,GinN);
                    localScheme = fdcfunc.partial_nu (i,h) + diag([0,0,fdcfunc.b_u(xx(i),tt(j),UnetGn(line)),0,0]);
                    admat (line,1:M*N) = utils.flatten(localScheme,j,i,M,N);
                    
                    advect (line) = advect (line) + fdcfunc.psi_u(xx(i),tt(j),UnetGn(line),UnetGn(gline)) ...
                        + (fdcfunc.psi_uu(xx(i),tt(j),UnetGn(line),UnetGn(gline)) - Wn(line) * fdcfunc.b_uu (xx(i),tt(j),UnetGn(line)))...
                        * (UetG(line)-UnetGn(line));
                end 
                for i=2:N-1 % intérieur
                    line = fdcfunc.ligne (i,j,N);
                    localScheme = fdcfunc.Euler(j,M,dt,0) - nu * fdcfunc.Lap (h,i,N) ...
                        + diag([0,0,fdcfunc.d_u(xx(i),tt(j),UnetGn(line)),0,0]);
                    admat(line,1:M*N) = utils.flatten(localScheme,j,i,M,N);
                    
                    advect (line) = advect (line) + fdcfunc.varphi_u(xx(i),tt(j),UnetGn(line)) ...
                        + (fdcfunc.varphi_uu(xx(i),tt(j),UnetGn(line)) - Wn(line) * fdcfunc.d_uu (xx(i),tt(j),UnetGn(line)))...
                        * (UetG(line)-UnetGn(line));
                end % for i
            end % for j
            % condition finale
            j=M;
            for i=1:N
                line = fdcfunc.ligne (i,j,N);
                admat(line,line) = 1.0;
                advect (line) = advect (line) + fdcfunc.phi_u(xx(i),UnetGn(line)) + fdcfunc.phi_uu(xx(i),UnetGn(line)) * (UetG(line)-UnetGn(line));
            end
        end % adjointSys2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % schémas pour SQP
        %   Entrées : 
        %       * k, K : int, indice du point courant et max
        %       * delta : double, pas (spatial/temporel) de la discrétisation
        %       * which : string, quel opérateur
        %       * sl : int, (stencil) rayon de la boule centrée en x pour une approx 
        %   Sorties : 
        %       * localScheme : petite matrice, filtre à appliquer 
        %       * scdm : terme à RETRANCHER au second membre de la ligne
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [localScheme, scdm] = getLocalScheme (k,K,delta,which,sl,Lb,Rb)
            switch (which)
                case "lap"
                    up = sl+min(k-1-sl,0)-min(K-k-sl,0); down = sl+min(K-k-sl,0)-min(k-1-sl,0);
                    scheme = utils.TaylorScheme (delta*(-up:down),2);
                    scdm = 0;
                    nsl = max(up,down); % new sl
                    localScheme = zeros(2*nsl+1);
                    localScheme (nsl+1,nsl+1-up:nsl+1+down) = scheme;
                case "partial_t"
                    up = sl+min(k-1-sl,0); down = max(0,sl-up);
                    scheme = utils.TaylorScheme (delta*(-up:down),1);
                    ls = max(0,up-(k-2)); rs = max(0,down-(K-k)); % left / right shifts
                    scdm = sum (scheme(1:ls).*Lb);
                    localScheme = zeros(2*sl+1);
                    localScheme (sl+1-up+ls:sl+1+down-rs,sl+1) = scheme(1+ls:end-rs)';                    
                case "-partial_t" % adjoint system : reverse time
                    down = sl+min(K-k-sl,0);
                    up = max(0,sl-down); 
                    scheme = utils.TaylorScheme (delta*(-up:down),1);
                    ls = max(0,up-(k-1)); rs = max(0,down-(K-1-k)); % left / right shifts
                    scdm = sum (scheme(end-rs+1:end).*Rb);
                    localScheme = zeros(2*sl+1);
                    localScheme (sl+1-up+ls:sl+1+down-rs,sl+1) = scheme(1+ls:end-rs)';
                case "partial_nu"
                    up=(sl+1)*(k==K); down = (sl+1)*(k==1);
                    scheme = utils.TaylorScheme (delta*(-up:down),1);
                    scdm = 0;
                    localScheme = zeros(2*(sl+1)+1);
                    localScheme (sl+2,sl+2-up:sl+2+down) = (1-2*(k==1)) * scheme;
                otherwise 
                    fprintf("erreur getLocalScheme, opérateur %s non reconnu\n", which);
            end % switch
        end % getLocalScheme
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % schéma local du laplacien pour SQP
        %   Entrées : 
        %       * h pas spatial du schéma
        %   Sorties : 
        %       * filtre à appliquer pour le Laplacien
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = Lap (h,i,N)
            res = 1/(h*h) * [zeros(2,5); [0,1,-2, 1,0]; zeros(2,5)];
%             if (i==2)
%                 res = 1/(h*h) * [zeros(2,5); [0, 13/12, -27/12, 15/12,-1/12]; zeros(2,5)];
%             else if (i==N-1)
%                 res = 1/(h*h) * [zeros(2,5); [-1/12, 15/12, -27/12, 13/12,0]; zeros(2,5)];
%                 else
%                 res = 1/(h*h) * [zeros(2,5); [-1/12, 4/3, -5/2, 4/3,-1/12]; zeros(2,5)];
%                 end
%             end
        end % Lap
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % schéma local d'Euler pour SQP
        %   Entrées : 
        %       * j point courant, M nb de points en temps
        %       * dt pas temporel
        %       * fwd logical : sens forward ?
        %   Sorties : 
        %       * filtre à appliquer 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = Euler (j,M,dt,fwd)
            % upwind
%             res = [zeros(5,2), [0; -1/dt*(j~=1)*fwd; 1/dt; -1/dt*(j~=M)*(1-fwd); 0], zeros(5,2)];

%             if (fwd)
%                 if (j==2) res = 1/(2*dt) .* [zeros(5,2), [0;-1;0;1;0], zeros(5,2)]; 
% %                 if (j==2) res = 1/(2*dt) .* [zeros(5,2), [0;-2;2;0;0], zeros(5,2)]; 
%                 else res = 1/(2*dt) * [zeros(5,2), [1;-4;3;0;0], zeros(5,2)]; end
%             else
%                 if (j==M-1) res = 1/(2*dt) .* [zeros(5,2), [0;1;0;-1;0], zeros(5,2)]; 
% %                 if (j==2) res = 1/(2*dt) .* [zeros(5,2), [0;0;2;-2;0], zeros(5,2)]; 
%                 else res = 1/(2*dt) * [zeros(5,2), [0;0;3;-4;1], zeros(5,2)]; end
%             end

            if (fwd)
                if (j==2) 
                    scheme = utils.TaylorScheme (dt*(-1:0),1); 
                    res = zeros(5,5);
                    res (2:3,3) = scheme;
                else res = [zeros(5,2), [0;utils.TaylorScheme(dt*(-1:0),1);0;0], zeros(5,2)]; end
            else
                if (j==M-1) 
                    scheme = utils.TaylorScheme (dt*(0:1),1); 
                    res = zeros(5,5);
                    res (3:4,3) = -scheme;
                else res = [zeros(5,2), [0;0;-utils.TaylorScheme(dt*(0:1),1);0], zeros(5,2)]; end
%                 res = [zeros(5,2), [0; -1/dt*(j~=1)*fwd; 1/dt; -1/dt*(j~=M)*(1-fwd); 0], zeros(5,2)];
            end
        end % Euler
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % schéma pour la dérivée normale en espace
        %   Entrées : 
        %       * i point courant en espace
        %       * h pas spatial du schéma
        %   Sorties : 
        %       * filtre à appliquer 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = partial_nu (i,h)
%             if (i==1)
%                 res = [zeros(2,5); [0,0,1/h,-1/h,0]; zeros(2,5)];
%             else
%                 res = [zeros(2,5); [0,-1/h, 1/h,0,0]; zeros(2,5)];
%             end

            if (i==1)
                res = 1/(2*h) * [zeros(2,5); [0,0,3,-4,1]; zeros(2,5)];
            else
                res = 1/(2*h) * [zeros(2,5); [1,-4,3,0,0]; zeros(2,5)];
            end
        end % partial_nu
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assemblage du coût à l'étape n
        %   Entrées : 
        %       * M,N int, dimensions (temps/espace)
        %       * h, dt real, pas du schéma
        %       * UetG vecteur de taille (M*N+M) approximation courante
        %       * W vecteur de taille (M*N), adjoint courant
        %       * Gin1, GinN logicals : a-t-on un contrôle en i=1, i=N ?
        %   Sorties : 
        %       * Jmat matrice M*N * dim(UetG)
        %       * Jvect vecteur même dim que UetG
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Jmat, Jvect] = cost (M,N,h,dt,UnetGn,W,Gin1,GinN)
            global xx tt 
            Jmat = sparse(length(UnetGn),length(UnetGn));
            Jvect = sparse(length(UnetGn),1);
            
            % code pour le vecteur z := (u-u^n, g-g^n)^t. Changement après
            for j=2:M
                for i=1:N
                    line = fdcfunc.ligne (i,j,N);
                    % termes sur tout le domaine U : quadratique z-z
                    Jmat(line,line) = Jmat(line,line) + h*dt/2 * (fdcfunc.varphi_uu(xx(i),tt(j),UnetGn(line))...
                        - W(line) * fdcfunc.d_uu(xx(i),tt(j),UnetGn(line)));
                    
                    % termes linéaires en z
                    Jvect (line) = Jvect(line) + h * dt * fdcfunc.varphi_u (xx(i),tt(j),UnetGn(line));
                end % for i
                for i=[1,N] % bords
                    line = fdcfunc.ligne (i,j,N);
                    gline = fdcfunc.gligne (i,j,N,M,Gin1,GinN);
                    % termes quadratiques de bord en u
                    Jmat(line,line) = Jmat(line,line) + dt/2 * (fdcfunc.psi_uu(xx(i),tt(j),UnetGn(line),UnetGn(gline))...
                        - W(line) * fdcfunc.b_uu(xx(i),tt(j),UnetGn(line)));
                    % termes incluant le contrôle
                    if ((i==1 && Gin1) || (i==N && GinN))
                        terme_croise = 0.5*(fdcfunc.psi_ug(xx(i),tt(j),UnetGn(line),UnetGn(gline)) ...
                                          + fdcfunc.psi_gu(xx(i),tt(j),UnetGn(line),UnetGn(gline)));
                        Jmat(line,gline) = Jmat(line,gline) + dt/2 * terme_croise;
                        Jmat(gline,line) = Jmat(gline,line) + dt/2 * terme_croise;
                        Jmat (gline,gline) = Jmat (gline,gline) + dt/2 * fdcfunc.psi_gg(xx(i),tt(j),UnetGn(line),UnetGn(gline));

                        % termes linéaire
                        Jvect (line)  = Jvect(line)  + dt * fdcfunc.psi_u (xx(i),tt(j),UnetGn(line),UnetGn(gline));
                        Jvect (gline) = Jvect(gline) + dt * fdcfunc.psi_g (xx(i),tt(j),UnetGn(line),UnetGn(gline));
                    end % if bord 1 ou bord N
                end % for i (bords)
            end % for j
            
            % termes à j=M sur Omega
            for i=1:N % bord j=M
                line = fdcfunc.ligne(i,j,N);
                Jmat(line,line) = Jmat(line,line) + 1/2 * h * fdcfunc.phi_uu (xx(i),UnetGn(line));
                Jvect (line)  = Jvect(line)  + h * fdcfunc.phi_u (xx(i),UnetGn(line));
            end % for i
            
            % change of variable z -> (u,g)
            Jvect = Jvect - (Jmat + Jmat') * UnetGn;
            % because quadprog solves for 1/2 X^t Jmat X + b . X (notice the 1/2)
            Jmat = 2 * Jmat;
        end % cost
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % récupère les box constraints en fonction de UnetGn et du temps
        %   Entrées : 
        %       * Gin1, GinN logicals : contrôle en x=0 ? x=l ?
        %       * shift le nombre de coords à ajouter avant 
        %   Sorties : 
        %       * lb, ub vecteurs de taille shift + taille Gn
        %   la variable de la minimisation est \tilde{g} = g_n+1 - g_n
        %   donc pour assurer que g_n+1 \in A, on impose \tilde{g} + g_n \in A
        %   soit \tilde{g} \in A - g_n
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [lb, ub] = boxes (Gin1,GinN,shift)
            global xx tt
            if (Gin1 && GinN)
                lb = [fdcfunc.g_min(xx(1),tt(2:end)'), fdcfunc.g_min(xx(end),tt(2:end)')];
                ub = [fdcfunc.g_max(xx(1),tt(2:end)'), fdcfunc.g_max(xx(end),tt(2:end)')];
                lb = reshape (lb', [], 1);
                ub = reshape (ub', [], 1);
            else if (Gin1 || GinN)
                    lb = fdcfunc.g_min (xx(1)*Gin1 + xx(end)*GinN, tt(2:end)');
                    ub = fdcfunc.g_max (xx(1)*Gin1 + xx(end)*GinN, tt(2:end)');
                else
                    lb = []; ub = []; % no control
                end
            end
            lb = [-Inf*ones(shift,1); lb];
            ub = [ Inf*ones(shift,1); ub];
        end % boxes        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % évaluation de la fonctionnelle de coût
        %   Entrées : 
        %       * UnetGn l'approx actuelle
        %       * dt,h pas du maillage
        %       * N,M dimensions du maillage
        %       * Gin1,GinN logicals
        %   Sorties : 
        %       * res la valeur du coût
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = Jvalue (UnetGn,dt,h,N,M,Gin1,GinN)
            global xx tt
            res = 0.0;
            for j=2:M
                for i=[1,N] % bord
                    line = fdcfunc.ligne(i,j,N);
                    gline = fdcfunc.gligne (i,j,N,M,Gin1,GinN);
                    res = res + dt * fdcfunc.psi (xx(i),tt(j),UnetGn(line),UnetGn(gline));
                end
                for i=1:N % intérieur
                    line = fdcfunc.ligne(i,j,N);
                    res = res + dt * h * fdcfunc.varphi (xx(i),tt(j),UnetGn(line));
                end % for i 
            end % for j
            % temps final
            for i=1:N                 
                line = fdcfunc.ligne(i,j,N);
                res = res + h * fdcfunc.phi (xx(i),UnetGn(line));
            end % for i
        end % Jvalue
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fonctions de la simulation
        %   Entrées : 
        %       * x,t le point courant
        %       * u (réel) la valeur (approchée) de u(x,t)
        %     < * g (réel) la valeur (approchée) de g(x,t) >
        %   Sorties : 
        %       * res un réel.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function res = phi (x,u)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.5 * (u - (exp(1) + exp(-1))*cos(x))^2;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0 * (u - cos(12*x))^2;
            end
        end % phi
        function res = phi_u (x,u)
            global EX            
            switch (EX)
                case "TRO"
                    res = u - (exp(1) + exp(-1))*cos(x);
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0;
            end
        end % phi_u
        function res = phi_uu (x,u)
            global EX
            switch (EX)
                case "TRO"
                    res = 1.0;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0;
            end
        end % phi_uu
        
        function res = obj (x,t)
            if (x < 0.3)
                ee = 1.1;
                res = sin(1/(ee - t)) - sin(1/ee);
            else 
               res = 0.0;
            end
        end % obj
        function res = varphi (x,t,u)
            global EX            
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0 * (u - fdcfunc.obj(x,t))^2;
                otherwise
                    res = 0.0;
            end
        end % varphi
        function res = varphi_u (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0*(u - fdcfunc.obj(x,t));
                otherwise
                    res = 0.0;
            end
        end % varphi_u
        function res = varphi_uu (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0;
            end
        end % varphi_uu
        
        function res = psi (x,t,u,g)
            global EX
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = 0.0;
                    else 
                        res = - exp(-2*t) * u + exp(1/3)/sqrt(2) * g ...
                            + ((exp(2/3)-exp(1/3))/sqrt(2))/2 * g^2; 
                    end
                case "LIN"
                    res = 0.5*(g - fdcfunc.bobj(x,t))^2;
                case "DISC"
                    res = 0.5*(u - fdcfunc.bobj(x,t))^2;
                otherwise
                    res = 0.0;
            end
        end % psi        
        function res = bobj (x,t)
            if (x<0.3)
                res = 0;
            else 
                res = - sqrt(0.01+abs(sin(5*t)));
            end
        end % bobj (boundary objective)
        function res = psi_u (x,t,u,g)
            global EX            
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = 0.0;
                    else 
                        res = - exp(-2*t); 
                    end
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = (u - fdcfunc.bobj(x,t));
                otherwise
                    res = 0.0;
            end
        end % psi_u
        function res = psi_g (x,t,u,g)
            global EX            
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = 0.0;
                    else 
                        res = exp(1/3)/sqrt(2) + (exp(2/3)-exp(1/3))/sqrt(2) * g;
                    end
                case "LIN"
                    res = g - fdcfunc.bobj(x,t);
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0;
            end
        end % psi_g
        function res = psi_uu (x,t,u,g)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 1.0;
                otherwise
                    res = 0.0;
            end
        end % psi_uu
        function res = psi_ug (x,t,u,g)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0;
            end
        end % psi_ug
        function res = psi_gu (x,t,u,g)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0;
            end
        end % psi_gu
        function res = psi_gg (x,t,u,g)
            global EX
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = 0.0;
                    else
                        res = (exp(2/3)-exp(1/3))/sqrt(2);
                    end
                case "LIN"
                    res = 1.0;
                case "DISC"
                    res = 0.0;
                otherwise
                    res = 0.0;
            end
        end % psi_gg
        
        function res = d (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = exp(x);
                case "DISC"
                    res = u;
                otherwise
                    res = 0.0;
            end
        end % d
        function res = d_u (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = exp(x);
                case "DISC"
                    res = 1;
                otherwise
                    res = 0.0;
            end
        end % d_u
        function res = d_uu (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    res = 0.0;
                case "LIN"
                    res = exp(x);
                case "DISC"
                    res = 0;
                otherwise
                    res = 0.0;
            end
        end % d_uu
                
        function res = b (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = 0.0;
                    else
                        res = u - fdcfunc.beta(t) + u * abs(u)^3;
                    end
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = exp(t*u);
                otherwise
                    res = 0.0;
            end
        end % b
        function res = beta (t) % pour l'exemple Tröltzsch
            res = (exp(-4*t)/4 - min(1,max(0,(exp(t)-exp(1/3))/(exp(2/3)-exp(1/3)))));
        end % beta
        function res = b_u (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = 0.0;
                    else
                        res = 1.0 + 4 * u * abs(u)^2;
                    end
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = t * exp(t*u);
                otherwise
                    res = 0.0;
            end
        end % b_u
        function res = b_uu (x,t,u)
            global EX
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = 0.0;
                    else
                        res = 12 * u * abs(u);
                    end
                case "LIN"
                    res = 0.0;
                case "DISC"
                    res = t^2 * exp(t*u);
                otherwise
                    res = 0.0;
            end
        end % b_uu
        
        function res = u_0 (x)
            global EX 
            switch (EX)
                case "TRO"
                    res = cos(x);
                case "LIN"
                    res = -1 + 2*x/0.78;
                case "DISC"
                    res = 0;
                otherwise
                    res = 0.0;
            end
        end % u_0
        
        function res = g_min (x,t)
            global EX
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = -Inf * ones(size(t));
                    else
                        res = 0.0 * ones (size(t));
                    end
                case "LIN"
                    res = - Inf * ones(size(t));
                case "DISC"
                    res = - Inf * ones(size(t));
                otherwise
                    res = - 90 * ones(size(t));
            end
        end % g_min
        function res = g_max (x,t)
            global EX
            switch (EX)
                case "TRO"
                    if (x < 0.3)
                        res = Inf * ones(size(t));
                    else
                        res = ones(size(t));
                    end
                case "LIN"
                    res = Inf * ones(size(t));
                case "DISC"
                    if (x < 0.3)
                        res = Inf * ones(size(t));
                    else 
                        res = 0.0 * ones(size(t));
                    end
                otherwise
                    res = 50 * ones(size(t));
            end
        end % g_max
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calcul de la ligne correspondante (évite les oublis)
        %   Entrées : 
        %       * i,j,N,M point courant & taille de maillage
        %       * Gin1, GinN logicals
        %   Sorties : 
        %       * indice de la ligne demandée
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = ligne (i,j,N)
            res = (j-1)*N + i;
        end % ligne
        function res = gligne (i,j,N,M,Gin1,GinN)
            res = M*N+(j-1)*(Gin1+GinN)-Gin1*GinN*(i==1);
        end % gligne        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fonctions pour le solveur de matlab pdepe
        %   Entrées : 
        %       * voir doc
        %   Sorties : 
        %       * voir doc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [c,f,s] = pdefun (x,t,u,dudx) % equation
            global nu
            c = 1; % coefficient in front of du/dt
            f = nu * dudx; % defines a classical laplacian * nu
            s = - fdcfunc.d(x,t,u); % uncontrolled equation in domain
        end % pdefun
        
        function u0 = icfun (x) % initial condition
            u0 = fdcfunc.u_0(x);
        end % icfun
        
        function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t) % boundary condition
            pL = fdcfunc.b(xL,t,uL) - fdcfunc.bobj(xL,t);
            qL = -1; % neumann exterior condition
            pR = fdcfunc.b(xR,t,uR) - fdcfunc.bobj(xR,t);
            qR = 1; % neumann exterior condition
        end % bcfun
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % affichage du contrôle
        %   Entrées : 
        %       * Gin1, GinN logicals
        %       * Gn le vecteur des contrôles à afficher
        %   Sorties : 
        %       * 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotControl (Gin1,GinN,Gn)
            global tt
            if (Gin1 || GinN)
                plot (tt(2:end), reshape (Gn, Gin1+GinN, []), marker='x');
                if (Gin1)
                    leg = ["contrôle en x=0"];
                end
                if (GinN)
                    if (Gin1) leg = ["contrôle en x=0","contrôle en x=l"]; else leg = ["contrôle en x=l"]; end
                end
                legend (leg);
                xlim padded 
                ylim padded
            end % if controls to plot
        end % plotControl
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        %   Entrées : 
        %       * 
        %   Sorties : 
        %       * 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function 
        % end % 
        
    end % methods static
end % class fdcfunc