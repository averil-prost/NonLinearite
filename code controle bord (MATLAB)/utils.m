%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tool box pour les différences finies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef utils
    methods (Static) % uniquement pour le rangement

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construction d'un opérateur discrétisé matriciel STATIONNAIRE 1D
        % adapté pour stencils arbitraires
        %   Entrées : 
        %       * localScheme le schéma local, de taille impaire
        %       * dim la taille du domaine 
        %       * left, right les localScheme du bord
        %   Sorties : 
        %       * Mat telle que Mat * V = S(V) le schéma
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Mat = getScheme1D (localScheme, dim, left, right)
            Mat = zeros(dim, dim);
            Mat (1:length(left(:,1)), 1:length(left(1,:))) = left;
            % pasting every diagonal (long live 1D)
            dep = 1+length(left(:,1)); arr = dim-length(right(:,1));
            for j=1:length(localScheme(1,:))
                Mat(dep:arr, j:arr-dep+j) = Mat(dep:arr, j:arr-dep+j) + localScheme(j) .* eye(arr-dep+1);
            end % for
            Mat (end-length(right(:,1))+1:end, end-length(right(1,:))+1:end) = right;
        end % getScheme1D
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construction d'un opérateur discrétisé matriciel STATIONNAIRE 2D
        % adapté pour stencil de 9 points
        %   Entrées : 
        %       * localScheme le schéma local, carré, de taille impaire
        %       * n,m la taille du domaine (x, y)
        %       * left, right, top, bottom les localScheme du bord
        %       * [tb][lr]corn schemas des coins (top,bottom,left,right)
        %   Sorties : 
        %       * Mat telle que Mat * V = S(V) le schéma
        %
        %    ---------------> x       Les variables sont numérotées
        %    | 1,   2 ..              x d'abord, i va de 1 à n
        %    | n+1, n+2               y en ligne, j va de 1 à m
        %    | ...                    Le schéma fait M (j,i) !!!
        %    | (m-1)*n+1, ...
        %  y v
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Mat = getScheme2D (localScheme,m,n,...
                left,right,top,bottom,tlcorn,trcorn,blcorn,brcorn)
            Mat = zeros(m*n, m*n); 
            %%% top
            start=0; % useless but easier for the reader
            Mat(start+1,:) = utils.flatten(tlcorn, 1,1, m,n); % top left corner
            for i=2:n-1 % top line but corners
                Mat(start+i,:) = utils.flatten(top, 1,i, m,n);
            end % for i
            Mat(start+n,:) = utils.flatten(trcorn, 1,n, m,n); % top right corner
            %%% body
            for j=2:m-1
                start = (j-1)*n;
                Mat(start+1,:) = utils.flatten (left, j,1, m,n); % right
                for i=2:n-1 % interior
                    Mat (start+i,:) = utils.flatten(localScheme, j,i,m,n);
                end % for i
                Mat (start+n,:) = utils.flatten (right, j,n, m,n); % left
            end % for j
            %%% bottom
            start=(m-1)*n;
            Mat (start+1,:) = utils.flatten (blcorn, m,1, m,n); % bottom left corner
            for i=2:n-1 % bottom line but corners
                Mat (start+i,:) = utils.flatten (bottom, m,i, m,n); 
            end % for i
            Mat (start+n,:) = utils.flatten (brcorn, m,n, m,n); % bottom right corner
        end % getScheme2D
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vector-equivalent of a 2D scheme on a m-n grid
        %   Entrées : 
        %       * localScheme a little matrix, size must be odd
        %       * j,i position of the center of localScheme in grid
        %           !! truncating if needed to stay into the grid
        %       * m,n size of aforementionned grid
        %   Sorties : 
        %       * flat a vector of m*n coordinates
        %
        %   Computes A by placing localScheme at position (j,i)
        %   then transforms the matrix A = (a_ji)_ji into
        %   the vector [a_11, ..., a_1n, a21,...,a2n, ... am1, ... amn]
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function flat = flatten (localScheme,j,i,m,n)
            thegrid = zeros(m,n);
            jshift = round((length(localScheme(:,1))-1)/2);
            ishift = round((length(localScheme(1,:))-1)/2);
            if (j-jshift >= 1 && j+jshift <= m && i-ishift >= 1 && i+ishift <= n)
                thegrid(j-jshift:j+jshift,i-ishift:i+ishift) = localScheme;
            else
                % corner-specific mode on
                
                % version where localScheme is shifted
%                 vshift = max(1+jshift-j,0) + min(m-jshift-j,0); % vertical
%                 hshift = max(1+ishift-i,0) + min(n-ishift-i,0); % horizontal
%                 jtshift = jshift-vshift; jbshift = jshift+vshift;
%                 ilshift = ishift-hshift; irshift = ishift+hshift;
%                 thegrid(j-jtshift:j+jbshift,i-ilshift:i+irshift) = localScheme;
                
                % version where localscheme is truncated
                jtshift = min(jshift,j-1); % top shift
                jbshift = min(jshift,m-j); % bottom shift
                ilshift = min(ishift,i-1); % left shift
                irshift = min(ishift,n-i); % right shift
%                  fprintf("%d, %d, %d, %d\n%d, %d, %d, %d\n", jtshift, jbshift, ilshift, irshift,...
%                      max(0,jshift-jtshift), max(0,jshift-jbshift), max(0,ishift-ilshift),max(0,ishift-irshift));
                thegrid(j-jtshift:j+jbshift,i-ilshift:i+irshift) = localScheme ...
                    (1+max(0,jshift-jtshift):end-max(0,jshift-jbshift), ...
                     1+max(0,ishift-ilshift):end-max(0,ishift-irshift)); 
            end % if
            flat = reshape (thegrid', 1, []);
        end % 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % même chose mais localScheme contient des fonctions anonymes
        %   Entrées : 
        %       * localScheme : (j,i) -> schéma local autour de (j,i)
        %   Sorties : 
        %       * Mat telle que Mat * V = S(V) le schéma
        %
        %    ---------------> x       Les variables sont numérotées
        %    | 1,   2 ..              x d'abord, i va de 1 à n
        %    | n+1, n+2               y en ligne, j va de 1 à m
        %    | ...                    Le schéma fait M (j,i) !!!
        %    | (m-1)*n+1, ...
        %  y v
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Mat = getScheme2DEval (localScheme,m,n)
            Mat = zeros(m*n, m*n); 
            for j=1:m
                for i=1:n
                    Mat((j-1)*n+i,:) = utils.flatten(localScheme(j,i), j,i, m,n);
                end % for i
            end % for j
        end % getScheme2DEval
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % renvoie la matrice d'un schéma espace (1D)-temps (1D)
        %   Entrées : 
        %       * spaceScheme le schéma 1D espace
        %       * timeScheme le schéma 1D temps
        %   Sorties : 
        %       * Mat la matrice complète, x = espace, y=temps
        %   -----> x
        %   |
        % y v
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Mat = spaceTimeScheme (spaceScheme, timeScheme)
            Mat = repmat(spaceScheme, length(timeScheme(1,:)),1) + repmat(timeScheme, 1, length(spaceScheme(1,:)));
        end % spaceTimeScheme
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % schéma numérique par développement de Taylor
        %   Entrées : 
        %       * eval : vecteur ligne des points d'évaluation
        %           du schéma autour de x=0 (ex. [-h, 0, h])
        %       * derivative : ordre de la dérivée (f'(x) = 1, f''(x)=2...)
        %   Sorties : 
        %       * scheme : vecteur tel que 
        %   [f(x+a) for a in eval] \cdot scheme ~ f^(derivative)(x)
        %   !!! aucun test : si derivative trop grand, renverra 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function scheme = TaylorScheme (eval,derivative)
            len = length(eval);
            Taylor = zeros(len,len);
            Taylor (1,:) = 1.0;
            for idev=2:len
                Taylor(idev,:) = eval/(idev-1) .* Taylor(idev-1,:);
            end % for dev
            scdm = zeros(len,1);
            scdm (derivative+1) = 1;
            scheme = Taylor \ scdm;
        end % TaylorScheme
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % adds a frame of value of size number around array
        % mimics padarray(array, [dim, dim], number, 'both')
        %   Entrées : 
        %       * array une matrice
        %       * number : int, rayon du pad
        %       * value : valeur à mettre au bord
        %   Sorties : 
        %       * res la matrice "encadrée"
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = mypad (array, number, value)
            res = [value*ones(number,length(array(1,:))+2*number);
                   value*ones(length(array(:,1)),number), array, value*ones(length(array(:,1)),number);
                   value*ones(number,length(array(1,:))+2*number)];
        end % mypad
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot in space-time domain
        %   Entrées : 
        %       * Approx l'approximation courante
        %               ordonnée dimension de xx d'abord
        %       * xx, tt les maillages
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ploot (Approx,xx,tt)
            hold on
            colormap cool
            YY = reshape(Approx,length(xx),length(tt))'; % x en colonne, t en ligne
            surfc(xx,tt,YY);
            view(-222,24);
            xlabel('x');
            ylabel('t');
        end % ploot
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % affichage
        %   Entrées : 
        %       * listeh : liste des pas spatiaux
        %       * erreur : le tableau à afficher. Chaque ligne = une série
        %       * nom : liste de string, correspondant aux séries
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotErreur (listeh, erreur, noms)
            hold on
            grid on
            mk = ["-o", "-x", "-d", "-p", "-^"]; 
            LEGEND = string(zeros(1,2*length(noms)));
            for i=1:length(erreur(:,1))            
                p = polyfit(log10(listeh),log10(erreur(i,:)),1);
                loglog(listeh,erreur(i,:),mk(i),listeh,10.^(p(1).*log10(listeh)+p(2)),"-"+mk(i));
                LEGEND(2*i-1) = noms(i);
                LEGEND(2*i) = sprintf("a = %.3f",p(1));
            end
            xlabel ("pas spatial"); 
            ylabel (sprintf("Erreur relative"));
            title (sprintf("Erreur, échelle log-log"));
            legend(LEGEND,'Location','southeast');
        end % plotErreur
        
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
end % class utils