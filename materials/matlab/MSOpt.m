% Class for computing several optimization criteria for a multistratum
% design with a full quadratic model, given the settings matrix
%
% Requires the distribution of factors across strata (cell array), the
% number of units in each stratum (array), the number of levels for each
% factor (array, or one number for all factors), the ratios of error
% variance between subsequent strata (array), the combination of
% criteria to be optimized (array) and the model type (string)
% in the constructor:
%
% mso = MSOpt(facts,units,levels,etas,criteria,model)
%
% criteria can be a cell array with any combination of:
% 'I'  - I-optimality
% 'Id' - Id-optimality
% 'D'  - D-optimality
% 'A'  - A-optimality
% 'Ds' - Ds-optimality
% 'As' - As-optimality
%
% model can by any of 'main', 'interaction' or 'quadratic'
%
% Scoring function:
%
% scores = mso.Score(settings,criteria)
%

% Francesco Sambo, Dept. of Information Engineering, Univ. of Padova, 2013

classdef MSOpt < handle

   properties
      facts    % distribution of factors over strata
      nfacts   % number of factors
      nstrat   % number of strata
      units    % number of units in each stratum
      runs     % number of runs
      etas     % ratios of error variance between subsequent strata
      model    % model type, among 'main', 'interaction', 'quadratic'
      Vinv     % inverse of the covariance matrix of the responses
      avlev    % cell array of available levels for each factor
      levs     % number of levels for each factor
      crit     % criteria to be optimized
      ncrit    % number of criteria
      % optional, depending on the criteria
      M        % matrix of moments of the cube (I-opt)
      M0       % matrix of moments of the cube (Id-opt)
      W        % diagonal matrix of weights (As-opt)
   end

   methods

      % Constructor
      function mso = MSOpt(facts,units,levels,etas,criteria,model)
         mso.facts = facts;
         mso.nfacts = length([facts{:}]);
         mso.nstrat = length(facts);
         mso.units = units;
         mso.runs = prod(units);
         mso.etas = etas;
         mso.avlev = cell(1,mso.nfacts);
         % uniformly distribute available levels in the interval [-1,1]
         if isscalar(levels)
            mso.levs = ones(1,mso.nfacts)*levels;
            % same number of levels for each factor
            for i = 1:mso.nfacts
               mso.avlev{i} = 2*(0:levels-1)'/(levels-1) - 1;
            end
         else
            mso.levs = levels;
            for i = 1:mso.nfacts
               mso.avlev{i} = 2*(0:levels{i}-1)'/(levels{i}-1) - 1;
            end
         end
         V = eye(mso.runs);
         for i = 1:mso.nstrat-1
            V = V + etas(i) * ...
               kron( eye(prod(units(1:i))), ones(prod(units(i+1:end))) );
         end
         mso.Vinv = inv(V);

         mso.model = model;

         mso.crit = criteria;
         mso.ncrit = length(criteria);

         % optional part, depending on criteria
         if ismember( 'I', criteria )
            k = mso.nfacts; k2 = k*(k-1)/2;
            switch model
               case 'main'
                  mso.M = [1      zeros(1,k);
                      zeros(k,1)  eye(k)/3];
               case 'interaction'
                  mso.M = [1      zeros(1,k)  zeros(1,k2);
                      zeros(k,1)  eye(k)/3    zeros(k,k2);
                      zeros(k2,1) zeros(k2,k) eye(k2)/9];
               case 'quadratic'
                  mso.M = [1      zeros(1,k)  ones(1,k)/3                zeros(1,k2);
                      zeros(k,1)  eye(k)/3    zeros(k,k)                 zeros(k,k2);
                      ones(k,1)/3 zeros(k,k)  (4*eye(k)+5*(ones(k)))/45  zeros(k,k2);
                      zeros(k2,1) zeros(k2,k) zeros(k2,k)                eye(k2)/9];
               otherwise
                  error('error:model','Model type not valid');
            end
         end

         if ismember( 'Id', criteria )
            k = mso.nfacts; k2 = k*(k-1)/2;
            switch model
               case 'main'
                  mso.M0 = [1      zeros(1,k);
                      zeros(k,1)  eye(k)/3];
               case 'interaction'
                  mso.M0 = [1      zeros(1,k)  zeros(1,k2);
                      zeros(k,1)  eye(k)/3    zeros(k,k2);
                      zeros(k2,1) zeros(k2,k) eye(k2)/9];
               case 'quadratic'
                  mso.M0 = [1      zeros(1,k)  ones(1,k)/3                zeros(1,k2);
                      zeros(k,1)  eye(k)/3    zeros(k,k)                 zeros(k,k2);
                      ones(k,1)/3 zeros(k,k)  (4*eye(k)+5*(ones(k)))/45  zeros(k,k2);
                      zeros(k2,1) zeros(k2,k) zeros(k2,k)                eye(k2)/9];
               otherwise
                  error('error:model','Model type not valid');
            end
            mso.M0(:,1) = 0;
            mso.M0(1,:) = 0;
         end

         if ismember( 'As', criteria )
            switch model
               case 'main'
                  w = ones( mso.nfacts, 1 );
               case 'interaction'
                  w = ones( mso.nfacts + mso.nfacts * (mso.nfacts - 1) / 2, 1);
               case 'quadratic'
                  w = [ ones( mso.nfacts, 1 ); ones( mso.nfacts, 1) / 4;
                     ones( mso.nfacts * (mso.nfacts - 1) / 2, 1) ];
               otherwise
                  error('error:model','Model type not valid');
            end
            mso.W = diag( w/sum(w) );
         end

      end

      % Scores
      function scores = Score(obj,settings)
         switch obj.model
            case 'main'
               X = [ones(obj.runs,1) settings];
            case 'interaction'
               X = [ones(obj.runs,1) settings colprod(settings)];
            case 'quadratic'
               X = [ones(obj.runs,1) settings settings.^2 colprod(settings)];
         end
         B = X'*obj.Vinv*X;
         determ = det(B);
         scores = Inf(size(obj.crit));
         if rcond(B) > 1e-5 && determ > 0

            ind = strcmp('D', obj.crit);
            if any(ind)
               scores(ind) = 1 / determ^(1/size(X,2));
            end

            if any(ismember({'I','Id','Ds','A','As'}, obj.crit))
               Binv = inv(B);
            end

            ind = strcmp('I', obj.crit);
            if any(ind)
               scores(ind) = sum(diag(Binv * obj.M));
            end

            ind = strcmp('Id', obj.crit);
            if any(ind)
               scores(ind) = sum(diag(Binv * obj.M0));
            end

            ind = strcmp('A', obj.crit);
            if any(ind)
               scores(ind) = sum(diag(Binv))/size(X,2);
            end

            ind = strcmp('Ds', obj.crit);
            if any(ind)
               scores(ind) = det(Binv(2:end,2:end)) ^ (1/(size(X,2)-1));
            end

            ind = strcmp('As', obj.crit);
            if any(ind)
               scores(ind) = sum(diag(obj.W * Binv(2:end,2:end)));
            end
         end
      end

   end
end


