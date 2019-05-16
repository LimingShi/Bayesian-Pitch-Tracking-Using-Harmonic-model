classdef BayesianfastF0NLS< handle

%     usuage

% obj= BayesianfastF0NLS(N, L, pitchBounds, A_var,voicingProb,varargin)
%  N is the frame length 
%  L is the maximum allowed harmonic order
%  pitchBounds denotes the pitch range e.g.,[f_min/fs, f_max/fs], where fs is
%  the sampling frequency
%  A_var denotes the frequncy transition variance.
%  voicingProb denotes the initial voicing probability
    
    
    
    % Can only be set in the constructor
    properties (SetAccess=immutable) 
        N
        L
        F
        pitchBoundsOuter = [0.0 0.5] % The outer bounds that are acceptable
        
        % default values
        epsilon = 0.0
        dcIsIncluded = false
        epsilon_ref = 0.0
        
        % Precomputed quantities
        crossCorrelationVectors
        Gamma1
        Gamma2
        fftShiftVector
    end

    properties (SetAccess=private)
        pitchBounds  % The current active bounds 
        fullPitchGrid
        validFftIndices
        defaultRefinementTol
        refinementTol
        gPriorParam = 3
        logPitchPdfs % the pitch pdf (column) for every candidate model order
        logModelPmf % an (L+1)x1 vector
         logpi
        A
        scaled_alpha_buffer
        scaled_alpha_buffer2
        unvoicing_scaled_alpha_buffer
        logModelPrior;
        varSpeech=5e-3;
        B
        C
        norml_factor=0;
        cnt=1;
        prew_state
    end

    methods

        
        function obj = BayesianfastF0NLS(N, L, pitchBounds, A_var,voicingProb,varargin)
            switch nargin
                case 0
                     N=0.03*16000;
                     L=10;
                     pitchBounds=[50/16000,400/16000];
                     A_var=[2/16000];
                     voicingProb=0.7;                    
                case 1
                     L=10;
                     pitchBounds=[50/16000,400/16000];
                     A_var=[2/16000];
                     voicingProb=0.7;
                case 2
                    pitchBounds=[50/16000,400/16000];
                     A_var=[2/16000];
                     voicingProb=0.7;
                case 3
                     A_var=[2/16000];
                     voicingProb=0.7;
                case 4
                     voicingProb=0.7;               
            end

            % validate input 
            if length(varargin) >= 1
                if islogical(varargin{1})
                    obj.dcIsIncluded = varargin{1};
                else
                    error('Argument 4 is not of type logical (true/false)')
                end
            end
            if length(varargin) >= 2
                if isscalar(varargin{2})
                    obj.F = varargin{2};
                else
                    error('Argument 5 is not a scalar')
                end
            else
                obj.F = 2^14;
            end
            obj.defaultRefinementTol = 1/obj.F;
            
            if length(varargin) >= 3
                if isscalar(varargin{3})
                    obj.epsilon = varargin{3};
                else
                    error('Argument 6 is not a scalar')
                end
            end
            
            if length(varargin) > 4
                error('Too many input arguments')
            end
            
            if ~isscalar(N)
                error('Input argument N is not a scalar')
            else
                obj.N = N;
            end
            if ~isscalar(L)
                error('Input argument L is not a scalar')
            else
                obj.L = L;
            end

            if ~ (isvector(pitchBounds) && length(pitchBounds) == 2) 
                error(['Input argument pitchBound must be a 2-' ...
                       'vector'])
            elseif pitchBounds(1) < obj.pitchBoundsOuter(1) || ...
                pitchBounds(2) > obj.pitchBoundsOuter(2)
                
                error(['Input argument pitchBounds must be within the ' ...
                       'bounds specified in the constructor (at ' ...
                       'least [0.0 0.5])'])
            else
                obj.pitchBounds = pitchBounds;
            end
            
            if pitchBounds(1) < 1/N
                warning(['The lower pitch bound is set lower than one '...
                    ' period/segment. Inaccurate results might be '...
                    'produced - especially if you do not set the'...
                    ' regularisation parameter.']);
            end

            % Init
            F = obj.F;
            minFftIndex = ceil(F*pitchBounds(1));
            maxFftIndex = floor(F*pitchBounds(2));
            obj.validFftIndices = (minFftIndex:maxFftIndex)';
            obj.fullPitchGrid = obj.validFftIndices/F;
            nPitches = length(obj.fullPitchGrid);
            obj.logModelPrior=[voicingProb,1-voicingProb];
            
            % cross-correlation vectors
            obj.crossCorrelationVectors = ...
                [N*ones(1, nPitches)/2 + N*obj.epsilon;...
                sin(pi*(1:2*L)'*obj.fullPitchGrid'*N)./...
                (2*sin(pi*(1:2*L)'*obj.fullPitchGrid'))];

            obj.fftShiftVector = ...
                exp(1i*2*pi*(0:ceil(F/2)-1)'*(N-1)/(2*F));

            % Compute Gamma for the T+H and T-H systems
            [obj.Gamma1, obj.Gamma2] = computeGamma(L, F, pitchBounds, ...
                obj.crossCorrelationVectors, nPitches, ...
                obj.validFftIndices, ...
                obj.dcIsIncluded);
            %%          added by Liming Shi
            obj.scaled_alpha_buffer=nan;
            grid_width=diff(obj.fullPitchGrid(1:2));
%             obj.logpi=(log(1/(nPitches-1))-log(grid_width))*ones(nPitches,1); 
            obj.logpi=log((voicingProb)/(nPitches*L))*ones(nPitches,1); 
            obj.A=zeros(nPitches,nPitches);
            
                   
            for jj=1:nPitches
%                 f0=obj.fullPitchGrid(jj);
%                  obj.A(jj,:)= abs(sin(pi*(obj.fullPitchGrid-f0)/f0)./pi./(obj.fullPitchGrid-f0)*f0)+.1;  
                 obj.A(jj,:)=(normpdf(obj.fullPitchGrid,obj.fullPitchGrid(jj),A_var));
%                  var_value=1e-1;
%                  obj.logA(jj,:)=log(lognpdf(obj.fullPitchGrid,log(obj.fullPitchGrid(jj)*exp(var_value)),var_value));
%                    obj.logA(jj,:)=log(normpdf(log(obj.fullPitchGrid),log(obj.fullPitchGrid(jj)),.1));
            end
            
            
            obj.A=obj.A./repmat(sum(obj.A,2),1,nPitches);
            
            obj.B=zeros(obj.L,obj.L);
            
                   
            for jj=1:obj.L
                 obj.B(jj,:)=normpdf(jj,[1:obj.L],1);
            end
            
            
            obj.B=obj.B./repmat(sum(obj.B,2),1,obj.L);
            
            obj.C=[.7 .3;.4,.6];
            
            obj.scaled_alpha_buffer2=1/obj.L/nPitches*ones(nPitches,obj.L);
            
        end
        
        function varargout = computeCostFunctions(obj, x)
            % Compute and returns the cost functions for l=1,...,L

            % validate input 
            if ~isvector(x) && length(x) == obj.N
                error(['First argument x must be vector of ' ...
                               'length N=', num2str(obj.N)]);
            end
            x = reshape(x, obj.N, 1); % make sure its a column
                                      % vector

                
            varargout{1} = computeAllCostFunctions(x, obj.L, ...
            	obj.fullPitchGrid, obj.fftShiftVector, ...
                obj.crossCorrelationVectors, obj.Gamma1, obj.Gamma2, ...
                obj.dcIsIncluded);
            
            if nargout == 2
                varargout{2} = obj.fullPitchGrid;
            end
            
        end
        
        function varargout = estimate(obj, x, flag,tinc,varargin)
            % Estimate fundamental frequency and order of the signal x
                
            % validate input 
            if ~isvector(x) && length(x) == obj.N
                error(['First argument x must be vector of ' ...
                               'length N=', num2str(obj.N)]);
            end
            x = reshape(x, obj.N, 1); % make sure its a column
                                      % vector
            
            if length(varargin) == 1
                if isscalar(varargin{1})
                    obj.refinementTol = varargin{1};
                else
                    error('Argument 2 is not a scalar')
                end
            else
                obj.refinementTol = obj.defaultRefinementTol;
            end
            if flag==1
                spec_pow_all=(abs(fft(x,2048)).^2)';
                if isempty(obj.prew_state)
                    [psdouput,obj.prew_state]=estnoiseg_orig(spec_pow_all,tinc);
                else
                    [psdouput,obj.prew_state]=estnoiseg_orig(spec_pow_all,obj.prew_state);
                end
                mmse_covariance_all = (real(ifft(psdouput')));
                mmse_lpc_all= levinson(mmse_covariance_all,30);
                [xp] = filter(mmse_lpc_all,1,...
                    x);
                  x=xp/norm(xp)*norm(x);
            end
            
            
            
            
            % if the data consists of all zeros then return null model 
            if  x'*x/obj.N<1e-10
                estimatedPitch = nan;
                estimatedOrder = 0;
                unvoicing_scaled_alpha=log(1);
            else
                % Step 2: compute the profile log-likelihood function 
                % over only the fundamental frequency (the linear
                % parameters, noise variance and g have been integrated
                % out)
                costs = obj.computeCostFunctions(x);
                cod = costs*(1/(x'*x+obj.varSpeech));
%                   cod = costs*(1/(x'*x));
                [~, pitchLogLikelihood] = ...
                    computePitchLogLikelihood(cod, obj.N, obj.gPriorParam);
                null_modellike=1;
                
                    %% added/modified by Liming Shi (implements the hidden markov chain)
                if isnan(obj.scaled_alpha_buffer)                
                    inx=isnan(pitchLogLikelihood);
                    pitchLogLikelihood(inx)=-inf;
                    bar_alpha=obj.logpi+pitchLogLikelihood';
                    
                    unvoicing_bar_alpha=log(obj.logModelPrior(2)*null_modellike);
                    log_scale=log_sumsum_exp_ls([[unvoicing_bar_alpha,-inf*ones(1,obj.L-1)];bar_alpha]);
                    scaled_alpha=bar_alpha-log_scale;
                    unvoicing_scaled_alpha=unvoicing_bar_alpha-log_scale;
 
                else
                    inx=isnan(pitchLogLikelihood);
                    pitchLogLikelihood(inx)=-inf;
                    state_prior=obj.C(1,1)*obj.A'*obj.scaled_alpha_buffer*obj.B;

                    state_prior=state_prior+obj.scaled_alpha_buffer2*obj.C(2,1)*obj.unvoicing_scaled_alpha_buffer;

                    bar_alpha=log(state_prior)+pitchLogLikelihood';
                    
                    temp=[1-obj.unvoicing_scaled_alpha_buffer;obj.unvoicing_scaled_alpha_buffer];

                    unvoicing_bar_alpha=log(obj.C(:,2)'*temp*null_modellike);
                    log_scale=log_sumsum_exp_ls([[unvoicing_bar_alpha,-inf*ones(1,obj.L-1)];bar_alpha]);
                    scaled_alpha=bar_alpha-log_scale;
                    unvoicing_scaled_alpha=unvoicing_bar_alpha-log_scale;
                    
                end

                    [value, inx] = ...
                        max(scaled_alpha);
                    [~,inx_2]=max(value);
                    pitchIndex=inx(inx_2);
                    estimatedOrder=inx_2;
                    estimatedPitch = obj.fullPitchGrid(pitchIndex(1));
             %% added by Liming Shi (buffering the previous meaningful estimates of the posterior states)


                if unvoicing_scaled_alpha < log(0.5)
                    obj.scaled_alpha_buffer2=exp(scaled_alpha-log_sumsum_exp_ls(scaled_alpha));
                else 
%                     estimatedOrder=nan;
                end
                obj.scaled_alpha_buffer=exp(scaled_alpha);
                obj.unvoicing_scaled_alpha_buffer=exp(unvoicing_scaled_alpha);
            end

            %%
                
            if nargout >= 1
                varargout{1} = estimatedPitch;
            end
            if nargout >= 2
                varargout{2} = estimatedOrder;
            end
            if nargout >= 3
%                 [~, alpha] = objFunction(estimatedPitch, x, ...
%                     estimatedOrder, obj.dcIsIncluded);
                varargout{3} = 1-exp(unvoicing_scaled_alpha);
            end
            
        end
        function obj=destroy(obj)
        end
    end
end
function y=log_sumsum_exp_ls(x)
max_temp=max(max(x));
inf_ind=isinf(max_temp);
y=log(sum(sum(exp(x-max_temp))))+max_temp;
y(inf_ind)=-inf;
end

