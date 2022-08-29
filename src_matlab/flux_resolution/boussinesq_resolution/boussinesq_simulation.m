classdef boussinesq_simulation
% class storing the different object in order to make run a boussinesq simulation for different conditions
    properties(Access=public)
        discretization          % space_discretization object contains also the spatial properties of your hillslopes (x_S, x_Q, width function w, soil_depth, angle, f, k) 
        source_terms            % source object containing the time properties and the recharge related to it
        boundary_cond           % boundary_conditions object containing the time of boundary on xmin & xmax (Q fixed, S fixed or free condtions)
        initial_conditions      % initial conditions for S, Q & QS
        ratio_P_R               % ratio Recharge over Precipitation
        sol_simulated           % contains all the information provided by ode15s
    end
    
    methods(Access=public)
        % constructor
        function obj=boussinesq_simulation
        end
        
        function obj=simulation_parametrization(obj,discretization,source_terms,boundary_cond,percentage_loaded,t_initial,ratio_P_R,Sinitial)
            if(nargin<7)
                obj.ratio_P_R=1;
            end
            if(nargin<8)
                Sinitial=nan;
            end
            obj.discretization=discretization;
            obj.source_terms=source_terms;
            obj.ratio_P_R=ratio_P_R;
            obj.boundary_cond=boundary_cond;
            obj.initial_conditions=obj.set_initial_conditions(percentage_loaded,t_initial,Sinitial);
        end
        
        function initial_conditions=set_initial_conditions(obj,percentage_loaded,t_initial,Sinitial)
            if(nargin<4)
                Sinitial=nan;
            end
            if(percentage_loaded==-1)
                 Sin_name=input('Enter the name where the initial soil moisture you want is stored: \n');
                 Sin_directory=which(Sin_name);
                 Sin=load(Sin_directory);
                 initial_conditions.Sin=Sin.Sin;
            elseif(percentage_loaded==-2 && ~isnan(sum(Sinitial)))
                % assign only the initial storage values
                if(length(Sinitial)==obj.discretization.Nx)
                    initial_conditions.Sin=Sinitial;
                 % possible to directy assign all the values of the simulation S, Q and QS to prevent error with initial values
                elseif(length(Sinitial)==3*obj.discretization.Nx+1)
                    initial_conditions.Sin=Sinitial(1:obj.discretization.Nx);
                    initial_conditions.Qin=Sinitial(obj.discretization.Nx+1:2*obj.discretization.Nx+1);
                    initial_conditions.QSin=Sinitial(2*obj.discretization.Nx+2:3*obj.discretization.Nx+1);
                else
                    fprintf('Error: length of the assigned initial state values of the simulation are not of the good length \n');
                end
            else
                [~,w,soil_depth,~,~,f]=obj.discretization.get_resampled_variables;
                Smax=f.*w.*soil_depth;
                if(percentage_loaded==-2 && isnan(Sinitial))
                    initial_conditions.Sin=0*Smax;
                else
                    initial_conditions.Sin=percentage_loaded*Smax;
                end
            end
            
            
            Edges=obj.boundary_cond.fixed_edge_matrix_values;
            Edges_bool=obj.boundary_cond.fixed_edge_matrix_boolean;
            initial_conditions.Sin(1)=Edges_bool(1)*Edges(1)+(1-Edges_bool(1))*initial_conditions.Sin(1); initial_conditions.Sin(end)=Edges_bool(2)*Edges(2)+(1-Edges_bool(2))*initial_conditions.Sin(end);
%             if(~isnan(sum(Sinitial)) && length(Sinitial)==3*obj.discretization.Nx+1)
%                 initial_conditions.Qin(1)=Edges_bool(3)*Edges(3)+(1-Edges_bool(3))*initial_conditions.Qin(1); initial_conditions.Qin(end)=Edges_bool(4)*Edges(4)+(1-Edges_bool(4))*initial_conditions.Qin(end);
%             else
%                 initial_conditions.Qin=-obj.compute_Q_from_S(initial_conditions.Sin)*initial_conditions.Sin;
%                 initial_conditions.Qin(1)=Edges_bool(3)*Edges(3)+(1-Edges_bool(3))*initial_conditions.Qin(1); initial_conditions.Qin(end)=Edges_bool(4)*Edges(4)+(1-Edges_bool(4))*initial_conditions.Qin(end);
%                 initial_conditions.QSin=obj.compute_QS_from_Q([initial_conditions.Sin;initial_conditions.Qin;zeros(size(initial_conditions.Qin))],t_initial)*initial_conditions.Qin;
%             end
        end

        % set the parameters for the resolution of the Boussinesq DAE
        function [t,x_center,x_edge,S,Q,QS,obj]=implicit_scheme_solver(obj,t_properties,solver_options)
            % initial vector for the DAE
            y0=[obj.initial_conditions.Sin];
           
            % Number of unknowns
%             Nun=3*obj.discretization.Nx+1;
            block_size=obj.discretization.Nx;
            
            % Function to integrate
            f=@(t,y)(obj.odefun(y,t));
            
             % time range to solve
            [time_range,unit]=t_properties.get_properties;
            
            % Initial slope of the DAE
            yp0=f(time_range(1),y0);
            
            options = odeset( 'InitialSlope', yp0);
            if(solver_options.Events==1)
                 event=@(t,y)obj.eventfun(t,y);
                 solver_options.Events=event;
            end
            if(~isempty(solver_options.AbsTol) && isnan(solver_options.AbsTol))
                solver_options.AbsTol=[];
            end
            if(~isempty(solver_options.RelTol) && isnan(solver_options.RelTol))
                solver_options.RelTol=[];
            end
            if(~isempty(solver_options.MaxStep) && isnan(solver_options.MaxStep))
                solver_options.MaxStep=[];
            end
            options=odeset(options,solver_options);
           
            % get consistent initial conditions thanks to decic
            [y0new,yp0new] = decic(@(t,y,yp) yp - (obj.odefun(y,t)), time_range(1), y0, [], yp0, [], options);
            change_init_slope=odeset('InitialSlope',yp0new);%,'Nonnegative',ones(length(y0new),1));
            options=odeset(options,change_init_slope);
            
            % DAE System to solve
%             time_range=[time_range(1) time_range(end)];
%             time_range=time_range(1):3600:time_range(end);
%             [t,S]=ode15s(f,time_range,y0,options);
            
            % for computation speed reason if not specified let matlab decide for the interior points computation in ode15s
            if(~isempty(options.Refine) && options.Refine==-1)
                options.Refine=[];
                obj.sol_simulated=ode15s(f,time_range,y0new,options);
            else
                options.Refine=[];
                obj.sol_simulated=ode15s(f,[time_range(1) time_range(end)],y0new,options);
            end
            
            t=time_range;
            if(time_range(end)>obj.sol_simulated.x(end))
                tmax_pos=find(time_range-obj.sol_simulated.x(end)<=0,1,'last');
                t=time_range(1:tmax_pos);
            end
            S=deval(obj.sol_simulated,t);
            size_S=size(S);
            Q=nan(size_S(1)+1,size_S(2));
            QS=nan(size_S);
            
            for i=1:length(t)
                Q(:,i)=-obj.compute_Q_from_S(S(:,i))*S(:,i);
                QS(:,i)=obj.compute_QS_from_Q(S(:,i),t(i))*Q(:,i)+obj.partition_source_terms_QS(S(:,i),t(i)); %2*block_size+2:3*block_size+1,:);
            end
            
            x_edge=obj.discretization.get_edges_coordinates;
            x_center=obj.discretization.get_center_coordinates;

%             sol=ode15s(f,time_range,y0,options);
%             t=(sol.x);
%             S=(sol.y);
%             Q=S(block_size+1:2*block_size+1,:);
%             QS=S(2*block_size+2:3*block_size+1,:);
%             S=S(1:block_size,:);
%             x_edge=obj.discretization.get_edges_coordinates;
%             x_center=obj.discretization.get_center_coordinates;
           
        end
        
        % Compute Direct Precipitations on Saturated Areas (DPSA) and Return Flow (RF) in m2/s aggregated 
        % at the hillslope scale and also their spatialized counterparts (at the scale of a discretized elements
        function [DPSA_spat,RF_spat]=partition_DPSA_RF(obj,S,Q,t) % unit m2/s
            block_size=obj.discretization.Nx;
            % return flow
            RF_spat=obj.compute_QS_from_Q(S,t)*Q;
            % direct precipitations onto saturated areas
%             Source_Partitioned=obj.partition_source_terms(S,t);
%             Recharge_rate_spatialized=obj.compute_source_term_spatialized(S,t);
%             DPSA_spat=Recharge_rate_spatialized-Source_Partitioned;
            DPSA_spat=obj.partition_source_terms_QS(S,t);
        end
    end
    
    methods(Access=private)
        
        % Computes array at a given y
        function dy=odefun(obj,y,t)
            %% some values necessary to compute
            block_size=obj.discretization.Nx; % size of the matrix (corresponds to the number of discretized elements)
            Edges=obj.boundary_cond.fixed_edge_matrix_boolean; % kind of boundary conditions
            [~,w,d,angle,~,f,k,f_edges]=obj.discretization.get_resampled_variables; % properties of the hillslope
            A=obj.discretization.A; % derivation matrix 1
            B=obj.discretization.B; % derivation matrix 2
            %% compute recharge rate spatialized
            [Recharge_rate_spatialized,Threshold]=obj.compute_source_term_spatialized(y,t,w,d,f);
            %% Compute darcy flux from one box to another with variable angle
            Q_from_S=obj.compute_Q_from_S(y,w,angle,f,k,f_edges);
            % compute Q spatial derivative
            dQ_from_S=A*Q_from_S;
            %% compute test if potential seepage
            Test_Deriv=dQ_from_S*y+Recharge_rate_spatialized;
            Test_Deriv=Test_Deriv>=0;
            alpha=Threshold.*Test_Deriv+(1-Test_Deriv); % regularization function : drives where goes variations of mass to S or to QS
            %% compute C
            C    = sparse(1:block_size,1:block_size,alpha)*dQ_from_S;
            C(1,:)=(1-Edges(1))*C(1,:);
            C(block_size,:)=(1-Edges(2))*C(block_size,:);
            %% compute D
            D=alpha.*Recharge_rate_spatialized;
            D(1,:)=(1-Edges(1))*D(1,:);
            D(block_size,:)=(1-Edges(2))*D(block_size,:);
            D=sparse(D); % ###check necessity
            %% compute dy
            dy=C*y+D;
        end
        
        function D=partition_source_terms_QS(obj,y,t)
            alpha_complementar=1-obj.alpha(y,t);
            Recharge_rate_spatialized=obj.compute_source_term_spatialized(y,t);
            D=alpha_complementar.*Recharge_rate_spatialized;
           
            D=sparse(D);
        end
        
        % Darcy equation
        function Q_from_S=compute_Q_from_S(obj,y,w,angle,f,k,f_edges)
            if(nargin<3)
                [~,w,~,angle,~,f,k,f_edges]=obj.discretization.get_resampled_variables;
            end
            block_size=obj.discretization.Nx;
            B=obj.discretization.B;
            Q_from_S=obj.discretization.Omega; % Omega is directly in Q_from_S to gain speed
            % Compute darcy flux from one box to another with variable angle
            Q1_from_S=k./f_edges.*cos(angle).*(B*(y./(f.*w)));
            Q_from_S=sparse(1:block_size+1,1:block_size+1,Q1_from_S)*Q_from_S;%sparse(diag(Q1_from_S))*Q_from_S;%=bsxfun(@times,Q1_from_S,Omega);%Q1_from_S.*Omega;
            if(sum(angle)~=0)
                Q1_from_S=k./f_edges.*sin(angle); Q1_from_S=Q1_from_S(1:end-1);%Q1_from_S(end)=0;
                Q1_from_S=sparse(1:block_size,1:block_size,Q1_from_S,block_size+1,block_size);%sparse(1:block_size+1,1:block_size+1,Q1_from_S)*D;%sparse(diag(Q1_from_S))*D;%bsxfun(@times,Q2_from_S,D);%Q2_from_S.*D;%
                Q_from_S=Q_from_S+Q1_from_S;
            end
            % put boundary conditions on Q in the matrix
            Edges_BC=obj.boundary_cond.fixed_edge_matrix_boolean;
            Q_from_S(1,:)=(1-Edges_BC(3))*Q_from_S(1,:);
            Q_from_S(block_size+1,:)=(1-Edges_BC(4))*Q_from_S(block_size+1,:);
        end
        
        function QS_from_Q=compute_QS_from_Q(obj,y,t)
            QS_from_Q=-sparse(diag(1-obj.alpha(y,t)))*obj.discretization.A;
        end
        
        function [Recharge_rate_spatialized,Threshold]=compute_source_term_spatialized(obj,y,t,w,d,f)
            if(nargin<4)
                [~,w,d,~,~,f]=obj.discretization.get_resampled_variables;
            end
            relative_occupancy_rate=y./(f.*w.*d);
            Threshold=threshold_function(relative_occupancy_rate);
            Recharge_rate=obj.source_terms.compute_recharge_rate(t);
            if(Recharge_rate<0)
                Recharge_rate=Recharge_rate.*[1;(1-exp(-10*((y(2:end)-f(2:end).*w(2:end).*d(1))./(f(2:end).*w(2:end).*(d(2:end)-d(1))))))];%                Recharge_rate=Recharge_rate.*(1-exp(-1000*(1-min(1,-Recharge_rate*3600./(y./(f.*w)-d(1))))));
            end
            
            Recharge_rate_spatialized=Recharge_rate.*w;
            Recharge_rate_spatialized=(Threshold*(obj.ratio_P_R-1)+1).*Recharge_rate_spatialized;
            % if there is an ETP time series given, compute ETR from ETP, Recharge and S/Smax
            if(logical((~isnan(obj.source_terms.ETP_chronicle)).*(~isempty(obj.source_terms.ETP_chronicle))))
                ETP_rate=obj.source_terms.compute_ETP_rate(t);
                ETP_rate_spatialized=ETP_rate.*w;
                ETR_rate_spatialized=(Recharge_rate_spatialized-ETP_rate_spatialized>0).*(ETP_rate_spatialized)+...
                    (Recharge_rate_spatialized-ETP_rate_spatialized<=0).*(Recharge_rate_spatialized+(ETP_rate_spatialized-Recharge_rate_spatialized).*(1-exp(-5*relative_occupancy_rate))); 
                Recharge_rate_spatialized=Recharge_rate_spatialized-ETR_rate_spatialized;
            end
        end
        
        function Recharge_deep_spatialized=compute_deep_recharge(obj,y,k,f)
            block_size=obj.discretization.Nx;
            Recharge_deep_spatialized=(k/100)*(y(1:block_size)/(100*f));
            Recharge_deep_spatialized(1)=0;
        end
        
        function OUT=Test_Derivative(obj,y,t)
            A=obj.discretization.A;
            Recharge_rate_spatialized=obj.compute_source_term_spatialized(y,t);
            OUT=-A*(-obj.compute_Q_from_S(y)*y)+Recharge_rate_spatialized;
            OUT=OUT>=0;
        end
        
        function OUT=alpha(obj,y,t)
            block_size=obj.discretization.Nx;
            [~,w,soil_depth,~,~,f]=obj.discretization.get_resampled_variables;
            
            Thresh=threshold_function(y(1:block_size)./(f.*soil_depth.*w));
            Test_Deriv=obj.Test_Derivative(y,t);
            OUT=Thresh.*Test_Deriv+(1-Test_Deriv);
        end
        
        % Mass matrix: Identity on S, O for Q & QS
        function M=compute_mass_matrix(obj)
            block_size=obj.discretization.Nx;
            Nunknowns=3*block_size+1;
            M=zeros(Nunknowns,Nunknowns);
            M(1:block_size,1:block_size)=eye(block_size);
            M=sparse(M);
        end
        
        function [x,isterm,dir] = eventfun(obj,t,y)
            block_size=obj.discretization.Nx;
           
            dy=obj.odefun(y,t); 
            
% %             [~,w]=obj.discretization.get_resampled_variables;
% %             D=obj.source_terms.recharge_mean*w;%obj.partition_source_terms(y,t);
% %             dy=abs(dy./D(1:block_size));
% %             dy=dy(D(1:block_size)>0);
% %             x=max(dy)-5e-11;
            
            [~,w,soil_depth,~,~,f]=obj.discretization.get_resampled_variables;
            dy=dy./(f.*w.*soil_depth);
            x = max(abs(dy)) - 5e-11;%5e-11; %5e-8; %#JM think to change it 5e-8
%             x = max(abs(dy)) - 5e-11;
            isterm = 1;
            dir = 0;  %or -1, doesn't matter
        end
        
        function dS=compute_storage_variation(obj,t,y)
            block_size=obj.discretization.Nx;
            dS=nan(block_size,length(t));
            for j=1:length(t)
                dy=obj.odefun(y(:,j),t(j));
                dS(:,j)=dy(1:block_size);
            end
        end
    end
    
    %% Test functions
    methods(Static)    
        function [obj,Mass_balance,Flux_in,Darcy_Flux_out,Seepage_Flux_out,Storage_Var]=synthetic_experiment(xcustom,Date,Pluvio)
            if(nargin<1) xcustom=-1; end
            if(nargin<2) Date=-1; end
            if(nargin<3) Pluvio=-1; end
            
%             obj=boussinesq_simulation(xcustom);
            obj=boussinesq_simulation;
            % Space discretization
            % Possible to load space discretization log or lin or custom 
            Nx=100; % Number of space steps
            obj=obj.set_space_discretization(Nx,'lin');
%             obj=obj.set_space_discretization(Nx,'custom',xcustom);


            if(Date==-1)
                % Time properties choice tmin, tmax
                Nt=4001;
                t=time_properties(0,40,Nt,'day');
                
                % First setting possibility: Drainage virtual experiment like in Troch paper
                recharge=0;
                obj=obj.set_source_terms(t,recharge);
                % For the initial conditions
                percentage_loaded=0.30;
            
%                 % Second Setting: Recharge virtual experiment like in Troch experiment
%                 recharge=1e-2/(24*3600);
%                 % Steady source terms (like Troch experiment)
%                 obj=obj.set_source_terms(t,recharge,'steady');
%                 % For the initial conditions
%                 percentage_loaded=0;
                
%                % Third setting: Recharge increased and source term variable (sinusoidal)
%                % Recharge Rate to get outflow
%                recharge=10e-2/(24*3600);
%                % Source varies periodically with a period of 10 days
%                obj=obj.set_source_terms(t,recharge,'periodical',10);
%                % For the initial conditions
%                percentage_loaded=0;
            else
                t=time_properties(Date(1),Date(end),length(Date),'day');
                % convert the pluvio from mm/d to m/s
                Pluvio=Pluvio*1e-3/(24*3600);
                obj=obj.set_data_based_source_terms('data_based',t,Pluvio);
                % For the initial conditions
                percentage_loaded=0;
            end
            
            % boundary conditions filled as a 2 value type array in a 2 value array for xmin & xmax
            % negative values mean Q directed minus x axis
            obj=obj.set_boundary_conditions({'S','Q'},[0,0]);
            obj=obj.set_initial_conditions(percentage_loaded,t_initial);
            
            % Run the DAE
            obj=obj.implicit_scheme_solver(t);
            
            % Plot some results
            obj.plot_results;
            [Mass_balance,Flux_in,Darcy_Flux_out,Seepage_Flux_out,Storage_Var]=obj.compute_mass_changes;
        end
    end
    
    methods(Access=public)     
        % infiltration (Flux_in_spat) in m2/s
        function Flux_in_spat=compute_source_term_spatialized_vectorized(obj,t,y)
            block_size=obj.discretization.Nx;
            Flux_in_spat=nan(block_size,length(t));
            for j=1:length(t)
                Flux_in_spat(:,j)=obj.compute_source_term_spatialized(y(:,j),t(j));
            end
        end
        
        function N_in=compute_recharge_total(obj,t)
            N_in=obj.source_terms.compute_recharge_rate(t);
            [~,w,~,~]=obj.discretization.get_resampled_variables;
            dx=obj.discretization.compute_dx;
            % (2:end) so we do not consider rain falling in the river: i.e. first point in w
            N_in=N_in*sum(w(2:end).*dx(2:end));
        end
        
        function N_in=compute_recharge(obj,t)
            N_in=obj.source_terms.compute_recharge_rate(t);
        end
        
        function ETR_OUT=compute_ETR_OUT(obj,t,S)
            if(logical((~isnan(obj.source_terms.ETP_chronicle)).*(~isempty(obj.source_terms.ETP_chronicle))))
                [~,w,d,~,~,f]=obj.discretization.get_resampled_variables;
                Recharge_rate=obj.source_terms.compute_recharge_rate(t);
                Recharge_rate_spatialized=Recharge_rate.*w;
                relative_occupancy_rate=bsxfun(@rdivide,S,f*w.*d);
                Thresh=threshold_function(relative_occupancy_rate);
                Recharge_rate_spatialized=(Thresh*(obj.ratio_P_R-1)+1).*Recharge_rate_spatialized;
                ETP_rate=obj.source_terms.compute_ETP_rate(t);
                ETP_rate_spatialized=ETP_rate.*w;
                ETR_OUT=(Recharge_rate_spatialized-ETP_rate_spatialized>0).*(ETP_rate_spatialized)+...
                    (Recharge_rate_spatialized-ETP_rate_spatialized<=0).*(Recharge_rate_spatialized+(ETP_rate_spatialized-Recharge_rate_spatialized).*(1-exp(-5*relative_occupancy_rate)));
            else
                ETR_OUT=nan;
            end
        end
        
        function dS=compute_raw_storage_variation(obj,t,y)
            if(nargin<2)
                dS=obj.compute_storage_variation(obj.sol_simulated.x,obj.sol_simulated.y);
            else
                dS=obj.compute_storage_variation(t,y);
            end
        end
        
        function [mass_balance,Flux_in_spatialized,Darcy_Flux_difference_spatialized,Seepage_Flux_spatialized,dS]=compute_raw_mass_balance(obj)
            block_size=obj.discretization.Nx;
            t=obj.sol_simulated.x;
            y=obj.sol_simulated.y;
            Q=y(block_size+1:2*block_size+1,:);
            Darcy_Flux_difference_spatialized=Q(1:end-1,:)-Q(2:end,:);
            qS=y(2*block_size+2:3*block_size+1,:);
            x_Q=obj.discretization.get_edges_coordinates;
            dx=x_Q(2:end)-x_Q(1:end-1);
            dx_QS=diag(dx);
            Seepage_Flux_spatialized=dx_QS*qS;
            dS=obj.compute_raw_storage_variation(t,y);
            Flux_in_spatialized=obj.compute_source_term_spatialized_vectorized(t,y);
            mass_balance=Flux_in_spatialized+Darcy_Flux_difference_spatialized-Seepage_Flux_spatialized-dS;
        end
        
        function error=save_error_file(obj,tmax,foldername)
            error=0;
            if(obj.sol_simulated.x(end)<tmax)
                error=1;
                filename_err=fullfile(foldername,'Error.err');
                fid = fopen(filename_err, 'w');
                fprintf(fid, '1');
                fclose(fid);
            end
        end
        
        function save_info_file(obj,foldername)
            filename_err=strcat(foldername,'\Stats_ode15s.info');
            fid = fopen(filename_err, 'w');
            Stats=obj.sol_simulated.stats;
            nsteps=Stats.nsteps;
            nfailed=Stats.nfailed;
            nfevals=Stats.nfevals;           
            npds=Stats.npds;
            ndecomps=Stats.ndecomps;
            nsolves=Stats.nsolves;
            fprintf(fid,[num2str(nsteps),' successful steps \n']);
            fprintf(fid,[num2str(nfailed),' failed attempts \n']);
            fprintf(fid,[num2str(nfevals),' function evaluations \n']);
            fprintf(fid,[num2str(npds),' partial derivatives \n']);
            fprintf(fid,[num2str(ndecomps),' LU decompositions \n']);
            fprintf(fid,[num2str(nsolves),' solutions of linear systems \n']);
            fclose(fid);
        end
    end
end