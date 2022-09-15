classdef boussinesq_simulation_unsat
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
        function obj=boussinesq_simulation_unsat
        end
        
        function obj=simulation_parametrization(obj,discretization,source_terms,boundary_cond,percentage_loaded,t_initial,ratio_P_R,Storages_initial)
            if(nargin<7)
                obj.ratio_P_R=1;
            end
            if(nargin<8)
                Storages_initial=nan;
            end
            obj.discretization=discretization;
            obj.source_terms=source_terms;
            obj.ratio_P_R=ratio_P_R;
            obj.boundary_cond=boundary_cond;
            obj.initial_conditions=obj.set_initial_conditions(percentage_loaded,t_initial,Storages_initial);
        end
        
        function initial_conditions=set_initial_conditions(obj,percentage_loaded,t_initial,Sinitial)
            if(nargin<4)
                Sinitial=nan;
            end
            if(percentage_loaded==-1)
                 Sin_name=input('Enter the name where the initial groundwater storage and soil moisture you want stored initially: \n');
                 Sin_directory=which(Sin_name);
                 Sin=load(Sin_directory);
                 initial_conditions.Sin=Sin.Sin;
            elseif(percentage_loaded==-2 && ~isnan(sum(Sinitial)))
                % assign only the initial storage values
                if(length(Sinitial)==2*obj.discretization.Nx)
                    initial_conditions.Sin=Sinitial;
                else
                    fprintf('Error: length of the assigned initial state values of the simulation are not of the good length \n');
                end
            else
                [~,w,soil_depth,~,~,f,~,~,phi]=obj.discretization.get_resampled_variables;
                Smax=f.*w.*soil_depth;
                if(percentage_loaded==-2 && isnan(Sinitial))
                    initial_conditions.Sin=0*Smax;
                else
                    initial_conditions.Sin=[percentage_loaded*Smax;percentage_loaded*(phi-f)./f.*Smax];
                end
            end
            
            
            Edges=obj.boundary_cond.fixed_edge_matrix_values;
            Edges_bool=obj.boundary_cond.fixed_edge_matrix_boolean;
            initial_conditions.Sin(1)=Edges_bool(1)*Edges(1)+(1-Edges_bool(1))*initial_conditions.Sin(1); initial_conditions.Sin(obj.discretization.Nx)=Edges_bool(2)*Edges(2)+(1-Edges_bool(2))*initial_conditions.Sin(obj.discretization.Nx);
%             if(~isnan(sum(Sinitial)) && length(Sinitial)==3*obj.discretization.Nx+1)
%                 initial_conditions.Qin(1)=Edges_bool(3)*Edges(3)+(1-Edges_bool(3))*initial_conditions.Qin(1); initial_conditions.Qin(end)=Edges_bool(4)*Edges(4)+(1-Edges_bool(4))*initial_conditions.Qin(end);
%             else
%                 initial_conditions.Qin=-obj.compute_Q_from_S(initial_conditions.Sin)*initial_conditions.Sin;
%                 initial_conditions.Qin(1)=Edges_bool(3)*Edges(3)+(1-Edges_bool(3))*initial_conditions.Qin(1); initial_conditions.Qin(end)=Edges_bool(4)*Edges(4)+(1-Edges_bool(4))*initial_conditions.Qin(end);
%                 initial_conditions.QSin=obj.compute_QS_from_Q([initial_conditions.Sin;initial_conditions.Qin;zeros(size(initial_conditions.Qin))],t_initial)*initial_conditions.Qin;
%             end
        end

        % set the parameters for the resolution of the Boussinesq DAE
        function [t,x_center,x_edge,S,Q,QS,obj,S_u,Su_max]=implicit_scheme_solver(obj,t_properties,solver_options)
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
            solver_options.NonNegative=1:2*block_size;
            options=odeset(options,solver_options);
           
            % get consistent initial conditions thanks to decic
            [y0new,yp0new] = decic(@(t,y,yp) yp - obj.odefun(y,t), time_range(1), y0, [], yp0, [], options);
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
            S_tot=deval(obj.sol_simulated,t);
            S_u=S_tot(block_size+1:end,:);
            S=S_tot(1:block_size,:);
            size_S=size(S);
            Q=nan(size_S(1)+1,size_S(2));
            QS=nan(size_S);
            
            for i=1:length(t)
                Q(:,i)=-obj.compute_Q_from_S(S_tot(:,i))*S(:,i);
                QS(:,i)=obj.compute_QS_from_Q(S_tot(:,i),t(i))*Q(:,i)+obj.partition_source_terms_QS(S_tot(:,i),t(i)); %2*block_size+2:3*block_size+1,:);    
            end
            Su_max=obj.get_max_unsaturated_storage;
            
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
            DPSA_spat=obj.partition_source_terms_QS(S,t);
        end
        
        function [Rainfall,Infiltration,Interception,ETRs,ETRu]=partition_rainfall(obj,y,t) % unit m2/s
            [Infiltration,~,ETRs,ETRu,Interception]=obj.compute_source_term_spatialized(y,t);
            Rainfall=Infiltration+Interception;
        end
    end
    
    methods(Access=private)
        
        % Computes array at a given y
        function dy=odefun(obj,y,t)
            % C = [-alpha*A    0 ]
            %     [   0        0 ]
            %% some values necessary to compute
            block_size=obj.discretization.Nx; % size of the matrix (corresponds to the number of discretized elements)
            Edges=obj.boundary_cond.fixed_edge_matrix_boolean; % kind of boundary conditions
            [~,w,d,angle,~,f,k,f_edges,phi]=obj.discretization.get_resampled_variables; % properties of the hillslope
            A=obj.discretization.A; % derivation matrix 1
            
            %% compute recharge rate spatialized
            [Recharge_rate_spatialized,Threshold,ETR_s,ETR_u]=obj.compute_source_term_spatialized(y,t,w,d,f,phi);
            %% Compute darcy flux from one box to another with variable angle
            Q_from_S=obj.compute_Q_from_S(y,w,angle,f,k,f_edges);
            % compute Q spatial derivative
            dQ_from_S=A*Q_from_S;
            %% compute test if potential seepage
            Test_Deriv=dQ_from_S*y(1:block_size)+Recharge_rate_spatialized-ETR_s;
            Test_Deriv=Test_Deriv>=0;
            alpha=Threshold.*Test_Deriv+(1-Test_Deriv); % regularization function : drives where goes variations of mass to S or to QS
            beta=obj.beta(y,t); % regularization function : drives where goes precip to saturated or unsaturated component
            Test_Deriv2=Recharge_rate_spatialized-ETR_u>=0;
            beta=beta.*Test_Deriv2;            
            gamma=1-beta;

            %% compute C
            C1=sparse(1:block_size,1:block_size,alpha)*dQ_from_S;
            C=[C1,sparse(block_size,block_size);sparse(block_size,block_size),sparse(block_size,block_size)];
            C(1,:)=(1-Edges(1))*C(1,:);
            C(block_size,:)=(1-Edges(2))*C(block_size,:);
            C(block_size+1,:)=0;
            %% compute D
            D1=beta.*alpha.*Recharge_rate_spatialized-ETR_s;
            D=[D1;gamma.*Recharge_rate_spatialized-ETR_u];           
            D(1,:)=(1-Edges(1))*D(1,:);
            D(block_size,:)=(1-Edges(2))*D(block_size,:);
            D(block_size+1,:)=0;
            D=sparse(D);
            %% compute dy
            dy=C*y+D;
        end
        
        function D=partition_source_terms_QS(obj,y,t)
            beta=obj.beta(y,t); % regularization function : drives where goes precip to saturated or unsaturated component 
            [Recharge_rate_spatialized,Threshold,~,ETR_u]=obj.compute_source_term_spatialized(y,t);
            Test_Deriv2=Recharge_rate_spatialized-ETR_u>=0;
            beta=beta.*Test_Deriv2;
            Test_Deriv=obj.Test_Derivative(y,t);
            alpha_complementar=1-(Threshold.*Test_Deriv+(1-Test_Deriv));
            D=beta.*alpha_complementar.*Recharge_rate_spatialized;
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
            Q1_from_S=k./f_edges.*cos(angle).*(B*(y(1:block_size)./(f.*w)));
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
            block_size=obj.discretization.Nx;
            [~,w,soil_depth,~,~,f]=obj.discretization.get_resampled_variables;
            
            Thresh=threshold_function(y(1:block_size)./(f.*soil_depth.*w));
            Test_Deriv=obj.Test_Derivative(y,t);            
            alpha=Thresh.*Test_Deriv+(1-Test_Deriv);
            QS_from_Q=sparse(diag(1-alpha))*obj.discretization.A;
        end
        
% %         function dSdt_from_Q=compute_dSdt_from_Q(obj,y,t)
% %             A=obj.discretization.A;
% %             alpha=sparse(diag(obj.alpha(y,t)));
% %             dSdt_from_Q=-alpha*A;
% % %             dSdt_from_Q=sparse(dSdt_from_Q);
% %         end
% %         
% %         function dSudt_from_S=compute_dSudt_from_S(obj,y,t)
% %             [~,~,~,~,~,f,~,~,phi]=obj.discretization.get_resampled_variables;
% %             dSudt_from_S=-(phi-f)./f.*obj.gamma(y,t).*obj.compute_dSdt_from_Q(y,t)*(-obj.compute_Q_from_S(y));
% %         end
% %         
% %         function dSdt_from_leakage=compute_deep_dSdt_from_leakage(obj,y)
% %             [~,~,~,~,~,f,k]=obj.discretization.get_resampled_variables;
% %             Recharge_deep_spatialized=obj.compute_deep_recharge(y,k,f);
% %             dx=obj.discretization.compute_dx;
% %             k2=1/3600/10;
% %             dSdt_from_leakage=sum(Recharge_deep_spatialized.*dx)-k2*y(end);
% %         end
        
        function [Recharge_rate_spatialized,Threshold,ETR_s,ETR_u,Interception_rate_spatialized]=compute_source_term_spatialized(obj,y,t,w,d,f,phi)
            if(nargin<4)
                [~,w,d,~,~,f,~,~,phi]=obj.discretization.get_resampled_variables;
            end
            block_size=obj.discretization.Nx;
            relative_occupancy_rate=y(1:block_size)./(f.*w.*d);
            Threshold=threshold_function(relative_occupancy_rate);
            [Recharge_rate,ETP_rate]=obj.source_terms.compute_recharge_rate(t);
            if(Recharge_rate<0)
                Recharge_rate=Recharge_rate.*[1;(1-exp(-10*((y(2:block_size)-f(2:block_size).*w(2:block_size).*d(1))./(f(2:block_size).*w(2:block_size).*(d(2:block_size)-d(1))))))];%                Recharge_rate=Recharge_rate.*(1-exp(-1000*(1-min(1,-Recharge_rate*3600./(y./(f.*w)-d(1))))));
            end
            
            Recharge_rate_spatialized=Recharge_rate.*w;
            % if there is an ETP time series given, compute ETR from ETP, Recharge and S/Smax
            if(~isnan(ETP_rate))
%                 ETP_rate=obj.source_terms.compute_ETP_rate(t);
                ETP_rate_spatialized=ETP_rate.*w;
                [Su_max,Smax]=obj.get_max_unsaturated_storage;
                relative_occupancy_rate_unsaturated_zone=y(1+block_size:end)./Su_max;
                relative_occupancy_rate_unsaturated_zone(Su_max<=0)=0; % to avoid division by zero leading to inf or nan values for relative_occupancy_rate_unsaturated_zone quantity
                r=5;
                r_u=5;
                f_Su=1-exp(-r_u*relative_occupancy_rate_unsaturated_zone);
                f_S=1-exp(-r*relative_occupancy_rate);
                % 1st option with interception
                Interception_rate_spatialized=(Recharge_rate_spatialized-ETP_rate_spatialized>0).*(ETP_rate_spatialized)+...
                    (Recharge_rate_spatialized-ETP_rate_spatialized<=0).*Recharge_rate_spatialized;
                ETR_u=(Recharge_rate_spatialized-ETP_rate_spatialized<=0).*(ETP_rate_spatialized-Recharge_rate_spatialized).*f_Su;
                ETR_s=(Recharge_rate_spatialized-ETP_rate_spatialized<=0).*(ETP_rate_spatialized-Recharge_rate_spatialized-ETR_u).*f_S;
                Recharge_rate_spatialized=Recharge_rate_spatialized-Interception_rate_spatialized;
                % 2nd option : former computation
%                 ETR_u=0;
%                 ETR_s=(Recharge_rate_spatialized-ETP_rate_spatialized<=0).*(ETP_rate_spatialized-Recharge_rate_spatialized).*(1-exp(-5*relative_occupancy_rate));
            else
                ETR_s=0;
                ETR_u=0;
                Interception_rate_spatialized=0;
            end
        end
        
        function Recharge_deep_spatialized=compute_deep_recharge(obj,y,k,f)
            block_size=obj.discretization.Nx;
            Recharge_deep_spatialized=(k/100)*(y(1:block_size)/(100*f));
            Recharge_deep_spatialized(1)=0;
        end
        
        function OUT=Test_Derivative(obj,y,t)
            block_size=obj.discretization.Nx;
            A=obj.discretization.A;
            [Recharge_rate_spatialized,~,ETR_s]=obj.compute_source_term_spatialized(y,t);
            OUT=A*(obj.compute_Q_from_S(y)*y(1:block_size))+Recharge_rate_spatialized-ETR_s;
            OUT=OUT>=0;
        end
        
        function OUT=alpha(obj,y,t)
            block_size=obj.discretization.Nx;
            [~,w,soil_depth,~,~,f]=obj.discretization.get_resampled_variables;
            
            Thresh=threshold_function(y(1:block_size)./(f.*soil_depth.*w));
            Test_Deriv=obj.Test_Derivative(y,t);
            OUT=Thresh.*Test_Deriv+(1-Test_Deriv);
        end
        
        function OUT=beta(obj,y,t)
            block_size=obj.discretization.Nx;
            Su_max=obj.get_max_unsaturated_storage;
            
%             Thresh=threshold_function2((0.05*Smax+y(block_size+1:end))./(0.05*Smax+Su_max));
            Thresh=threshold_function2(y(block_size+1:end)./Su_max);
            Thresh(Su_max<=0)=0;
%             Test_Deriv=obj.Test_Derivative(y,t);
%             OUT=Thresh.*Test_Deriv+(1-Test_Deriv);
            OUT=1-Thresh;
        end
        
        function [Su_max,Smax]=get_max_unsaturated_storage(obj)
            [~,w,soil_depth,~,~,f,~,~,phi]=obj.discretization.get_resampled_variables;
            Su_max=(phi-f).*soil_depth.*w;
            Smax=f.*soil_depth.*w;
        end
        
        % Mass matrix: Identity on S, O for Q & QS
        function M=compute_mass_matrix(obj)
            block_size=obj.discretization.Nx;
            Nunknowns=2*block_size;
            M=zeros(Nunknowns,Nunknowns);
            M(1:2*block_size,1:2*block_size)=eye(block_size);
            M=sparse(M);
        end
        
        function [x,isterm,dir] = eventfun(obj,t,y)
            block_size=obj.discretization.Nx;
            dy=obj.odefun(y,t); 
            dy1=dy(1:block_size); 
            dy2=dy(1+block_size:end);
            [~,w,soil_depth,~,~,f,~,~,phi]=obj.discretization.get_resampled_variables;
            dy1=dy1./(f.*w.*soil_depth);
            dy2=dy2./((phi-f).*w.*soil_depth);
            x = max(abs([dy1;dy2])) - 1e-5;%5e-11; %5e-8; %#JM think to change it 5e-8
            if(x<0)
                AAA=1;
            end
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
            
%             obj=boussinesq_simulation_unsat(xcustom);
            obj=boussinesq_simulation_unsat;
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
            [~,w,d,~,~,f]=obj.discretization.get_resampled_variables;
            ETP_rate=obj.source_terms.compute_ETP_rate(t);
            if(~isnan(ETP_rate))
                ETP_rate_spatialized=bsxfun(@times,ETP_rate,w);
                relative_occupancy_rate=bsxfun(@rdivide,S,f*w.*d);
                r=0.15;
                ETR_OUT=ETP_rate_spatialized.*(1-exp(-relative_occupancy_rate*1/r));
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
                filename_err=strcat(foldername,'\Error.err');
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