classdef transport_1D_par
    properties(Access=public)
        t           % [1 x Nt] array [s]
        N           % [1 x Nt] array containing recharge time series corresponding to the time array [m/s]
        x_traj      % [(size(t_inj)*N_x) x Nt]  matrix containing the trajectories on synchronised timesteps [m]
        t_inj       % [1 x size(t_inj)] array containing the timesteps when injection is done (when N~=0)
        t_inj_pos   % [1 x size(t_inj)] array containing the position of the injection time steps compared to the obj.t array
        N_inj       % [(size(t_inj)*N_x) x 1] array containing the recharge corresponding to individual trajectory
        x           % [N_x x 1] array containing the river distance aray x [m]
        weight      % [(size(t_inj)*N_x) x 1] array containing the weight to apply to the trajectories
        ET          % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via ET
        DGW         % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via deep groundwater
        RF          % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via seepage
        DPSA        % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via direct precipitations
        NA
    end
    
    methods(Access=public)
        % precipitation time series (t is the time series of the simulation)
        function obj=transport_1D_par(t,N)        
            obj.t=t;
            obj.N=N;
        end
        
        % load the injection time series (time at which particles will be injected)
        function obj=compute_t_inj(obj,spacing_,threshold)
            if(nargin<2)
                spacing_=-1;
            end
            if(nargin<3 || spacing_~=-1)
                threshold=0;
            end

            if(spacing_==-1)
                size_N=size(obj.N);
                t_inj_pos=1:1:length(obj.t);
                if(size_N(1)>1)
                    N_max=obj.N(1,:);%max(obj.N);
                    obj.t_inj=obj.t(N_max>threshold);
                    obj.N_inj=obj.N(:,N_max>threshold);
                    obj.N_inj=reshape(obj.N_inj,size(obj.N_inj,1)*size(obj.N_inj,2),1);
                    obj.t_inj_pos=t_inj_pos(N_max>threshold);
                    if(obj.t_inj(end)==obj.t(end))
                        obj.t_inj=obj.t_inj(1:end-1);
                        obj.t_inj_pos=obj.t_inj_pos(1:end-1);
                        obj.N_inj=obj.N_inj(1:end-length(obj.x));
                    end
                else
                    obj.t_inj=obj.t(obj.N>threshold);
                    obj.N_inj=obj.N(obj.N>threshold);
                    obj.t_inj_pos=t_inj_pos(obj.N>threshold);
                    if(obj.t_inj(end)==obj.t(end))
                        obj.t_inj=obj.t_inj(1:end-1);
                        obj.t_inj_pos=obj.t_inj_pos(1:end-1);
                        obj.N_inj=obj.N_inj(1:end-1);
                    end
                end
            else
%                 spacing_=30*24*3600;
                t_min=spacing_/2; t_max=max(obj.t)-spacing_/2;
                % one per month
                t_inj_temp=t_min:spacing_:t_max;
                %             t_inj_pos=1:1:length(obj.t);
                Distance_=pdist2((obj.t)',(t_inj_temp)');
                [~,Idx_]=min(Distance_,[],1);
                obj.t_inj_pos=Idx_;
                obj.t_inj=obj.t(Idx_);
                window_width=floor(spacing_/(obj.t(2)-obj.t(1)));
                %             N_sum_month=movsum(obj.N,window_width);
                N_sum_month=conv(obj.N, ones(1, window_width),'same');
                obj.N_inj=N_sum_month(Idx_);
                obj.t_inj_pos=obj.t_inj_pos(obj.N_inj>0);
                obj.t_inj=obj.t_inj(obj.N_inj>0);
                obj.N_inj=obj.N_inj(obj.N_inj>0);
            end
        end
        
        % compute the weighting factor of the particles path
        function obj=compute_weight(obj,discretized_area) % discretized area is the area of the element dx where transport is computed (equal to dx.*w where w is the width function of the hillslope)
            size_N=size(obj.N);
            if(size_N(1)>1)
                dt=obj.t(2:end)-obj.t(1:end-1);
                t_edge1=obj.t-[dt(1),dt]/2;
                t_edge2=obj.t+[dt,dt(end)]/2;
                dt=t_edge2-t_edge1;
                inmass_weight=(discretized_area*dt)*1000; % incoming mass in m2/s
                inmass_weight=inmass_weight(:,obj.t_inj_pos);
                size_inmass_weight=size(inmass_weight);
                inmass_weight=reshape(inmass_weight,size_inmass_weight(1)*size_inmass_weight(2),1);
                weight_temp=obj.N_inj.*inmass_weight;
%                 inmass=obj.N.*(discretized_area*dt)*1000; % incoming mass in L
% %                 inmass=bsxfun(@times,discretized_area,dt)*1000;
%                 weight_temp=inmass(:,obj.t_inj_pos);
%                 size_=size(weight_temp);
%                 weight_temp=reshape(weight_temp,size_(1)*size_(2),1);
                obj.weight=weight_temp;
            else
%                 dt_inj=obj.t_inj(2:end)-obj.t_inj(1:end-1);
%                 t_inj_edge1=obj.t_inj-[dt_inj(1),dt_inj]/2;
%                 t_inj_edge2=obj.t_inj+[dt_inj,dt_inj(end)]/2;
%                 dt_inj=t_inj_edge2-t_inj_edge1;
                dt=obj.t(2:end)-obj.t(1:end-1);
                t_edge1=obj.t-[dt(1),dt]/2;
                t_edge2=obj.t+[dt,dt(end)]/2;
                dt=t_edge2-t_edge1;
                dt_inj=dt(obj.t_inj_pos);
                dflux=dt_inj.*obj.N_inj*1000;
                weight_temp=discretized_area*dflux;
                obj.weight=reshape(weight_temp,length(dflux)*length(discretized_area),1);
            end
        end
        
        % properties of transport object that informs about where particles exit and / or where a certain proportion of it exits
        function obj=instantiate_exit_tags(obj,DPSA_prop,GW_prop)
            if(nargin<2)
                obj.DGW=ones(size(obj.weight));
                obj.DPSA=zeros(size(obj.weight));
            else
                obj.DGW=GW_prop;
                obj.DPSA=DPSA_prop;
            end
            obj.RF=zeros(size(obj.weight));
            obj.ET=zeros(size(obj.weight));
        end
        
        % compute the trajectories of the flowpaths
        function obj=compute_trajectories(obj,velocity,block_size,x_S,x_Q,recharge,hydraulic_head,hydraulic_head_gradient,folder_mex_option)
            %#JM to change in a future release threslhold should be defined outside the method
            threshold=0;
            obj.x=x_S;
            if(nargin>8)
                %% C code in mex file
                load(fullfile(folder_mex_option,'velocity_field.mat'));
                load(fullfile(folder_mex_option,'t_inj.mat'));
                load(fullfile(folder_mex_option,'x_inj.mat'));
                obj.x_traj=solve_ode_trajectory(Velocity_2,t_inj,x_inj);
            else
                
                % define the velocity vector V = [(Vx)_i] with i between 1 and N_x, the number of discretized elements
                Velocity_1D= @(t,y) transport_1D_par.velocity_field_1D(obj.t,x_Q,velocity,t,y);
                
                % define a stop event for ode solving when particles cross the river boundary
                S_mat=sparse(diag(ones(block_size,1))); %#JMIPG check if S_mat is still necessary
                events=@(t,y) obj.eventfun(t,y,x_Q(2));
                options_reg = odeset('Vectorized','on','Jpattern',S_mat,'Events',events); %#JMIPG check if the vectorized option improves efficiency
                
                size_row=length(obj.t_inj)*length(obj.x);
                size_column=length(obj.t);
                % choose the format of x_traj if not problem for matlab for allocating memory
                if(size_row*size_column<15e9)
                    mat_pos_allocate=cell(length(obj.t_inj_pos),1);
                    for i=1:length(obj.t_inj_pos)
                        pos_temp=obj.t_inj_pos(i);
                        size_N=size(obj.N);
                        if(size_N(1)>1)
                            Bool_inj_spat=obj.N(:,pos_temp)>threshold;
                            % if there is a non-river cell where injection has to happen
                            if(max(x_S(Bool_inj_spat))>x_Q(2))
                                % compute trajectories
                                [~,traj_temp] = ode45(Velocity_1D,obj.t(pos_temp:end),x_S(Bool_inj_spat),options_reg);
                                % else there only injection at the river cell the trajectory is constituted of 1 point as the particle leave the river right away
                            else
                                traj_temp=[x_S(Bool_inj_spat)]';
                            end
                            x_traj_temp=zeros(length(Bool_inj_spat),size(traj_temp,2));
                            x_traj_temp(Bool_inj_spat,:)=(traj_temp(:,1:sum(Bool_inj_spat)))';
                        else
                            % compute trajectories
                            [~,traj_temp] = ode45(Velocity_1D,obj.t(pos_temp:end),x_S,options_reg);
                            x_traj_temp=(traj_temp(:,1:length(x_S)))';
                        end
                        if(length(obj.t(pos_temp:end))==2)
                            x_traj_temp=x_traj_temp(:,[1,end]);
                        end
                        matrix_positions=(combvec(block_size*(i-1)+1:block_size*i,pos_temp:size(x_traj_temp,2)+pos_temp-1))';
                        % build a boolean to delete particles once they have arrived in the river
                        bool_delete=x_traj_temp<x_Q(2);
                        bool_delete=logical(cumsum(bool_delete,2));
                        bool_delete=[false(length(obj.x),1),bool_delete(:,1:end-1)];
                        bool_delete=bool_delete(:);
                        % store x trajectories
                        x_traj_temp=x_traj_temp(:);
                        mat_pos_allocate{i}=[matrix_positions(~bool_delete,:),x_traj_temp(~bool_delete)];
                        % #JM comment to gain speed
                        % fprintf(strcat(num2str(i),'/',num2str(length(obj.t_inj)),'\n'));
                    end
                    % rebuild the (x,z) trajectories in the result matrix
                    mat_pos_allocate=vertcat(mat_pos_allocate{:});
                    obj.x_traj=sparse(mat_pos_allocate(:,1),mat_pos_allocate(:,2),mat_pos_allocate(:,3),size_row,size_column);
                else
                    fprintf('Warning due to memory allocation problems, x_traj will be in cell format \n');
                    %%%%%%% This would be to code if there is an interest in developing cell format trajectories
                end
            end
        end
        
        function obj=compute_trajectories_old(obj,velocity,block_size,x_S,x_Q,Bool_sat,folder_mex_option)
            %#JM to change in a future release threslhold should be defined outside the method
            threshold=0;
            obj.x=x_S;
            if(nargin>6)
                %% C code in mex file
                load(fullfile(folder_mex_option,'velocity_field.mat'));
                load(fullfile(folder_mex_option,'t_inj.mat'));
                load(fullfile(folder_mex_option,'x_inj.mat'));
                obj.x_traj=solve_ode_trajectory(Velocity_2,t_inj,x_inj);
            else
                
                % define the velocity vector V = [(Vx)_i , (Vz)_i] with i between 1 and N_x, the number of discretized elements
                Velocity_1D= @(t,y) transport_1D_par.velocity_field_1D(obj.t,x_Q,velocity,t,y);
                
                % define a stop event for ode solving when particles cross the river boundary
                S_mat=sparse(diag(ones(block_size,1))); %#JMIPG check if S_mat is still necessary
                events=@(t,y) obj.eventfun(t,y,x_Q(2));
                options_reg = odeset('Vectorized','on','Jpattern',S_mat,'Events',events); %#JMIPG check if the vectorized option improves efficiency
                
                size_row=length(obj.t_inj)*length(obj.x);
                size_column=length(obj.t);
                % choose the format of x_traj if not problem for matlab for allocating memory
                if(size_row*size_column<15e9)
                    mat_pos_allocate_x=cell(length(obj.t)-1,1);%zeros(size_row*size_column,3);%[];%
                    mat_pos_tags=cell(length(obj.t)-1,1);
                    for i=1:(length(obj.t)-1)
                        % %                         i=obj.t_inj_pos(i);
                        % %                         size_N=size(obj.N);
                        % %                         if(size_N(1)>1)
                        Bool_sat_temp=Bool_sat(:,i);
                        Bool_inj_spat=logical((obj.N(:,i)>threshold).*(ones(size(Bool_sat_temp))));
                        Bool_inj_nonsat=logical(Bool_inj_spat.*(~Bool_sat_temp));
                        Bool_inj_sat=logical(Bool_inj_spat.*(Bool_sat_temp));
                        % if there is a particle injection
                        if(sum(Bool_inj_spat)~=0)
                            % if there is a non-river cell where injection has to happen
                            if(sum(Bool_inj_nonsat)~=0)
                                % compute trajectories on non saturated areas
                                [~,traj_temp_nonsat] = ode45(Velocity_1D,obj.t(i:end),x_S(Bool_inj_nonsat),options_reg);
                            end
                            % compute trajectories on saturated areas
                            traj_temp_sat=x_S(Bool_inj_sat)';
                            % concatenate in single matrices for x traj and z traj
                            x_traj_temp=zeros(length(Bool_inj_spat),size(traj_temp_nonsat,1));
                            x_traj_temp(Bool_inj_nonsat,:)=traj_temp_nonsat';
                            x_traj_temp(Bool_inj_sat,1)=traj_temp_sat';
                            
                            
                            % %                         else
                            % %                             % compute trajectories
                            % %                             [~,x_traj_temp] = ode45(Velocity_2D,obj.t(pos_temp:end),[x_S;hydraulic_head2(:,pos_temp)],options_reg);
                            % %                             x_traj_temp=x_traj_temp';
                            % %                         end
                            % not to take into accounts internal timesteps if there are only 2 timesteps required
                            if(length(obj.t(i:end))==2)
                                x_traj_temp=x_traj_temp(:,[1,end]);
                            end
                            matrix_positions=(combvec(block_size*(i-1)+1:block_size*i,i:size(x_traj_temp,2)+i-1))';
                            % build a boolean to delete particles once they have arrived in the river
                            bool_delete=x_traj_temp<x_Q(2);
                            bool_delete=logical(cumsum(bool_delete,2));
                            bool_delete=[false(length(obj.x),1),bool_delete(:,1:end-1)];
                            bool_delete=bool_delete(:);
                            % store x trajectories
                            x_traj_temp=x_traj_temp(:);
                            mat_pos_allocate_x{i}=[matrix_positions(~bool_delete,:),x_traj_temp(~bool_delete)];%mat_pos_allocate(compt:compt+sum(~bool_delete)-1,:)=[matrix_positions(~bool_delete,:),x_traj_temp(~bool_delete)];%=[mat_pos_allocate_x;[matrix_positions(~bool_delete,:),x_traj_temp(~bool_delete)]];%
                            mat_pos_tags{i}=[Bool_inj_nonsat,Bool_inj_sat];
                        end
                        % #JM comment to gain speed
                         fprintf(strcat(num2str(i),'/',num2str(length(obj.t_inj)),'\n'));
                        
                    end
                    % rebuild the (x,z) trajectories in the trajectory matrix
                    mat_pos_allocate_x=vertcat(mat_pos_allocate_x{:});
                    obj.x_traj=sparse(mat_pos_allocate_x(:,1),mat_pos_allocate_x(:,2),mat_pos_allocate_x(:,3),size_row,size_column);
                    
                    mat_pos_tags=vertcat(mat_pos_tags{:});
                    obj.DGW=mat_pos_tags(:,1);
                    obj.DPSA=mat_pos_tags(:,2);
                else
                    fprintf('Warning due to memory allocation problems, x_traj will be in cell format \n');
                    %%%%%%% This would be where stands the code if there is an interest in developing cell format trajectories
                    fprintf('Not yet coded ... \n');
                end
            end
        end
        
        function [obj,mat_pos_allocate_x]=compute_trajectories3(obj,velocity,block_size,x_S,x_Q,speed_option,folder_mex_option)
            %#JM to change in a future release threslhold should be defined outside the method
            threshold=0;
            obj.x=x_S;
            if(nargin<6 | strcmp(speed_option,'fast'))
                stop_conditions=x_Q(2);
            elseif(strcmp(speed_option,'slow'))
                stop_conditions=x_S(1)/1e12;%100;%000;%1e-9;%1e12;%
            end
            if(nargin>6)
                %% C code in mex file
                load(fullfile(folder_mex_option,'velocity_field.mat'));
                load(fullfile(folder_mex_option,'t_inj.mat'));
                load(fullfile(folder_mex_option,'x_inj.mat'));
                obj.x_traj=solve_ode_trajectory(Velocity_2,t_inj,x_inj);
            else
                
                % define the velocity vector V = [(Vx)_i , (Vz)_i] with i between 1 and N_x, the number of discretized elements
                Velocity_1D= @(t,y) transport_1D_par.velocity_field_1D(obj.t,x_Q,velocity,t,y);
                
                % define a stop event for ode solving when particles cross the river boundary
                S_mat=sparse(diag(ones(block_size,1))); %#JMIPG check if S_mat is still necessary
                events=@(t,y) obj.eventfun(t,y,stop_conditions);
                options_reg = odeset('Vectorized','on','Events',events);%,'Jpattern',S_mat,'Events',events); %#JMIPG check if the vectorized option improves efficiency
                
                size_row=length(obj.t_inj)*length(obj.x);
                size_column=length(obj.t);
                % choose the format of x_traj if not problem for matlab for allocating memory
                if(size_row*size_column<1e11)%15e9)
%                     numCores = feature('numcores');
%                     p = parpool(numCores);
                    p = parpool(24);
                    mat_pos_allocate_x=cell(length(obj.t_inj)-1,1);%zeros(size_row*size_column,3);%[];%
                    N_b=obj.N;
                    t_inj_pos_b=obj.t_inj_pos;
                    t_inj_b=obj.t_inj;
                    DPSA_b=obj.DPSA;
                    t_b=obj.t;
                    x_b=obj.x;
                    parfor i=1:(length(t_inj_b)-1)
                        Bool_inj_spat=logical((N_b(:,t_inj_pos_b(i))>threshold).*(ones(size(x_S))));
                        Bool_inj_sat=logical((DPSA_b(1+(i-1)*block_size:i*block_size)==1).*(Bool_inj_spat));
                        Bool_inj_nonsat=logical(Bool_inj_spat-Bool_inj_sat);
                        % if there is a particle injection
                        if(sum(Bool_inj_spat)~=0) 
                            if(sum(Bool_inj_nonsat)~=0)
                                % compute trajectories on non saturated areas
%                                 [~,traj_temp_nonsat] = ode45(Velocity_1D,t_b(t_inj_pos_b(i):end),x_S(Bool_inj_nonsat),options_reg);
                                [~,traj_temp_nonsat] = ode45(Velocity_1D,t_b(t_inj_pos_b(i):end),x_S(Bool_inj_nonsat),options_reg);
                            else
                                traj_temp_nonsat=[];
                            end
                            % compute trajectories on saturated areas
                            traj_temp_sat=x_S(Bool_inj_sat)';
                            % concatenate in single matrices for x traj and z traj
                            x_traj_temp=zeros(length(Bool_inj_spat),size(traj_temp_nonsat,1));
                            x_traj_temp(Bool_inj_nonsat,:)=(traj_temp_nonsat(:,1:sum(Bool_inj_nonsat)))';
                            x_traj_temp(Bool_inj_sat,1)=(traj_temp_sat(:,1:sum(Bool_inj_sat)))';
                            % not to take into accounts internal timesteps if there are only 2 timesteps required
                            if(length(t_b(t_inj_pos_b(i):end))==2)
                                x_traj_temp=x_traj_temp(:,[1,end]);
                            end
                            matrix_positions=(combvec(block_size*(i-1)+1:block_size*i,t_inj_pos_b(i):size(x_traj_temp,2)+t_inj_pos_b(i)-1))';
                            % build a boolean to delete particles once they have arrived in the river
                            if(strcmp(speed_option,'fast'))
                                bool_delete=x_traj_temp<stop_conditions;
                            else
                                bool_delete=x_traj_temp<(stop_conditions/1000);
                            end
                            bool_delete=logical(cumsum(bool_delete,2));
                            bool_delete=[false(length(x_b),1),bool_delete(:,1:end-1)];
                            bool_delete=bool_delete(:);
                            % store x trajectories
                            x_traj_temp=x_traj_temp(:);
                            mat_pos_allocate_x{i}=[matrix_positions(~bool_delete,:),x_traj_temp(~bool_delete)];%mat_pos_allocate(compt:compt+sum(~bool_delete)-1,:)=[matrix_positions(~bool_delete,:),x_traj_temp(~bool_delete)];%=[mat_pos_allocate_x;[matrix_positions(~bool_delete,:),x_traj_temp(~bool_delete)]];%
                            mat_pos_allocate_x{i}=mat_pos_allocate_x{i}(mat_pos_allocate_x{i}(:,3)~=0,:);
                        end
                        % to know at what injection we are
                         fprintf(strcat(num2str(i),'/',num2str(length(t_inj_b)-1),'\n'));
                        
                    end
                    poolobj = gcp('nocreate');
                    delete(poolobj);
                    % rebuild the (x,z) trajectories in the trajectory matrix
                    mat_pos_allocate_x=vertcat(mat_pos_allocate_x{:});
%                     obj.x_traj=sparse(mat_pos_allocate_x(:,1),mat_pos_allocate_x(:,2),mat_pos_allocate_x(:,3),size_row,size_column);
                else
                    fprintf('Warning due to memory allocation problems, x_traj will be in cell format \n');
                    %%%%%%% This would be where stands the code if there is an interest in developing cell format trajectories
                    fprintf('Not yet coded ... \n');
                end
            end
        end
        
        function x_traj_soil=compute_soil_trajectories(obj,velocity_soil,block_size,x_S,x_Q)
            Velocity_reg=@(t,y) nakeinterp1(x_S,-velocity_soil,y);
            events=@(t,y)obj.eventfun(t,y,x_Q(2));
            options_reg = odeset('Vectorized','on','Events',events);
            [~,traj_temp_soil] = ode45(Velocity_reg,obj.t,x_S,options_reg);
            traj_temp_soil=traj_temp_soil';
            bool_delete=traj_temp_soil<x_Q(2);
            bool_delete=logical(cumsum(bool_delete,2));
            bool_delete=[false(length(obj.x),1),bool_delete(:,1:end-1)];
            traj_temp_soil(bool_delete)=0;
            traj_temp_soil(logical((traj_temp_soil<x_Q(2)).*(traj_temp_soil~=0)))=x_S(1);
            
            % build a trajectory matrix for the soil part
            size_trajectories=size(obj.x_traj);
            x_traj_soil=zeros(size_trajectories);
            AA=(1:1:block_size)';
            size_=size(traj_temp_soil);
            size_2=size_(2);
            for i=1:length(obj.t_inj_pos)
                if((obj.t_inj_pos(i)+size_2-1)<=size_trajectories(2))
                    x_traj_soil(AA+block_size*(i-1),obj.t_inj_pos(i):(obj.t_inj_pos(i)+size_2-1))=traj_temp_soil;
                else
                    x_traj_soil(AA+block_size*(i-1),obj.t_inj_pos(i):(size_trajectories(2)))=traj_temp_soil(:,1:(end-(obj.t_inj_pos(i)+size_2-1-size_trajectories(2))));
                end
            end
            x_traj_soil=sparse(x_traj_soil);
        end
        
        function [value,isterminal,direction] = events(obj,t,y,ylim)
            % Locate the time when height passes through zero in a decreasing
            direction=0;
            % and stop integration.
            value = max(y-ylim); % detect y-1/2 = 0
            isterminal = 1; % stop the integration
            direction = -1; % negative direction
        end
        
        function x_traj_single=compute_single_trajectory(obj,t_inj_pos,x_inj,Flux_reg,options_reg)
            [~,x_traj_single] = ode15s(Flux_reg,obj.t(t_inj_pos:end),x_inj,options_reg);
            x_traj_single=x_traj_single';
            %             [t_traj,x_traj] = ode15s(Flux,test.obj.t(1:end),x_S,options);
        end
        
        function obj=aggregate_trajectories_computation(obj,sol_simulated,block_size,x_S,x_Q)
            obj.x=x_S;
            x_t_inj=combvec((x_S)',(obj.t_inj));
            compt=1;
            obj.x_traj=nan(length(x_t_inj),length(obj.t));
            
            Flux=@(t,y) nakeinterp1([-1;x_Q],[0;deval(sol_simulated,t,block_size+1:2*block_size+1)],y);
            Storage=@(t,y) nakeinterp1(x_S,deval(sol_simulated,t,1:block_size),y);
            MassMatrix=@(t,y) Storage(t,y);
            options = odeset('Mass', MassMatrix);
            Flux_reg=@(t,y) Flux(t,y).*double(Storage(t,y)>0);
            Storage_reg=@(t,y) Storage(t,y).*double(Storage(t,y)>0)+5e-5.*double(Storage(t,y)<=0);
            MassMatrix_reg=@(t,y) Storage_reg(t,y);
            options_reg = odeset('Mass', MassMatrix_reg);
            for i=1:length(x_t_inj)
                t_inj_pos=obj.t_inj_pos(floor(i/(0.01+block_size))+1);%find(obj.t==x_t_inj(2,i));
                x_inj_pos=mod(i,block_size);
                x_inj=x_S(x_inj_pos+100*(x_inj_pos==0));
                x_traj_single=obj.compute_single_trajectory(t_inj_pos,x_inj,Flux_reg,options_reg);
                obj.x_traj(i,t_inj_pos:size(x_traj_single,2)+t_inj_pos-1)=x_traj_single;
                fprintf(strcat(num2str(compt),'/',num2str(length(x_t_inj)),'\n'));
                compt=compt+1;
            end
        end
        
        % if a particle enters the river element then it goes out of the medium the next time step
        function obj=cut_trajectory_in_river(obj,x_Q)
            x_bool=sparse((obj.x_traj<x_Q(2)) & (obj.x_traj>0));
%             x_bool2=obj.x_traj>0;
%             x_bool=x_bool.*x_bool2;
            size_=size(x_bool);
            x_bool=[sparse(size_(1),1),x_bool(:,1:end-1)];
            x_bool=logical(x_bool);
% % %#20180116                         obj.x_traj(x_bool)=nan;
            obj.x_traj(x_bool)=0;
            obj.x_traj=sparse(obj.x_traj);
% % %#20180116                         x_traj=obj.x_traj;
        end
        
        function seepage_proportion=evaluate_seepage_proportion(obj,sol_simulated,block_size,x_S,x_Q)
            Subsurface_Flux=@(t,y) nakeinterp1([-1;x_Q],[0;deval(sol_simulated,t,block_size+1:2*block_size+1)],y);
            Seepage_Flux=@(t,y) nakeinterp1(x_S,deval(sol_simulated,t,2*block_size+2:3*block_size+1),y);
            Seepage_proportion=@(t,y) Seepage_Flux(t,y)./(Seepage_Flux(t,y)+abs(Subsurface_Flux(t,y)));
            seepage_proportion=nan(size(obj.x_traj));
            for i=1:length(obj.t)
                seepage_proportion(:,i)=Seepage_proportion(obj.t(i),obj.x_traj(:,i));
            end
        end
        
        % analyze particle paths present in the same element of the hillslope at a given time. If seepage is occuring in this element
        % then a certain number of particle is exiting the medium according to the ratio of seepage vs total flow the element is experiencing
%1         function obj=cut_trajectory_seepage(obj,Discretized_Aquifer_Volume,block_size,x_S,x_Q,w,RF_spat,Subsurface_flux,x_traj_soil)
        function obj=cut_trajectory_seepage(obj,block_size,x_S,x_Q,w,RF_spat,x_traj_soil)
            Seepage=RF_spat;
            Seepage(Seepage<0)=0;
            Seepage(1,:)=0;
            NetPrecip=obj.N;
            dx_Q=x_Q(2:end)-x_Q(1:end-1);
            area_spatialized=w.*dx_Q;
            if(size(obj.N,1)==1)
                NetPrecip=area_spatialized*NetPrecip;
            else
                NetPrecip=bsxfun(@times,NetPrecip,area_spatialized);
            end
            NetPrecip(1,:)=0;

            %#JM_IPGP #JM test if the normalization should happen to the volume
            dt=obj.t(2:end)-obj.t(1:end-1);
            t_edge1=obj.t-[dt(1),dt]/2;
            t_edge2=obj.t+[dt,dt(end)]/2;
            dt=t_edge2-t_edge1;

            for i=1:length(obj.t)
                % #JM #JM_IPGP "delete" the particles that are retrieved by ET
                ET_pos=find(-NetPrecip(:,i)>0);
                Particle_to_tag=[];
                ET_m3s=0;
                for j=1:length(ET_pos)
                    if(~isempty(Particle_to_tag))
                        ET_m3s=ET_m3s-Weight_cum_ET(Weight_cum_ET<=ET_m3s);
                        ET_m3s=ET_m3s(end);
                    end
                    pos_ET=ET_pos(end-j+1); %pos_ET=ET_pos(j); % 
                    ET_m3s=ET_m3s-NetPrecip(pos_ET,i);
                    
                    particle_subject_to_ET=find((obj.x_traj(:,i)<=x_Q(pos_ET+1)) & (obj.x_traj(:,i)>x_Q(pos_ET)));
                    Weight_partial_ET=obj.weight(particle_subject_to_ET)/(1000*dt(i));
                    pos_2_ET=mod(particle_subject_to_ET,block_size); pos_2_ET(pos_2_ET==0)=block_size;
                    Initial_infiltration_point_ET=x_S(pos_2_ET);
                    [~,Index_]=sort(Initial_infiltration_point_ET);
                    Weight_cum_ET=cumsum(Weight_partial_ET(Index_));
                    Particle_to_tag=particle_subject_to_ET(Index_(Weight_cum_ET<=ET_m3s));
%                     obj.x_traj(Particle_to_tag,1:end)=0;
                     obj.ET(Particle_to_tag)=obj.ET(Particle_to_tag)+1;
                     obj.DGW(Particle_to_tag)=0;
                     obj.x_traj(Particle_to_tag,i+1:end)=0;
                end
                
                % not considering the "seepage" coming out the river as it is not real seepage
                
                Seepage_pos=find(Seepage(:,i)>0);
                Particle_Position_to_delete=[];
                Seep_m3s=0;
                for j=1:length(Seepage_pos)
                    if(~isempty(Particle_Position_to_delete))
%                         Seep_prop=Seep_prop+Seepage_proportion(pos_,i);
%                     else
                        Seep_m3s=Seep_m3s-Weight_cum(Weight_cum<=Seep_m3s);
                        Seep_m3s=Seep_m3s(end);
                    end
                    pos_=Seepage_pos(end-j+1); %pos_=Seepage_pos(j); %
                    Seep_m3s=Seep_m3s+Seepage(pos_,i);
                    
                    particle_subject_to_seep=find((obj.x_traj(:,i)<=x_Q(pos_+1)) & (obj.x_traj(:,i)>x_Q(pos_)));
                    Weight_partial=obj.weight(particle_subject_to_seep)/(1000*dt(i));
                    pos_2=mod(particle_subject_to_seep,block_size); pos_2(pos_2==0)=block_size;
                    Initial_infiltration_point=x_S(pos_2);
                    [~,Index_]=sort(Initial_infiltration_point);
                    Weight_cum=cumsum(Weight_partial(Index_));
                    Particle_Position_to_delete=particle_subject_to_seep(Index_(Weight_cum<=Seep_m3s));
                    obj.RF(Particle_Position_to_delete)=obj.RF(Particle_Position_to_delete)+1;
                    obj.DGW(Particle_Position_to_delete)=0;
                    
                    if(nargin>6)% && ~isnan(x_traj_soil))
                        obj.x_traj(Particle_Position_to_delete,i+1:end)=x_traj_soil(Particle_Position_to_delete,i+1:end);
                    else
                        obj.x_traj(Particle_Position_to_delete,i+1:end)=0;
                    end
                end
                fprintf(strcat(num2str(i),'/',num2str(length(obj.t)),'\n'));
                
%                 Distance=pdist2(obj.x_traj(:,i),x_S);
%                 [Value_,Distance_pos]=nanmin(Distance,[],2);
%                 Distance_pos(isnan(Value_))=nan;
%                 Seepage_proportion_x_traj=nan(size(Distance_pos));
%                 Seepage_proportion_x_traj(~isnan(Distance_pos))=Seepage_proportion(Distance_pos(~isnan(Distance_pos)));
%                 seepage_proportion=deval(sol_simulated,t,2*block_size+2:3*block_size+1)
            end
        end
        
        % analyze particle paths present in the same element of the hillslope at a given time. If seepage is occuring in this element
        % then a certain number of particle is exiting the medium according to the ratio of seepage vs total flow the element is experiencing
        function [obj,Error_RF_DGW,Error_ET,mat_pos_allocate_x_sorted]=cut_trajectory_ET_RF(obj,block_size,x_S,x_Q,Flux_in_spat,RF_spat,speed_option,mat_pos_allocate_x_sorted)
            if(nargin<7 | strcmp(speed_option,'fast'))
                % not considering the "seepage" coming out the river and then deep groundwater flow is computed as soon as the particle cross the aquifer river border
                speed_option='fast';
                RF_spat(1,:)=0;
            end
                
            % not considering ET in the river ? To check #JM_IPGP
%             Flux_in_spat(1,:)=0;

            %#JM_IPGP #JM test if the normalization should happen to the volume
            dt=obj.t(2:end)-obj.t(1:end-1);
            t_edge1=obj.t-[dt(1),dt]/2;
            t_edge2=obj.t+[dt,dt(end)]/2;
            dt=t_edge2-t_edge1;
            
            % Tags position in the matrix trajectories to delete
            Position_to_tag_before_deletion=[0,0,0]; % last column enables to track for which flux the particle exits: 0=RF exit, 1=DGW exit, 2=ET exit
            
            % Array registering error on fluxes made
            Error_ET=zeros(1,length(obj.t)-1);
            Error_RF_DGW=zeros(1,length(obj.t)-1);
            [size_row,~]=size(obj.x_traj);
%             Position_Number=(1:size_row)';
%             Position_Number2=ones(size_row,1);
            
            for i=1:(length(obj.t)-1)
                if(i==1215)
                    AA=1;
                end
                % "delete" the particles that are retrieved by ET
                ET_pos=find(-Flux_in_spat(:,i)>0);
                ET_m3s=0;
                for j=1:length(ET_pos)
                    pos_ET=ET_pos(end-j+1); %pos_ET=ET_pos(j); % 
                    ET_m3s=ET_m3s-Flux_in_spat(pos_ET,i);
                    
                    particle_subject_to_ET=find((obj.x_traj(:,i)<=x_Q(pos_ET+1)) & (obj.x_traj(:,i)>x_Q(pos_ET)) & (obj.DPSA~=1));
                    particle_subject_to_ET = setdiff(particle_subject_to_ET,Position_to_tag_before_deletion(:,1));
                                        
                    Weight_partial_ET=obj.weight(particle_subject_to_ET)/(1000*dt(i)).*obj.DGW(particle_subject_to_ET);
                    pos_2_ET=mod(particle_subject_to_ET,block_size); pos_2_ET(pos_2_ET==0)=block_size;
                    Initial_infiltration_point_ET=x_S(pos_2_ET);
                    [~,Index_]=sort(Initial_infiltration_point_ET);
                    Weight_cum_ET=cumsum(Weight_partial_ET(Index_));
                    Particle_to_tag=particle_subject_to_ET(Index_(Weight_cum_ET<=ET_m3s));
                    Position_to_tag_before_deletion=[Position_to_tag_before_deletion;[Particle_to_tag,ones(size(Particle_to_tag))*i,ones(size(Particle_to_tag))*2]];
                    if(~isempty(Particle_to_tag))
                        ET_m3s=ET_m3s-Weight_cum_ET(Weight_cum_ET<=ET_m3s);
                        ET_m3s=ET_m3s(end);
                    end
                end
                Error_ET(i)=ET_m3s;
                
                
                Seepage_pos=find(RF_spat(:,i)>0);
                Seep_m3s=0;
                %%
                Dist_Trajectories_MeshCentres = pdist2(mat_pos_allocate_x_sorted{i}(:,3),obj.x(Seepage_pos));
                [DistMin,Pos_Seepage] = min(Dist_Trajectories_MeshCentres,[],2);
                Delta_X = (x_Q(2:end)-x_Q(1:end-1))/2;
                Delta_X = Delta_X(Seepage_pos); 
                % assume X_Q are uniformly spaced #JM to change
                Delta_X = Delta_X(1);
                Pos_Seepage=Pos_Seepage(DistMin<=Delta_X);
                particle_subject_to_seep=mat_pos_allocate_x_sorted{i}(DistMin<=Delta_X,1);
                final_location=mat_pos_allocate_x_sorted{i}(DistMin<=Delta_X,3);
                Bool_1=obj.DPSA(particle_subject_to_seep);
                particle_subject_to_seep=particle_subject_to_seep(Bool_1<1);
                Pos_Seepage=Pos_Seepage(Bool_1<1);
                final_location=final_location(Bool_1<1);

                [particle_subject_to_seep,Indexes_delete] = setdiff(particle_subject_to_seep,Position_to_tag_before_deletion(:,1));
                Pos_Seepage=Pos_Seepage(Indexes_delete);
                final_location=final_location(Indexes_delete);
                Weight_partial=obj.weight(particle_subject_to_seep)/(1000*dt(i)).*obj.DGW(particle_subject_to_seep);
                
                pos_2=mod(particle_subject_to_seep,block_size); pos_2(pos_2==0)=block_size;
                Initial_infiltration_point=x_S(pos_2);
                
                [Initial_infiltration_point,Index_]=sort(Initial_infiltration_point);
                Weight_partial=Weight_partial(Index_);
                Pos_Seepage=Pos_Seepage(Index_);
                final_location=final_location(Index_);
                particle_subject_to_seep=particle_subject_to_seep(Index_);
                
                Seep_m3s=RF_spat(Seepage_pos,i);
                
%                 Weight_cum=accumarray(Pos_Seepage,Weight_partial,@cumsum);
                [Pos_Seepage_sorted, I_sorted] = sort(Pos_Seepage);
                Weight_cum = accumarray(Pos_Seepage_sorted, Weight_partial(I_sorted), [], @(r){cumsum(r)});
                particle_subject_to_seep=accumarray(Pos_Seepage_sorted,particle_subject_to_seep(I_sorted),[],@(r){r});
                final_location=accumarray(Pos_Seepage_sorted,final_location(I_sorted),[],@(r){r});
                Initial_infiltration_point=accumarray(Pos_Seepage_sorted,Initial_infiltration_point(I_sorted),[],@(r){r});
% %                 Weight_cum=accumarray(Pos_Seepage,Weight_partial,[],@(r){cumsum(r)});
% %                 particle_subject_to_seep=accumarray(Pos_Seepage,particle_subject_to_seep,[],@(r){r});
                if(~isempty(Weight_cum) && ~isempty(Weight_cum{1}))
%                     [final_location{1},Index_location_1]=sort(final_location{1});
%                     Weight_cum{1}=[Weight_cum{1}(2:end)-Weight_cum{1}(1:end-1);Weight_cum{1}(1)];
%                     Weight_cum{1}=Weight_cum{1}(Index_location_1);
%                     Weight_cum{1}=cumsum(Weight_cum{1});
%                     particle_subject_to_seep{1}=particle_subject_to_seep{1}(Index_location_1);
                    
                    Weight_cum{1}=[Weight_cum{1}(1);Weight_cum{1}(2:end)-Weight_cum{1}(1:end-1)];
                    Weight_cum{1}=flip(Weight_cum{1});
% %                     Weight_cum{1}=[Weight_cum{1}(1:end-1)-Weight_cum{1}(2:end);Weight_cum{1}(1)];
% % %                     Weight_cum{1}=cumsum(Weight_cum{1});
                    particle_subject_to_seep{1}=flip(particle_subject_to_seep{1});
                    
%                     Weight_cum{1}=[Weight_cum{1}(1);Weight_cum{1}(2:end)-Weight_cum{1}(1:end-1)];
                    [particle_certain_to_seep,Indexes_particle_certain]=setdiff(particle_subject_to_seep{1},mat_pos_allocate_x_sorted{i+1}(mat_pos_allocate_x_sorted{i+1}(:,3)<x_Q(2),1));
                    [particle_potential_to_seep,Indexes_potential]=setdiff(particle_subject_to_seep{1},particle_certain_to_seep);
                    particle_subject_to_seep{1}=[particle_certain_to_seep;particle_potential_to_seep];
                    Weight_cum{1}=[Weight_cum{1}(Indexes_particle_certain);Weight_cum{1}(Indexes_potential)];
                    Weight_cum{1}=cumsum(Weight_cum{1});
                end
                Particle_Position_to_delete=cell(length(Weight_cum),1);
                Seepage_value=0;
                for k=length(Weight_cum):-1:1
                    Seepage_value=Seepage_value+Seep_m3s(k);
                    pos_=Seepage_pos(k);
                    [ ~, ix ] = min( abs( [flip(Weight_cum{k});0]-Seepage_value ) );
                    ix=length(Weight_cum{k})-ix+2;
                    if(ix>1)
                        Particle_Position_to_delete{k}=[particle_subject_to_seep{k}(1:1:(ix-1)),i*ones(ix-1,1),double(pos_==1)*ones(ix-1,1)];%ones(size(Particle_Position_to_delete))*]
                        Seepage_value=Seepage_value-Weight_cum{k}(ix-1);
                    end
                end
                Error_RF_DGW(i)=Seepage_value;
                if(i>2000 & Error_RF_DGW(i)/(sum(RF_spat(:,i)))>0.05)
                    AA=1;
                end
                Particle_Position_to_delete=vertcat(Particle_Position_to_delete{:});
                Position_to_tag_before_deletion=[Position_to_tag_before_deletion;Particle_Position_to_delete];

                fprintf(strcat(num2str(i),'/',num2str(length(obj.t)),'\n'));

                
                %%
% % % %                 for j=1:length(Seepage_pos)
% % % %                     pos_=Seepage_pos(end-j+1); %pos_=Seepage_pos(j); %
% % % %                     Seep_m3s=Seep_m3s+RF_spat(pos_,i);
% % % %                     
% % % % % %                     particle_subject_to_seep=find((obj.x_traj(:,i)<=x_Q(pos_+1)) & (obj.x_traj(:,i)>x_Q(pos_)) & (obj.DPSA<1));%find((obj.x_traj(:,i)<=x_Q(pos_+1)) & (obj.x_traj(:,i)>x_Q(pos_)));%
% % % %                     Bool_1=(mat_pos_allocate_x_sorted{i}(:,3)<=x_Q(pos_+1)) & (mat_pos_allocate_x_sorted{i}(:,3)>x_Q(pos_));% & (obj.DPSA<1);
% % % %                     particle_subject_to_seep = mat_pos_allocate_x_sorted{i}(Bool_1,1);
% % % %                     Bool_1=obj.DPSA(particle_subject_to_seep);
% % % %                     particle_subject_to_seep=particle_subject_to_seep(Bool_1<1);
% % % % 
% % % %                     particle_subject_to_seep = setdiff(particle_subject_to_seep,Position_to_tag_before_deletion(:,1));
% % % %                     Weight_partial=obj.weight(particle_subject_to_seep)/(1000*dt(i)).*obj.DGW(particle_subject_to_seep);
% % % % 
% % % %                     pos_2=mod(particle_subject_to_seep,block_size); pos_2(pos_2==0)=block_size;
% % % %                     Initial_infiltration_point=x_S(pos_2);
% % % %                     [~,Index_]=sort(Initial_infiltration_point);
% % % %                     if(pos_==1)
% % % %                         Index_=flip(Index_);
% % % %                     end
% % % %                     Weight_cum=cumsum(Weight_partial(Index_));
% % % %                     % option 1
% % % % %                     Particle_Position_to_delete=particle_subject_to_seep(Index_(Weight_cum<=Seep_m3s));
% % % %                     % option 2
% % % %                     [ ~, ix ] = min( abs( [0;Weight_cum]-Seep_m3s ) );
% % % %                     if(ix>1)
% % % %                         Particle_Position_to_delete=particle_subject_to_seep(Index_(1:1:(ix-1)));
% % % %                     else
% % % %                         Particle_Position_to_delete=[];
% % % %                     end
% % % %                    
% % % % %                     if(pos_==1)
% % % % %                         Index_surplus=find(Weight_cum>Seep_m3s,1);
% % % % %                         Particle_Position_to_delete=[Particle_Position_to_delete;particle_subject_to_seep(Index_(Index_surplus))];
% % % % %                         Weight_division_particle=(Seep_m3s-Weight_cum(Index_surplus-1))/(Weight_cum(Index_surplus)-Weight_cum(Index_surplus-1));
% % % % %                         Position_to_tag_before_deletion=[Position_to_tag_before_deletion;[Particle_Position_to_delete,ones(size(Particle_Position_to_delete))*i,ones(size(Particle_Position_to_delete))*double(pos_==1),[ones(length(Particle_Position_to_delete)-1,1);Weight_division_particle]]];
% % % % %                     else
% % % %                         Position_to_tag_before_deletion=[Position_to_tag_before_deletion;[Particle_Position_to_delete,ones(size(Particle_Position_to_delete))*i,ones(size(Particle_Position_to_delete))*double(pos_==1)]];
% % % % %                     end
% % % %                     if(~isempty(Particle_Position_to_delete))
% % % %                         % option 1
% % % % %                         Seep_m3s=Seep_m3s-Weight_cum(Weight_cum<=Seep_m3s);
% % % % %                         Seep_m3s=Seep_m3s(end);
% % % %                         % option 2
% % % %                         Seep_m3s=Seep_m3s-Weight_cum(ix-1);
% % % %                     end
% % % %                 end
% % % %                 
% % % % %                 last_particle_subject_to_seep = setdiff(particle_subject_to_seep,Particle_Position_to_delete);
% % % % %                 if(~isempty(last_particle_subject_to_seep))
% % % % %                     last_weight_partial=obj.weight(last_particle_subject_to_seep)/(1000*dt(i)).*obj.DGW(last_particle_subject_to_seep);
% % % % %                     [LastValue,Index_last]=min(abs(last_weight_partial-Seep_m3s));
% % % % %                     if(LastValue<Seep_m3s)
% % % % %                         Position_to_tag_before_deletion=[Position_to_tag_before_deletion;[last_particle_subject_to_seep(Index_last),i,double(pos_==1)]];
% % % % %                         Seep_m3s=Seep_m3s-last_weight_partial(Index_last);
% % % % %                     end
% % % % %                 end
% % % %                 Error_RF_DGW(i)=Seep_m3s;
% % % %                 fprintf(strcat(num2str(i),'/',num2str(length(obj.t)),'\n'));
            end
            
            % effectively delete particles in ET, RF (and DGW if slow option)
            Position_to_tag_before_deletion=Position_to_tag_before_deletion(2:end,:);
%             [size_row,size_column]=size(obj.x_traj);
%             Position_to_tag_before_deletion2=[(1:1:size_row)',ones(size_row,1)*size_column];
%             Position_to_tag_before_deletion2(Position_to_tag_before_deletion(:,1),2)=Position_to_tag_before_deletion(:,2);
%             [I,J,S]=find(obj.x_traj);
%             keep=J<=Position_to_tag_before_deletion2(I,2);
%             I=I(keep); J=J(keep); S=S(keep);

            size_row=length(obj.x)*length(obj.t_inj);
            size_column=length(obj.t);
            Position_to_tag_before_deletion2=[(1:1:size_row)',ones(size_row,1)*size_column];
            Position_to_tag_before_deletion2(Position_to_tag_before_deletion(:,1),2)=Position_to_tag_before_deletion(:,2);
            mat_pos_allocate_x_sorted=vertcat(mat_pos_allocate_x_sorted{:});
            keep=mat_pos_allocate_x_sorted(:,2)<=Position_to_tag_before_deletion2(mat_pos_allocate_x_sorted(:,1),2);
            mat_pos_allocate_x_sorted=mat_pos_allocate_x_sorted(keep,:);
            

            clear obj.x_traj Position_to_tag_before_deletion2 keep
%             obj.x_traj=sparse(I, J, S ,size_row,size_column);
            obj.RF(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==0,1))=obj.RF(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==0,1))...
                                                                                                    +obj.DGW(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==0,1));
            obj.ET(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==2,1))=obj.ET(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==2,1))...
                                                                                                    +obj.DGW(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==2,1));
            if(strcmp(speed_option,'fast'))
                obj.DGW(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==0,1))=0;
                obj.DGW(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==2,1))=0;
            else
                obj.NA=obj.DGW;
                obj.DGW=zeros(size(obj.DGW));
                obj.DGW(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==1,1))=obj.NA(Position_to_tag_before_deletion(Position_to_tag_before_deletion(:,3)==1,1));
                obj.NA(Position_to_tag_before_deletion(:,1))=0;
            end
        end
        
        function obj=cut_trajectory_groundwater(obj,sol_simulated,x_Q)
            dt=obj.t(2:end)-obj.t(1:end-1);
            t_edge1=obj.t-[dt(1),dt]/2;
            t_edge2=obj.t+[dt,dt(end)]/2;
            dt=t_edge2-t_edge1;
            
            Q_GW=-deval(sol_simulated,obj.t,block_size+1);
            
            for i=1:length(obj.t)
                particle_subject_to_gw=find((obj.x_traj(:,i)<=x_Q(3)));
                particle_position=obj.x_traj(particle_subject_to_gw,i);
                [~,Index_]=sort(particle_position);
                Weight_partial=obj.weight(particle_subject_to_gw)/(1000*dt(i));
                Weight_cum=cumsum(Weight_partial(Index_));
                Particle_Position_to_gw=particle_subject_to_gw(Index_(Weight_cum<=Q_GW(i)));
                obj.x_traj(Particle_Position_to_gw,i+1:end)=0;
            end
        end
        
        function obj=update_DGW(obj,x_Q)
            Bool_fin=obj.x_traj(:,end)>x_Q(2);
            obj.DGW(Bool_fin)=0;
        end
        
        % if the particle injected on a saturated area -> goes out of the medium right away
        % if travel_time_soil is given 1D rootage depending on soil permeability
        function obj=cut_trajectory_saturated_areas(obj,Bool_sat,x_traj_soil)
            %% add to take into account gravity driven flow in soil/saprolite profile
            Bool_sat=Bool_sat(:,obj.t_inj_pos);
            size_=size(Bool_sat);
            Bool_sat=reshape(Bool_sat,size_(1)*size_(2),1);
            if(nargin>2)% && ~isnan(x_traj_soil))
                obj.x_traj(Bool_sat,:)=x_traj_soil(Bool_sat,:);
            else
                %%
                obj.x_traj=obj.x_traj(Bool_sat,:);
                % % %#20180116             detect_no_nan_values=~isnan(x_traj_temp);
                detect_no_nan_values=(obj.x_traj~=0);
                size_2=size(detect_no_nan_values);
                %#20180124 possibility to comment next line not to take into account direct precipitations onto saturated areas
                detect_no_nan_values=[false(size_2(1),1),detect_no_nan_values(:,1:end-1)];
                % % %#20180116                         x_traj_temp(detect_no_nan_values)=nan;
                obj.x_traj(detect_no_nan_values)=0;
% % %                 x_traj_temp=obj.x_traj(Bool_sat,:);
% % %                 % % %#20180116             detect_no_nan_values=~isnan(x_traj_temp);
% % %                 detect_no_nan_values=(x_traj_temp~=0);
% % %                 size_2=size(detect_no_nan_values);
% % %                 %#20180124 possibility to comment next line not to take into account direct precipitations onto saturated areas
% % %                 detect_no_nan_values=[false(size_2(1),1),detect_no_nan_values(:,1:end-1)];
% % %                 % % %#20180116                         x_traj_temp(detect_no_nan_values)=nan;
% % %                 x_traj_temp(detect_no_nan_values)=0;
% % %                 obj.x_traj(Bool_sat,:)=x_traj_temp;
                obj.DPSA(Bool_sat)=obj.DPSA(Bool_sat)+1;
                obj.DGW(Bool_sat)=0;
            end
        end
        
        % retrieve all the particles that are exiting just after sample_time 
        function [ttd_array,transit_times,weights_,transit_distance,mean_velocity,number_of_particle]=get_ttd(obj,sample_time)
            [~,Index_time]=min(abs(obj.t-sample_time));
% % %#20180116                       Particle_constituting_ttd=find(isnan(obj.x_traj(:,Index_time+1)) & ~isnan(obj.x_traj(:,Index_time)));
            Particle_constituting_ttd=find((obj.x_traj(:,Index_time+1)==0) & (obj.x_traj(:,Index_time)~=0));
            if(~isempty(Particle_constituting_ttd))
                t_inj=obj.t_inj(floor((Particle_constituting_ttd-0.1)/length(obj.x))+1);
                transit_times=(sample_time-t_inj)';
                weights_=obj.weight(Particle_constituting_ttd);
                % #JM comment to gain speed
% %                 if(length(transit_times)==1)
% %                     ttd_array=(repelem(transit_times,floor(weights_)))';
% %                 else
% %                     ttd_array=(repelem(transit_times,floor(weights_)));
% %                 end
                ttd_array=nan;
                x_traj_temp=obj.x_traj(Particle_constituting_ttd,:);
% % %#20180116                         [r,c] = find(~isnan(x_traj_temp));
                [r,c] = find((x_traj_temp)~=0);
                if(length(unique(r))==1)
                    x_init=x_traj_temp(r(1),c(1));
                    x_fin=obj.x(1);%x_fin=x_traj_temp(r(1),c(end));
                else
                    firstIndex = accumarray(r,c,[size(x_traj_temp,1),1],@min,NaN);
                    lastIndex = accumarray(r,c,[size(x_traj_temp,1),1],@max,NaN);
                    LinearIndices_init=sub2ind(size(x_traj_temp),(1:1:size(x_traj_temp,1))',firstIndex);
                    LinearIndices_fin=sub2ind(size(x_traj_temp),(1:1:size(x_traj_temp,1))',lastIndex);
                    x_init=x_traj_temp(LinearIndices_init);
                    x_fin=obj.x(1);%x_fin=x_traj_temp(LinearIndices_fin);
                end
                transit_distance=x_init-x_fin;
                mean_velocity=transit_distance./transit_times;
                number_of_particle=length(Particle_constituting_ttd);
            else
                ttd_array=nan;
                transit_times=nan;
                weights_=nan;
                transit_distance=nan;
                mean_velocity=nan;
                number_of_particle=0;
            end
        end
        
        function [t_out,delta_t,x_fin,delta_x_groundwater,t_in,weights]=get_trajectory_properties(obj,mat_pos_allocate_x)
            if(nargin<2)
                [row,col] = find(obj.x_traj);
                size_trajectories=size(obj.x_traj);
            else
                row = mat_pos_allocate_x(:,1);
                col = mat_pos_allocate_x(:,2);
                size_trajectories=[length(obj.t_inj)*length(obj.x),length(obj.t)];
            end
            min_=accumarray(row,col,[],@min);
            max_=accumarray(row,col,[],@max);
            max_=[max_;zeros(length(obj.weight)-length(max_),1)];
            min_=[min_;zeros(length(obj.weight)-length(min_),1)];
            
            t_out=nan(size(max_));
            t_in=nan(size(min_));
            
            t_out(max_>0)=obj.t(max_(max_>0));
            t_out(max_==length(obj.t))=nan;
            t_in(min_>0)=obj.t(min_(min_>0));
            delta_t=t_out-t_in;
            
            
            AA=(1:1:size_trajectories(1))';
            linearInd_fin = sub2ind(size_trajectories,AA(max_>0),max_(max_>0));
            linearInd_init = sub2ind(size_trajectories,AA(min_>0),min_(min_>0));
            
            if(nargin<2)
                x_fin=nan(size(max_));
                x_fin(max_>0)=obj.x_traj(linearInd_fin);
                x_fin(max_(max_>0)==length(obj.t))=nan;
                x_init=nan(length(min_),1);
                x_init(min_>0)=obj.x_traj(linearInd_init);
                delta_x_groundwater=x_init-x_fin;
                if(length(unique(obj.DPSA))>2)
                    t_out=[t_out;t_in(obj.DPSA>0 & obj.DPSA<1)];
                    delta_t=[delta_t;zeros(sum(obj.DPSA>0 & obj.DPSA<1),1)];
                    x_fin=[x_fin;x_init(obj.DPSA>0 & obj.DPSA<1)];
                    delta_x_groundwater=[delta_x_groundwater;zeros(sum(obj.DPSA>0 & obj.DPSA<1),1)];
                    t_in=[t_in;t_in(obj.DPSA>0 & obj.DPSA<1)];
                    % append the weights depending on the partitioning between DPSA and Infiltration
                    weights=obj.weight;
                    weights(obj.DGW>0)=weights(obj.DGW>0).*obj.DGW(obj.DGW>0);
                    weights(obj.NA>0)=weights(obj.NA>0).*obj.NA(obj.NA>0);
                    weights(obj.RF>0)=weights(obj.RF>0).*obj.RF(obj.RF>0);
                    weights(obj.DPSA==1)=weights(obj.DPSA==1);
                    weights=[weights;obj.weight(obj.DPSA>0 & obj.DPSA<1).*obj.DPSA(obj.DPSA>0 & obj.DPSA<1)];
                end                
            else
                weights=obj.weight;
                x_init=[];
                x_fin=[];
                delta_x_groundwater=[];
                if(length(unique(obj.DPSA))>2)
                    t_out=[t_out;t_in(obj.DPSA>0 & obj.DPSA<1)];
                    delta_t=[delta_t;zeros(sum(obj.DPSA>0 & obj.DPSA<1),1)];
                    t_in=[t_in;t_in(obj.DPSA>0 & obj.DPSA<1)];
                    % append the weights depending on the partitioning between DPSA and Infiltration
                    weights=obj.weight;
                    weights(obj.DGW>0)=weights(obj.DGW>0).*obj.DGW(obj.DGW>0);
                    weights(obj.NA>0)=weights(obj.NA>0).*obj.NA(obj.NA>0);
                    weights(obj.RF>0)=weights(obj.RF>0).*obj.RF(obj.RF>0);
                    weights(obj.DPSA==1)=weights(obj.DPSA==1);
                    weights=[weights;obj.weight(obj.DPSA>0 & obj.DPSA<1).*obj.DPSA(obj.DPSA>0 & obj.DPSA<1)];
                else
                    weights=obj.weight;
                end
            end
        end
        
        function [transit_times,weights_,number_of_particle,delta_t_unsat]=get_ttds_bis(obj,sample_time,t_out,delta_t,weights,delta_t_unsat)
            %#JM without weights provided seems to be an outdated methods chose the second instead
            if(nargin<5)
                if(nargin<3)
                    [t_out,delta_t]=get_trajectory_properties(obj);
                end
                [~,Index_time]=min(abs(obj.t-sample_time));
                transit_times=delta_t(t_out==obj.t(Index_time));
                weights_=obj.weight(t_out==obj.t(Index_time));
                number_of_particle=sum(t_out==obj.t(Index_time));
            else
                [~,Index_time]=min(abs(obj.t-sample_time));
                transit_times=delta_t(t_out==obj.t(Index_time));
                weights_=weights(t_out==obj.t(Index_time));
                number_of_particle=sum(t_out==obj.t(Index_time));
                if(nargin==6)
                    delta_t_unsat=delta_t_unsat(t_out==obj.t(Index_time));
                end
            end
        end
        
        function [transit_times,weights_]=get_rtds_bis(obj,sample_time,t_out,delta_t)
            if(nargin<3)
                [t_out,delta_t]=get_trajectory_properties(obj);
            end
            
        end
        
        function  [transit_times,weights_,number_of_particle,x_fin]=get_rtd(obj,sample_time)
            [~,Index_time]=min(abs(obj.t-sample_time));
            Particle_constituting_ttd=find((obj.x_traj(:,Index_time)~=0)); %Particle_constituting_ttd=find((obj.x_traj(:,Index_time)~=0) & (obj.x_traj(:,Index_time+1)~=0));%
            
            x_traj_temp=obj.x_traj(Particle_constituting_ttd,:);
            test=sum(x_traj_temp~=0,2);
            Particle_constituting_ttd=Particle_constituting_ttd(test>=3);
            
            if(~isempty(Particle_constituting_ttd))
                t_inj=obj.t_inj(floor((Particle_constituting_ttd-0.1)/length(obj.x))+1);
                transit_times=(sample_time-t_inj)';
                weights_=obj.weight(Particle_constituting_ttd);
                number_of_particle=length(Particle_constituting_ttd);
                x_fin=obj.x_traj(Particle_constituting_ttd,Index_time);
%                 idx = sub2ind(size(obj.x_traj), Particle_constituting_ttd, floor((Particle_constituting_ttd-0.1)/length(obj.x))+1);
%                 x_init=obj.x_traj(idx);
%                 transit_distance=x_init-x_fin;
            else
                transit_times=nan;
                weights_=nan;
                number_of_particle=0;
            end
        end
        
        function [ttd_array,transit_times,weights_,transit_distance,mean_velocity,number_of_particle]=retrieve_particle_characteristics(obj,Particle_constituting_ttd,sample_time)
            if(~isempty(Particle_constituting_ttd))
                t_inj=obj.t_inj(floor(Particle_constituting_ttd/length(obj.x))+1);
                transit_times=(sample_time-t_inj)';
                weights_=obj.weight(Particle_constituting_ttd);
                if(length(transit_times)==1)
                    ttd_array=(repelem(transit_times,floor(weights_)))';
                else
                    ttd_array=(repelem(transit_times,floor(weights_)));
                end
                x_traj_temp=obj.x_traj(Particle_constituting_ttd,:);
                [r,c] = find(~isnan(x_traj_temp));
                if(length(unique(r))==1)
                    x_init=x_traj_temp(r(1),c(1));
                    x_fin=x_traj_temp(r(1),c(end));
                else
                    firstIndex = accumarray(r,c,[size(x_traj_temp,1),1],@min,NaN);
                    lastIndex = accumarray(r,c,[size(x_traj_temp,1),1],@max,NaN);
                    LinearIndices_init=sub2ind(size(x_traj_temp),(1:1:size(x_traj_temp,1))',firstIndex);
                    LinearIndices_fin=sub2ind(size(x_traj_temp),(1:1:size(x_traj_temp,1))',lastIndex);
                    x_init=x_traj_temp(LinearIndices_init);
                    x_fin=x_traj_temp(LinearIndices_fin);
                end
                transit_distance=x_init-x_fin;
                mean_velocity=transit_distance./transit_times;
                number_of_particle=length(Particle_constituting_ttd);
            else
                ttd_array=nan;
                transit_times=nan;
                weights_=nan;
                transit_distance=nan;
                mean_velocity=nan;
                number_of_particle=0;
            end
        end
        
        function [time_support,weighted_pdf]=compute_full_distribution_from_sample_time(obj,sample_time)
            if(nargin<2)
                [~,transit_times,weights_]=obj.get_ttd(obj.t(end-2));
            else
                [~,transit_times,weights_]=obj.get_ttd(sample_time);
            end
            [time_support,weighted_pdf]=compute_full_distribution(obj,transit_times,weights_);
        end
        
        function [time_support,weighted_pdf,mean_,std_]=compute_full_distribution(obj,transit_times,weights_,time_support,transit_distance,travel_time_soil)
            if(nargin<4)
                time_support=(0:0.1:1000)*24*365*3600;
            end
            if(nargin>4)
                travel_time_particle_soil=interpn(obj.x-obj.x(1),travel_time_soil,transit_distance);
                transit_times(transit_times==0)=travel_time_particle_soil(transit_times==0);
                time_support=(0:0.01:1000)*24*365*3600;
            end

                % #JM20180617 Test add a one day delay for particle injected on the river (which have a zero transit time)
%                 transit_times(transit_times==0)=0.1*24*3600; % maybe we need to refine time_support discretization to see the difference
                %             pdf=@(t)(0.25.*(transit_distance/alpha_disp).*(transit_time./t).^3).^0.5.*exp(-0.25*transit_distance/alpha_disp.*(t-transit_time).^2./(t.*transit_time));
                % Gelhar et al. 1992 WRR : transit_distance/alpha_disp = 10
%#JM20180617                 pdf=@(t)(10/(2*pi).*bsxfun(@rdivide,transit_times,t.^3)).^0.5.*exp(-10/2.*(bsxfun(@minus,transit_times,t)).^2./((bsxfun(@times,transit_times,t))));
                pdf=@(t)(25/(4*pi).*bsxfun(@rdivide,transit_times,t.^3)).^0.5.*exp(-25/4.*(bsxfun(@minus,transit_times,t)).^2./((bsxfun(@times,transit_times,t))));
                %             k=0.5*ones(size(transit_times)); theta=2*transit_times;
                %             fun1=@(k,t) t.^(k-1)./gamma(k); fun2=@(theta,t) exp(-t./theta); fun3=@(a,b) a./b;
                %             pdf=@(t) bsxfun(fun3,bsxfun(fun1,k,t).*bsxfun(fun2,theta,t),theta.^k);
                
                %             time_support=(0:1:365)*24*3600;
                if(time_support(1)~=0)
                    transit_times(transit_times==0)=time_support(1);
                end
                compute_pdf=pdf(time_support);
                if(time_support(1)==0)
                    compute_pdf(:,1)=0;
                    compute_pdf(transit_times==0,1)=2/(time_support(2)-time_support(1));
                    compute_pdf(transit_times==0,2:end)=0;
                end
                if(length(weights_)~=1)
                    weight_total=sum(weights_);
                    weighted_pdf=1/weight_total*nansum(bsxfun(@times,compute_pdf,weights_));
                else
                    weighted_pdf=compute_pdf;
                end
                %             figure;
                weighted_pdf_area=trapz(time_support,weighted_pdf);
                weighted_pdf=weighted_pdf/weighted_pdf_area;
                %             plot(time_support/(24*365*3600),weighted_pdf*24*365*3600/weighted_pdf_area);
                % #JM comment to gain speed
                mean_=trapz(time_support/(24*365*3600),time_support/(24*365*3600).*weighted_pdf*24*365*3600/weighted_pdf_area);
                std_=sqrt(trapz(time_support/(24*365*3600),(time_support/(24*365*3600)).^2.*weighted_pdf*24*365*3600/weighted_pdf_area)-mean_^2);
                % % %             fprintf(strcat('mean: ',num2str(mean_),'years \n','std: ',num2str(std_),'years \n'));
                %             plot(time_support/(24*3600),weighted_pdf*24*3600/weighted_pdf_area);
                %             pdf=@(t)(0.25.*10.*(transit_times./t).^3).^0.5.*exp(-0.25*10.*(t-transit_times).^2./(t.*transit_times));
        end
        
        function [time_support,weighted_pdf,mean_,std_]=compute_full_distribution_unsat(obj,transit_times,weights_,time_support,transit_times_unsat)
%             pdf=@(t)(25/(4*pi).*bsxfun(@rdivide,transit_times,t.^3)).^0.5.*exp(-25/4.*(bsxfun(@minus,transit_times,t)).^2./((bsxfun(@times,transit_times,t))));
            transit_times_unsat(isinf(transit_times_unsat))=2171700;
            
            pdf=@(t) (bsxfun(@times,exp(-bsxfun(@rdivide,bsxfun(@minus,t,transit_times_unsat),transit_times)),1./transit_times)).*(t>=transit_times_unsat);
            if(time_support(1)~=0)
                transit_times(transit_times==0)=time_support(1);
            end
            compute_pdf=pdf(time_support);
            if(time_support(1)==0)
                compute_pdf(:,1)=0;
                compute_pdf(transit_times==0,1)=2/(time_support(2)-time_support(1));
                compute_pdf(transit_times==0,2:end)=0;
            end
            weighted_pdf_area=trapz(time_support,compute_pdf,2);
            weighted_pdf=compute_pdf./weighted_pdf_area;
            if(length(weights_)~=1)
                weight_total=sum(weights_);
                weighted_pdf=1/weight_total*nansum(bsxfun(@times,weighted_pdf,weights_));
            end
            %             figure;
           
            %             plot(time_support/(24*365*3600),weighted_pdf*24*365*3600/weighted_pdf_area);
        end
        
        
        function [ttd_arrays,mean_,transit_times,weights_,number_of_particle,time_of_sampling]=get_ttds(obj,spacing_)
            if(nargin<2)
                spacing_=-1;
            end
            
            mean_=nan(1,length(obj.t)-1);
            number_of_particle=nan(1,length(obj.t)-1);
            time_of_sampling=nan(1,length(obj.t)-1);
            for i=1:(length(obj.t)-1)
                [ttd_arrays{i},transit_times{i},weights_{i},~,~,number_of_particle(i)]=obj.get_ttd(obj.t(i));
%                 [ttd_arrays{i},transit_times{i},weights_{i},~,~,number_of_particle(i)]=obj.get_rtd(obj.t(i));
                mean_(i)=nanmean(ttd_arrays{i});
                time_of_sampling(i)=obj.t(i);
                fprintf(strcat(num2str(i),'\n'));
            end
            
            if(spacing_~=-1)
                window_width=floor(spacing_/(obj.t(2)-obj.t(1)));
                J=floor(length(obj.t)/window_width); %length(obj.t_inj); 
                ttd_arrays2=cell(1,J); transit_times2=cell(1,J); weights_2=cell(1,J); number_of_particle2=nan(1,J); time_of_sampling2=nan(1,J);
                for j=1:J
%                     k=floor(obj.t_inj(j)/window_width);
                    time_of_sampling2(j)=(window_width*(j-1)+1+window_width*(j))/2;
                    ttd_arrays2{j}=vertcat(ttd_arrays{window_width*(j-1)+1:window_width*(j)});
                    ttd_arrays2{j}=ttd_arrays2{j}(~isnan(ttd_arrays2{j}));
                    transit_times2{j}=vertcat(transit_times{window_width*(j-1)+1:window_width*(j)});%,transit_times{i},weights_{i},~,~,number_of_particle(i)]=
                    transit_times2{j}=transit_times2{j}(~isnan(transit_times2{j}));
                    weights_2{j}=vertcat(weights_{window_width*(j-1)+1:window_width*(j)});
                    weights_2{j}=weights_2{j}(~isnan(weights_2{j}));
                    number_of_particle2(j)=nansum(number_of_particle(window_width*(j-1)+1:window_width*(j)));
                end
                ttd_arrays=ttd_arrays2;
                transit_times=transit_times2;
                weights_=weights_2;
                number_of_particle=number_of_particle2;
                time_of_sampling=time_of_sampling2;
            end

%             else
%                 mean_=nan(1,length(obj.t_inj));
%                 number_of_particle=nan(1,length(obj.t_inj));
%                 for i=1:(length(obj.t)-1)
%                     [ttd_arrays{i},transit_times{i},weights_{i},~,~,number_of_particle(i)]=obj.get_ttd(obj.t(i));
%                     mean_(i)=nanmean(ttd_arrays{i});
%                     fprintf(strcat(num2str(i),'\n'));
%                 end
%             end
        end
        
        function [time_support_cell,weighted_pdf_cell,time_of_sampling,mean_,std_,number_of_particle]=compute_evoluting_ttds(obj,spacing_)
            [~,~,transit_times,weights_,number_of_particle,time_of_sampling]=get_ttds(obj,spacing_);
            time_support_cell=cell(size(transit_times));
            weighted_pdf_cell=cell(size(transit_times));
            mean_=nan(size(transit_times));
            std_=nan(size(transit_times));
            for i=1:length(transit_times)
                [time_support_cell{i},weighted_pdf_cell{i},mean_(i),std_(i)]=obj.compute_full_distribution(transit_times{i},weights_{i});
            end
        end
        
        function DSi_time_series=compute_DSi_conc(obj)
            DSi_time_series=nan(1,length(1704-1531+1));
            for i=1531:1704
                [~,transit_times{i},weights_{i},~,~,~]=obj.get_ttd(obj.t(i));
%                 [ttd_arrays{i},transit_times{i},weights_{i},~,~,number_of_particle(i)]=obj.get_rtd(obj.t(i));
                [time_support_cell,weighted_pdf_cell{i},~,~]=obj.compute_full_distribution(transit_times{i},weights_{i});
                DSi_support=0.2*time_support_cell/(365*2600*24)+4;
                DSi_support(time_support_cell==0)=0;
                DSi_time_series(i-1531+1)=trapz(time_support_cell,DSi_support.*weighted_pdf_cell{i});
            end
            
        end
        
        function save_injection_var(obj,filename)
            t_inj=obj.t_inj;
            filename_t_inj=fullfile(filename,'t_inj.txt');
            dlmwrite(filename_t_inj,t_inj, '-append', 'precision', '%E','delimiter','\t');
            save(fullfile(filename,'t_inj.mat'),'t_inj');
            x_inj=obj.x;
            filename_x_inj=fullfile(filename,'x_inj.txt');
            dlmwrite(filename_x_inj,x_inj, '-append', 'precision', '%E','delimiter','\t');
            save(fullfile(filename,'x_inj.mat'),'x_inj');
        end
        
        function [x,isterm,dir] = eventfun(obj,t,y,xlim)
            x = y(end)-xlim;
            isterm = 1;
            dir = 0;  %or -1, doesn't matter
        end
        
       function [DPSA_prop,DGW_prop,RF_spat,Flux_in_spat]=compute_DPSA_GW_proportion(obj,runs)
            [DPSA_,~,DPSA_spat,RF_spat]=compute_DPSA_RF(runs.simulation_results,runs.boussinesq_simulation);
            [~,Flux_in_spat]=compute_infiltration(runs.simulation_results,runs.boussinesq_simulation);
            size_mat=size(DPSA_spat);
            DPSA_spat=reshape(DPSA_spat(:,obj.t_inj_pos),size_mat(1)*(length(obj.t_inj)),1); %DPSA_spat=reshape(DPSA_spat(:,1:end-1),size_mat(1)*(size_mat(2)-1),1);
            Flux_in_spat_reshaped=reshape(Flux_in_spat(:,obj.t_inj_pos),size_mat(1)*(length(obj.t_inj)),1); %Flux_in_spat_reshaped=reshape(Flux_in_spat(:,1:end-1),size_mat(1)*(size_mat(2)-1),1);
            DPSA_prop=DPSA_spat./Flux_in_spat_reshaped;
%             DPSA_prop(obj.N_inj==0)=0;
            DPSA_prop(1,:)=1;
            if(sum(DPSA_prop<0)>0)
                fprintf(strcat('WARNING: \n',num2str(sum(DPSA_prop<0)),' values in the spatialized matrix representing Direct Precipitations on Saturated Areas \n','were negative (minimal values found:',num2str(DPSA_spat(DPSA_prop==min(DPSA_prop))),'). They have been replaced by 0.\n'));
                DPSA_prop(DPSA_prop<0)=0;
            end
            if(sum(DPSA_prop>1)>0)
                percent_above=(max(DPSA_prop)-1)*100;
                fprintf(strcat('WARNING: \n',num2str(sum(DPSA_prop>1)),' values in the spatialized matrix representing Direct Precipitations on Saturated Areas \n','were above the Flux_in infiltration matrix values \n','(maximal value was ',num2str(percent_above),'%% above the corresponding infiltration value). \n','They have been replaced by the corresponding infiltration value.\n'));
                DPSA_prop(DPSA_prop>1)=1;
            end
            DGW_prop=1-DPSA_prop;
%             DGW_prop(obj.N_inj==0)=0;
       end
    end
    methods(Static)
        function [obj,DSi]=test(runs,folder_root,velocity_soil)
            sim_res=runs.simulation_results;
            bouss_sim=runs.boussinesq_simulation;
%             f=bouss_sim.hydraulic_properties.f;
            t=sim_res.t;
            N_in=bouss_sim.compute_recharge(t);
            [~,width,soil_depth,~,~,f]=bouss_sim.discretization.get_resampled_variables;
            Smax=f.*width.*soil_depth;
            relstor=bsxfun(@rdivide,sim_res.S,Smax);
            % first mesh is always considered as saturated (it is the river)
            relstor(1,:)=1;
            Bool_sat=relstor>=1;
            x_S=sim_res.x_S;
            x_Q=sim_res.x_Q;
            dx=x_Q(2:end)-x_Q(1:end-1);
            block_size=length(x_S);
            obj=transport_1D_par(t,N_in);
            obj.x=x_S;
%             spacing_=1*24*3600;
            spacing_=-1;
            obj=obj.compute_t_inj(spacing_);
            obj=obj.compute_weight(dx.*width);
%             if(nargin>1)
%                 obj.save_injection_var(folder_root);
%                 obj=obj.compute_trajectories(bouss_sim.sol_simulated,block_size,x_S,x_Q,folder_root);
%             else
% % % % %                 obj=obj.compute_trajectories(bouss_sim.sol_simulated,block_size,x_S,x_Q);
                S_edges=bouss_sim.discretization.Omega*sim_res.S;
                S_edges(1,:)=1;
                S_edges(end,:)=1;
                velocity=sim_res.Q./S_edges;
                
                [~,c]=find(isnan(velocity));
                if(max(c)>2)
                    fprintf('Warning: Nan values in velocity vector \n');
                end
                velocity(isnan(velocity))=0;
                obj=obj.compute_trajectories(velocity,block_size,x_S,x_Q);
%             end
            if(nargin>2)
                % compute the trajectory for typical soil particles
                x_traj_soil=obj.compute_soil_trajectories(velocity_soil,block_size,x_S,x_Q);
                
                obj=obj.cut_trajectory_saturated_areas(Bool_sat,x_traj_soil);
                obj=obj.cut_trajectory_seepage(bouss_sim.sol_simulated,block_size,x_S,x_Q,width,x_traj_soil);
            else
                obj=obj.cut_trajectory_saturated_areas(Bool_sat);%,travel_time_soil);
                %             obj=obj.cut_trajectory_in_river(sim_res.x_Q);
                obj=obj.cut_trajectory_seepage(bouss_sim.sol_simulated,block_size,x_S,x_Q,width);
            end
            
            
            time_support=(0:0.1:1000)*24*365*3600;
            DSi_conc_dis=0.2*time_support/(24*365*3600)+4;
%             DSi_conc_dis(1)=0;
            for i=1531:1704
                [~,transit_times,weights_,~,~,number_of_particle(i-1530)]=obj.get_ttd(obj.t(i));
                [~,weighted_pdf,mean_(i),std_(i)]=obj.compute_full_distribution(transit_times,weights_);
                DSi(i-1530)=trapz(time_support,DSi_conc_dis.*weighted_pdf);
            end
%             obj=obj.aggregate_trajectories_computation(bouss_sim.sol_simulated,block_size,x_S,x_Q);
        end
        
%         function [obj,DSi]=compute_DSi(run_obj,folder_root)
%             
%         end

        function [obj,t_out_groundwater,transit_times_groundwater,x_fin_groundwater,weights_groundwater]=transport_with_rooting(runs,x,slope_angle,k_soil)
            % instantiate transport_1D_par object
            [obj,x_S,x_Q,width,velocity,RF_spat,Flux_in_spat]=transport_1D_par.instantiate_transport_and_velocity_field(runs);
            block_size=length(x_S);
            %             dx=x_Q(2:end)-x_Q(1:end-1);
            
            % compute the trajectories inside the aquifer
            % option 1 all integrated
            tic
            speed_option='slow';
            [obj,mat_pos_allocate_x]=obj.compute_trajectories3(velocity,block_size,x_S,x_Q,speed_option);%obj=obj.compute_trajectories3(velocity,block_size,x_S,x_Q,speed_option);%compute_trajectories2(velocity,block_size,x_S,x_Q,Bool_sat);
            if(sum(RF_spat(:)<0)>0)
                min_index=find(RF_spat(:)==min(RF_spat(:)));
                fprintf(strcat('WARNING: \n',num2str(sum(RF_spat(:)<0)),' values in the spatialized matrix representing Return Flow \n','were negative (minimal values found:',num2str(RF_spat(min_index)),'). They have been replaced by 0.\n'));
                RF_spat(RF_spat<0)=0;
            end
            toc
%             RF_spat(1,:)=0;
            [~,~,X] = unique(mat_pos_allocate_x(:,2));
            mat_pos_allocate_x = accumarray(X,1:size(mat_pos_allocate_x,1),[],@(r){mat_pos_allocate_x(r,:)});
            [obj,Error_RF_DGW,Error_ET,mat_pos_allocate_x]=obj.cut_trajectory_ET_RF(block_size,x_S,x_Q,Flux_in_spat,RF_spat,speed_option,mat_pos_allocate_x);%cut_trajectory_seepage(block_size,x_S,x_Q,width,RF_spat);%ob
            
            % option 2 all separated
%             obj=obj.compute_trajectories(velocity,block_size,x_S,x_Q,recharge,hydraulic_head,hydraulic_head_gradient);
%             obj=obj.cut_trajectory_saturated_areas(Bool_sat);
%             
%             obj=obj.cut_trajectory_seepage(block_size,x_S,x_Q,width,RF_spat);%obj=obj.cut_trajectory_seepage(Discretized_Aquifer_Volume,block_size,x_S,x_Q,width,RF_spat,Subsurface_flux);

%             obj=obj.cut_trajectory_groundwater(sol_simulated,x_Q);
% % %             obj=obj.update_DGW(x_Q);
            
            % retrieve the particles essential properties for their travel inside the aquifer
%             [t_out_groundwater,transit_times_groundwater,x_fin_groundwater]=obj.get_trajectory_properties;
            [t_out_groundwater,transit_times_groundwater,x_fin_groundwater,delta_x_groundwater,t_in,weights_groundwater]=obj.get_trajectory_properties(mat_pos_allocate_x);
            toc
            dt=obj.t(2:end)-obj.t(1:end-1);
            t_edge1=obj.t-[dt(1),dt]/2;
            t_edge2=obj.t+[dt,dt(end)]/2;
            dt=t_edge2-t_edge1;
            
            DPSA_select=[obj.DPSA==1;obj.DPSA(obj.DPSA>0 & obj.DPSA<1)];
            DGW_select=[obj.DGW;zeros(sum(obj.DPSA>0 & obj.DPSA<1),1)];
            RF_select=[obj.RF;zeros(sum(obj.DPSA>0 & obj.DPSA<1),1)];
            DGW_select2=[obj.NA;zeros(sum(obj.DPSA>0  & obj.DPSA<1),1)];
            Q_Seep=nan(4505,1);
            Q_DP_out=nan(4505,1);
            Q_DGW_out=nan(4505,1);
            Q_NA_out=nan(4505,1);
            for i=1:4505
                Q_Seep(i)=(sum(weights(obj.t(i)==t_out_groundwater & RF_select>0)/1000)/(dt(i)))';
                Q_DP_out(i)=(sum(weights(obj.t(i)==t_out_groundwater & DPSA_select>0)/1000)/(dt(i)))';
                Q_DGW_out(i)=(sum(weights(obj.t(i)==t_out_groundwater & DGW_select>0)/1000)/(dt(i)))';
                Q_NA_out(i)=(sum(weights(obj.t(i)==t_out_groundwater & DGW_select2>0)/1000)/(dt(i)))';
            end
            
            RF_=sum(RF_spat(2:end,:));
            DGW_=RF_spat(1,:);
            DPSA_=compute_DPSA_RF(runs.simulation_results,runs.boussinesq_simulation);
            figure; hold on
            plot(temp.t_real,DPSA_(4332:4505))
            plot(temp.t_real,Q_DP_out)
            figure; hold on
            
            % compute the typical travel time inside the subsurface (soil layer)
            
            %% check DPSA quanti
%             DPSA_m3s_part=obj.DPSA.*obj.weight/(1000);
%             t_out_pos=(obj.t_inj_pos(ceil((1:1:349600)'/80)))';
%             DPSA_time=accumarray(t_out_pos,DPSA_m3s_part,[],@nansum);
%             dt=obj.t(2:end)-obj.t(1:end-1);
%             dt=obj.t(2:end)-obj.t(1:end-1);
%             t_edge1=obj.t-[dt(1),dt]/2;
%             t_edge2=obj.t+[dt,dt(end)]/2;
%             dt=t_edge2-t_edge1;
%             DPSA_time=[DPSA_time',0]./dt;
            
% % %             slope_angle1_top=interpn(x,slope_angle,x_S);
% % %             k_soil2=ones(size(x_S))*k_soil;
% % %             velocity_soil=k_soil2.*sin(slope_angle1_top);
% % %             travel_time_soil=cumtrapz(dx./velocity_soil);
% % %             travel_time_soil2=[0;travel_time_soil];
% % %             
% % %             % travel time of the particles inside the soil layer
% % %             transit_times_soil=interpn(x_Q,travel_time_soil2,full(x_fin_groundwater),'linear',0);
% % %             
% % %             % t_out_river : time at which the particles arrives finally inside the river 
% % %             t_out_river=cat(2,t_out_groundwater,transit_times_soil); 
% % %             t_out_river= nansum(t_out_river,2); 
% % %             t_out_river=interp1([obj.t,2*obj.t(end)-obj.t(end-1)],[obj.t,2*obj.t(end)-obj.t(end-1)],t_out_river,'nearest',nan);
% % %             transit_times_total=cat(2,transit_times_groundwater,transit_times_soil); transit_times_total= nansum(transit_times_total,2); 
% % %             
% % %             transit_times_total(obj.ET>0)=nan;
% % %             t_out_river(obj.ET>0)=nan;
% % %             t_out_groundwater(obj.ET>0)=nan;
% % %             transit_times_groundwater(obj.ET>0)=nan;
% % %             x_fin_groundwater(obj.ET>0)=nan;
% % %             transit_times_soil(obj.ET>0)=nan;
        end
        
        function transport_2compartments(run_shallow,run_deep)
            obj=transport_1D_par.run_transport_simu(run_shallow);
            % retrieve the particles essential properties for their travel inside the aquifer
            [t_out_groundwater,transit_times_groundwater,x_fin_groundwater]=obj.get_trajectory_properties;
            
            % exponential approximation for run_deep
            [~,~,~,RF_spat_deep]=compute_DPSA_RF(run_deep.simulation_results,run_deep.boussinesq_simulation);
            Q_out_deep=sum(RF_spat_deep);
            Volume_deep=trapz(x_S,run_deep.simulation_results.S);
            tau_deep=Volume_deep./Q_out_deep;
            
            time_support=logspace(-4,3,10001)*24*365*3600;%(0:0.001:10)*24*365*3600;
            weighted_pdfs_shallow=nan(length(obj.t)-1,length(time_support));
            weighted_pdfs_deep=nan(length(obj.t)-1,length(time_support));
            weighted_pdfs=nan(length(obj.t)-1,length(time_support));
            mean_shallow=nan(1,length(obj.t)-1);
            std_shallow=nan(1,length(obj.t)-1);
            mean_tot=nan(1,length(obj.t)-1);
            std_tot=nan(1,length(obj.t)-1);
            median_tot=nan(1,length(obj.t)-1);
            median_shallow=nan(1,length(obj.t)-1);

            number_of_particle_shallow=nan(1,length(obj.t)-1);
            Q_shallow=nan(1,length(obj.t)-1);
            dt=obj.t(2:end)-obj.t(1:end-1);
            t_edge1=obj.t-[dt(1),dt]/2;
            t_edge2=obj.t+[dt,dt(end)]/2;
            dt=t_edge2-t_edge1;
            for i=1:(length(obj.t)-1)
                [transit_times,weights_,number_of_particle_shallow(i)]=obj.get_ttds_bis(obj.t(i),t_out_groundwater,transit_times_groundwater);
                [~,weighted_pdfs_shallow(i,:),mean_shallow(i),std_shallow(i)]=obj.compute_full_distribution(transit_times,weights_,time_support);
                weighted_pdfs_deep(i,:)=1/tau_deep(i)*exp(-time_support/tau_deep(i));
                weighted_pdfs(i,:)=(Q_out_deep(i)*weighted_pdfs_deep(i,:)+(DPSA_shallow(i)+RF_shallow(i))*weighted_pdfs_shallow(i,:))/(Q_out_deep(i)+DPSA_shallow(i)+RF_shallow(i));
                Q_shallow(i)=sum(weights_/1000)/(dt(i));
                mean_tot(i)=trapz(time_support/(24*365*3600),weighted_pdfs(i,:)*(24*365*3600).*time_support/(24*365*3600));
                std_tot(i)=sqrt(trapz(time_support/(24*365*3600),(time_support/(24*365*3600)).^2.*weighted_pdfs(i,:)*24*365*3600)-mean_tot(i)^2);
                [~,Index_]=min(abs(0.5-cumtrapz(time_support/(24*365*3600),weighted_pdfs_shallow(i,:)*(24*365*3600))));
                median_shallow(i)=time_support(Index_)/(3600*24);
                [~,Index_]=min(abs(0.5-cumtrapz(time_support/(24*365*3600),weighted_pdfs(i,:)*(24*365*3600))));
                median_tot(i)=time_support(Index_)/(3600*24);
            end
        end
        
        function [obj_shallow,obj_reg,obj_soil,time_support,weighted_pdfs,weighted_pdfs_shallow,weighted_pdfs_deep_approx,...
                weighted_pdfs_reg,weighted_pdfs_soil]=transport_3compartments(run_soil,run_reg,run_bed)
            tic
            obj_soil=transport_1D_par.run_transport_simu(run_soil);
            toc
            [obj_reg,DPSA_shallow,RF_shallow,x_S]=transport_1D_par.run_transport_simu(run_reg);
            toc
            
            
            obj_shallow=obj_reg;
            obj_shallow.weight(obj_reg.DPSA==1)=obj_soil.weight(obj_reg.DPSA==1);
            obj_shallow.x_traj(obj_reg.DPSA==1,:)=obj_soil.x_traj(obj_reg.DPSA==1,:);
            

            % transport for run_deep
%             [obj_deep,DPSA_deep,RF_deep]=transport_1D_par.run_transport_simu(run_bed);
%             toc
%             % transport total
%             obj_tot=obj_deep;
%             obj_tot.weight(obj_deep.DPSA==1)=obj_shallow.weight(obj_deep.DPSA==1);
%             obj_tot.x_traj(obj_deep.DPSA==1,:)=obj_shallow.x_traj(obj_deep.DPSA==1,:);
%             [t_out_gw,transit_times_gw]=obj_deep.get_trajectory_properties;
%             [t_out_tot,transit_times_tot]=obj_tot.get_trajectory_properties;
            
            % retrieve the particles essential properties for their travel inside the aquifer
            [t_out_shallow,transit_times_shallow,x_fin_shallow]=obj_shallow.get_trajectory_properties;
            [t_out_reg,transit_times_reg,delta_x_reg,t_in_reg,weights_reg]=obj_reg.get_trajectory_properties;
            [t_out_soil,transit_times_soil,delta_x_soil,t_in_soil,weights_soil]=obj_soil.get_trajectory_properties;
            
            % exponential approximation for run_deep
            [~,~,~,RF_spat_deep]=compute_DPSA_RF(run_bed.simulation_results,run_bed.boussinesq_simulation);
            Q_out_deep=sum(RF_spat_deep);
            Q_out_deep(Q_out_deep<0)=0;
            Volume_deep=trapz(x_S,run_bed.simulation_results.S);
            tau_deep=Volume_deep./Q_out_deep;
            
            time_support=logspace(-4,3,10001)*24*365*3600;%(0:0.001:10)*24*365*3600;
            weighted_pdfs_shallow=nan(length(obj_shallow.t)-1,length(time_support));
            weighted_pdfs_deep=nan(length(obj_shallow.t)-1,length(time_support));
            weighted_pdfs_deep_approx=nan(length(obj_shallow.t)-1,length(time_support));
            weighted_pdfs=nan(length(obj_shallow.t)-1,length(time_support));
            weighted_pdfs2=nan(length(obj_shallow.t)-1,length(time_support));
            
            
            weighted_pdfs_soil=nan(length(obj_shallow.t)-1,length(time_support));
            weighted_pdfs_reg=nan(length(obj_shallow.t)-1,length(time_support));
% %             mean_tot=nan(1,length(obj_shallow.t)-1);
% %             std_tot=nan(1,length(obj_shallow.t)-1);
% %             median_tot=nan(1,length(obj_shallow.t)-1);
% %             median_shallow=nan(1,length(obj_shallow.t)-1);
% %             delta_18O_synth=nan(1,length(obj_shallow.t)-1);

            number_of_particle_shallow=nan(1,length(obj_shallow.t)-1);
            Q_shallow=nan(1,length(obj_shallow.t)-1);
            dt=obj_shallow.t(2:end)-obj_shallow.t(1:end-1);
            t_edge1=obj_shallow.t-[dt(1),dt]/2;
            t_edge2=obj_shallow.t+[dt,dt(end)]/2;
            dt=t_edge2-t_edge1;
            for i=1:(length(obj_shallow.t)-1)
                [transit_times,weights_,number_of_particle_shallow(i)]=obj_shallow.get_ttds_bis(obj_shallow.t(i),t_out_shallow,transit_times_shallow);
                [~,weighted_pdfs_shallow(i,:)]=obj_shallow.compute_full_distribution(transit_times,weights_,time_support);
                
                [transit_times,weights_]=obj_reg.get_ttds_bis(obj_reg.t(i),t_out_reg,transit_times_reg);
                [~,weighted_pdfs_reg(i,:)]=obj_reg.compute_full_distribution(transit_times,weights_,time_support);
                
                 [transit_times,weights_]=obj_soil.get_ttds_bis(obj_soil.t(i),t_out_soil,transit_times_soil);
                [~,weighted_pdfs_soil(i,:)]=obj_soil.compute_full_distribution(transit_times,weights_,time_support);
                
%                 [transit_times,weights_]=obj_deep.get_ttds_bis(obj_deep.t(i),t_out_gw,transit_times_gw);
%                 [~,weighted_pdfs_deep(i,:)]=obj_deep.compute_full_distribution(transit_times,weights_,time_support);
                
                weighted_pdfs_deep_approx(i,:)=1/tau_deep(i)*exp(-time_support/tau_deep(i));
                
                
%                 [transit_times,weights_]=obj_tot.get_ttds_bis(obj_tot.t(i),t_out_tot,transit_times_tot);
%                 [~,weighted_pdfs(i,:)]=obj_tot.compute_full_distribution(transit_times,weights_,time_support);
                weighted_pdfs2(i,:)=(Q_out_deep(i)*weighted_pdfs_deep_approx(i,:)+(DPSA_shallow(i)+RF_shallow(i))*weighted_pdfs_shallow(i,:))/(Q_out_deep(i)+DPSA_shallow(i)+RF_shallow(i));
                Q_shallow(i)=sum(weights_/1000)/(dt(i));
                
% %                 mean_tot(i)=trapz(time_support/(24*365*3600),weighted_pdfs(i,:)*(24*365*3600).*time_support/(24*365*3600));
% %                 std_tot(i)=sqrt(trapz(time_support/(24*365*3600),(time_support/(24*365*3600)).^2.*weighted_pdfs(i,:)*24*365*3600)-mean_tot(i)^2);
% %                 [~,Index_]=min(abs(0.5-cumtrapz(time_support/(24*365*3600),weighted_pdfs_shallow(i,:)*(24*365*3600))));
% %                 median_shallow(i)=time_support(Index_)/(3600*24);
% %                 [~,Index_]=min(abs(0.5-cumtrapz(time_support/(24*365*3600),weighted_pdfs(i,:)*(24*365*3600))));
% %                 median_tot(i)=time_support(Index_)/(3600*24);
            end
        end
        
        function [obj,x_S,x_Q,width,velocity,RF_spat,Flux_in_spat]=instantiate_transport_and_velocity_field(runs,threshold)
            if(nargin<2)
                threshold=0;
            end
            sim_res=runs.simulation_results;
            bouss_sim=runs.boussinesq_simulation;
            %             f=bouss_sim.hydraulic_properties.f;
            t=sim_res.t;
            N_in=bouss_sim.compute_recharge(t);
            [~,width]=bouss_sim.discretization.get_resampled_variables;
            x_S=sim_res.x_S;
            x_Q=sim_res.x_Q;
            dx=x_Q(2:end)-x_Q(1:end-1);

            
            obj=transport_1D_par(t,N_in);
            obj.x=x_S;
            
            %             spacing_=1*24*3600;
            spacing_=-1;
            obj=obj.compute_t_inj(spacing_,threshold);
            [DPSA_prop,GW_prop,RF_spat,Flux_in_spat]=obj.compute_DPSA_GW_proportion(runs);
            obj=obj.compute_weight(dx.*width);%(Flux_in_spat);%
%             if(nargin<2)
%                 obj=obj.instantiate_exit_tags;
%             else
                obj=obj.instantiate_exit_tags(DPSA_prop,GW_prop);
%             end
            
            S_edges=bouss_sim.discretization.Omega*sim_res.S;
            S_edges(1,:)=1;
            S_edges(end,:)=1;
            velocity=sim_res.Q./S_edges;
            
            [~,c]=find(isnan(velocity));
            if(max(c)>2)
                fprintf('Warning: Nan values in velocity vector \n');
            end
            velocity(isnan(velocity))=0;
        end
        
        function [obj,DPSA_,RF_,x_S]=run_transport_simu(run_boussinesq)
            threshold=0;
            
            [DPSA_,RF_,DPSA_spat,RF_spat]=compute_DPSA_RF(run_boussinesq.simulation_results,run_boussinesq.boussinesq_simulation);
            [~,Flux_in_spat]=compute_infiltration(run_boussinesq.simulation_results,run_boussinesq.boussinesq_simulation);
            Bool_spat_rest=max(Flux_in_spat)>threshold;
            Flux_in_spat=Flux_in_spat(:,Bool_spat_rest);
            DPSA_spat=DPSA_spat(:,Bool_spat_rest);
            size_mat=size(DPSA_spat);
            DPSA_spat=reshape(DPSA_spat(:,1:end-1),size_mat(1)*(size_mat(2)-1),1);
            Flux_in_spat_reshaped=reshape(Flux_in_spat(:,1:end-1),size_mat(1)*(size_mat(2)-1),1);
            DPSA_prop=DPSA_spat./Flux_in_spat_reshaped;
            if(sum(DPSA_prop<0)>0)
                fprintf(strcat('WARNING: \n',num2str(sum(DPSA_prop<0)),' values in the spatialized matrix representing Direct Precipitations on Saturated Areas \n','were negative (minimal values found:',num2str(DPSA_spat(DPSA_prop==min(DPSA_prop))),'). They have been replaced by 0.\n'));
                DPSA_prop(DPSA_prop<0)=0;
            end
            if(sum(DPSA_prop>1)>0)
                percent_above=(max(DPSA_prop)-1)*100;
                fprintf(strcat('WARNING: \n',num2str(sum(DPSA_prop>1)),' values in the spatialized matrix representing Direct Precipitations on Saturated Areas \n','were above the Flux_in infiltration matrix values \n','(maximal value was ',num2str(percent_above),'%% above the corresponding infiltration value). \n','They have been replaced by the corresponding infiltration value.\n'));
                DPSA_prop(DPSA_prop>1)=1;
            end
            DGW_prop=1-DPSA_prop;
            
            % instantiate transport_1D_par object
            [obj,x_S,x_Q,width,velocity]=transport_1D_par.instantiate_transport_and_velocity_field(run_boussinesq,DPSA_prop,DGW_prop,threshold);
            block_size=length(x_S);
            
                        %             dx=x_Q(2:end)-x_Q(1:end-1);
            
            % compute the trajectories inside the aquifer
            % option 1 all integrated
            speed_option='slow';
            obj=obj.compute_trajectories3(velocity,block_size,x_S,x_Q,speed_option);%compute_trajectories2(velocity,block_size,x_S,x_Q,Bool_sat);
            if(sum(RF_spat(:)<0)>0)
                min_index=find(RF_spat(:)==min(RF_spat(:)));
                fprintf(strcat('WARNING: \n',num2str(sum(RF_spat(:)<0)),' values in the spatialized matrix representing Return Flow \n','were negative (minimal values found:',num2str(RF_spat(min_index)),'). They have been replaced by 0.\n'));
                RF_spat(RF_spat<0)=0;
            end
            [~,Flux_in_spat]=compute_infiltration(run_boussinesq.simulation_results,run_boussinesq.boussinesq_simulation);
%             RF_spat(1,:)=0;
            [obj,Error_RF_DGW,Error_ET]=obj.cut_trajectory_ET_RF(block_size,x_S,x_Q,Flux_in_spat,RF_spat,speed_option);%cut_trajectory_seepage(block_size,x_S,x_Q,width,RF_spat);%ob
            
        end
        
        function obj=post_transport(folder_dir)
            directory=[folder_dir,'\Simulations'];
            directory_list=list_folder_of(directory);
            compt=1;
            for i=1:length(directory_list)
                A=dir(fullfile(directory,directory_list{i},'*.mat'));
                if(~isempty(A))
%                     output_list_directory{compt}=fullfile(directory,directory_list{i},A.name);
                    bouss_sim_dir=fullfile(directory,directory_list{i},'boussinesq_simulation.mat');
                    sim_dir=fullfile(directory,directory_list{i},'simulation_results.mat');
                    load(bouss_sim_dir); load(sim_dir);
                    runs1=runs; runs1.boussinesq_simulation=bouss_sim; runs1.simulation_results=obj;
                    clear bouss_sim obj
                    transport_1D_par=transport_1D_par.test(runs1);
                    save(fullfile(directory,directory_list{i},'transport_1D_par.mat'),'-v7.3');
                    clear transport_1D_par runs1
                end
            end
        end
        
        function Velocity=velocity_field_1D(T,XQ,velocity,t,y) % integrated_parsec = ( - recharge + velocity_horiz .* grad_h)./ h;
            % preparation of the coefficients for the 2D interpolation
            idx_time=findidxmex(T,t); %             idx_time=find_idx(t,T);
            idx_space=findidxmex(XQ,y);%             idx_space=find_idx(x,XQ);

            % C++ bilinear interpolation with a mex file
            Velocity=splinterp2(velocity,idx_time*ones(size(idx_space)),idx_space);
        end
        
        function [mean_,std_]=compute_mean(time_support,weighted_pdfs)
            mean_=trapz(time_support/(24*365*3600),weighted_pdfs*(24*365*3600).*time_support/(24*365*3600),2);
            std_=sqrt(trapz(time_support/(24*365*3600),(time_support/(24*365*3600)).^2.*weighted_pdfs*24*365*3600,2)-mean_.^2);
        end
        
        function median_=compute_median(time_support,weighted_pdfs)
            [~,Index_]=min(abs(0.5-cumtrapz(time_support/(24*365*3600),weighted_pdfs*(24*365*3600),2)),[],2);
            median_=time_support(Index_)/(3600*24);
        end
        
        function prop_less_3_months=compute_Fyw(time_support,weighted_pdfs)
            idx_time=find_idx(3*30,time_support/(3600*24));
            prop_time=(idx_time-floor(idx_time));
            prop_less_3_months=prop_time*trapz(time_support(1:ceil(idx_time)),weighted_pdfs(:,1:ceil(idx_time)),2)+(1-prop_time)*trapz(time_support(1:floor(idx_time)),weighted_pdfs(:,1:floor(idx_time)),2);
        end
        
        function Solute=compute_solute(time_support,weighted_pdfs,k,Ceq,C1)
            Solute_discretized=C1+Ceq*(1-exp(-time_support*k/Ceq));
            Solute=trapz(time_support/(24*365*3600),weighted_pdfs*(24*365*3600).*Solute_discretized,2);
        end
    end
end