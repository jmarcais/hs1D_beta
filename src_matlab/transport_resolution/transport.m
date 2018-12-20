classdef transport
    properties(Access=public)
        t           % [1 x Nt] array [s]
        N           % [1 x Nt] array containing recharge time series corresponding to the time array [m/s]
        x_traj      % [(size(t_inj)*N_x) x Nt]  matrix containing the trajectories on synchronised timesteps [m]
        t_inj       % [1 x size(t_inj)] array containing the timesteps when injection is done (when N~=0)
        t_inj_pos
        N_inj       % [(size(t_inj)*N_x) x 1] array containing the recharge corresponding to individual trajectory
        x           % [N_x x 1] array containing the river distance aray x [m]
        weight      % [(size(t_inj)*N_x) x 1] array containing the weight to apply to the trajectories
        ET_out      % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via ET
        DGW_out     % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via deep groundwater
        Seep_out    % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via seepage
        DP_out      % [(size(t_inj)*N_x) x 1] boolean array tagging the particles that are retrieved from the model via direct precipitations
    end
    
    methods(Access=public)
        % precipitation time series (t is the time series of the simulation)
        function obj=transport(t,N)        
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
                    N_max=max(obj.N);
                    obj.t_inj=obj.t(N_max>threshold);
                    obj.N_inj=obj.N(:,N_max>threshold);
                    obj.N_inj=reshape(obj.N_inj,size(obj.N_inj,1)*size(obj.N_inj,2),1);
                    obj.t_inj_pos=t_inj_pos(N_max>threshold);
                else
                    obj.t_inj=obj.t(obj.N>threshold);
                    obj.N_inj=obj.N(obj.N>threshold);
                    obj.t_inj_pos=t_inj_pos(obj.N>threshold);
                end
                if(obj.t_inj(end)==obj.t(end))
                    obj.t_inj=obj.t_inj(1:end-1);
                    obj.t_inj_pos=obj.t_inj_pos(1:end-1);
                    obj.N_inj=obj.N_inj(1:end-1);
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
        
%         % load the injection time series (time at which particles will be injected)
%         function obj=compute_t_inj2(obj,threshold)
%             if(nargin<2)
%                 threshold=0;
%             end
%             t_min=15*3600*24; t_max=max(obj.t)-15*3600*24;
%             % one per month
%             spacing_=30*24*3600;
%             obj.t_inj=t_min:spacing_:t_max;
% %             t_inj_pos=1:1:length(obj.t);
%             Distance_=pdist2((obj.t)',(obj.t_inj)');
%             [min_ Idx_]=min(Distance_,[],1);
%             obj.t_inj_pos=Idx_;
%             window_width=floor(spacing_/(obj.t(2)-obj.t(1)));
% %             N_sum_month=movsum(obj.N,window_width);
%             N_sum_month=conv(obj.N, ones(1, window_width),'same');
%             obj.N_inj=N_sum_month(Idx_);
%         end
        
        % compute the weighting factor of the particles path
        function obj=compute_weight(obj,discretized_area) % discretized area is the area of the element dx where transport is computed (equal to dx.*w where w is the width function of the hillslope)
            % convert N_inj in mm/d (~1mm/d in average)
% % %             dt=obj.t(2:end)-obj.t(1:end-1);
% % %             dt=(dt(2:end)+dt(1:end-1))/2;
% % %             dt=[dt(1),dt];
% % %             N_inj=obj.N_inj.*dt*(1000);
            size_N=size(obj.N);
            if(size_N(1)>1)
                dt=obj.t(2:end)-obj.t(1:end-1);
                t_edge1=obj.t-[dt(1),dt]/2;
                t_edge2=obj.t+[dt,dt(end)]/2;
                dt=t_edge2-t_edge1;
                dt_inj=dt(obj.t_inj_pos);
                dt_inj=dt_inj';
                dt_inj=repelem(dt_inj,length(obj.x));
                dflux=dt_inj.*obj.N_inj*1000;
                % % %             weight_temp=discretized_area*N_inj;
                weight_temp=repelem(discretized_area,length(obj.t_inj)).*dflux;
                obj.weight=weight_temp;
            else
                dt_inj=obj.t_inj(2:end)-obj.t_inj(1:end-1);
                t_inj_edge1=obj.t_inj-[dt_inj(1),dt_inj]/2;
                t_inj_edge2=obj.t_inj+[dt_inj,dt_inj(end)]/2;
                dt_inj=t_inj_edge2-t_inj_edge1;
                dflux=dt_inj.*obj.N_inj*1000;
                % % %             weight_temp=discretized_area*N_inj;
                weight_temp=discretized_area*dflux;
                obj.weight=reshape(weight_temp,length(dflux)*length(discretized_area),1);
            end
        end
        
        function obj=instantiate_exit_tags(obj)
            obj.DP_out=zeros(size(obj.weight));
            obj.Seep_out=zeros(size(obj.weight));
            obj.DGW_out=ones(size(obj.weight));
            obj.ET_out=zeros(size(obj.weight));
        end
        
        % compute the trajectories of the flowpaths
%         function obj=compute_trajectories(obj,sol_simulated,block_size,x_S,x_Q,folder_mex_option)
        function obj=compute_trajectories(obj,velocity,block_size,x_S,x_Q,folder_mex_option)
            %#JM to change in a future release threslhold should be defined outside the method
            threshold=0;
            obj.x=x_S;
            if(nargin>5)
                %% C code in mex file
                load(fullfile(folder_mex_option,'velocity_field.mat'));
                load(fullfile(folder_mex_option,'t_inj.mat'));
                load(fullfile(folder_mex_option,'x_inj.mat'));
                obj.x_traj=solve_ode_trajectory(Velocity_2,t_inj,x_inj);
            else
% % %                 Storage=@(t,y) nakeinterp1(x_S,deval(sol_simulated,t,1:block_size),y);
% % %                 
% % %                 Flux_reg=@(t,y) nakeinterp1([-1;x_Q],[0;deval(sol_simulated,t,block_size+1:2*block_size+1)],y);
% % %                 Storage_reg=@(t,y) Storage(t,y)+1e-3.*double(Storage(t,y)<=0);
% % %                 Velocity_reg=@(t,y) Flux_reg(t,y)./Storage_reg(t,y);
                [T,XQ]=meshgrid(obj.t,x_Q);
%                 Velocity_reg = @(t,y) qinterp2(T,XQ,velocity,t,y);
%                 Velocity_reg = @(t,y) interp2(T,XQ,velocity,t,y,'linear',0);
                Velocity_reg = @(t,y) qinterp2(T,XQ,velocity,t*ones(size(y)),y,2);
                S_mat=sparse(diag(ones(block_size,1)));
                events=@(t,y)obj.eventfun(t,y,x_Q(2));
                options_reg = odeset('Vectorized','on','Jpattern',S_mat,'Events',events);

                size_row=length(obj.t_inj)*length(obj.x);
                size_column=length(obj.t);
                % choose the format of x_traj if not problem for matlab for allocating memory
                if(size_row*size_column<5e9)
% % %#20180116                                             obj.x_traj=nan(size_row,size_column);
%                     obj.x_traj=sparse(size_row,size_column);
                    mat_pos_allocate=cell(length(obj.t_inj_pos),1);
                    for i=1:length(obj.t_inj_pos)
                        pos_temp=obj.t_inj_pos(i);
                        size_N=size(obj.N);
                        if(size_N(1)>1)
                            Bool_inj_spat=obj.N(:,pos_temp)>threshold;
                            [~,traj_temp2] = ode45(Velocity_reg,obj.t(pos_temp:end),x_S(Bool_inj_spat),options_reg);
                            traj_temp2=traj_temp2';
                            traj_temp=zeros(length(Bool_inj_spat),size(traj_temp2,2));
                            traj_temp(Bool_inj_spat,:)=traj_temp2;
                        else
                            [~,traj_temp] = ode45(Velocity_reg,obj.t(pos_temp:end),x_S,options_reg);
                            traj_temp=traj_temp';
                        end
                        if(length(obj.t(pos_temp:end))==2)
                            traj_temp=traj_temp(:,[1,end]);
                        end
% % %#20180116                                             obj.x_traj(block_size*(i-1)+1:block_size*i,pos_temp:size(traj_temp,1)+pos_temp-1)=traj_temp';
                        AA=(combvec(block_size*(i-1)+1:block_size*i,pos_temp:size(traj_temp,2)+pos_temp-1))';
                        bool_delete=traj_temp<x_Q(2);
                        bool_delete=logical(cumsum(bool_delete,2));
                        bool_delete=[false(length(obj.x),1),bool_delete(:,1:end-1)];
                        traj_temp(bool_delete)=0;
                        traj_temp=traj_temp(:);
                        mat_pos_allocate{i}=[AA,traj_temp];
                        % delete particles arriving in the river
%                         bool_delete=traj_temp<x_Q(2);
                        bool_delete=traj_temp==0;
                        mat_pos_allocate{i}=mat_pos_allocate{i}(~bool_delete,:);
                        % #JM comment to gain speed
% % %                         fprintf(strcat(num2str(i),'/',num2str(length(obj.t_inj)),'\n'));
                    end
                      mat_pos_allocate=vertcat(mat_pos_allocate{:});
% %                     bool_delete=mat_pos_allocate(:,3)<1e-3;
% %                     mat_pos_allocate=mat_pos_allocate(~bool_delete,:);
                     obj.x_traj=sparse(mat_pos_allocate(:,1),mat_pos_allocate(:,2),mat_pos_allocate(:,3),size_row,size_column);
% % % % % % % % % % % % %                     obj.x_traj=sparse(mat_pos_allocate(1:90084959,1),mat_pos_allocate(1:90084959,2),mat_pos_allocate(1:90084959,3),size_row,size_column);
% % % % % % % % % % % % %                     x_traj2=sparse(mat_pos_allocate((1+90084959):end,1),mat_pos_allocate((1+90084959):end,2),mat_pos_allocate((1+90084959):end,3),size_row,size_column);
% % % % % % % % % % % % %                     obj.x_traj=obj.x_traj+x_traj2;
%                     obj.x_traj=sparse(obj.x_traj);
                else
                    fprintf('Warning due to memory allocation problems, x_traj will be in cell format \n');
                    obj.x_traj=cell(size_row,1);
                    for i=1:100
                        %                 for i=1:length(obj.t_inj_pos)
                        pos_temp=obj.t_inj_pos(i);
                        %                   [~,traj_temp] = ode15s(Flux,obj.t(pos_temp:end),x_S,options);
                        %                   [~,traj_temp] = ode15s(Flux_reg,obj.t(pos_temp:end),x_S,options_reg);
                        [~,traj_temp] = ode45(Velocity_reg,obj.t(pos_temp:end),x_S,options_reg);
                        temp_mat=nan(length(x_S),pos_temp-1);
                        traj_temp2=[temp_mat,traj_temp'];
                        row_cut=ones(1,length(x_S)); column_cut=[length(obj.t)];
                        traj_temp_cell=mat2cell(traj_temp2,row_cut,column_cut);
                        obj.x_traj(block_size*(i-1)+1:block_size*i)=traj_temp_cell;
                        % #JM comment to gain speed
% % %                         fprintf(strcat(num2str(i),'/',num2str(length(obj.t_inj)),'\n'));
                    end
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
        function obj=cut_trajectory_seepage(obj,Discretized_Aquifer_Volume,block_size,x_S,x_Q,w,RF_spat,Subsurface_flux,x_traj_soil)
            Seepage=RF_spat;
            Seepage(Seepage<0)=0;
            Subsurface_flux2=Subsurface_flux;
            NetPrecip=obj.N;
            dx_Q=x_Q(2:end)-x_Q(1:end-1);
            area_spatialized=w.*dx_Q;
            if(size(obj.N,1)==1)
                NetPrecip=area_spatialized*NetPrecip;
            else
                NetPrecip=bsxfun(@times,NetPrecip,area_spatialized);
            end
%             
%             area_spatialized=w;
%             area_spatialized= repmat(area_spatialized,1,size(obj.N,2));
%             precip=obj.N;
%             if(size(obj.N,1)==1)
%                 precip=repmat(precip,length(w),1);
%             end
%             infilt_spat=precip.*area_spatialized;
%             Seepage=deval(sol_simulated,obj.t,2*block_size+2:3*block_size+1)-infilt_spat;
%             Seepage(Seepage<0)=0;
%             
%             Subsurface_flux2=deval(sol_simulated,obj.t,block_size+1:2*block_size+1);
            Subsurface_flux2=[zeros(1,size(obj.N,2));Subsurface_flux2];
            x_Q2=[-1;x_Q];
            Subsurface_flux=zeros(length(x_S),size(obj.N,2));
            for i=1:size(obj.N,2)
                Subsurface_flux(:,i)=nakeinterp1(x_Q2,Subsurface_flux2(:,i),x_S);
            end
            %#JM_IPGP #JM test if the normalization should happen to the volume
            dt=obj.t(2:end)-obj.t(1:end-1);
            t_edge1=obj.t-[dt(1),dt]/2;
            t_edge2=obj.t+[dt,dt(end)]/2;
            dt=t_edge2-t_edge1;
            Discretized_Aquifer_Volume_m3s=bsxfun(@rdivide,Discretized_Aquifer_Volume,dt);
            Negative_NetPrecip_proportion=zeros(size(NetPrecip));
            Negative_NetPrecip_proportion(NetPrecip<0)=-NetPrecip(NetPrecip<0)./(Discretized_Aquifer_Volume_m3s(NetPrecip<0));
            Negative_NetPrecip_proportion(isnan(Negative_NetPrecip_proportion))=0;
            Negative_NetPrecip_proportion(1,:)=0;
            Seepage_proportion=Seepage./(Discretized_Aquifer_Volume_m3s);
            Seepage_proportion(isnan(Seepage_proportion))=0;
            Seepage_proportion(1,:)=0;
            
%             Negative_NetPrecip_proportion=zeros(size(NetPrecip));
%             Negative_NetPrecip_proportion(NetPrecip<0)=-NetPrecip(NetPrecip<0)./(-NetPrecip(NetPrecip<0)+Seepage(NetPrecip<0)+abs(Subsurface_flux(NetPrecip<0)));
%             Seepage_proportion=Seepage./(Seepage+abs(Subsurface_flux)); 
%             Seepage_proportion(NetPrecip<0)=Seepage(NetPrecip<0)./(-NetPrecip(NetPrecip<0)+Seepage(NetPrecip<0)+abs(Subsurface_flux(NetPrecip<0))); 
%             Seepage_proportion(isnan(Seepage_proportion))=0;
%             Seepage_proportion(1,:)=0;
            BB=zeros(1,length(obj.t));
            for i=1:length(obj.t)
% %#JM20180616ChangeTest                 Surface_flux=deval(sol_simulated,obj.t(i),2*block_size+2:3*block_size+1).*dx;
% %                 Seepage=Surface_flux; Seepage(Surface_flux>0)=Seepage(Surface_flux>0)-obj.N(i)*area_spatialized(Surface_flux>0);
% %                 Seepage(Seepage<0)=0;
% %                 Subsurface_flux=nakeinterp1([-1;x_Q],[0;deval(sol_simulated,obj.t(i),block_size+1:2*block_size+1)],x_S);
% %                 Seepage_proportion=Seepage./(Seepage+abs(Subsurface_flux));
% %                 Seepage_proportion(isnan(Seepage_proportion))=0;
                % #JM #JM_IPGP "delete" the particles that are retrieved by ET
                % #JM check if the fact to put the entire line at 0 has for effect not to consider it in the TTD computations
                ET_pos=find(Negative_NetPrecip_proportion(:,i)>0);
                Particle_to_tag=[];
                ET_prop=0;
                for j=1:length(ET_pos)
                    if(~isempty(Particle_to_tag))
                        ET_prop=ET_prop-Weight_cum_ET(Weight_cum_ET<=ET_prop);
                        ET_prop=ET_prop(end);
                        ET_m3s=ET_prop*(-NetPrecip(pos_ET,i))/Negative_NetPrecip_proportion(pos_ET,i);
                    end
                    pos_ET=ET_pos(end-j+1); %pos_ET=ET_pos(j); % 
                    if(~isempty(Particle_to_tag))
                        ET_prop=ET_m3s./Discretized_Aquifer_Volume_m3s(pos_ET,i);
                    end
                    ET_prop=ET_prop+Negative_NetPrecip_proportion(pos_ET,i);
                    ET_prop=Negative_NetPrecip_proportion(pos_ET,i);
                    
                    particle_subject_to_ET=find((obj.x_traj(:,i)<=x_Q(pos_ET+1)) & (obj.x_traj(:,i)>x_Q(pos_ET)));
                    Weight_partial_ET=obj.weight(particle_subject_to_ET);
                    Weight_total_ET=sum(Weight_partial_ET);
                    Weight_partial_ET=Weight_partial_ET./Weight_total_ET;
                    pos_2_ET=mod(particle_subject_to_ET,block_size); pos_2_ET(pos_2_ET==0)=block_size;
                    Initial_infiltration_point_ET=x_S(pos_2_ET);
                    [~,Index_]=sort(Initial_infiltration_point_ET);
                    Weight_cum_ET=cumsum(Weight_partial_ET(Index_));
                    Particle_to_tag=particle_subject_to_ET(Index_(Weight_cum_ET<=ET_prop));
%                     obj.x_traj(Particle_to_tag,1:end)=0;
                     obj.ET_out(Particle_to_tag)=obj.ET_out(Particle_to_tag)+1;
                     obj.DGW_out(Particle_to_tag)=0;
                     obj.x_traj(Particle_to_tag,i+1:end)=0;
% %                      if(j==length(ET_pos))
% %                          AA=ET_prop-Weight_cum_ET(Weight_cum_ET<=ET_prop);
% %                          if(isempty(AA))
% %                              AA=ET_prop;
% %                          end
% %                          if(~isempty(AA))
% %                              AA=AA(end);
% %                              BB(i)=AA*(-NetPrecip(pos_ET,i))/Negative_NetPrecip_proportion(pos_ET,i);
% %                              fprintf(strcat(num2str(BB(i)),'\n'));
% %                          end
% %                      end
                end
                
                % not considering the "seepage" coming out the river as it is not real seepage
                
                Seepage_pos=find(Seepage_proportion(:,i)>0);
                Particle_Position_to_delete=[];
                Seep_prop=0;
                for j=1:length(Seepage_pos)
                    if(~isempty(Particle_Position_to_delete))
%                         Seep_prop=Seep_prop+Seepage_proportion(pos_,i);
%                     else
                        Seep_prop=Seep_prop-Weight_cum(Weight_cum<=Seep_prop);
                        Seep_prop=Seep_prop(end);
                        Seep_m3s=Seep_prop*(Seepage(pos_,i))/Seepage_proportion(pos_,i);
                    end
                    pos_=Seepage_pos(end-j+1); %pos_=Seepage_pos(j); %
                    if(~isempty(Particle_Position_to_delete))
                        Seep_prop=Seep_m3s./Discretized_Aquifer_Volume_m3s(pos_,i);
                    end
                    Seep_prop=Seep_prop+Seepage_proportion(pos_,i);
                    Seep_prop=Seepage_proportion(pos_,i);
                    
%                     x_S_seep=x_S(pos_);
                    particle_subject_to_seep=find((obj.x_traj(:,i)<=x_Q(pos_+1)) & (obj.x_traj(:,i)>x_Q(pos_)));
                    Weight_partial=obj.weight(particle_subject_to_seep);
                    Weight_total=sum(Weight_partial);
                    Weight_partial=Weight_partial./Weight_total;
                    pos_2=mod(particle_subject_to_seep,block_size); pos_2(pos_2==0)=block_size;
                    Initial_infiltration_point=x_S(pos_2);
                    [~,Index_]=sort(Initial_infiltration_point);
                    Weight_cum=cumsum(Weight_partial(Index_));
                    Particle_Position_to_delete=particle_subject_to_seep(Index_(Weight_cum<=Seep_prop));
                    obj.Seep_out(Particle_Position_to_delete)=obj.Seep_out(Particle_Position_to_delete)+1;
                    obj.DGW_out(Particle_Position_to_delete)=0;
                    
% % %#20180116                                             obj.x_traj(Particle_Position_to_delete,i+1:end)=nan;
                    if(nargin>8)% && ~isnan(x_traj_soil))
                        obj.x_traj(Particle_Position_to_delete,i+1:end)=x_traj_soil(Particle_Position_to_delete,i+1:end);
                    else
                        obj.x_traj(Particle_Position_to_delete,i+1:end)=0;
                    end
                end
                
%                 Distance=pdist2(obj.x_traj(:,i),x_S);
%                 [Value_,Distance_pos]=nanmin(Distance,[],2);
%                 Distance_pos(isnan(Value_))=nan;
%                 Seepage_proportion_x_traj=nan(size(Distance_pos));
%                 Seepage_proportion_x_traj(~isnan(Distance_pos))=Seepage_proportion(Distance_pos(~isnan(Distance_pos)));
%                 seepage_proportion=deval(sol_simulated,t,2*block_size+2:3*block_size+1)
            end
            % % %#20180116
            obj.x_traj=sparse(obj.x_traj);
        end
        
        function obj=update_DGW(obj,x_Q)
            Bool_fin=obj.x_traj(:,end)>x_Q(2);
            obj.DGW_out(Bool_fin)=0;
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
                x_traj_temp=obj.x_traj(Bool_sat,:);
                % % %#20180116             detect_no_nan_values=~isnan(x_traj_temp);
                detect_no_nan_values=(x_traj_temp~=0);
                size_2=size(detect_no_nan_values);
                %#20180124 possibility to comment next line not to take into account direct precipitations onto saturated areas
                detect_no_nan_values=[false(size_2(1),1),detect_no_nan_values(:,1:end-1)];
                % % %#20180116                         x_traj_temp(detect_no_nan_values)=nan;
                x_traj_temp(detect_no_nan_values)=0;
                obj.x_traj(Bool_sat,:)=x_traj_temp;
                obj.DP_out(Bool_sat)=obj.DP_out(Bool_sat)+1;
                obj.DGW_out(Bool_sat)=0;
                %%%#20180116
                obj.x_traj=sparse(obj.x_traj);
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
        
        function [t_out,delta_t,x_fin]=get_trajectory_properties(obj)
            [row,col] = find(obj.x_traj);
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
            
            size_trajectories=size(obj.x_traj);
            AA=(1:1:size_trajectories(1))';
            linearInd_fin = sub2ind(size_trajectories,AA(max_>0),max_(max_>0));
            linearInd_init = sub2ind(size_trajectories,AA(min_>0),min_(min_>0));
%             linearInd_fin = sub2ind(size_trajectories,(1:1:size_trajectories(1))',max_);
%             linearInd_init = sub2ind(size_trajectories,(1:1:size_trajectories(1))',min_);
            x_fin=nan(size(max_));
            x_fin(max_>0)=obj.x_traj(linearInd_fin);
            x_fin(max_(max_>0)==length(obj.t))=nan;
            x_init=nan(length(min_),1);
            x_init(min_>0)=obj.x_traj(linearInd_init);
            delta_x_groundwater=x_fin-x_init;
        end
        
        function [transit_times,weights_]=get_ttds_bis(obj,sample_time,t_out,delta_t)
            if(nargin<3)
                [t_out,delta_t]=get_trajectory_properties(obj);
            end
            [~,Index_time]=min(abs(obj.t-sample_time));
            transit_times=delta_t(t_out==obj.t(Index_time));
            weights_=obj.weight(t_out==obj.t(Index_time));
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
        
        function [time_support,weighted_pdf,mean_,std_]=compute_full_distribution(obj,transit_times,weights_,transit_distance,travel_time_soil)
            if(nargin>3)
                travel_time_particle_soil=interpn(obj.x-obj.x(1),travel_time_soil,transit_distance);
                transit_times(transit_times==0)=travel_time_particle_soil(transit_times==0);
                time_support=(0:0.01:1000)*24*365*3600;
            else
                time_support=(0:0.1:1000)*24*365*3600;
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
                compute_pdf=pdf(time_support);
                compute_pdf(:,1)=0;
                compute_pdf(transit_times==0,1)=2/(time_support(2)-time_support(1));
                compute_pdf(transit_times==0,2:end)=0;
                if(length(weights_)~=1)
                    weight_total=sum(weights_);
                    weighted_pdf=1/weight_total*nansum(bsxfun(@times,compute_pdf,weights_));
                else
                    weighted_pdf=compute_pdf;
                end
                %             figure;
                weighted_pdf_area=trapz(time_support,weighted_pdf);
% % % %                 weighted_pdf=weighted_pdf/weighted_pdf_area;
                %             plot(time_support/(24*365*3600),weighted_pdf*24*365*3600/weighted_pdf_area);
                % #JM comment to gain speed
                mean_=trapz(time_support/(24*365*3600),time_support/(24*365*3600).*weighted_pdf*24*365*3600/weighted_pdf_area);
                std_=sqrt(trapz(time_support/(24*365*3600),(time_support/(24*365*3600)).^2.*weighted_pdf*24*365*3600/weighted_pdf_area)-mean_^2);
                % % %             fprintf(strcat('mean: ',num2str(mean_),'years \n','std: ',num2str(std_),'years \n'));
                %             plot(time_support/(24*3600),weighted_pdf*24*3600/weighted_pdf_area);
                %             pdf=@(t)(0.25.*10.*(transit_times./t).^3).^0.5.*exp(-0.25*10.*(t-transit_times).^2./(t.*transit_times));
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
            obj=transport(t,N_in);
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

        function obj=transport_with_rooting(runs,x,slope_angle,k_soil)
            [obj,x_S,x_Q,width,velocity,Bool_sat,sol_simulated]=transport.instantiate_transport_and_velocity_field(runs);
            block_size=length(x_S);
            
            [~,~,DPSA_spat,RF_spat]=compute_DPSA_RF(runs.simulation_results,runs.boussinesq_simulation);
            Subsurface_flux=runs.simulation_results.Q;
            
            % compute the trajectories inside the aquifer
            obj=obj.compute_trajectories(velocity,block_size,x_S,x_Q);
            obj=obj.cut_trajectory_saturated_areas(Bool_sat);
            dx=x_Q(2:end)-x_Q(1:end-1);
            Discretized_Aquifer_Volume=bsxfun(@times,runs.simulation_results.S,dx);
            obj=obj.cut_trajectory_seepage(Discretized_Aquifer_Volume,block_size,x_S,x_Q,width,RF_spat,Subsurface_flux);
            obj=obj.update_DGW(x_Q);
            
            % retrieve the particles essential properties for their travel inside the aquifer
            [t_out_groundwater,transit_times_groundwater,x_fin_groundwater]=obj.get_trajectory_properties;
            
            % compute the typical travel time inside the subsurface (soil layer)
            
            slope_angle1_top=interpn(x,slope_angle,x_S);
            k_soil2=ones(size(x_S))*k_soil;
            velocity_soil=k_soil2.*sin(slope_angle1_top);
            travel_time_soil=cumtrapz(dx./velocity_soil);
            travel_time_soil2=[0;travel_time_soil];
            
            % travel time of the particles inside the soil layer
            transit_times_soil=interpn(x_Q,travel_time_soil2,full(x_fin_groundwater),'linear',0);
            
            % t_out_river : time at which the particles arrives finally inside the river 
            t_out_river=cat(2,t_out_groundwater,transit_times_soil); 
            t_out_river= nansum(t_out_river,2); 
            t_out_river=interp1([obj.t,2*obj.t(end)-obj.t(end-1)],[obj.t,2*obj.t(end)-obj.t(end-1)],t_out_river,'nearest',nan);
            transit_times_total=cat(2,transit_times_groundwater,transit_times_soil); transit_times_total= nansum(transit_times_total,2); 
            
            transit_times_total(obj.ET_out>0)=nan;
            t_out_river(obj.ET_out>0)=nan;
            t_out_groundwater(obj.ET_out>0)=nan;
            transit_times_groundwater(obj.ET_out>0)=nan;
            x_fin_groundwater(obj.ET_out>0)=nan;
            transit_times_soil(obj.ET_out>0)=nan;
        end
        
        function [obj,x_S,x_Q,width,velocity,Bool_sat,sol_simulated]=instantiate_transport_and_velocity_field(runs)
            sim_res=runs.simulation_results;
            bouss_sim=runs.boussinesq_simulation;
            %             f=bouss_sim.hydraulic_properties.f;
            t=sim_res.t;
            N_in=bouss_sim.compute_recharge(t);
            [~,width,soil_depth,~,~,f]=bouss_sim.discretization.get_resampled_variables;
            Smax=f.*width.*soil_depth;
            x_S=sim_res.x_S;
            x_Q=sim_res.x_Q;
            dx=x_Q(2:end)-x_Q(1:end-1);

            
            obj=transport(t,N_in);
            obj.x=x_S;
            
            % Bool_sat computation
            relstor=bsxfun(@rdivide,sim_res.S,Smax);
            % first mesh is always considered as saturated (it is the river)
            relstor(:,1)=1;
            Bool_sat=relstor>=0.999925;%0.99995;%0.99999;%0.9999;%0.999;
            
            %             spacing_=1*24*3600;
            spacing_=-1;
            obj=obj.compute_t_inj(spacing_);
            obj=obj.compute_weight(dx.*width);
            obj=obj.instantiate_exit_tags;
            
            S_edges=bouss_sim.discretization.Omega*sim_res.S;
            S_edges(1,:)=1;
            S_edges(end,:)=1;
            velocity=sim_res.Q./S_edges;
            
            [~,c]=find(isnan(velocity));
            if(max(c)>2)
                fprintf('Warning: Nan values in velocity vector \n');
            end
            velocity(isnan(velocity))=0;
            sol_simulated=bouss_sim.sol_simulated;
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
                    transport=transport.test(runs1);
                    save(fullfile(directory,directory_list{i},'transport.mat'),'-v7.3');
                    clear transport runs1
                end
            end
        end
    end
end