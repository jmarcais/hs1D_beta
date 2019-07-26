classdef ttds
    properties(Access=public)
        time_support     % [1xN_t] time support of the pdf [s]
        pdfs             % [N_tsampledxN_t] transit_time distributions matrix storing for each sampled time (each row) the  corresponding pdf [s-1]
        means            % [N_tsampledx1] mean transit time of the distribution [s]
        stds             % [N_tsampledx1] mean transit time of the distribution [s]
        fyw              % [N_tsampledx1] youg water proportion (less than 3 months old) [-]
        sampling_time    % [N_tsampledx1] sampling time where 
    end

    properties(Access=private)
        transit_times    % {N_t_sampledx1} cell matrix storing the corresponding particles individual transit times
        travel_distances % {N_t_sampledx1} cell matrix storing the corresponding particles travel_distances
        weights          % {N_t_sampledx1} cell matrix storing the corresponding particles weights
    end
    
    methods(Access=public)
        function obj=ttds
            obj.time_support=[];
            obj.pdfs=[];
        end
        
        function obj=instantiate_ttds(obj,time_support)
            if(nargin<2)
                obj.time_support=[linspace(0,95,95*4+1)*3600,linspace(4,364,361)*24*3600,linspace(1,100,100)*24*3600*365];
            else
                obj.time_support=time_support;
            end
        end
        
        function obj=store_cell_matrices_data(obj,times_out,transit_times,weights,travel_distances)
            [~,~,X] = unique(times_out);
            obj.transit_times=accumarray(X,1:size(transit_times,1),[],@(r){transit_times(r)});
            obj.travel_distances=accumarray(X,1:size(travel_distances,1),[],@(r){travel_distances(r)});
            obj.weights=accumarray(X,1:size(weights,1),[],@(r){weights(r)});
        end
        
        function obj=compute_ttds(obj,times_out,transit_times,weights,travel_distance)
            % Gelhar et al. 1992 WRR : alpha_disp/transit_distance = 0.04
            alpha_disp=0.04*travel_distance;
            pdf=@(t)(bsxfun(@rdivide,(transit_times./(4*pi*alpha_disp./travel_distance)),t.^3)).^0.5.*...
                exp(-(bsxfun(@minus,transit_times,t)).^2./((bsxfun(@times,4*alpha_disp.*transit_times./travel_distance,t))));
            
            if(obj.time_support(1)~=0)
                transit_times(transit_times==0)=obj.time_support(1);
            end
            weighted_pdf=pdf(obj.time_support);   
            % #JM a better way to handle this rather than a Dirac ??
            if(obj.time_support(1)==0)
                weighted_pdf(:,1)=0;
                weighted_pdf(transit_times==0,1)=2/(obj.time_support(2)-obj.time_support(1));
                weighted_pdf(transit_times==0,2:end)=0;
            end
            weighted_pdf=bsxfun(@rdivide,weighted_pdf,trapz(obj.time_support,weighted_pdf,2));
            weighted_pdf=bsxfun(@times,weighted_pdf,weights);
            
            [obj.sampling_time,~,X] = unique(times_out);
            
            weights_summed=accumarray(X,weights,[],@sum);
            %             obj.pdfs=accumarray(X,weighted_pdf,[],@sum);
            [xx, yy] = ndgrid(X,1:size(weighted_pdf,2));
            obj.pdfs=accumarray([xx(:) yy(:)],weighted_pdf(:));
            obj.pdfs=bsxfun(@rdivide,obj.pdfs,weights_summed);
        end
        
        function obj=compute_ttds_moments(obj)
            obj.means=trapz(obj.time_support,obj.pdfs.*obj.time_support,2);
            obj.stds=sqrt(trapz(obj.time_support,obj.time_support.^2.*obj.pdfs,2)-obj.means.^2);
        end
        
        function obj=compute_youngwaterfraction(obj)
            threemonths_in_sec=365.25/12*3*24*3600;
            k=find_idx(threemonths_in_sec,obj.time_support);
            fyw1=trapz(obj.time_support(1:round(k)),obj.pdfs(:,1:round(k)),2);
            obj.fyw=trapz(obj.time_support(1:round(k)+1),obj.pdfs(:,1:round(k)+1),2);
            obj.fyw=(k-round(k))*fyw1+(round(k)+1-k)*obj.fyw;
        end
    end
    
    methods(Static)
        function obj=retrieve_ttds(times_out,transit_times,weights,travel_distances,time_support)
            transit_times=transit_times(~isnan(times_out));
            weights=weights(~isnan(times_out));
            travel_distances=travel_distances(~isnan(times_out));
            times_out=times_out(~isnan(times_out));
            if(nargin<5)
                time_support=[linspace(0,95,95*4+1)*3600,linspace(4,364,361)*24*3600,linspace(1,100,100)*24*3600*365];
            end
            obj=ttds;
            obj=obj.instantiate_ttds(time_support);
            obj=obj.compute_ttds(times_out,transit_times,weights,travel_distances);
            obj=compute_ttds_moments(obj);
            obj=compute_youngwaterfraction(obj);
            obj=obj.store_cell_matrices_data(times_out,transit_times,weights,travel_distances);
        end
    end
end