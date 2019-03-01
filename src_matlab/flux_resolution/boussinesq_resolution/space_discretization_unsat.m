classdef space_discretization_unsat < space_discretization
% class defining the discretization in space. Discretization has not to be regular (linearly spaced)

    properties(Access=private)
        phi                       % (Nx X 1) total porosity in the unsaturated portion of the hillslope
                                  % by opposition to f which is the drainable porosity
    end
    
    methods(Access=public)
        function obj=space_discretization_unsat
        end
        
        function obj=set_space_discretization_properties(obj,xmin,xmax,Nx,discretization_type,xcustom)
            if(nargin<5) discretization_type='lin'; end
            if(nargin<6) xcustom=nan; end 
            obj.xmin=xmin;
            obj.xmax=xmax;
            obj.Nx=Nx;
            obj.discretization_type=discretization_type;
            switch obj.discretization_type
                case 'lin'
                    obj.x=obj.xmin:(obj.xmax-obj.xmin)/obj.Nx:obj.xmax;
                case 'log'
                    if(obj.xmin==0) 
                        log_x_min=-2;
                        obj.xmin=10^log_x_min;
                    else
                        log_x_min=log10(obj.xmin);
                    end
                    log_x_max=log10(obj.xmax);
                    obj.x=logspace(log_x_min,log_x_max,obj.Nx+1);
                case 'custom'
                    obj.x=xcustom;
                    obj.Nx=length(obj.x)-1;
                    obj.xmin=min(obj.x);
                    obj.xmax=max(obj.x);
                case 'square'
                    obj.x=sqrt(obj.xmin):(sqrt(obj.xmax)-sqrt(obj.xmin))/obj.Nx:sqrt(obj.xmax);
                    obj.x=(obj.x).^2;
                otherwise
                    fprintf(strcat('no discretization type corresponding to ',obj.discretization_type,'  \n'));
            end
            obj.x=obj.x';
%             obj.dx=obj.compute_dx;
            obj=obj.compute_x_centered;
        end
           
        function [obj,Matrix_link]=resample_hs1D_spatial_variables(obj,x,w,soil_depth,angle,f,k,phi)
            Distance=pdist2(x,obj.x_S);
            [~,MinDistPos]=min(Distance,[],2);
            Matrix_link=zeros(length(x),length(obj.x_S));
            
            if(~isempty(MinDistPos))
                Position_distance=(1:1:length(x))';
                Ind_matrix_link=sub2ind(size(Matrix_link),Position_distance,MinDistPos);
                Matrix_link(Ind_matrix_link)=1;
                w_resampled_temp=accumarray(MinDistPos,w,[],@nanmean,double(NaN));
            end
%             obj.w_resampled=interpn(x,w,obj.x_S);
%             obj.soil_depth_resampled=interpn(x,soil_depth,obj.x_S);
%             obj.angle_resampled=interpn(x,angle,obj.x);
%             smooth_width_function = fit(x, w,  'smoothingspline', 'SmoothingParam', 0.0001);
% #JM to change
            smooth_width_function = fit(x, w,  'smoothingspline', 'SmoothingParam', 0.9);
            obj.w_resampled=smooth_width_function(obj.x_S);
%             smooth_slope_function = fit(x, angle,  'smoothingspline', 'SmoothingParam', 0.0001);
% #JM to change
            smooth_slope_function = fit(x, angle,  'smoothingspline', 'SmoothingParam', 0.9);
            obj.angle_resampled=smooth_slope_function(obj.x);
            obj.soil_depth_resampled=interpn(x,soil_depth,obj.x_S);
            
            if(length(k)==1)
                k=k*ones(size(x));
            end
            if(length(f)==1)
                f=f*ones(size(x));
            end
            if(length(phi)==1)
                phi=phi*ones(size(x));
            end
            obj.k=interpn(x,k,obj.x);
            obj.f=interpn(x,f,obj.x_S);
            obj.f_edges=interpn(x,f,obj.x);
            obj.phi=interpn(x,phi,obj.x_S);
        end
        
        function [x_S,w_resampled,soil_depth_resampled,angle_resampled,x,f,k,f_edges,phi]=get_resampled_variables(obj)
            x_S=obj.x_S;
            w_resampled=obj.w_resampled;
            soil_depth_resampled=obj.soil_depth_resampled;
            angle_resampled=obj.angle_resampled;
            x=obj.x;
            f=obj.f;
            k=obj.k;
            f_edges=obj.f_edges;
            phi=obj.phi;
        end
    end
end