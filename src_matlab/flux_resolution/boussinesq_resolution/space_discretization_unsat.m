classdef space_discretization_unsat < space_discretization
% class defining the discretization in space. Discretization has not to be regular (linearly spaced)

    properties(Access=private)
        
    end
    
    methods(Access=public)
        function obj=space_discretization_unsat
        end
        
        function obj=set_space_discretization_properties(obj,xmin,xmax,Nx,discretization_type,xcustom)
            
        end
           
        function [obj_new,hs1D]=resample_hs1D_spatial_variables_unsat(obj,x,w,soil_depth,angle,f,k,phi)
           
        end
        
        function [x_S,w_resampled,soil_depth_resampled,angle_resampled,x,f,k,f_edges,phi]=get_resampled_variables(obj)
           
        end
    end
end