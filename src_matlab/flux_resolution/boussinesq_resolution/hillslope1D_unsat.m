classdef hillslope1D_unsat < hillslope1D
    properties(Access=public)
        phi % unsaturated total porosity
    end
    
    methods(Access=public)
        % constructor
        function obj=hillslope1D_unsat
        end
        
        function obj=set_properties(obj,Id,f,k,phi)
            if(nargin<5)
                phi=0.4;
                if(nargin<4)
                    k=1/3600; %fprintf('default hydraulic conductivity 1 m/h \n');
                    if(nargin<3)
                        f=0.3; %fprintf('default porosity: f=0.3 \n'); %if(nargin<2) i=0; fprintf('flat hillslopes choosen by default: i=0 \n'); end
                        if(nargin<2)
                            Id=-1;
                        end
                    end
                end
            end
            obj.Id=Id;
            obj.f=f;
            obj.k=k;
            obj.phi=phi;
        end
        
        function obj=set_hydraulic_parameters(obj,k,f,phi)
            obj.k=k/(3600);
            obj.f=f;
            obj.phi=phi;
        end
        
        function [f,k,phi]=get_hydraulic_properties(obj)
            f=obj.f;
            k=obj.k;
            phi=obj.phi;
        end
    end

    methods(Static)
        function [w]=generate_troch_hillslope(x,n,omega,ymax)
            L=x(end); Lmin=x(1);
            H=5.01;
            if(omega>0)
                if(n==2)
                    K=ymax/(L^(omega*L^n/H));
                    w=2*K*x.^(omega*L^n/H);
                else
                    K=ymax*exp(-2*omega*L^2/(n*(2-n)*H));
                    w=2*K*exp(x.^(2-n)*2*omega*L^n/(n*(2-n)*H));
                end
            elseif(omega<0)
                if(n==2)
                    K=ymax/(Lmin^(omega*L^n/H));
                    w=2*K*x.^(omega*L^n/H);
                else
                    K=ymax*exp(-2*omega*L^n*Lmin^(2-n)/(n*(2-n)*H));
                    w=2*K*exp(x.^(2-n)*2*omega*L^n/(n*(2-n)*H));
                end
            else
                
            end
        end
        
        function [x,w,i,z_true,z_mod]=read_hs1D_data(hillslope_struct_file_path)
            fid=fopen(hillslope_struct_file_path);
            if(fid>0)
                data_hillslope=dlmread(hillslope_struct_file_path,'\t',3,0);
                x=data_hillslope(:,1);
                w=data_hillslope(:,2);
                i=data_hillslope(:,3);
                z_true=data_hillslope(:,4);
                z_mod=data_hillslope(:,5);
            else
                fprintf('Cannot open hillslope structure data file where lies morphologic information \n');
                x=nan;
                w=nan;
                i=nan;
                z_true=nan;
                z_mod=nan;
            end
            fclose(fid);
        end
        
        function [f,k,d]=read_hillslope_parametrization_data(hillslope_param_file_path)
            fid=fopen(hillslope_param_file_path);
            if(fid>0)
                param_hillslope=dlmread(hillslope_param_file_path,'\t',3,0);
                f=param_hillslope(1,1);
                k=param_hillslope(1,2);
                d=param_hillslope(1,3);
            else
                fprintf('Cannot open hillslope structure data file where lies morphologic information \n');
                f=nan;
                k=nan;
                d=nan;
            end
        end
        
        function obj=load_hs1D_class_from_txt_file(hillslope_struct_file_path,hillslope_param_file_path)
            [x,w,i,~,~]=hillslope1D.read_hs1D_data(hillslope_struct_file_path);
            obj=hillslope1D;
            obj=obj.set_properties(-10,-10,-10);
            if(nargin>1)
                [f,k,d]=hillslope1D.read_hillslope_parametrization_data(hillslope_param_file_path);
                d=d*ones(size(x));
                obj=obj.set_spatial_parameters(x,w,i,d);
                obj=obj.set_hydraulic_parameters(k,f);
            else
                obj=obj.set_spatial_parameters(x,w,i);
            end
        end
        
        % #JM Obsolete?
        function generate_Troch_hs1D(number_hillslopes,number_parameters)
            omega={[-1e-3,1e-3],'lin'};
            n={[0.1,4],'lin'};
%             folder_initial='C:\Users\Jean\Documents\Données\SyntheticDataHillslope';
            folder_initial='C:\Users\Jean\Documents\Données\SyntheticData';
            n_param=(number_hillslopes)^0.5;
            p{1}=parameter_grid(omega{1}(1),omega{1}(2),n_param,strcmp(omega{2},'log'));
            p{2}=parameter_grid(n{1}(1),n{1}(2),n_param,strcmp(n{2},'log'));
            grid_list_hillslope_parameters=parameter_grid_list(p);
            val=grid_list_hillslope_parameters.cross_product;
            for i=1:length(val)
%                 if(val(i,1)>0)
%                     folder_directory=strcat(folder_initial,'\ConvergentHillslope');
%                 elseif(val(i,1)<0)
%                     folder_directory=strcat(folder_initial,'\DivergentHillslope');
%                 else
%                     folder_directory=strcat(folder_initial,'\StraightHillslope');
%                 end
%                 folder_directory=strcat(folder_directory,'\','omega_',num2str(val(i,1)),'_n_',num2str(val(i,2)));
                folder_directory=folder_initial;
                hs1D=hillslope1D;
                hs1D=hs1D.set_properties(i,-1,-1);
                x=0:1:100; ymax=25;
                w=hs1D.generate_troch_hillslope(x,val(i,2),val(i,1),ymax);
                hs1D=hs1D.set_spatial_parameters(x,w,-1,-1);
                geologic_input_set.run_simulations(hs1D,folder_directory,number_parameters);
            end
        end
    end
end