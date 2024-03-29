classdef hillslope1D
    properties(Access=public)
        Id           % same as hillslope from where it derives Id=-1 is for hillslope with an unsaturated component
        x            % 1D (Nix1) array distance from stream (m)
        z            % 1D (Nix1) array mean elevation of the hillslope according to the distance (m)
        w            % 1D (Nix1) array width of the hillslope according to the distance x (m)
        z_predicted  % linear elevation increment according to x in the watershed (m)   
        R2           % Quality of the linear fitting (ie approximating z by z_predicted)
        i            % angle (rad) that corresponds to the mean slope of the hillslope 
        soil_depth   % 1D array mean soil depth according to the distance to the stream x (m)
        f            % drainable porosity f constant for the all hillslope (-)
        k            % soil hydraulic conductivity (m/s)
        phi          % unsaturated total porosity
        saturation   % 1D (Nix1) array saturation proportion according to real data proportion of saturation of each storage for a given lateral distance x
        beven_ratio  % A/(w0.sin(i0))
        saturated_area
        saturated_length
        link_hs1D
        volume % m3 of porous medium above the stream base level
    end
    
    methods(Access=public)
        % constructor
        function obj=hillslope1D
        end
        
        function obj=set_properties(obj,Id,f,k,phi)
            if(nargin<4)
                k=1/3600; %fprintf('default hydraulic conductivity 1 m/h \n');
                if(nargin<3)
                    f=0.3; %fprintf('default porosity: f=0.3 \n'); %if(nargin<2) i=0; fprintf('flat hillslopes choosen by default: i=0 \n'); end
                    if(nargin<2)
                        Id=1;
                    end
                end
            end
            obj.Id=Id;
            obj.f=f;
            obj.k=k;
            if(nargin>4)
                obj.phi=phi;
            end
        end
        
        function obj=set_space_default_parameters(obj)
            obj.x=0:1:100;
%             obj.w=ones(size(obj.x));
            obj.w=0.1*(100-obj.x)+0.1;
            obj.soil_depth=2*ones(size(obj.x));
            obj.i=atan(0.05)*ones(size(obj.x));
        end
        
        function obj=set_spatial_parameters(obj,x,w,angle,soil_depth,z)
            obj.x=x;
            obj.w=w;
            obj.i=angle;
            if(nargin>4)
                obj.soil_depth=soil_depth;
            end
            if(nargin>5)
                obj.z=z;
            end
        end
        
        function obj=set_hydraulic_parameters(obj,k,f,phi)
            obj.k=k/(3600);
            obj.f=f;
            if(nargin>3)
                obj.phi=phi;
            end
        end

        function [obj,Matrix_link]=set_parameters(obj,distance,z,DEM_resolution,volume,z0)
            if(nargin<4)
                DEM_resolution=5;
            end
            obj.x=(min(distance):DEM_resolution:max(distance))';
            Distance=pdist2(distance,obj.x);
            [~,MinDistPos]=min(Distance,[],2);
            Matrix_link=sparse(length(distance),length(obj.x));
            if(~isempty(MinDistPos))
                Position_distance=(1:1:length(distance))';
                Ind_matrix_link=sub2ind(size(Matrix_link),Position_distance,MinDistPos);
                Matrix_link(Ind_matrix_link)=1;
                % Number of tile at the same distance from the stream and each tile is 5 m large
                obj.w=accumarray(MinDistPos,1,[],@nansum,double(NaN))*DEM_resolution;
                if(nargin<5)
                    obj.z=accumarray(double(MinDistPos),double(z),[],@nanmean,double(NaN));
                else
                    obj.z=accumarray(double(MinDistPos),double(volume),[],@nansum,double(NaN))./(double(obj.w)*double(DEM_resolution))+double(z0);
                    obj.volume=accumarray(double(MinDistPos),double(volume),[],@nansum,double(NaN));
                end
                obj.x=double(obj.x);
                obj.z=double(obj.z);
                ztemp=obj.z(~isnan(obj.z));
                xtemp=obj.x(~isnan(obj.z));
                if(~isempty(xtemp) && length(xtemp)>1 && license('test', 'curve_fitting_toolbox'))
                    f3 = fit(xtemp, ztemp,  'smoothingspline', 'SmoothingParam', 0.0001);
                    obj.z(isnan(obj.z))=f3(obj.x(isnan(obj.z)));
                end
                wtemp=obj.w(~isnan(obj.w));
                xtemp=obj.x(~isnan(obj.w));
                if(~isempty(xtemp) && length(xtemp)>1 && license('test', 'curve_fitting_toolbox'))
                    f3 = fit(xtemp, wtemp,  'smoothingspline', 'SmoothingParam', 0.0001);
                    obj.w(isnan(obj.w))=f3(obj.x(isnan(obj.w)));
                end
                obj.w(obj.w<5)=5;
%                 obj.z=fillmissing(obj.z,'spline');
%                 obj.w=fillmissing(obj.w,'spline');
            else
                Idstring=num2str(obj.Id);
                fprintf(['cannot set parameters for hillslope ',Idstring ,' \n']);
            end
            % #JM temporary will be extracted from soil data shpfile GIP "sol de bretagne"
            obj.soil_depth=ones(length(obj.x),1);
        end
        
        function obj=set_saturation(obj,distance,sat)
            Distance=pdist2(distance,obj.x);
            [~,MinDistPos]=min(Distance,[],2);
            if(~isempty(MinDistPos))
                obj.saturation=accumarray(MinDistPos,sat,[],@mean);
            else
                Idstring=num2str(obj.Id);
                fprintf(['cannot set saturation for hillslope ',Idstring ,' \n']);
            end
        end
        
        function obj=transform_to_constant_slope(obj)
            if(~isempty(obj.x))
                parameters=polyfit(obj.x,obj.z,1);
                parameters=double(parameters);
                obj.z_predicted=parameters(1)*obj.x+parameters(2);
                obj.R2=sum((obj.z_predicted-mean(obj.z)).^2)/sum((obj.z-mean(obj.z)).^2);
                obj.i=atan((obj.z_predicted(end)-obj.z_predicted(1))/(obj.x(end)-obj.x(1)));
                obj.i=obj.i*ones(size(obj.x));
            else
                Idstring=num2str(obj.Id);
                fprintf(['cannot transform to constant slope for hillslope ',Idstring ,' \n']);
            end
        end
        
        function obj=transform_to_spline_slope(obj,SmoothingParam)
            if(length(obj.x)>2)
                z=double(obj.z);
                x=double(obj.x);
                if(license('test', 'curve_fitting_toolbox'))
                    f3 = fit(x, z,  'smoothingspline', 'SmoothingParam', SmoothingParam);
                    obj.z_predicted=f3(obj.x);
                else
                    fprintf('WARNING: transformation to spline slope for z failed because you do not have curve fitting toolbox \n');
                end
                x2=obj.x+1; z2=f3(x2); 
                x3=obj.x-1; z3=f3(x3);
%                 z2(end)=2*obj.z(end)-z3(end);
%                 z3(1)=2*obj.z(1)-z2(1);
                obj.R2=sum((obj.z_predicted-mean(obj.z)).^2)/sum((obj.z-mean(obj.z)).^2);
%                 x=[2*obj.x(1)-obj.x(2);obj.x;2*obj.x(end)-obj.x(end-1)];
%                 z=[2*obj.z_predicted(1)-obj.z_predicted(2);obj.z_predicted;2*obj.z_predicted(end)-obj.z_predicted(end-1)];
                obj.i=atan((z2-z3)/2);
            else
                Idstring=num2str(obj.Id);
                fprintf(['cannot transform to splined slope for hillslope ',Idstring ,' \n']);
            end
        end
        
        function obj=compute_beven_ratio(obj)
            obj.beven_ratio=trapz(obj.x,obj.w)/(mean(obj.w(obj.x<50))*sin(mean(obj.i(obj.x<50))));
        end
        
        function obj=compute_saturation_ind(obj)
            obj.saturated_area=trapz(obj.x,obj.w.*obj.saturation)/trapz(obj.x,obj.w);
            x_max_pos=find(obj.saturation<0.67,1);
            obj.saturated_length=obj.x(x_max_pos)/max(obj.x);
        end
        
        function [x,w,soil_depth,i,z]=get_spatial_properties(obj)
            x=obj.x;
            w=obj.w;
            soil_depth=obj.soil_depth;
            i=obj.i;
            z=obj.z;
        end
        
        function [f,k,phi]=get_hydraulic_properties(obj)
%             i=obj.i;
            f=obj.f;
            k=obj.k;
            phi=obj.phi;
        end
        
        function [obj,x_mean,z,s]=compute_elevation(obj)
            x_mean=(obj.x(2:end)+obj.x(1:end-1))/2;
            angle_mean=(obj.i(2:end)+obj.i(1:end-1))/2;
            dx=(obj.x(2:end)-obj.x(1:end-1));
            z=cumsum(dx.*tan(angle_mean));
            obj.z=interpn(x_mean,z,obj.x,'linear');
            % curvilinear absciss
            s=cumsum(dx./cos(angle_mean));
        end
        
        function folder_directories=save_hillslope(obj,file_output,name_watershed)
            X_pos=strfind(file_output,'X'); 
            Y_pos=strfind(file_output,'Y'); 
            X_coord=file_output(X_pos(end)+2:Y_pos(end)-2);
            Y_coord=file_output(Y_pos(end)+2:length(file_output));
            folder_directories=[];
% % %             % save one hillslope with constant slope (approximation with a linear fitting)
% % %             file_output1=[file_output,'_slopcst'];
% % %             folder_create(file_output1);
% % %             filename=strcat(file_output1,'/morphologic.input');
% % % %             obj=obj.transform_to_constant_slope;
% % %             M=nan(length(obj.x),5); M(:,1)=obj.x; M(:,2)=obj.w; M(:,3)=obj.i; M(:,4)=obj.z; M(:,5)=obj.z_predicted;
% % %             fid = fopen(filename, 'w');
% % %             if(nargin>=3)
% % %                 string_char=sprintf(['Real morphologic data taken in ',name_watershed ,' watershed. Slope is assumed constant \n']);
% % %             else
% % %                 string_char=sprintf('Real morphologic data taken. Slope is assumed constant \n');
% % %             end
% % %             fprintf(fid, string_char);
% % %             string_char=sprintf(['Hillslope coordinates: X= ',X_coord ,' Y= ',Y_coord,'\n']);
% % %             fprintf(fid, string_char);
% % %             string_char=sprintf('x\tw\ti\tz_true\tz_mod\n');
% % %             fprintf(fid, string_char);
% % %             fclose(fid);
% % %             dlmwrite(filename,M, '-append', 'precision', '%E','delimiter','\t');
% % %             save([file_output1,'/hs1D.mat'],'obj');
% % %             obj.plot_save_width_function(file_output1);
% % %             obj.plot_save_elevation_function(file_output1);
% % %             obj.plot_save_slope_angle_function(file_output1);
% % %             close all;
% % %             folder_directories{2}=file_output1; 
            
            % save one hillslope with non constant slope (approximation with spline curves)
            file_output2=[file_output,'_slopvar'];
            folder_create(file_output2);
            filename=strcat(file_output2,'/morphologic.input');
            SmoothingParam=1e-6;
%             SmoothingParam=0.5;
%             obj=obj.transform_to_spline_slope(SmoothingParam);
            M=nan(length(obj.x),5); M(:,1)=obj.x; M(:,2)=obj.w; M(:,3)=obj.i; M(:,4)=obj.z; M(:,5)=obj.z_predicted;
            fid = fopen(filename, 'w');
            if(nargin>=3)
                string_char=sprintf(['Real morphologic data taken in ',name_watershed ,' watershed. Slope is variable \n']);
            else
                string_char=sprintf('Real morphologic data taken. Slope is variable \n');
            end
            fprintf(fid, string_char);
            string_char=sprintf(['Hillslope coordinates: X= ',X_coord ,' Y= ',Y_coord,'\n']);
            fprintf(fid, string_char);
            string_char=sprintf('x\tw\ti\tz_true\tz_mod\n');
            fprintf(fid, string_char);
            fclose(fid);
            dlmwrite(filename,M, '-append', 'precision', '%E','delimiter','\t');
            save([file_output2,'/hs1D.mat'],'obj');
            obj.plot_save_width_function(file_output2);
            obj.plot_save_elevation_function(file_output2);
            obj.plot_save_slope_angle_function(file_output2);
            close all;
            folder_directories{1}=file_output2;
        end
       
        function plot_save_width_function(obj,filename)
            figure; hold on;
            plot(obj.x,obj.w,'linewidth',3,'Color',[0 0.4470 0.7410]);
            plot(obj.x,obj.w,'o','linestyle','none','Color',[0 0.4470 0.7410]);
            title('width function');
            xlabel('distance to the channel [m]');
            ylabel('width [m]');
            box on;
            savefig([filename,'/width_function.fig']); 
            print([filename,'/width_function.png'],'-dpng'); 
        end
        
        function plot_save_elevation_function(obj,filename)
            figure; hold on;
            plot(obj.x,obj.z,'linewidth',2,'Color',[0 0.4470 0.7410]);
            if(~isempty(obj.z_predicted))
            plot(obj.x,obj.z_predicted,'linewidth',2,'Color',[0.8500 0.3250 0.0980]);
            end
            plot(obj.x,obj.z,'o','linestyle','none','Color',[0 0.4470 0.7410]);
            title('elevation vs distance');
            xlabel('distance to the channel [m]');
            ylabel('elevation [m]');
            box on;
            if(~isempty(obj.z_predicted))
                legend('obtained from the dem','after the transformation');
            else
                legend('obtained from the dem');
            end
            savefig([filename,'/elevation_function.fig']);
            print([filename,'/elevation_function.png'],'-dpng');
            close all
        end
        
        function plot_save_slope_angle_function(obj,filename)
            figure; hold on;
            plot(obj.x,obj.i,'linewidth',2,'Color',[0 0.4470 0.7410]);
            plot(obj.x,obj.i,'o','linestyle','none','Color',[0 0.4470 0.7410]);
            title('slope angle vs distance');
            xlabel('distance to the channel [m]');
            ylabel('angle [%]');
            if(length(obj.i)>2)
                imin=min(0,min(obj.i)*0.98);
                imax=max(obj.i)*1.02;
                if(imin>=imax)
                    fprintf(['WARNING: potentially a hillslope not realistic \n hillslope number ID: ',num2str(obj.Id),'\n']);
                    ylim([imin imin+1]);
                else
                    ylim([imin imax]);
                end
            end
            savefig([filename,'/slope_angle_function.fig']);
            print([filename,'/slope_angle_function.png'],'-dpng');
            close all
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
%             folder_initial='C:\Users\Jean\Documents\Donn�es\SyntheticDataHillslope';
            folder_initial='C:\Users\Jean\Documents\Donn�es\SyntheticData';
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