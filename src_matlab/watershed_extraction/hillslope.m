classdef hillslope
    properties(Access=public)
        Id          % Hillslope identifier number
        coordinates_id % Hillslope coordinates identifier
        x           % coord of the hillslope
        y           % coord of the hillslope
        z           % elevation
        distance    % distance to the first node in the river following the drainage network
        Hilltype    % hillslope type 1: Head Hillslope 2 Channel Hillslope
        hsB         % hillslope1D (from Troch et al. 2003) ie class hillslope1D
        link_hillslope_hs1D % link hillslope hs1D
    end
    
    methods(Access=public)
        % constructor
        function obj=hillslope(Id)
            obj.Id=Id;
        end
        
        function obj=set_hillslope_prop(obj,x,y,z,distance)
            obj.x=x;
            obj.y=y;
            obj.z=z;
            obj.distance=distance;
        end
        
        function obj=extract_hillslopes_from_channel(obj,DEM,Channel,FD,Distances_to_stream)
            [X,Y]=getcoordinates(DEM);
            hillslope=drainagebasins(FD,Channel(end,1),Channel(end,2));
            if(length(Channel(:,1))>1)
                drainage_temp=drainagebasins(FD,Channel(1,1),Channel(1,2));
                hillslope.Z=hillslope.Z-drainage_temp.Z;
            end
            [r,c]=find(hillslope.Z==1);
            obj.x=X(c)'; obj.y=Y(r);
            obj.z=DEM.Z(hillslope.Z==1);
            obj.distance=Distances_to_stream.Z(hillslope.Z==1);
            obj.Hilltype=1;
%             [~,pos_coord_mean]=min((obj.distance-mean(obj.distance)).^2);
%             obj.coordinates_id=[obj.x(pos_coord_mean),obj.y(pos_coord_mean)];
            obj.coordinates_id=[Channel(end,1),Channel(end,2)];
            % remove channel's points
%             obj=obj.remove_points_in_stream(Channel);
        end
        
        % extract equivalent hillslope from a channel 
        function obj=extract_equivalent_hillslope_from_channel(obj,DEM,Channel,FD,Distances_to_stream,Outlet,vertdistance2stream)
            [X,Y]=getcoordinates(DEM);
            if(nargin<6)
                hillslope=drainagebasins(FD,Channel(end,3));
            else
                hillslope=drainagebasins(FD,Outlet(3));
            end
            [r,c]=find(hillslope.Z==1);
            obj.x=X(c)'; obj.y=Y(r);
%             obj.z=DEM.Z(hillslope.Z==1);
            if(nargin>6)
                obj.z=vertdistance2stream.Z(hillslope.Z==1);
            else
                obj.z=DEM.Z(hillslope.Z==1);
            end
            obj.distance=Distances_to_stream.Z(hillslope.Z==1);
            obj.Hilltype=3;
        end
        
        function obj=remove_points_in_stream(obj,Channel)
            Distance_stream = pdist2([obj.x,obj.y],Channel,'euclidean');
            Distance_stream = min(Distance_stream,[],2);
            In_stream= Distance_stream==0;
            obj.x(In_stream)=[]; obj.y(In_stream)=[];
            obj.z(In_stream)=[]; obj.distance(In_stream)=[];
%             for i=1:length(Channel(:,1))
%                 In_channel=(Channel(i,1)==obj.x & Channel(i,2)==obj.y);
%                 obj.x(In_channel)=[]; obj.y(In_channel)=[];
%                 obj.z(In_channel)=[]; obj.distance(In_channel)=[];
%             end
        end
        
        function [obj,obj2]=extract_unique_hillslope_from_channel(obj,Channel,DEM_resolution)
            if(nargin<3)
                DEM_resolution=2;
            end
            [Bank1,Bank2]=obj.generate_banks(Channel,DEM_resolution);
            coord_hillslopes=[obj.x,obj.y];
            D_Bank1 = pdist2(coord_hillslopes,Bank1,'euclidean');
            D_Bank2 = pdist2(coord_hillslopes,Bank2,'euclidean');
            D_Bank1=min(D_Bank1,[],2);
            D_Bank2=min(D_Bank2,[],2);
            % create second hillslope
            obj2=hillslope(obj.Id+1);
            obj2.x=obj.x(D_Bank1>=D_Bank2,:);
            obj2.y=obj.y(D_Bank1>=D_Bank2,:);
            obj2.z=obj.z(D_Bank1>=D_Bank2,:);
            obj2.distance=obj.distance(D_Bank1>=D_Bank2,:);
            obj2.Hilltype=2;
%             [~,pos_coord_mean]=min((obj2.distance-mean(obj2.distance)).^2);
%             obj2.coordinates_id=[obj2.x(pos_coord_mean),obj2.y(pos_coord_mean)];
            obj2.coordinates_id=[Channel(end,1),Channel(end,2)];
            [~,pos_coord_mean]=min(pdist2(Bank2,obj2.coordinates_id,'euclidean'));
            obj2.coordinates_id=[Bank2(pos_coord_mean,1),Bank2(pos_coord_mean,2)];
            % append first hillslope
            obj.x=obj.x(D_Bank1<=D_Bank2,:);
            obj.y=obj.y(D_Bank1<=D_Bank2,:);
            obj.z=obj.z(D_Bank1<=D_Bank2,:);
            obj.distance=obj.distance(D_Bank1<=D_Bank2,:);
            obj.Hilltype=2;
%             [~,pos_coord_mean]=min((obj.distance-mean(obj.distance)).^2);
%             obj.coordinates_id=[obj.x(pos_coord_mean),obj.y(pos_coord_mean)];
            obj.coordinates_id=[Channel(end,1),Channel(end,2)];
            [~,pos_coord_mean]=min(pdist2(Bank1,obj.coordinates_id,'euclidean'));
            obj.coordinates_id=[Bank1(pos_coord_mean,1),Bank1(pos_coord_mean,2)];
        end
        
        function obj=compute_hillslope1D(obj,DEM_resolution)
            if(nargin<2)
                DEM_resolution=5;
            end
            obj.hsB=hillslope1D;
            obj.hsB=obj.hsB.set_properties(obj,obj.Id);
            [obj.hsB,obj.link_hillslope_hs1D]=obj.hsB.set_parameters(obj.distance,obj.z,DEM_resolution);
%             obj.hsB=obj.hsB.transform_to_constant_slope;
            obj.hsB=obj.hsB.transform_to_spline_slope(0.0001);
%             obj.hsB=obj.hsB.compute_beven_ratio;
        end
        
        function save_hillslope(obj,folder_directory,name_watershed)
%             c=clock; time_string_folder=strcat(num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)));
%             coordinate_string_folder=['_X_',num2str(round(obj.coordinates_id(1))),'_Y_',num2str(round(obj.coordinates_id(2)))];            
%             folder_name=[time_string_folder,coordinate_string_folder];
%             coordinate_string_folder=['X_',num2str(round(obj.coordinates_id(1))),'_Y_',num2str(round(obj.coordinates_id(2)))];
            [x_bound,y_bound]=obj.compute_hillslope_contour;
            x_barycentre=mean(x_bound);
            y_barycentre=mean(y_bound);
            % in case barycentre is not in the contour
            %compute Euclidean distances:
            A=[x_barycentre,y_barycentre];
            B=[obj.x,obj.y];
            distances = sqrt(sum(bsxfun(@minus, B, A).^2,2));
            %find the smallest distance and use that as an index into B:
            closest = B(find(distances==min(distances)),:);
            if(isempty(closest))
                x_coord_hillslope=x_barycentre;
                y_coord_hillslope=y_barycentre;
                fprintf(['WARNING: Potential problem with hillslope',num2str(obj.Id),'X= ',num2str(x_coord_hillslope),'Y=',num2str(y_coord_hillslope),'\n']);
            else
                x_coord_hillslope=closest(1,1);
                y_coord_hillslope=closest(1,2);
            end

            coordinate_string_folder=['X_',num2str(round(x_coord_hillslope)),'_Y_',num2str(round(y_coord_hillslope))];
            folder_name=coordinate_string_folder;
            file_output=strcat(folder_directory,'\',folder_name);

            folder_directories=obj.hsB.save_hillslope(file_output,name_watershed);
            for j=1:length(folder_directories)
                obj.save_hillslope_object(folder_directories{j});
                obj.save_hillslope_contour(x_bound,y_bound,folder_directories{j});
                obj.plot_hillslope_contour(x_bound,y_bound,x_coord_hillslope,y_coord_hillslope,folder_directories{j});
            end
        end
        
        function save_hillslope_object(obj,folder)
            save(strcat(folder,'\hillslope.mat'),'obj');
        end
        function [x_bound,y_bound]=compute_hillslope_contour(obj)
            k=boundary(obj.x,obj.y);
            x_bound=obj.x(k);
            y_bound=obj.y(k);
        end
            
        function save_hillslope_contour(obj,x_bound,y_bound,folder)
            filename=strcat(folder,'\contour.struct');
            M=nan(length(x_bound),2); M(:,1)=x_bound; M(:,2)=y_bound;
            fid = fopen(filename, 'w');
            string_char=sprintf('Contour data of the hillslope \n');
            fprintf(fid, string_char);
            string_char=sprintf('First column x coordinates and second column y coordinates (of the hillslopes boundaries) \n');
            fprintf(fid, string_char);
            string_char=sprintf('x\ty\n');
            fprintf(fid, string_char);
            fclose(fid);
            dlmwrite(filename,M, '-append', 'precision', '%E','delimiter','\t');
        end
        
        function plot_hillslope_contour(obj,x_bound,y_bound,x_coord_hillslope,y_coord_hillslope,folder)
            x_riv=obj.x(obj.distance<5);
            y_riv=obj.y(obj.distance<5);
            figure; hold on;
            plot(obj.x,obj.y,'Marker','o','Color',[0.9290 0.6940 0.1250],'LineStyle','none');
            plot(x_bound,y_bound,'linewidth',3,'Color',[0.8500 0.3250 0.0980]);
            plot(x_riv,y_riv,'linewidth',3,'Color',[0 0.447 0.7410]);%,'LineStyle','none','Marker','.');
            plot(x_coord_hillslope,y_coord_hillslope,'Color',[0.4940 0.1840 0.5560],'Marker','+','MarkerSize',12);
            title('Hillslope 2D structure');
            xlabel('X coordinates [m]');
            ylabel('Y coordinates [m]');
            legend('hillslope','contour','river');
            if(nargin>=4)
                savefig([folder,'\contour.fig']);
                print([folder,'\contour.png'],'-dpng');
            end
            close all;
        end
        
        function [x,y,k]=extract_hillslope_contour(obj)
            k=boundary(obj.x,obj.y,1);
            x=obj.x(k);
            y=obj.y(k);
        end
        
        function [obj]=extract_oriented_hillslope_from_channel(obj,Channel,DEM_resolution)
            if(nargin<3)
                DEM_resolution=2;
            end
            [Bank1,Bank2]=obj.generate_banks(Channel,DEM_resolution);
            coord_hillslopes=[obj.x,obj.y];
            size_hillslopes_points=length(obj.x);
            size_bank=size(Bank1,1);
            if(size_bank*size_hillslopes_points<1.6e9)
                D_Bank1 = pdist2(coord_hillslopes,Bank1,'euclidean');
                D_Bank2 = pdist2(coord_hillslopes,Bank2,'euclidean');
                D_Bank1=min(D_Bank1,[],2);
                D_Bank2=min(D_Bank2,[],2);
            else
                max_size_bank=floor(1.6e9/size_hillslopes_points);
                N_step=ceil(size_bank/max_size_bank);
                N_iter=floor((N_step).^(1/3));
                D_Bank1=[]; D_Bank2=[];
                for i=1:N_iter
                    D_Bank1 = [D_Bank1,pdist2(coord_hillslopes,Bank1((1+(i-1)*N_iter):N_step:end,:),'euclidean')];
                    D_Bank2 = [D_Bank2,pdist2(coord_hillslopes,Bank2((1+(i-1)*N_iter):N_step:end,:),'euclidean')];
                    D_Bank1=min(D_Bank1,[],2);
                    D_Bank2=min(D_Bank2,[],2);
                end
%                 D_Bank1 = pdist2(coord_hillslopes,Bank1,'euclidean');
%                 D_Bank2 = pdist2(coord_hillslopes,Bank2,'euclidean');
            end
            % 
            obj.distance(D_Bank1>=D_Bank2)=-obj.distance(D_Bank1>=D_Bank2,:);
        end
        
        function plot_spatialized_data(obj,z)
            if(nargin<2)
                z=obj.z;
            end         
            X1=unique(obj.x);
            Y1=unique(obj.y);
            Z1=nan(length(Y1),length(X1));

            Distance_X=pdist2(obj.x,X1,'euclidean');
            Distance_Y=pdist2(obj.y,Y1,'euclidean');

            [~,Pos_X]=min(Distance_X,[],2);
            [~,Pos_Y]=min(Distance_Y,[],2);

            Ind_matrix_link=sub2ind(size(Z1),Pos_Y,Pos_X);
            Z1(Ind_matrix_link)=z;
%             figure; hold on
            h=imagesc(X1,Y1,Z1);
            set(gca,'Ydir','normal');
            set(h,'alphadata',~isnan(Z1));
            axis equal
            colorbar
        end

    end
    
    methods(Static)
        function [Bank1,Bank2]=generate_banks(Channel,resolution)
            Vector=Channel(2:end,1:2)-Channel(1:end-1,1:2);
            Midpoints=(Channel(2:end,1:2)+Channel(1:end-1,1:2))/2;
            Distance_vec=(Vector(:,1).^2+Vector(:,2).^2).^0.5;
            Distance_crit=1.5*sqrt(2*resolution^2);
            Midpoints=Midpoints(Distance_vec<Distance_crit,:);
            Vector=Vector(Distance_vec<Distance_crit,:);
            Channel_bis=Channel(1:end-1,1:2);
            Channel_bis=Channel_bis(Distance_vec<Distance_crit,:);
            Bank1=Midpoints+resolution/2*[Vector(:,2),-Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
            Bank1_temp=Channel_bis+resolution/2*[Vector(:,2),-Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
            Bank1=[Bank1;Bank1_temp];
            Bank2=Midpoints+resolution/2*[-Vector(:,2),Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
            Bank2_temp=Channel_bis+resolution/2*[-Vector(:,2),Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
            Bank2=[Bank2;Bank2_temp];
            %
%             Bank1=Midpoints+resolution/2*[Vector(:,2),-Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
%             Bank1_temp=Channel_bis(1:end-1,1:2)+resolution/2*[Vector(:,2),-Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
%             Bank1=[Bank1;Bank1_temp];
%             Bank2=Midpoints+resolution/2*[-Vector(:,2),Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
%             Bank2_temp=Channel_bis(1:end-1,1:2)+resolution/2*[-Vector(:,2),Vector(:,1)]./[(Vector(:,2).^2+Vector(:,1).^2).^0.5,(Vector(:,2).^2+Vector(:,1).^2).^0.5];
%             Bank2=[Bank2;Bank2_temp];
        end
    end
    

    
end