classdef watershed
    % Set of hillslope
    properties(Access=public)
        DEM            % GRIDobj topotoolbox
        Outlet         % x y coordinates and ID (from topotoolbox) (1x3 array)
        Confluences    % x y coordinates and ID (from topotoolbox) (nx3 matrix)
        Channelheads   % x y coordinates and ID (from topotoolbox)  (nx3 matrix)
        S              % STREAMobj topotoolbox
        Channels       % Stream divided into channels from one singular points (channelheads, confluences or outlet) to another
        hillslopes     % set of hillslopes (see class hillslope)
        StreamOrderChannels
        FD 
        geomorphologic_prop
    end
    
    methods(Static)
        function [obj,hillslope_equiv]=test(DEM_filepath,outlet_coordinates,critic_drainage_area,wgs84arg)
            if(nargin<3)
                critic_drainage_area=40000;
            end
            if(nargin<2)
                outlet_coordinates=[364778.7,6834822.7];
            end
            if(nargin<1)
                DEM_filepath='MNT_PF_5m.tif';
            end
            obj=watershed;
            if(nargin<4)
                wgs84arg=0;
            end
            obj=obj.analyze_hillslopes(DEM_filepath,outlet_coordinates,critic_drainage_area,wgs84arg);
            hillslope_equiv=extract_one_equivalent_hillslope(obj);
% % %             outlet_coord_BV6=[265549.0,6783308.0];%[214629.641,2346729.811];% (265549 6783308)
% % %             critic_drainage_area=40000;
% % %             obj=obj.analyze_hillslopes('MNT_5m_Naizin.tif',outlet_coord_BV6,critic_drainage_area);
% %             outlet_coord=[181036.9995,6864336.5019];
% %             critic_drainage_area=400;
% %             obj=obj.analyze_hillslopes('Guillec50.tif',outlet_coord,critic_drainage_area);
% % %             outlet_coord=[175981.46,6819377.47];
% % %             critic_drainage_area=400;
% % %             obj=obj.analyze_hillslopes('Douffine50.tif',outlet_coord,critic_drainage_area);
%             folder_directory='C:\Users\Jean\Documents\Donn�es\Synthetic\temp';
%             name_watershed='Pleine Fougeres'; 
%             obj.save_hillslopes(folder_directory,name_watershed);
%             obj=watershed;
%             outlet_coord_BV6=[364778.7,6834822.7];
%             critic_drainage_area=40000;
%             obj=obj.analyze_hillslopes('MNT_PF_5m.tif',outlet_coord_BV6,critic_drainage_area);
            
%             obj=watershed;
%             outlet_coord1=[365738.54,6823382.22];
%             outlet_coord2=[365110.2962,6824211.47];
%             outlet_coord3=[366204.899,6820675.00];
%             outlet_coord4=[370252.401,6825459.003];
%             outlet_coord5=[387853.90,6812521.00];
%             critic_drainage_area=40000;
%             obj=obj.analyze_hillslopes('MNT_Couesnon3_5m.tif',outlet_coord5,critic_drainage_area);
% % %                       obj=watershed;
% % %             outlet_coord_BV6=[351525,6790117];
% % %             critic_drainage_area=40000;
% % %             obj=obj.analyze_hillslopes('MNT_Vilaine_5m.tif',outlet_coord_BV6,critic_drainage_area);
% % %             folder_directory='C:\Users\Jean\Documents\Donn�es\Synthetic\temp';
% % %             name_watershed='Upper Vilaine'; 
% % %             obj.save_hillslopes(folder_directory,name_watershed);
%             obj=watershed;
%             outlet_coord_Binic=[268863.41 ,6849415.51];
%             critic_drainage_area=40000;
%             obj=obj.analyze_hillslopes('Ic_BV2.tif',outlet_coord_Binic,critic_drainage_area);
%             folder_directory='C:\Users\Jean\Documents\Donn�es\Synthetic\MorphologicInputs';
%             name_watershed='Ic at Binic'; 
%             obj.save_hillslopes(folder_directory,name_watershed);
%             name_watershed2 = name_watershed(~isspace(name_watershed));
%             save(strcat(folder_directory,'\watershed',name_watershed2,'.mat'),'obj');
% % % % % % % % % %             obj=watershed;
% % % % % % % % % %             outlet_coord_Kerb=[169091.50,6784372.86388];
% % % % % % % % % % %             outlet_coord_Kerb=[169565.30,6784657.50];%[168965.1255,6784200.191937];
% % % % % % % % % %             critic_drainage_area=10000;
% % % % % % % % % %             obj=obj.analyze_hillslopes('MNT_Kerb_5m.tif',outlet_coord_Kerb,critic_drainage_area);
% %             folder_directory='C:\Users\Jean\Documents\Donn�es\Synthetic\MorphologicInputs';
%             folder_directory='C:\Users\Jean\Documents\Donn�es\Synthetic\temp';
%             name_watershed='Kerbernez extended'; 
%             obj.save_hillslopes(folder_directory,name_watershed);
%             name_watershed2 = name_watershed(~isspace(name_watershed));
%             save(strcat(folder_directory,'\watershed',name_watershed2,'.mat'),'obj');
            
%             obj=obj.load_DEM('Naizin_surrounding.tif');
%             outlet_coord_naizin=[266520.40,6779439.00];
%             obj=obj.get_watershed_DEM(outlet_coord_naizin);
% 
%             obj=obj.load_DEM('Blavet.tif');
%             outlet_coord_blavet=[255371.81,6774575.91];
%             obj=obj.get_watershed_DEM(outlet_coord_blavet);
%             obj=obj.load_DEM('Ic_BV2.tif');
%             outlet_coord_Binic=[268863.41 ,6849415.51];
%             obj=obj.get_watershed_DEM(outlet_coord_Binic);
%             obj=obj.load_DEM('MNT_Kerb_5m.tif');
%             outlet_coord_E30=[169091.7,6784372.8]; %E6=[118189.57,2347135.81];
%             obj=obj.get_watershed_DEM(outlet_coord_E30);
%             obj=obj.extract_stream_and_singular_points;
%             obj=obj.extract_channels_from_stream;
%             obj=obj.set_hillslopes_from_channels;
%             obj=obj.create_hillslopes1D;
%             folder_directory='C:\Users\Jean\Documents\Donn�es\Synthetic\MorphologicInputs';
%             name_watershed='Pleine Fougeres'; 
%             name_watershed='Nazin extended'; 
%             folder_directory='C:\Users\Jean\Documents\Donn�es\HillslopesData\Parameters';
%             folder_directory='C:\Users\Jean\Documents\Donn�es\SyntheticData\MorphologicInputs';
%             obj.save_hillslopes(folder_directory,name_watershed);
% % %             outlet_coord_Svartberget=[731404.4794,7133661.9441]; % 731406	7133634
%             critic_drainage_area=2000;%12800;%floor(469000/25);%4000;%1000;%floor(469000/25);%
% % %             obj=obj.analyze_hillslopes('C:\Users\Jean\Documents\ProjetTamara\SHP_Raster\Svartberget_DEM5m_2.tif',outlet_coord_Svartberget,critic_drainage_area);
% % %             outlet_coord_Icacos=[-65.7854572,18.2752528];
% % %             critic_drainage_area=2000;
% % % % %             obj=obj.analyze_hillslopes('C:\Users\Jean\Documents\ProjectLou\DEMPuertoRico\Icacos.tif',outlet_coord_Icacos,critic_drainage_area);
% %             outlet_coord_Icacos2=[205450.48,2022982.24];
% %             critic_drainage_area=50000/5;%500;
% % %             obj=obj.analyze_hillslopes('C:\Users\Jean\Documents\ProjectLou\GisDataLuquilloCZO_\Icacos2.tif',outlet_coord_Icacos2,critic_drainage_area);
% %             obj=obj.analyze_hillslopes('C:\Users\Jean\Documents\ProjectLou\DEM_1m\outputmean.tif',outlet_coord_Icacos2,critic_drainage_area);
            
%                         outlet_coord_Hawai=[239536.33,2125675.80];
%             critic_drainage_area=1e6/100;%500;
% %             obj=obj.analyze_hillslopes('C:\Users\Jean\Documents\ProjectLou\GisDataLuquilloCZO_\Icacos2.tif',outlet_coord_Icacos2,critic_drainage_area);
%             obj=obj.analyze_hillslopes('C:\Users\Jean\Documents\ProjectLouHawai\AlidaKestimates\Hawaii_DEM\hawaii_dem\DEMStream1.tif',outlet_coord_Hawai,critic_drainage_area);
% % % % % % % % % % % % % % % % % outlet_coord_Icacos=[2.0557e+05,2.0186e+06];%[205793,2017411];
% % % % % % % % % % % % % % % % % critic_drainage_area=1e5/25;%500;
% % % % % % % % % % % % % % % % % obj=obj.analyze_hillslopes('C:\Users\Jean\Documents\ProjectLou\DEM_1m\USGS_LiDAR_DEM\PuertoRico_Extended_UTM20N_5m.tif',outlet_coord_Icacos,critic_drainage_area);
         
%             outlet_coord_Quentin=[217874.42,6790718.674];
%             critic_drainage_area=400;
%             obj=obj.analyze_hillslopes('MNT2_temp4.tif',outlet_coord_Quentin,critic_drainage_area);
        end
        
%         function [x,y]=contour(DEM_filepath,outlet_coordinates,critic_drainage_area,Output_folder,wgs84arg)
        function Err=contour(DEM_filepath,outlet_coordinates,critic_drainage_area,Output_folder,Station_Name,accumulation_vectorize,wgs84arg)
             if(nargin<5)
                Station_Name='PeineFougeres.shp';
            end
            if(nargin<4)
                Output_folder='/home/jean.marcais/Bureau/';
            end

            if(nargin<2)
                outlet_coordinates=[364778.7,6834822.7];
            end
            if(nargin<1)
                DEM_filepath='MNT_PF_5m.tif';
            end
            obj=watershed;
            if(nargin<7)
                wgs84arg=0;
            end
            
            if(ischar(DEM_filepath) || isstr(DEM_filepath))
                obj=obj.load_DEM(DEM_filepath);
            else
                obj = DEM_filepath;
            end
            if(nargin<3)
                critic_drainage_area=1e6;
            end
            if(nargin>6 && wgs84arg>0)
                obj.DEM=reproject2utm(obj.DEM,30,'zone','30U');
            end
            if(outlet_coordinates==-1)
                [obj,xriv,yriv]=obj.get_watershed_DEM;
            elseif(nargin>5 && ~isnumeric(accumulation_vectorize))
                [obj,xriv,yriv]=obj.get_watershed_DEM(outlet_coordinates,floor(critic_drainage_area/obj.DEM.cellsize^2),accumulation_vectorize);
            else
                [obj,xriv,yriv]=obj.get_watershed_DEM(outlet_coordinates,floor(critic_drainage_area/obj.DEM.cellsize^2));
            end
            
            % first alternative
%             if(critic_drainage_area==-1)
%                 obj=obj.extract_stream_and_singular_points(xriv,yriv);
%             else
%                 obj=obj.extract_stream_and_singular_points(xriv,yriv,critic_drainage_area);
%             end
%             
%             hillslope_equiv=obj.extract_one_equivalent_hillslope;
%             
%             k=boundary(hillslope_equiv.x,hillslope_equiv.y,1);
%             x=hillslope_equiv.x(k);
%             y=hillslope_equiv.y(k);
            
            % second alternative
            [x,y] = getcoordinates(obj.DEM,'matrix');
            size_=size(obj.DEM.Z);

            x= reshape(x,size_(1)*size_(2),1);
            y= reshape(y,size_(1)*size_(2),1);
            z= reshape(obj.DEM.Z,size_(1)*size_(2),1);
            x=x(~isnan(z));
            y=y(~isnan(z));
            
%             if(length(x)<1.1e7)
                k1=boundary(x,y,1);
                x=x(k1);
                y=y(k1);
                polyout = polyshape(x,y);
%             elseif(length(x)>=1.1e7 && length(x)<2e7)
%                 x1=x(1:floor(5*length(x)/9));
%                 y1=y(1:floor(5*length(x)/9));
%                 x2=x(end-floor(5*length(x)/9):end);
%                 y2=y(end-floor(5*length(x)/9):end);
%                 
%                 k1=boundary(x1,y1,1);
%                 k2=boundary(x2,y2,1);
%                 x1=x1(k1);
%                 y1=y1(k1);
%                 x2=x2(k2);
%                 y2=y2(k2);
%                 polyout1 = polyshape(x1,y1);
%                 polyout2 = polyshape(x2,y2);
%                 polyout = union(polyout1,polyout2);
%             elseif(length(x)>=2e7 && length(x)<1e8)
%                 x1=x(1:floor(length(x)/3));
%                 y1=y(1:floor(length(x)/3));
%                 x2=x(floor(length(x)/3)+1:2*floor(length(x)/3));
%                 y2=y(floor(length(x)/3)+1:2*floor(length(x)/3));
%                 x3=x(2*floor(length(x)/3)+1:end);
%                 y3=y(2*floor(length(x)/3)+1:end);
%                 
%                 k1=boundary(x1,y1,1);
%                 k2=boundary(x2,y2,1);
%                 k3=boundary(x3,y3,1);
%                 
%                 x1=x1(k1);
%                 y1=y1(k1);
%                 x2=x2(k2);
%                 y2=y2(k2);
%                 x3=x3(k3);
%                 y3=y3(k3);
%                 polyout1 = polyshape(x1,y1);
%                 polyout2 = polyshape(x2,y2);
%                 polyout3 = polyshape(x3,y3);
%                 polyout = union(polyout1,polyout2);
%                 polyout = union(polyout,polyout3);
%             else
%             end

            
% % %             [X_1,~,index]=unique(x);
% % %             y_1=accumarray(index,y,[],@max)+obj.DEM.cellsize/2;
% % %             y_2=accumarray(index,y,[],@min)-obj.DEM.cellsize/2;
% % %             [Y_1,~,index_2]=unique(y);
% % %             x_1=accumarray(index_2,x,[],@max)+obj.DEM.cellsize/2;
% % %             x_2=accumarray(index_2,x,[],@min)-obj.DEM.cellsize/2;
% % %             X_11=[X_1;flip(X_1)];
% % %             Y_11=[y_1;flip(y_2)];
% % %             X_12=[x_1;flip(x_2)];
% % %             Y_12=[Y_1;flip(Y_1)];
% % %             
% % %             poly1 = polyshape(X_11,Y_11);
% % %             poly2 = polyshape(X_12,Y_12);
% % %             polyout = intersect(poly1,poly2);
            
            
            Data.Geometry = 'Polygon' ;
            Data.X = polyout.Vertices(:,1) ;  % latitude
            Data.Y = polyout.Vertices(:,2)  ;  % longitude
            Data.Name = Station_Name ;   % some random attribute/ name
            Data.AreaTopoTB = area(polyout); %sum(~isnan(obj.DEM.Z(:)))*obj.DEM.cellsize^2;
            Data.AreaBH = critic_drainage_area/0.7;
            Err = (critic_drainage_area/0.7 - Data.AreaTopoTB)/(critic_drainage_area/0.7);
            Data.DiffRelArea = Err;
            Data.X_Out_TopoTB = xriv;
            Data.Y_Out_TopoTB = yriv;
            Data.X_Out_BH = outlet_coordinates(1);
            Data.Y_Out_BH = outlet_coordinates(2);
            Data.Outlet_Dist_Err = sqrt(sum((outlet_coordinates(1)-xriv)^2+(outlet_coordinates(2)-yriv)^2));
            shapewrite(Data, strcat(Output_folder,'/',Station_Name,'.shp'));
        end
        
        function generate_hillslopes_from_dems
            dem_filenames={'MNT_PF_5m.tif','MNT_Kerb_5m.tif','Naizin_surrounding.tif','Blavet.tif','Ic_BV2.tif'};
            name_watersheds={'Pleine Fougeres','Kerbernez extended','Naizin extended','Blavet','Ic at Binic'};
            outlet_coords=[364778.7,6834822.7;169565.30,6784657.50;266520.40,6779439.00;255371.81,6774575.91;268863.41 ,6849415.51];
            critic_drainage_areas=[40000,10000,40000,20000,40000];
            folder_directory='C:\Users\Jean\Documents\Donn�es\Synthetic\MorphologicInputs';
            
            for j=1:length(dem_filenames)
                obj=watershed;
                obj=obj.analyze_hillslopes(dem_filenames{j},outlet_coords(j,:),critic_drainage_areas(j));
                obj.save_hillslopes(folder_directory,name_watersheds{j});
                name_watershed =name_watersheds{j};
                name_watershed = name_watershed(~isspace(name_watershed));
                save(strcat(folder_directory,'\watershed',name_watershed,'.mat'),'obj');
            end
        end
    end
    
    methods(Access=public)
        % Constructor
        function obj=watershed
        end
        
        % Load DEM
        function obj=load_DEM(obj,filename)
%             file_directory=which(filename);
            obj.DEM=GRIDobj(filename);
            % remove aberrant values into nan
%             obj.DEM.Z(obj.DEM.Z<-10)=nan;
            obj.DEM=fillsinks(obj.DEM);
        end
        
        function [obj,xriv,yriv]=get_watershed_DEM(obj,manual_outlet_coord,critical_drainage_area,accumulation_option)
            if(nargin<2) manual_outlet_coord=[];    end
            if(nargin<3) critical_drainage_area=floor(1e6/obj.DEM.cellsize^2);    end
            
            if(nargin<4 || ( isnumeric(accumulation_option) && accumulation_option==0))
                FD=FLOWobj(obj.DEM);
                obj.FD=FD;
                A=flowacc(FD);
            else 
                A=accumulation_option;
            end
            
            if(isempty(manual_outlet_coord))
                DB = drainagebasins(FD);
                DB = shufflelabel(DB);
                
                STATS = regionprops(DB.Z,'PixelIdxList','Area','Centroid');
                [~,IX] = max([STATS.Area]);
                obj.DEM.Z(DB.Z~=IX)=NaN;
            else
                W=A>critical_drainage_area;%>1111;%>40000;%
                S=STREAMobj(obj.FD,W);
                [xriv,yriv] = snap2stream(S,manual_outlet_coord(1),manual_outlet_coord(2));
                DB = drainagebasins(obj.FD,xriv,yriv);
%                 Outlets_coord=streampoi(FD,W,'outlets','xy');
%                 Outlets_coord=Outlets_coord(find(min(((Outlets_coord(:,1)-manual_outlet_coord(1)).^2+(Outlets_coord(:,2)-manual_outlet_coord(2)).^2).^0.5)));
%                 DB = drainagebasins(FD,Outlets_coord(1),Outlets_coord(2));
                obj.DEM.Z(DB.Z==0)=NaN;
            end
%             obj.FD=FLOWobj(obj.DEM);
        end
        

        function obj=extract_stream_and_singular_points(obj,xriv,yriv,critic_drainage_area)
            if(nargin<4) 
                % if not specified chose 40000 cells so for a 5m DEM it is 1km2
                critic_drainage_area=1111;%40000; 
            end
            FD=obj.FD;
            A=flowacc(FD);
            W = A>critic_drainage_area;
            obj.S = STREAMobj(FD,W);
            
            % Extract singular points
%             obj.Outlet=streampoi(FD,W,'outlets','xy');
            obj.Outlet=[xriv,yriv];
            Channelheads=streampoi(FD,W,'channelheads','xy');
            Confluences=streampoi(FD,W,'confluences','xy');
            % Get Topotoolbox ID
%             obj.Outlet=[obj.Outlet,streampoi(FD,W,'outlets','ix')];
            obj.Outlet=[obj.Outlet,coord2ind(obj.DEM,xriv,yriv)];
            Channelheads=[Channelheads,streampoi(FD,W,'channelheads','ix')];
            Confluences=[Confluences,streampoi(FD,W,'confluences','ix')];
            
            % retrieve singular points outside catchments
            [dem,x,y]=GRIDobj2mat(obj.DEM);
            [x,y]=meshgrid(x,y);
            x=x(~isnan(dem)); y=y(~isnan(dem));
            ix=coord2ind(obj.DEM,x,y);
            
            
            
            obj.Channelheads=[];
            obj.Confluences=[];
            for i=1:size(Channelheads,1 )
                Test=sum((ix==Channelheads(i,3)));
                if(Test~=0)
                    obj.Channelheads=[obj.Channelheads;Channelheads(i,:)];
                end
            end
            for i=1:size(Confluences,1)
                Test=sum((ix==Confluences(i,3)));
                if(Test~=0)
                    obj.Confluences=[obj.Confluences;Confluences(i,:)];
                end
            end
        end
        
        function obj=extract_channels_from_stream(obj)
            % Store all singular points in a vector
            StreamOrder=streamorder(obj.S);
            if(~isempty(obj.Confluences))
                Singular_points=[obj.Channelheads(:,3);obj.Confluences(:,3)];
                Singular_points_coord=[obj.Channelheads(:,:);obj.Confluences(:,:)];
                % Find first point upstream of the confluences
                Confluences_coord_prec=[];
                Confluences_prec=[];
                for k=1:length(obj.Confluences(:,1))
                    Id=find(obj.S.IXgrid==obj.Confluences(k,3));
                    Precedent_nodes=obj.S.ix(obj.S.ixc==Id);
                    Confluences_coord_prec=[Confluences_coord_prec; obj.S.x(Precedent_nodes),obj.S.y(Precedent_nodes)];
                    Confluences_prec=[Confluences_prec; obj.S.IXgrid(Precedent_nodes)];
                end
            else
                Singular_points=[obj.Channelheads(:,3)];
                Singular_points_coord=[obj.Channelheads(:,:)];
                Confluences_coord_prec=[];
                Confluences_prec=[];
            end
            
            % %             Singular_points=[obj.Confluences(:,3)];
            % %             Singular_points_coord=[obj.Confluences(:,:)];
            

            
            Singular_points_prec=[Confluences_prec;obj.Outlet(:,3)];
            %             Singular_points_coord_prec=[Confluences_coord_prec;obj.Outlet(:,1:2)];
            
            obj.Channels=cell(1,length(Singular_points)+length(obj.Channelheads(:,1)));
            obj.Channels=cell(1,length(Singular_points));
            for i=1:length(Singular_points)
                obj.StreamOrderChannels(i)=StreamOrder(obj.S.IXgrid==Singular_points(i));
                obj.Channels{i}=Singular_points_coord(i,:);
                Test=[];
                A=Singular_points(i);
                while(isempty(Test))
                    Id=find(obj.S.IXgrid==A);
                    Next_node=obj.S.ixc(obj.S.ix==Id);
                    A=obj.S.IXgrid(Next_node);
                    if(length(A)>1 || isempty(A))
                        AAA=1;
                    end
                    Test=find(Singular_points_prec==A);
                    obj.Channels{i}=[obj.Channels{i};obj.S.x(Next_node),obj.S.y(Next_node),obj.S.IXgrid(Next_node)];
                end
            end
            
            %% lines to comment if you want to retrieve only side channels hillslopes
            % % %             for i=(length(Singular_points)+1):(length(Singular_points)+length(obj.Channelheads(:,1)))
            % % %                 obj.Channels{i}=obj.Channelheads(i-length(Singular_points),1:3);
            % % %             end
        end
        
        function obj=set_hillslopes_from_channels(obj)
            FD=obj.FD;
            Distances_to_stream=flowdistance(FD,obj.S);
            compt=1;
            DEM_resolution=obj.DEM.cellsize;
            for i=1:length(obj.Channels)
                obj.hillslopes{i}=hillslope(compt);
                obj.hillslopes{i}=obj.hillslopes{i}.extract_hillslopes_from_channel(obj.DEM,obj.Channels{i},FD,Distances_to_stream);
                compt=compt+1;
                if(length(obj.Channels{i}(:,1))>1)
                    [obj.hillslopes{i}(1),obj.hillslopes{i}(2)]=obj.hillslopes{i}.extract_unique_hillslope_from_channel(obj.Channels{i},DEM_resolution);
                    compt=compt+1;
                end
            end
        end
        
        function hillslope_equiv=extract_one_equivalent_hillslope(obj)
            Distances_to_stream=flowdistance(obj.FD,obj.S);
            DZ = vertdistance2stream(obj.FD,obj.S,obj.DEM);
            hillslope_equiv=hillslope(1);
            Channel_equiv=[];
            for i=1:length(obj.Channels)
                Channel_equiv=[Channel_equiv;obj.Channels{i}];
            end
            FD=obj.FD;
            hillslope_equiv=hillslope_equiv.extract_equivalent_hillslope_from_channel(obj.DEM,Channel_equiv,FD,Distances_to_stream,obj.Outlet,DZ);
            hillslope_equiv=hillslope_equiv.compute_hillslope1D(obj.DEM.cellsize);
        end
        
        function obj=create_hillslopes1D(obj)
            DEM_resolution=obj.DEM.cellsize;
            for i=1:length(obj.hillslopes)
                if(length(obj.hillslopes{i})==1)
                    obj.hillslopes{i}=obj.hillslopes{i}.compute_hillslope1D(DEM_resolution);
                else
                    obj.hillslopes{i}(1)=obj.hillslopes{i}(1).compute_hillslope1D(DEM_resolution);
                    obj.hillslopes{i}(2)=obj.hillslopes{i}(2).compute_hillslope1D(DEM_resolution);
                end
            end
        end
        
        function obj=create_hillslopes1D_customed(obj,Id)
            DEM_resolution=obj.DEM.cellsize;
            if(nargin<2) Id=1; end
            if(Id>0)
                for i=1:length(obj.hillslopes)
                    if(obj.StreamOrderChannels(i)==1)
                        obj.hillslopes{i}=obj.hillslopes{i}.compute_hillslope1D(DEM_resolution);
                    else
                        obj.hillslopes{i}=obj.hillslopes{i}.extract_oriented_hillslope_from_channel(obj.Channels{i},DEM_resolution);
                        obj.hillslopes{i}=obj.hillslopes{i}.compute_hillslope1D(DEM_resolution);
                    end
                end
            else
                for i=1:length(obj.hillslopes)
                    if(floor(obj.StreamOrderChannels(i))==obj.StreamOrderChannels(i))
                        obj.hillslopes{i}=obj.hillslopes{i}.compute_hillslope1D(DEM_resolution);
                    else
                        obj.hillslopes{i}=obj.hillslopes{i}.extract_oriented_hillslope_from_channel(obj.Channels{i},DEM_resolution);
                        obj.hillslopes{i}=obj.hillslopes{i}.compute_hillslope1D(DEM_resolution);
                    end
                end
            end
        end
        
        function obj=set_hillslopes_from_channels_2(obj)
            FD=obj.FD;
            Distances_to_stream=flowdistance(FD,obj.S);
            compt=1;
            DEM_resolution=obj.DEM.cellsize;
            DZ = vertdistance2stream(FD,obj.S,obj.DEM);
            for i=1:length(obj.Channels)
                if(obj.StreamOrderChannels(i)>1)
                    obj.hillslopes{i}=hillslope(compt);
                    obj.hillslopes{i}=obj.hillslopes{i}.extract_hillslopes_from_channel(obj.DEM,obj.Channels{i},FD,Distances_to_stream);
                    compt=compt+1;
                    [obj.hillslopes{i}(1),obj.hillslopes{i}(2)]=obj.hillslopes{i}.extract_unique_hillslope_from_channel(obj.Channels{i},DEM_resolution);
                    compt=compt+1;
                else
                    obj.hillslopes{i}=hillslope(compt);
                    obj.hillslopes{i}=obj.hillslopes{i}.extract_equivalent_hillslope_from_channel(obj.DEM,obj.Channels{i},FD,Distances_to_stream,obj.Channels{i}(end,:),DZ);
                    compt=compt+1;
                end
            end
        end
        
        function obj=set_hillslopes_from_channels_3(obj)
            FD=obj.FD;
            Distances_to_stream=flowdistance(FD,obj.S);
            compt=1;
            DEM_resolution=obj.DEM.cellsize;
            DZ = vertdistance2stream(FD,obj.S,obj.DEM);
            for i=1:length(obj.Channels)
                if(obj.StreamOrderChannels(i)>1)
                    obj.hillslopes{i}=hillslope(compt);
                    obj.hillslopes{i}=obj.hillslopes{i}.extract_hillslopes_from_channel(obj.DEM,obj.Channels{i},FD,Distances_to_stream);
                    compt=compt+1;
%                     [obj.hillslopes{i}(1),obj.hillslopes{i}(2)]=obj.hillslopes{i}.extract_unique_hillslope_from_channel(obj.Channels{i},DEM_resolution);
%                     compt=compt+1;
                else
                    obj.hillslopes{i}=hillslope(compt);
                    obj.hillslopes{i}=obj.hillslopes{i}.extract_equivalent_hillslope_from_channel(obj.DEM,obj.Channels{i},FD,Distances_to_stream,DZ);
                    compt=compt+1;
                end
            end
        end
        
        function save_hillslopes(obj,folder_directory,name_watershed)
            for i=1:length(obj.hillslopes)
                if(length(obj.hillslopes{i})==1)
                    obj.hillslopes{i}.save_hillslope(folder_directory,name_watershed);
                else
                    obj.hillslopes{i}(1).save_hillslope(folder_directory,name_watershed);
                    obj.hillslopes{i}(2).save_hillslope(folder_directory,name_watershed);
                end
            end
        end
        
        function obj=analyze_hillslopes(obj,dem_filename,outlet_coord,critic_drainage_area,wgs84arg)
            obj=obj.load_DEM(dem_filename);
            if(nargin>4 && wgs84arg>0)
                obj.DEM=reproject2utm(obj.DEM,30,'zone','30U');
            end
            if(outlet_coord==-1)
                [obj,xriv,yriv]=obj.get_watershed_DEM;
            else
                [obj,xriv,yriv]=obj.get_watershed_DEM(outlet_coord,critic_drainage_area);
            end
            if(critic_drainage_area==-1)
                obj=obj.extract_stream_and_singular_points(xriv,yriv);
            else
                obj=obj.extract_stream_and_singular_points(xriv,yriv,critic_drainage_area);
            end
            obj=obj.extract_channels_from_stream;
            obj=obj.set_hillslopes_from_channels_2;
            obj=obj.create_hillslopes1D;
        end
        
        function obj=customed_analysis1_hillslopes(obj,dem_filename,outlet_coord,critic_drainage_area)
            obj=obj.load_DEM(dem_filename);
            if(outlet_coord==-1)
                [obj,xriv,yriv]=obj.get_watershed_DEM;
            else
                [obj,xriv,yriv]=obj.get_watershed_DEM(outlet_coord);
            end
            if(critic_drainage_area==-1)
                obj=obj.extract_stream_and_singular_points(xriv,yriv);
            else
                obj=obj.extract_stream_and_singular_points(xriv,yriv,critic_drainage_area);
            end
        end
        
        function obj=customed_analysis2_hillslopes(obj,Channelheads_ID,Confluences_ID)
            obj.Channelheads=obj.Channelheads(Channelheads_ID,:);
            if(isvector(Confluences_ID)==1)
                obj.Confluences=obj.Confluences(Confluences_ID,:);
                obj=obj.extract_channels_from_stream;
                obj=obj.set_hillslopes_from_channels_3;
                obj=obj.create_hillslopes1D_customed;
            else
                [Conf_temp_x,Conf_temp_y,Conf_temp_id]=snap2stream(obj.S,Confluences_ID(:,1),Confluences_ID(:,2));
                obj.Confluences=[Conf_temp_x,Conf_temp_y,Conf_temp_id];
                obj=obj.extract_channels_from_stream_and_hillslopes;
                obj=obj.create_hillslopes1D_customed(-1);
            end
        end
        
        function obj=extract_channels_from_stream_and_hillslopes(obj)
              % Store all singular points in a vector
            Distances_to_stream=flowdistance(obj.FD,obj.S);
            StreamOrder=streamorder(obj.S);
            [X,Y]=getcoordinates(obj.DEM);
            for i=1:(length(obj.Confluences)+1)
                obj.StreamOrderChannels(i)=0;
                if(i==length(obj.Confluences)+1) 
                     critical_point=obj.Outlet(1,3);
                else
                    critical_point=obj.Confluences(i,3);
                end
                hillslope_1=drainagebasins(obj.FD,critical_point);
                [r,c]=find(hillslope_1.Z==1);
                x_temp=X(c)'; y_temp=Y(r);
                ix_temp=coord2ind(obj.DEM,x_temp,y_temp);
                for k=1:length(obj.Confluences)
                    if(critical_point~=obj.Confluences(k,3) && sum((obj.Confluences(k,3)==ix_temp))~=0)
                        hillslope_temp=drainagebasins(obj.FD,obj.Confluences(k,3));
                        hillslope_1.Z=hillslope_1.Z-hillslope_temp.Z;
                        obj.StreamOrderChannels(i)=0.5;
                    end
                end
                [r,c]=find(hillslope_1.Z==1);
                x_temp=X(c)'; y_temp=Y(r);
                ix_temp=coord2ind(obj.DEM,x_temp,y_temp);
                Channel_id=intersect(obj.S.IXgrid,ix_temp,'stable');
                [x_s,y_s]=STREAMobj2XY(obj.S);
                ind_s=coord2ind(obj.DEM,x_s,y_s);
                Channel_id=intersect(ind_s,Channel_id,'stable');
                [Channel_x,Channel_y]=ind2coord(obj.DEM,Channel_id);
                
%                 [Channel_x,Channel_y]=ind2coord(obj.DEM,Channel_id);
                obj.Channels{i}=[Channel_x,Channel_y,Channel_id];
                obj.StreamOrderChannels(i)=obj.StreamOrderChannels(i)+StreamOrder(obj.S.IXgrid==critical_point);
                % create hillslope
                obj.hillslopes{i}=hillslope(i);
                obj.hillslopes{i}=obj.hillslopes{i}.set_hillslope_prop(x_temp,y_temp,obj.DEM.Z(hillslope_1.Z==1),Distances_to_stream.Z(hillslope_1.Z==1));
            end
        end
        
        function obj=set_geomorphologic_properties(obj,x_max,w_mean,w_riv,delta_z_top,Area,Area_width,Area_width_abs,Area_elevation,Area_elevation_abs,slope,convergence,StreamOrder)
            obj.geomorphologic_prop.x_max=x_max;
            obj.geomorphologic_prop.w_mean=w_mean;
            obj.geomorphologic_prop.w_riv=w_riv;
            obj.geomorphologic_prop.delta_z_top=delta_z_top;
            obj.geomorphologic_prop.Area=Area;
            obj.geomorphologic_prop.Area_width=Area_width;
            obj.geomorphologic_prop.Area_width_abs=Area_width_abs;
            obj.geomorphologic_prop.Area_elevation=Area_elevation;
            obj.geomorphologic_prop.Area_elevation_abs=Area_elevation_abs;
            obj.geomorphologic_prop.slope=slope;
            obj.geomorphologic_prop.convergence=convergence;
            obj.geomorphologic_prop.StreamOrder=StreamOrder;
        end
    end
    
    methods(Static)
        function [Catchment_error_list,Catchment_error_values]=extract_BH_stations
            warning('off','all');
            Stations_Shp = shaperead('~/Donnees/BanqueHydro/Shapefiles/BanqueHydro_StationsSignature.shp');
%             RegHy =
%             {'A','B','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','U','V','W','X','Y'};
%             Pb on H & L & M & Q & V
%             RegHy =
%             {'I','J','K','L','M','N','O','P','Q','R','S','U','V','W','X','Y'};
%             Pb on L & M & Q & V
%             RegHy = { 'N','O','P','Q','R','S','U','V','W','X','Y'};
%             RegHy = {'R','S','U','V','W','X','Y'};
%             RegHy = {'A','B','D','E','F','G','I','J','K'};%'W','X',
            RegHy = {'O'};
            save_folder='/home/jean.marcais/Donnees/Carthage/tmp/ContourWatershedFrance/';
            Catchment_error_list={};
            Catchment_error_values=[];
            count = 1;
            for k=1:length(RegHy)
                DEM_filepath=strcat('~/Donnees/Carthage/tmp/ZoneHydro_UniqueCdRegionHy_2/CdRegionHy_',RegHy{k},'.tif');
                obj=watershed;
                obj=obj.load_DEM(DEM_filepath);
%                 obj.DEM.Z(obj.DEM.Z<=0)=nan;
                
                obj.FD=FLOWobj(obj.DEM);
                A=flowacc(obj.FD);
                save_folder_RegHy=strcat(save_folder,'/RegHy_',RegHy{k},'_4/');
                test_ = mkdir(save_folder_RegHy);
                for i=1070:length(Stations_Shp)
                    fprintf(strcat('Hydrologic Region ',RegHy{k},' ---- Station number',num2str(i),' over ',num2str(length(Stations_Shp)),'\n'));
                    if(Stations_Shp(i).CdStationH(1)==RegHy{k})
                        fprintf(strcat('Station BELONGS to the Hydrologic Region', RegHy{k},'. \n'));
                        Err = watershed.contour(obj,[Stations_Shp(i).X,Stations_Shp(i).Y],0.7*Stations_Shp(i).Area_BH*1e6,save_folder_RegHy,Stations_Shp(i).CdStationH,A); 
                        if(abs(Err)>0.1)
                            fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
                            fprintf(strcat('Potential error on watershed ', Stations_Shp(i).CdStationH,' \n'));
                            fprintf(strcat('Error is ', num2str(Err*100),' %%. \n'));
                            fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
                            Catchment_error_values(count) = Err;
                            Catchment_error_list{count} = Stations_Shp(i).CdStationH;
                            count = count+1;
                        else
                            fprintf(strcat('Error is ', num2str(Err*100),' %%. SEEMS OK \n'));
                        end
                    else
                        fprintf(strcat('Station NOT in the Hydrologic Region', RegHy{k},'. \n'));
                    end
                end
                fprintf('==================  ====================\n');
                fprintf('==================  ====================\n');
                fprintf(strcat('================== Hydrologic Region ',RegHy{k},' is completed !! ====================\n'));
                fprintf('==================  ====================\n');
                fprintf('==================  ====================\n');

            end
        end
    end
end