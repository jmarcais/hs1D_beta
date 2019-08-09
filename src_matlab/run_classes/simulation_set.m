classdef simulation_set
    properties(Access=public)
        mother_folder_directory
        geologic_inputs_directory
        morphologic_inputs_directory
        hydrologic_inputs_directory
        combination_inputs
    end
    
    methods(Access=public)
        function obj=simulation_set(mother_folder_directory)
            obj.mother_folder_directory=mother_folder_directory;
        end
        
        function obj=instantiate_all_inputs_directory(obj)
            obj.geologic_inputs_directory=obj.get_inputs(obj.mother_folder_directory,'GeologicInputs');
            obj.morphologic_inputs_directory=obj.get_inputs(obj.mother_folder_directory,'MorphologicInputs');
            obj.hydrologic_inputs_directory=obj.get_inputs(obj.mother_folder_directory,'HydrologicInputs');
            obj=obj.compute_all_possible_combination_between_inputs;
        end
        
        % find directories in simulations results that encountered errors
        function output_list_directory=get_error_simulation_directory(obj)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            output_list_directory=[];
            compt=1;
            for i=1:length(directory_list)
                A=dir(fullfile(directory,directory_list{i},'*.err'));
                if(~isempty(A))
                    output_list_directory{compt}=fullfile(directory,directory_list{i},A.name);
                    compt=compt+1;
                    if(length(A)>1)
                        fprintf(['WARNING: Potential bad use of .err file chosen: more than one error file found in',fullfile(directory_list{i})]);
                    end
                end
            end
        end
        
        % find directories in simulations results that are incomplete
        function output_list_directory=detect_incomplete_simulations(obj)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            for i=1:length(directory_list)
                A=dir(fullfile(directory,directory_list{i},'mass_balance.png'));
                if(isempty(A))
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                    if(length(A)>1)
                        fprintf(['WARNING: Potential bad use of file chosen: more than one massbalance file found in',fullfile(directory_list{i})]);
                    end
                end
            end
        end
        
        % find directories in simulation_set that share the same infiltration chronicle
        function output_list_directory=find_simulation_infiltration_chronicle(obj,infiltration_identifier_string)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            length_string=length(infiltration_identifier_string);
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                Bool=strcmp(directory_list{i}(end-length_string+1:end),infiltration_identifier_string);
                if(Bool==1)
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                end
            end
        end
        
        % find directories in simulation_set that share the same hillslope shape
        function output_list_directory=find_simulation_hillslope_coordinates(obj,Xcoordinates,Ycoordinates)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                [X_coordinates_file,Y_coordinates_file]=obj.extract_coordinates(directory_list{i});
                
                Bool1=strcmp(X_coordinates_file,num2str(Xcoordinates));
                Bool2=strcmp(Y_coordinates_file,num2str(Ycoordinates));
                Bool3=Bool1*Bool2;
                if(Bool3==1)
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                end
            end
        end
        
        function output_list_directory=find_simulation_hillslope_sloptype(obj,sloptype_str)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                [sloptype]=obj.extract_slope_type(directory_list{i});
                
                Bool=strcmp(sloptype,sloptype_str);
                if(Bool==1)
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                end
            end
        end
        
        function [X_coordinates_file,Y_coordinates_file]=extract_coordinates(obj,directory_path_name)
            Id1=strfind(directory_path_name,'X');
            Id2=strfind(directory_path_name,'Y');
            Id3=strfind(directory_path_name,'slop');
            X_coordinates_file=directory_path_name(Id1+2:Id2-2);
            Y_coordinates_file=directory_path_name(Id2+2:Id3-2);
        end
        
        function [sloptype]=extract_slope_type(obj,directory_path_name)
            Id3=strfind(directory_path_name,'slop');
            sloptype=directory_path_name(Id3+4:Id3+6);
        end
                
        % find directories in simulation_set that share the same hillslope shape
        function output_list_directory=find_simulation_geometric_properties(obj,f,k,d)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                directory_total=fullfile(directory,directory_list{i},'geologic.input');
                fid=fopen(directory_total,'r');
                if(fid>0)
                    C = textscan(fid, '%s','delimiter', '\t');
                    f_file=str2num(C{1}{end-2});
                    k_file=str2num(C{1}{end-1});
                    d_file=str2num(C{1}{end});
                    Bool=(f_file-f)^2+(k_file-k)^2+(d_file-d)^2;
                    
                    if(Bool==0)
                        output_list_directory{compt}=fullfile(directory,directory_list{i});
                        compt=compt+1;
                    end
                    fclose(fid);
                end
            end
        end
      
        function output_list_directory=get_inputs(obj,directory,type_of_inputs)
            if(nargin>=3)
                directory=[directory,'\',type_of_inputs];
            end
            list=list_folder_of(directory);
            output_list_directory=[];
            compt=1;
            for i=1:length(list)
                A=dir(fullfile(directory,list{i},'*.input'));
                if(~isempty(A))
                    output_list_directory{compt}=fullfile(directory,list{i},A.name);
                    compt=compt+1;
                    if(length(A)>1)
                        fprintf(['WARNING: Potential bad use of .input file chosen: more than one file in',fullfile(directory,list{i})]);
                    end
                end
            end
        end
        
        function simulation_list_directory=get_output_simulation(obj)
            directory=obj.mother_folder_directory;
            directory=[directory,'\Simulations'];
            simulation_list_directory=list_folder_of(directory);
        end
        
        function obj=compute_all_possible_combination_between_inputs(obj)
            obj.combination_inputs=allcomb(obj.geologic_inputs_directory,obj.hydrologic_inputs_directory,obj.morphologic_inputs_directory);
        end
        
        function obj=run_simulation_set(obj)
            % create a Simulation folder in the mother_folder_directory
            simulation_folder_root=fullfile(obj.mother_folder_directory,'Simulations');
            folder_create(simulation_folder_root);
%             obj.initialize_summary_file(simulation_folder_root);
            % launch in a for loop all the simulations
            max_iter=size(obj.combination_inputs);
            max_iter=max_iter(1);
            
%             perm_pos=randperm(max_iter);
            
            % generate some pdfs with wich the runs will update f k d
            % ##JM maybe to change - to indroduce in new versions
%             pd_f=makedist('Normal','mu',0.3,'sigma',0.1);
%             t_f=truncate(pd_f,0.05,0.5);
%             pd_k=makedist('Lognormal','mu',1,'sigma',1);
%             t_k=truncate(pd_k,0.05,15);
%             pd_d=makedist('LogNormal','mu',1,'sigma',1);
%             t_d=truncate(pd_d,0.2,11);
%             r_f=random(t_f,max_iter,1);
%             r_k=random(t_k,max_iter,1);
%             r_d=random(t_d,max_iter,1);
            
            for i=1:max_iter               
                % locations of different inputs file
                geo_loc=obj.combination_inputs(i,1);
                hydro_loc=obj.combination_inputs(i,2);
                morpho_loc=obj.combination_inputs(i,3);
                 % ##JM maybe to change - to indroduce in new versions
%                 geo_loc=obj.combination_inputs(perm_pos(i),1);
%                 hydro_loc=obj.combination_inputs(perm_pos(i),2);
%                 morpho_loc=obj.combination_inputs(perm_pos(i),3);
%                 AA=obj.combination_inputs(i,3);
%                 M2=obj.read_input_file(AA{1});
%                 w_test=M2(:,2);
%                 if(w_test(1)<25)
                % 1/ create a specific folder to store parameters and results files of the run
                %    copy paste the inputs file and some figures in the simulation results folder
%                 for j=1:61
                [destination_geol_file,destination_hydro_file,destination_morpho_file,folder_output]=obj.create_specific_results_folder(simulation_folder_root,morpho_loc,hydro_loc,geo_loc);

                
                % 2/ read input files and "recreate" objects (hs1D and source objects)
                % ##JM maybe to change - to indroduce in new versions
                M=obj.read_input_file(destination_geol_file);
                %#JM change after
                f=M(2,1); k=M(2,2); d=M(2,3);
                % ##JM maybe to change - to indroduce in new versions
%                 f=r_f(i); k=r_k(i); d=r_d(i);
%                 dlmwrite(destination_geol_file,[f,k,d],'-append','delimiter','\t','precision','%E');
                M=obj.read_input_file(destination_morpho_file);
                %#JM change after test
                x=M(:,1); w=M(:,2); slope_angle=M(:,3); z_top=M(:,4); %d=M(:,6);
% % %                 z_bottom_init=min(z_top)-2;
% % %                 z_bottom=z_bottom_init+0.0*x;
% % %                 slope_angle=0.0*ones(size(x));
% % %                 d=z_top-z_bottom;
                d=d*ones(size(x));
                f=f*ones(size(x));
                k=k*ones(size(x));
% % %                 k2=k/100;
%                 f2=f/100;
% % %                 typical_depth=25;
% % % %                 f(d>typical_depth)=(f(d>typical_depth)*typical_depth+f2(d>typical_depth).*(d(d>typical_depth)-typical_depth))./d(d>typical_depth);
% % %                 k(d>typical_depth)=(k(d>typical_depth)*typical_depth+k2(d>typical_depth).*(d(d>typical_depth)-typical_depth))./d(d>typical_depth);
%                 if(j==31)
%                     slope_angle=0.15*ones(size(x));
%                 else
%                     slope_angle=((((j-1)/2):(62-2*j)/(2*(length(x)-1)):((61-j)/2))/100)';
%                 end
%                 slope_angle(x>=100)=0.05;
%                 slope_angle(x>150)=0.1;
%                 slope_angle=0.01*j*ones(size(slope_angle));
%                     w(w<25)=30;
                    % ##JM maybe to change - to indroduce in new versions
%                                    hs1D=hillslope1D;
%                     hs1D=hs1D.set_properties(perm_pos(i),f,k);
                 hs1D=hillslope1D;
                 hs1D=hs1D.set_properties(1,f,k);
%                 hs1D=hs1D.set_spatial_parameters(x,w,slope_angle,d*ones(size(x)));
%                 d=(0.5:1/(length(x)-1):1.5)';
%                 d=ones(size(x));
%                 d=0.5*j*ones(size(x));
%                 d(x<50)=0.5;
%                 d=0.5*ones(size(x));
                hs1D=hs1D.set_spatial_parameters(x,w,slope_angle,d);
                hs1D=hs1D.compute_elevation;
                % #JM maybe to change after
%                 hs1D.plot_save_width_function(folder_output);
                hs1D.plot_save_elevation_function(folder_output);
                hs1D.plot_save_slope_angle_function(folder_output);
                [M,input_type]=obj.read_input_file(destination_hydro_file);
                t=M(:,1)*3600*24;
                recharge_chronicle=(M(:,2))';
                if(length(M(1,:))>2)
                    ETP_chronicle=(M(:,3))';
                    ratio_P_R=1;
                    source_terms=source('data_based');
                    [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s',ETP_chronicle);
%                     [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle-ETP_chronicle,'m/s');
                    ratio_P_R=1;
                else
                    ratio_P_R=1;%/0.38;
                    source_terms=source('data_based');
                    [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                    ratio_P_R=1;%0.38;
                end
                
                % 3/ create a runs object and run the simulation
                run_obj=runs;
                % set the solver options default or assigned in parameters via an odeset structure
                % specify Refine options for real infiltrations chronicle because for accuracy you need 
                % to force matlab ode15s to compute where you know sthg is happening
                if(strcmp(input_type,'Real1') || strcmp(input_type,'Real2') || strcmp(input_type,'Real3'))% || strcmp(input_type,'RealNbis') || strcmp(input_type,'RealN1bis') || strcmp(input_type,'RealN') || strcmp(input_type,'RealN1'))
%                     odeset_struct=odeset('RelTol',2.5e-14,'AbsTol',1e-14,'Refine',-1);
                    odeset_struct=odeset('RelTol',2.5e-14,'Refine',-1);
                % otherwise let matlab define its own time range for efficiency reasons (speed)
                else
%                     odeset_struct=odeset('RelTol',2.5e-14,'AbsTol',1e-14);
                    odeset_struct=odeset('RelTol',2.5e-14);
                end
                solver_options=run_obj.set_solver_options(odeset_struct);

                % if the simulation is not already a steady state in itself...
                if(~strcmp(input_type,'Synthetic1'))
                    % ... preprocess simulations to reach a steady state with averaged forcing conditions
                    % ##JM better change after
%                     recharge_averaged=0;
                    recharge_averaged=1e3*24*3600*source_terms.recharge_mean; % recharge averaged in mm/d
                    state_values_initial=obj.prerun_steady_state(hs1D,recharge_averaged,ratio_P_R);
                    presteadystate_percentage_loaded=-2; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
                    % run the simulation starting from the steady state condition 
                    run_obj=run_obj.run_simulation(hs1D,source_terms,presteadystate_percentage_loaded,solver_options,ratio_P_R,state_values_initial);
                else
                    % run the simulation starting from empty hillslope
                    percentage_loaded=0;
                    run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
                end
                tmax=t(end);
                error=run_obj.boussinesq_simulation.save_error_file(tmax,folder_output);
                run_obj.boussinesq_simulation.save_info_file(folder_output);
                
                % save boussinesq_simulation & hs1D objects
                space_limited_code=2; % if no code, assumed that enough space to store directly .mat objects, if code=1, enough space to convert it in .txt, if code=2, ode15s details are not saved
%                 run_obj.save_sim_run(folder_output,space_limited_code);
                run_obj.save_sim_run(folder_output);
                
                % save simulation_results objects and simulation_results output files
                t_sim_results=run_obj.simulation_results.t;
                
                [S_max,Watershed_area,w,Flux_in_total,Watershed_area_spatialized]=run_obj.save_key_parameters(t_sim_results,folder_output);
                run_obj.simulation_results.save_results(S_max,Watershed_area,folder_output);
                run_obj.simulation_results.plot_results(S_max,Watershed_area,w,Flux_in_total,folder_output);
%                 run_obj.simulation_results.plot_mass_changes(Watershed_area,Flux_in_total',Watershed_area_spatialized,folder_output);
% % % % % % % %                 dS=run_obj.boussinesq_simulation.compute_raw_storage_variation;
% % % % % % % %                 Flux_in=run_obj.compute_flux_in(t_sim_results);
% % % % % % % %                 size_dS=size(dS);
% % % % % % % %                 dS_rough=nan(size_dS(1),length(t_sim_results));
% % % % % % % %                 for k=1:size_dS(1)
% % % % % % % %                     dS_rough(k,:)=interpn(run_obj.boussinesq_simulation.sol_simulated.x,dS(k,:),t_sim_results);
% % % % % % % %                 end
% % % % % % % %                 [Mass_balance_tot_mmd,Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd]= ...
% % % % % % % %                     run_obj.simulation_results.mass_changes_study(Watershed_area_spatialized,Flux_in,dS_rough);
% % % % % % % %                 run_obj.simulation_results.plot_mass_balance(Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd,Mass_balance_tot_mmd,folder_output);
                 % ##JM maybe to change - to indroduce in new versions
%                 obj.append_summary_file(simulation_folder_root,perm_pos(i),folder_output,error);
                obj.append_summary_file(simulation_folder_root,i,folder_output,error);
                clear run_obj;
%                 end
            end
        end
        
        function state_values_initial=prerun_steady_state(obj,hs1D,recharge_averaged,ratio_P_R,bound_river)
            if(nargin<4)
                ratio_P_R=1;
            end
            percentage_loaded=0;
            source_steady=source('steady');
            tsteady=0:1:36500;
            tmin=tsteady(1); tmax=tsteady(end); Nt=length(tsteady); time_unity_type='day';
            t=time_properties(tmin,tmax,Nt,time_unity_type);
            source_steady=source_steady.set_recharge_chronicle(t,recharge_averaged);
            prerun_steady=runs;
            % add an event equal to 1 to detect steady_state
            detect_steady_state=1;
%             odeset_struct=odeset('RelTol',2.5e-14,'AbsTol',1e-17,'Events',detect_steady_state);
            odeset_struct=odeset('RelTol',2.5e-14,  'Events',detect_steady_state);
            solver_options=prerun_steady.set_solver_options(odeset_struct);
            Sinitial=nan;
            if(nargin<5)
                prerun_steady=prerun_steady.run_simulation(hs1D,source_steady,percentage_loaded,solver_options,ratio_P_R,Sinitial);
            else
                prerun_steady=prerun_steady.run_simulation(hs1D,source_steady,percentage_loaded,solver_options,ratio_P_R,Sinitial,bound_river);
            end
%             Sinitial=prerun_steady.get_final_storage;
            state_values_initial=prerun_steady.get_final_state_values;
        end
        
        function [M,input_type]=read_input_file(obj,filename)
            fid = fopen(filename, 'r');
            if(fid>0)
                input_type= textscan(fid,'%s', 1,'headerlines',1);
                input_type=input_type{1};
                M=dlmread(filename,'\t',3,0);
                fclose(fid);
            else
                fprintf('');
            end
        end
        
        function copypaste_png_figures(obj,file_location,folder_sep_pos,folder_output)
            png_directories=dir(fullfile(file_location{1}(1:folder_sep_pos{1}(end)-1),'*.png'));
            for k=1:length(png_directories)
                filename_png=fullfile(file_location{1}(1:folder_sep_pos{1}(end)-1),png_directories(k).name);
                copypaste_file_customed(filename_png,folder_output);
            end
        end
        
        function [destination_geol_file,destination_hydro_file,destination_morpho_file,folder_output]=create_specific_results_folder(obj,simulation_folder_root,morpho_loc,hydro_loc,geo_loc)
            % 1/ create a specific folder to store parameters and results files of the run
            c=clock; time_string_folder=strcat(num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)),'_',num2str(round(c(6))));
            folder_sep_pos=strfind(morpho_loc,'\');
            coordinate_string=morpho_loc{1}(folder_sep_pos{1}(end-1)+1:folder_sep_pos{1}(end)-1);
            folder_sep_pos_hydro=strfind(hydro_loc,'\');
            hydro_type_string=hydro_loc{1}(folder_sep_pos_hydro{1}(end-1)+1:folder_sep_pos_hydro{1}(end)-1);
            folder_sep_pos_geo=strfind(geo_loc,'\');
            geo_type_string=geo_loc{1}(folder_sep_pos_geo{1}(end-1)+1:folder_sep_pos_geo{1}(end)-1);
            % only the hydro and the morpho identifiers in folders name
            name_folder=[time_string_folder,'_',coordinate_string,'_',hydro_type_string];
            folder_output=fullfile(simulation_folder_root,name_folder);
            folder_create(folder_output);
            
            % 2/ write an parameters file with the different input properties
            obj.write_input_file_parameters(folder_output,coordinate_string,hydro_type_string,geo_type_string);
            
            % 3/ copy paste the inputs file and some figures and the simulation results folder
            % get png images in simulations folder (previously in morphologic & hydrologic inputs directory)
            obj.copypaste_png_figures(morpho_loc,folder_sep_pos,folder_output);
            obj.copypaste_png_figures(hydro_loc,folder_sep_pos_hydro,folder_output);
            % copy paste the parameter folders in the simulation file folder file
            destination_geol_file=copypaste_file_customed(geo_loc{1},folder_output);
            destination_hydro_file=copypaste_file_customed(hydro_loc{1},folder_output);
            destination_morpho_file=copypaste_file_customed(morpho_loc{1},folder_output);
            % copypaste hillslope objects
            folder_sep_pos=strfind(morpho_loc,'\');
            dir_hillslope=[morpho_loc{1}(1:folder_sep_pos{1}(end)),'hillslope.mat'];
            copypaste_file_customed(dir_hillslope,folder_output);
        end
        
        function write_input_file_parameters(obj,foldername,coordinate_string,hydro_type_string,geo_type_string)
            folder_sep_pos_coord=strfind(coordinate_string,'_');
            X_coordinates=coordinate_string(folder_sep_pos_coord(1)+1:folder_sep_pos_coord(2)-1);
            Y_coordinates=coordinate_string(folder_sep_pos_coord(3)+1:folder_sep_pos_coord(4)-1);
            folder_sep_pos=strfind(geo_type_string,'_');
            f_values=geo_type_string(folder_sep_pos(1)+1:folder_sep_pos(2)-1);
            k_values=geo_type_string(folder_sep_pos(3)+1:folder_sep_pos(4)-1);
            d_values=geo_type_string(folder_sep_pos(5)+1:end);

            filename_input=strcat(foldername,'\input.param');
            fid = fopen(filename_input, 'w');
            
            fprintf(fid,'Hydrologic Forcing Type \n');
            fprintf(fid,[hydro_type_string,'\n']);
            fprintf(fid,'X\tY\n');
            fprintf(fid,[X_coordinates,'\t',Y_coordinates,'\n']);
            fprintf(fid,'f\tk\td\n');
            fprintf(fid,[f_values,'\t',k_values,'\t',d_values,'\n']);
            fclose(fid);
        end
        
        function initialize_summary_file(obj,simulation_folder_root)
            filename=strcat(simulation_folder_root,'\summary_file.txt');
            fid = fopen(filename, 'w');
            fprintf(fid, 'SimId\tFolderPath\tErr\n');
            fprintf(fid, '\n');
            fclose(fid);
        end
        
        function append_summary_file(obj,simulation_folder_root,i,simulation_folder_path,error)
            filename=strcat(simulation_folder_root,'\summary_file.txt');
            fid = fopen(filename, 'a');
            simulation_folder_path=strrep(simulation_folder_path, '\', '/');
            fprintf(fid, [num2str(i),'\t',simulation_folder_path,'\t',num2str(error),'\n']);
            fclose(fid);
        end
        
        function runs_=run_simulation_from_struct(obj,x,f,k,width,slope,d,runs_below)
             k=k*ones(size(x));
             f=f*ones(size(x));
                
                if(sum(d<0)==0)
                    hs1D=hillslope1D;
                    hs1D=hs1D.set_properties(1,f,k);
                    hs1D=hs1D.set_spatial_parameters(x,width,slope,d);%hs1D=hs1D.set_spatial_parameters(x./cos(slope),width,slope,d.*cos(slope));%
                    if(nargin<8)
                        hydro_loc=obj.hydrologic_inputs_directory{1};
                        [M,input_type]=obj.read_input_file(hydro_loc);
                        t=M(:,1);
%                         recharge_chronicle=(M(:,2:end))';
                        recharge_chronicle=(M(1:1000,2))';
                        % change timestep every 6 hours
%                         t=t(17:end-24);
%                         recharge_chronicle=recharge_chronicle(17:end-24);
%                         t=(reshape(t,6,740))';
%                         recharge_chronicle=(reshape(recharge_chronicle,6,740))';
%                         recharge_chronicle=mean(recharge_chronicle,2);
%                         t=mean(t,2);
                        
% %                         t=[t-t(end)+2*t(1)-t(2);t];
% %                         recharge_chronicle=[recharge_chronicle,recharge_chronicle];
                        
                        ratio_P_R=1;%0.83;%1.1588;%0.875;%.33;%/0.38;
                        source_terms=source('data_based');
                        [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                        time_2006_bis=time_properties(t(1),t(1000),1000,'sec');%time_2006_bis=time_properties(t(1),t(end),length(t),'sec');
                        time_2006_1=time_properties(t(1),t(1000),(1000-1)*4+1,'sec');%time_2006_1=time_properties(t(1),t(end),(length(t)-1)*4+1,'sec');
                        source_terms.time=time_2006_1;
                        source_terms.recharge_chronicle=interp1((time_2006_bis.get_properties)',recharge_chronicle*ratio_P_R,(time_2006_1.get_properties));
                        source_terms.recharge_mean=mean(source_terms.recharge_chronicle,2);
                        
                        odeset_struct=odeset('RelTol',1e-5,'AbsTol',1e-6,'MaxStep',3600*24);%30*3600*24);
                        ratio_P_R=1;
                        recharge_averaged=1e3*24*3600*source_terms.recharge_mean; % recharge averaged in mm/d
                        state_values_initial=obj.prerun_steady_state(hs1D,recharge_averaged,ratio_P_R,'empty');
                        presteadystate_percentage_loaded=-2; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
                        % run transient simulation
                        runs_=runs;
                        solver_options=runs_.set_solver_options(odeset_struct);
                        runs_=runs_.run_simulation(hs1D,source_terms,presteadystate_percentage_loaded,solver_options,ratio_P_R,state_values_initial,'empty');
                    else
                        [~,w_1]=get_resampled_variables(runs_below.boussinesq_simulation.discretization);
                        [~,~,DPSA_spat_bed]=compute_DPSA_RF(runs_below.simulation_results,runs_below.boussinesq_simulation);
                        qs1=DPSA_spat_bed;
                        
                        %option 1
%                         x_Q=x./cos(slope);
%                         dx_Q=x_Q(2:end)-x_Q(1:end-1);
%                         qs1=bsxfun(@rdivide,qs1,w_1.*dx_Q);
                        % option 2
                        dx_Q=runs_below.simulation_results.x_Q(2:end)-runs_below.simulation_results.x_Q(1:end-1);
                        qs1=bsxfun(@rdivide,qs1,w_1.*dx_Q);
                        
                        time_2006_2=time_properties(runs_below.simulation_results.t(1),runs_below.simulation_results.t(end),length(runs_below.simulation_results.t),'sec');
                        source_terms=source('data_based');
                        source_terms.time=time_2006_2;
                        source_terms.recharge_chronicle=qs1;
                        source_terms.recharge_mean=mean(source_terms.recharge_chronicle,2);
                        dt_mean=mean(diff(runs_below.simulation_results.t));
                        odeset_struct=odeset('RelTol',1e-5,'MaxStep',dt_mean);%3600*24);
                        ratio_P_R=1;
                        presteadystate_percentage_loaded=-2; %presteadystate_percentage_loaded=0; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
                        recharge_averaged=1e3*24*3600*source_terms.recharge_mean; % recharge averaged in mm/d
                        state_values_initial=obj.prerun_steady_state(hs1D,recharge_averaged,ratio_P_R,'empty');
                        % run transient simulation
                        runs_=runs;
                        solver_options=runs_.set_solver_options(odeset_struct);
                        runs_=runs_.run_simulation(hs1D,source_terms,presteadystate_percentage_loaded,solver_options,ratio_P_R,state_values_initial,'empty');
                    end
                    
                    
                end
        end
    end
    
    methods(Static)
        function run_simulations(folder_root)
            sim_set=simulation_set(folder_root);
            sim_set=sim_set.instantiate_all_inputs_directory;
            sim_set=sim_set.run_simulation_set;
        end
        
        function rerun_simulation(SimulationPath)
            sim_set=simulation_set(SimulationPath);
            output_list_directory=dir(fullfile(SimulationPath,'*.input'));
            sim_set.geologic_inputs_directory{1}=fullfile(SimulationPath,output_list_directory(1).name);
            sim_set.morphologic_inputs_directory{1}=fullfile(SimulationPath,output_list_directory(3).name);
            sim_set.hydrologic_inputs_directory{1}=fullfile(SimulationPath,output_list_directory(2).name);
            sim_set=sim_set.compute_all_possible_combination_between_inputs;
            sim_set.run_simulation_set;
        end
        
        function run_simulations2(hs1D,folder_root,n)
            % generate ad hoc parametrization defined in geologic_input_set static methods
            [geol_input_set,val]=geologic_input_set.generate_customed_hillslope_parametrization;
            % generate sistematic parametrization defined in generate_systematic_hillslope_parametrization methods
            %             [geol_input_set,val]=geologic_input_set.generate_systematic_hillslope_parametrization(n);
            folder_root=strcat(folder_root,'\Simulations\');
            for i=1:length(val)
                % create folder output
                c=clock; time_string_folder=strcat(num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)),'_',num2str(c(6)));
                folder_output=strcat(folder_root,time_string_folder);
                folder_create(folder_output);
                simulation_run=run_set;
                [simulation_run,discretization,percentage_loaded]=simulation_run.choose_structure_parameters(hs1D);
                if(size(val,2)==3)
                    simulation_run=simulation_run.choose_hydraulic_parameters(val(i,1),val(i,2),val(i,3));
                elseif(size(val,2)==4)
                    simulation_run=simulation_run.choose_hydraulic_parameters(val(i,1),val(i,2),val(i,3),val(i,4));
                end
                [t,source_term]=simulation_run.set_source_options;
                boundary_cond=simulation_run.set_boundary_conditions(discretization);
                solver_options=simulation_run.set_solver_options;
                simulation_run=simulation_run.initialize_boussinesq_simulation(discretization,source_term,boundary_cond,percentage_loaded,t.get_tmin);
                simulation_run.save_parametrization(folder_output);
                % Run the DAE
                [t,x_center,x_edge,S,Q,QS,obj.boussinesq_simulation]=simulation_run.boussinesq_simulation.implicit_scheme_solver(t,solver_options);
                simulation_run.simulation_results=simulation_results(t,x_center,x_edge,S,Q,QS);
                simulation_run.save_results(folder_output);
            end
        end
        
        function watershed_simulation=reconstruct_watershed_simulation(folder_dir,type_of_recharge)
            % get different simulations 
            directory=[folder_dir,'\Simulations'];
            directory_list=list_folder_of(directory);
            compt=1;
            for i=1:length(directory_list)
                A=dir(fullfile(directory,directory_list{i},'*.mat'));
                if(~isempty(A) && ~isempty(strfind(fullfile(directory,directory_list{i}),type_of_recharge)))
                    bouss_sim_dir{compt}=fullfile(directory,directory_list{i},'boussinesq_simulation.mat');
                    sim_dir{compt}=fullfile(directory,directory_list{i},'simulation_results.mat');
                    hillslope_dir{compt}=fullfile(directory,directory_list{i},'hillslope.mat');
                    hs1D_dir{compt}=fullfile(directory,directory_list{i},'hs1D.mat');
                    transport_dir{compt}=fullfile(directory,directory_list{i},'transport.mat');
                    compt=compt+1;
                end
            end
            
            Q={}; prof_WT_spat={}; number_of_particle={}; transit_distance={}; weights_={}; transit_time_particles={}; weighted_pdf={};
            for j=1:length(bouss_sim_dir)
                load(bouss_sim_dir{j}); load(sim_dir{j}); hill=load(hillslope_dir{j}); load(hs1D_dir{j});
                %                     runs1=runs; runs1.boussinesq_simulation=bouss_sim; runs1.simulation_results=obj;
                [x_S,width,depth,angle,x_Q,f,k]=get_resampled_variables(bouss_sim.discretization);
                dx=x_Q(2:end)-x_Q(1:end-1);
                z_top=cumtrapz(x_Q,angle);
                z_top2=interpn(x_Q,z_top,x_S);
                z_bottom=z_top-depth(1);
                z_bottom2=interpn(x_Q,z_bottom,x_S);
                Flux_in=bouss_sim.source_terms.recharge_chronicle;
                % compute flow time series
                river_outflow=obj.compute_seepage_total; DPSA=obj.compute_DPES(dx.*width,f.*width.*depth,Flux_in);
                Q{hill.obj.Id}=[river_outflow;DPSA];
                
                % compute depth of the Water Table
                h=bsxfun(@rdivide,obj.S,(f.*width));
                prof_WT=bsxfun(@minus,-h,-z_top2+z_bottom2); %z_top2-(z_bottom2+h);
                prof_WT_spat{hill.obj.Id}=hill.obj.link_hillslope_hs1D*hs1D.link_hs1D*prof_WT;
                clear bouss_sim obj hs1D
                % retrieve particles
                trans=load(transport_dir{j});
                L=length(trans.transport.t);
                for k=1:(L-1)
                    [~,transit_time_particles{hill.obj.Id}{k},weights_{hill.obj.Id}{k},transit_distance{hill.obj.Id}{k},~,number_of_particle{hill.obj.Id}(k)]=get_ttd(trans.transport,trans.transport.t(k));
                    [time_support{hill.obj.Id},weighted_pdf_temp]=trans.transport.compute_full_distribution(transit_time_particles{hill.obj.Id}{k},weights_{hill.obj.Id}{k});
                    weighted_pdf{hill.obj.Id}(k,:)=weighted_pdf_temp;
                end
                tt=trans.transport.t;
                clear trans
            end
            
            watershed_simulation.t=tt;
            watershed_simulation.Q=Q;
            watershed_simulation.prof_WT_spat=prof_WT_spat;
            watershed_simulation.transit_time_particles=transit_time_particles;
            watershed_simulation.weights_=weights_;
            watershed_simulation.number_of_particle=number_of_particle;
            watershed_simulation.transit_distance=transit_distance;
            watershed_simulation.time_support=time_support;
            watershed_simulation.weighted_pdf=weighted_pdf;
        end
        
        function Q_out=run_simulation(f1,k1)
            folder_root='C:\Users\Jean\Documents\ProjectDSi\GuillecEquiv';
            obj=simulation_set(folder_root);
            obj=obj.instantiate_all_inputs_directory;
            
            for i=1:length(obj.morphologic_inputs_directory)               
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{i};
                
                % 2/ read input files
                M=obj.read_input_file(morpho_loc);
                %#JM change after test
                x=M(:,1); w=M(:,2); slope_angle=M(:,3);
                d=25*ones(size(x));
                f=f1*ones(size(x));
                k=k1*ones(size(x));
                hs1D=hillslope1D;
                hs1D=hs1D.set_properties(i,f,k);

                hs1D=hs1D.set_spatial_parameters(x,w,slope_angle,d);
                
                [M,input_type]=obj.read_input_file(hydro_loc);
                t=M(:,1);
                recharge_chronicle=(M(:,2))';
                ratio_P_R=1.33;%/0.38;
                source_terms=source('data_based');
                [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                ratio_P_R=1;%0.38;
                
                % 3/ create a runs object and run the simulation
                run_obj=runs;
                % set the solver options default or assigned in parameters via an odeset structure
                % specify Refine options for real infiltrations chronicle because for accuracy you need 
                % to force matlab ode15s to compute where you know sthg is happening
                odeset_struct=odeset('RelTol',1e-10);%2.5e-14);%,'Refine',-1);
                solver_options=run_obj.set_solver_options(odeset_struct);

                % run the simulation starting from half empty hillslope
                percentage_loaded=0;
                run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
                Q_temp=run_obj.simulation_results.compute_seepage_total;
                Q(i,:)=Q_temp(1531:1704);
                
            end
            Q_out=sum(Q,1);
        end
        
        function Q_out=run_simulation2(f1,k1,d1)
            folder_root='C:\Users\Jean\Documents\ProjectDSi\GuillecEquiv';
            obj=simulation_set(folder_root);
            obj=obj.instantiate_all_inputs_directory;
            
            for i=1:length(obj.morphologic_inputs_directory)               
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{i};
                
                % 2/ read input files
                M=obj.read_input_file(morpho_loc);
                %#JM change after test
                x=M(:,1); w=M(:,2); slope_angle=M(:,3);
                d=(linspace(d1,25,length(x)))';
                f=f1*ones(size(x));
                k=k1*ones(size(x));
                hs1D=hillslope1D;
                hs1D=hs1D.set_properties(i,f,k);

                hs1D=hs1D.set_spatial_parameters(x,w,slope_angle,d);
                
                [M,input_type]=obj.read_input_file(hydro_loc);
                t=M(:,1);
                recharge_chronicle=(M(:,2))';
                ratio_P_R=1.33;%/0.38;
                source_terms=source('data_based');
                [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                ratio_P_R=1;%0.38;
                
                % 3/ create a runs object and run the simulation
                run_obj=runs;
                % set the solver options default or assigned in parameters via an odeset structure
                % specify Refine options for real infiltrations chronicle because for accuracy you need 
                % to force matlab ode15s to compute where you know sthg is happening
                odeset_struct=odeset('RelTol',1e-10);%2.5e-14);%,'Refine',-1);
                solver_options=run_obj.set_solver_options(odeset_struct);

                % run the simulation starting from half empty hillslope
                percentage_loaded=0;
                run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
                Q_temp=run_obj.simulation_results.compute_seepage_total;
                Q(i,:)=Q_temp(1531:1704);
                
            end
            Q_out=sum(Q,1);
        end
        
% %         function [Q_out,time_aqui,DPSA_prop,residual]=run_simulation3(k1,d1,file_path,Q_mean)
        function [Q_out,time_aqui,DPSA_prop,residual,Q_out2,residual2]=run_simulation3(k1,d1,file_path,Q_mean)
% %             if(nargin<3)
% %                 f1=0.2;
% %             end
            tic
            range_=  1531:1704; %4332:4505; %
            if(nargin<3)
%                 folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Guillec2';
%                 d_init_add=2.7707;
                folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Synthetic';
                d_init_add=0;
            else
                occurence_slash=strfind(file_path,'\');
                folder_root=strcat(file_path(1:occurence_slash(end-1)),file_path(occurence_slash(end)+1:end-4));
                if(strcmp(file_path(occurence_slash(end)+1:end-4),'Douffine2'))
                    d_init_add=0;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Guillec2'))
                    d_init_add=2.7707;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Dossen'))
                    d_init_add=0.9249;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Dourduff'))
                    d_init_add=2.2063;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Penze'))
                    d_init_add=0.9392;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Ris'))
                    d_init_add=1.7655;
                else
                    d_init_add=0;
                end
            end
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\GuillecEquiv';
            
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Douffine2';
            obj=simulation_set(folder_root);
            obj=obj.instantiate_all_inputs_directory;
            
            for i=1:length(obj.morphologic_inputs_directory)               
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{i};
                
                % 2/ read input files
                M=obj.read_input_file(morpho_loc);
                %#JM change after test
                x=M(:,1); w=M(:,2); slope_angle=M(:,3);
                z_top=cumtrapz(x,slope_angle);
%                 d=(linspace(d1,z_top+d1,length(x)))';
                
                %#JM 20180223 change permeability homogeneity
                k=k1*ones(size(x));
                f=0.2*ones(size(x));
%                 angle_reduc=(d(2:end)-d(1:end-1))./(x(2:end)-x(1:end-1));
%                 angle_reduc=[angle_reduc;angle_reduc(end)];
%                 slope_angle2=slope_angle-angle_reduc;
                slope_angle2=d1*slope_angle;%atan((d1)/(x(end)-x(1)))*ones(size(slope_angle));%zeros(size(slope_angle));
                z_bottom=cumtrapz(x,slope_angle2)-(1-d1)*d_init_add;
                d=z_top-z_bottom;
                if(sum(d<0)==0)
%                 k=k1*(exp(0.1*(d(1)-d)));%k=k1*(d(1)./d).^2;
%                 f=0.4*d(1)./d;
%                 k=k1*1./(d-d(1)+1);
%                 f=0.4*1./(d-d(1)+1);
                hs1D=hillslope1D;
                hs1D=hs1D.set_properties(i,f,k);

                hs1D=hs1D.set_spatial_parameters(x,w,slope_angle2,d);
                
                [M,input_type]=obj.read_input_file(hydro_loc);
                t=M(:,1);
                recharge_chronicle=(M(:,2))';
                if(nargin>3)
                    ratio_P_R=Q_mean/(trapz(x,w)*mean(M(range_,2)));%*M(1200,2));%
                else
                    ratio_P_R=1;%0.875;%1;%.33;%/0.38;
                end
                source_terms=source('data_based');
                [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                ratio_P_R=1;%0.38;
                
                % 3/ create a runs object and run the simulation
                run_obj=runs;
                % set the solver options default or assigned in parameters via an odeset structure
                % specify Refine options for real infiltrations chronicle because for accuracy you need 
                % to force matlab ode15s to compute where you know sthg is happening
                odeset_struct=odeset('RelTol',1e-10);%2.5e-14);%,'Refine',-1);
                solver_options=run_obj.set_solver_options(odeset_struct);

                % run the simulation starting from half empty hillslope
                percentage_loaded=0;
                run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
                Q_temp=run_obj.simulation_results.compute_river_flow;%compute_seepage_total;
                [x_S1,w_1,d1_2,angle1,x_Q1,f1,k1_2]=get_resampled_variables(run_obj.boussinesq_simulation.discretization);
                slope_angle_top=interpn(x,slope_angle,x_Q1);
                Q_temp2=run_obj.simulation_results.compute_river_flow_with_rooting(k1_2,slope_angle_top);
                Q(i,:)=Q_temp(range_);
                Q2(i,:)=Q_temp2(range_);
%                 [~,DSi(i,:)]=transport.test(run_obj,folder_root);
                end
            end
            if(sum(d<0)==0)
            Q_out=sum(Q,1);
            Q_out2=sum(Q2,1);
% %             time_aqui=nan;
% %             DPSA_prop=nan;
            %#JM
            
%             h_sim=bsxfun(@rdivide,run_obj.simulation_results.S,f1.*w_1);
%             z_bottom1=cumtrapz(x_Q1,angle1)-10;
%             z_bottom1=interpn(x_Q1,z_bottom1,x_S1);
%             z_top1=z_bottom1+d1;
%             WT_height=bsxfun(@plus,h_sim,z_bottom1);
            Q_out=Q_out+run_obj.boussinesq_simulation.source_terms.recharge_chronicle(range_)*x_S1(1)*w_1(1);
            Q_out2=Q_out2+run_obj.boussinesq_simulation.source_terms.recharge_chronicle(range_)*x_S1(1)*w_1(1);
            Vol_aqui=trapz(x_S1(:),run_obj.simulation_results.S(:,1000));
            Area=trapz(x_S1,w_1);
            time_aqui=Vol_aqui/(Area*mean(run_obj.boussinesq_simulation.source_terms.recharge_chronicle(range_))*365*3600*24);
            [DPSA_tot,DPSA_spat]=run_obj.simulation_results.compute_DPES((x_Q1(2:end)-x_Q1(1:end-1)).*w_1,f1.*w_1.*d1_2,run_obj.boussinesq_simulation.source_terms.recharge_chronicle);
            DPSA_prop=nanmean(DPSA_tot(range_)./Q_out);
            toc
%  #JM           [~,DSi(i,:)]=transport.test(run_obj,folder_root);
%             save(strcat(folder_root,'\k_',num2str(k1),'_alpha_',num2str(d1),'.mat'),'-v7.3');
            load(file_path);
            residual=Q_out-Q_real;
            residual=nansum(residual.^2)/nansum((Q_real-nanmean(Q_real)).^2);
            residual2=Q_out2-Q_real;
            residual2=nansum(residual2.^2)/nansum((Q_real-nanmean(Q_real)).^2);
            else
                Q_out=nan;
                time_aqui=nan;
                DPSA_prop=nan;
            end
%             DSi_out=sum(DSi,1);      
        end
        
        function [Q_out2,residual2,Q_out,residual,run_obj,obj,x,w,slope_angle]=run_simulation_rooting(k1,soil_coef,file_path,f1)
% %             if(nargin<3)
% %                 f1=0.2;
% %             end
            tic
            range_=  1531:1704; %1:1465;%1:8759; %1:1500;%5332:5505;%4332:4505; %
            if(nargin<4)
                f1=0.2;
            end
            if(nargin<3)
%                 folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Guillec2';
%                 d_init_add=2.7707;
                folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Synthetic';
                d_init_add=0;
            else
                occurence_slash=strfind(file_path,'\');
                folder_root=strcat(file_path(1:occurence_slash(end-1)),file_path(occurence_slash(end)+1:end-4));
                if(strcmp(file_path(occurence_slash(end)+1:end-4),'Douffine2'))
                    d_init_add=0;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Guillec'))
                    d_init_add=2.7707;%24;%6.5;%
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Dossen'))
                    d_init_add=0.9249;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Dourduff'))
                    d_init_add=2.2063;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Penze'))
                    d_init_add=0.9392;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Ris'))
                    d_init_add=1.7655;
                else
                    d_init_add=2.7707;%20;
                end
            end
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\GuillecEquiv';
            
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Douffine2';
            obj=simulation_set(folder_root);
            obj=obj.instantiate_all_inputs_directory;
            
                           
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{1};
                
                % 2/ read input files
                M=obj.read_input_file(morpho_loc);
                x=M(:,1); w=M(:,2); slope_angle=M(:,3); z=M(:,4);
% %                 w=mean(w)*ones(size(x));
                slope_angle=0.4*ones(size(x));
                
                % first option
                z_top=cumtrapz(x,slope_angle);
                slope_angle2=0*slope_angle;%atan((d1)/(x(end)-x(1)))*ones(size(slope_angle));%zeros(size(slope_angle));
                z_bottom=cumtrapz(x,slope_angle2)-d_init_add;
                d=z_top-z_bottom;
                
                % second option
% %                 z_landsurface=z(1)+cumtrapz(x,slope_angle);
% %                 z_top=z_landsurface-2;%z_landsurface-[(linspace(2.5,0,10))';zeros(length(z_landsurface)-10,1)]%-[(linspace(2.5,0,35))';zeros(length(z_landsurface)-35,1)];%-(linspace(1.5,0,length(z_landsurface)))';
% % % % % % % %                 z_top(x<=210)=z_top(x<=210)-1.5;
% % % % % % % %                 z_top(x==200.0433)=277;
% % % % % % % %                 z_top(x==205.0443)=277;
% % % % % %                 z_top=z(1)+cumtrapz(x,slope_angle)-1.5; % take the peat layer into account assume it has a thickness of 1.5m
% %                 z_bottom=z_top-d_init_add;
% %                 d=z_top-z_bottom;
% %                 slope_angle2=slope_angle;
% % %                 slope_angle2=0.0*ones(size(slope_angle));%mean(slope_angle)
% % %                 z_bottom1=z_top(1)-d_init_add(1);
% % %                 z_bottom=z_bottom1+cumtrapz(x,slope_angle2);
% % %                 d=z_top-z_bottom;
                

                % third option
%                 z_top=cumtrapz(x,slope_angle);
%                 d=(linspace(d_init_add,8,length(z_top)))';
%                 z_bottom=z_top-d;
%                 slope_angle2=slope_angle;

                % fourth option
% %                 z_top=z;
% %                 z_bottom=z_top-d_init_add;
% %                 d=z_top-z_bottom;
% %                 slope_angle2=slope_angle;
                
                % hydraulic parameters
                k=k1*ones(size(x));
                f=f1*ones(size(x));
                
                
                if(sum(d<0)==0)
                   
                    hs1D=hillslope1D;
                    hs1D=hs1D.set_properties(1,f,k);
                    
                    hs1D=hs1D.set_spatial_parameters(x,w,slope_angle2,d);
                    
                    [M,input_type]=obj.read_input_file(hydro_loc);
                    t=M(:,1);
                    recharge_chronicle=(M(:,2:end))';
                    
% %                     time_1=time_properties(t(1),t(end),(length(t)-1)*4+1,'sec');
% %                     source_terms=source('data_based');
% %                     source_terms.time=time_1;
% %                     source_terms.recharge_chronicle=interp1((t)',(M(:,2:end))',(time_1.get_properties));
% %                     source_terms.recharge_mean=mean(source_terms.recharge_chronicle,2);
                    
                    ratio_P_R=1;%0.875;%.33;%/0.38;
                    source_terms=source('data_based');
                    [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                    ratio_P_R=1;%0.38;
                    
                    % 3/ create a runs object and run the simulation
                    run_obj=runs;
                    % set the solver options default or assigned in parameters via an odeset structure
                    % specify Refine options for real infiltrations chronicle because for accuracy you need
                    % to force matlab ode15s to compute where you know sthg is happening
                    dt=diff(t);
                    dt_mean=mean(dt);
                    odeset_struct=odeset('RelTol',1e-5,'MaxStep',dt_mean);%30*2.5e-14);%,'Refine',-1);%odeset('RelTol',1e-3);%,'AbsTol',1e-7);%
                    solver_options=run_obj.set_solver_options(odeset_struct);
                    
% % %                     % run the simulation starting from half empty hillslope
% % %                     percentage_loaded=0.5;    
% % %                     % run transient simulation 
% % %                     run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
                    % run the simulation starting from the steady state condition
                    percentage_loaded=0;
                    recharge_averaged=1e3*24*3600*source_terms.recharge_mean; % recharge averaged in mm/d
                    state_values_initial=obj.prerun_steady_state(hs1D,recharge_averaged,ratio_P_R,'empty');
                    presteadystate_percentage_loaded=-2; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
                    % run transient simulation 
                    run_obj=run_obj.run_simulation(hs1D,source_terms,presteadystate_percentage_loaded,solver_options,ratio_P_R,state_values_initial,'empty');
                    
                    [x_S1,w_1,d1_2,angle1,x_Q1,f1,k1_2]=get_resampled_variables(run_obj.boussinesq_simulation.discretization);
                    slope_angle_top=interpn(x,slope_angle,x_Q1);
                    toc
                    
                    Q_temp=run_obj.simulation_results.compute_river_flow;%compute_seepage_total;
                    Q_out=Q_temp(range_);
                    
                    Q_temp2=run_obj.simulation_results.compute_river_flow_with_rooting(k1_2,slope_angle_top,soil_coef);
                    Q_out2=Q_temp2(range_);

                    
                    Q_out=Q_out+run_obj.boussinesq_simulation.source_terms.recharge_chronicle(range_)*x_S1(1)*w_1(1);
                    Q_out2=Q_out2+run_obj.boussinesq_simulation.source_terms.recharge_chronicle(range_)*x_S1(1)*w_1(1);
%                     toc
                    t1=datetime(1998,01,15);
                    t2=datetime(2012,06,15);
                    tt=t1:calmonths(1):t2;
                    
                    residual=nan;
                    residual2=nan;
%                     load(file_path);
%                     Q_real=(Ecoflux_data.Q(1:174))';
%                     t_real=Ecoflux_data.t(1:174);
%                     
%                     if(length(Q_real)==length(Q_out))
%                         residual2=Q_out2-Q_real;
%                         residual2=nansum(residual2.^2)/nansum((Q_real-nanmean(Q_real)).^2);
%                         
%                         residual=Q_out-Q_real;
%                         residual=nansum(residual.^2)/nansum((Q_real-nanmean(Q_real)).^2);
%                     else
%                         residual2=nan;
%                         residual=nan;
%                     end
                else
                    Q_out2=nan;
                    
                end
                %% test
% % %                 f_soil=0.02;
% % %                 k_soil=2;
% % %                 d_soil=3;
% % %                 run_soil=obj.run_simulation_from_struct(x,f_soil,k_soil,w,slope_angle,d_soil*ones(size(x)),run_obj);
%             DSi_out=sum(DSi,1);      
        end
        
        function [Q_out2,residual2,Q_out,residual,run_obj]=run_simulation_unsat(k1,f1,phi1,file_path)
% %             if(nargin<3)
% %                 f1=0.2;
% %             end
            tic
            range_= 4332:4505; % 1531:1704; %1:1465;%1:8759; %1:1500;%
            if(nargin<4)
                f1=0.2;
            end
            if(nargin<3)
%                 folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Guillec2';
%                 d_init_add=2.7707;
                folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Synthetic';
                d_init_add=0;
            else
                occurence_slash=strfind(file_path,'\');
                folder_root=strcat(file_path(1:occurence_slash(end-1)),file_path(occurence_slash(end)+1:end-4));
                if(strcmp(file_path(occurence_slash(end)+1:end-4),'Douffine2'))
                    d_init_add=0;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Guillec2'))
                    d_init_add=2.7707;%6.5;%
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Dossen'))
                    d_init_add=0.9249;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Dourduff'))
                    d_init_add=2.2063;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Penze'))
                    d_init_add=0.9392;
                elseif(strcmp(file_path(occurence_slash(end)+1:end-4),'Ris'))
                    d_init_add=1.7655;
                else
                    d_init_add=2.7707;%20;
                end
            end
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\GuillecEquiv';
            
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Douffine2';
            obj=simulation_set(folder_root);
            obj=obj.instantiate_all_inputs_directory;
            
                           
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{1};
                
                % 2/ read input files
                M=obj.read_input_file(morpho_loc);
                x=M(:,1); w=M(:,2); slope_angle=M(:,3); z=M(:,4);
                
                z_top=z(1)+cumtrapz(x,slope_angle);
                slope_angle2=([(linspace(0,0.1,6)),(linspace(0.1,1,length(x)-6)).^0.2])'.*slope_angle;%(linspace(0.7,1,length(x)))'.*slope_angle;%slope_angle;
                z_bottom=z_top(1)+cumtrapz(x,slope_angle2);
                d=z_top-z_bottom;
               
                figure; hold on
                plot(x,z_bottom)
                plot(x,z_top)
                
                
%                 M=obj.read_input_file(morpho_loc);
%                 x=M(:,1); w=M(:,2); slope_angle=M(:,3); z=M(:,4);
%                 
%                 % first option
%                 z_top=cumtrapz(x,slope_angle);
%                 slope_angle2=0*slope_angle;%atan((d1)/(x(end)-x(1)))*ones(size(slope_angle));%zeros(size(slope_angle));
%                 z_bottom=cumtrapz(x,slope_angle2)-d_init_add;
%                 d=z_top-z_bottom;
                
                % hydraulic parameters
                k=k1*ones(size(x));
                f=f1*ones(size(x));
                phi=phi1*ones(size(x));
                
                if(sum(d<0)==0)
                   
                    hs1D=hillslope1D_unsat;
                    hs1D=hs1D.set_properties(-1,f,k,phi);
                    
                    hs1D=hs1D.set_spatial_parameters(x,w,slope_angle2,d);
                    
                    [M,input_type]=obj.read_input_file(hydro_loc);
                    t=M(:,1);
                    recharge_chronicle=(M(:,2:end))';
                    
% %                     time_1=time_properties(t(1),t(end),(length(t)-1)*4+1,'sec');
% %                     source_terms=source('data_based');
% %                     source_terms.time=time_1;
% %                     source_terms.recharge_chronicle=interp1((t)',(M(:,2:end))',(time_1.get_properties));
% %                     source_terms.recharge_mean=mean(source_terms.recharge_chronicle,2);
                    
                    ratio_P_R=1;%0.875;%.33;%/0.38;
                    source_terms=source('data_based');
                    [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                    ratio_P_R=1;%0.38;
                    
                    % 3/ create a runs object and run the simulation
                    run_obj=runs;
                    % set the solver options default or assigned in parameters via an odeset structure
                    % specify Refine options for real infiltrations chronicle because for accuracy you need
                    % to force matlab ode15s to compute where you know sthg is happening
                    odeset_struct=odeset('RelTol',1e-5,'MaxStep',30*3600*24);%2.5e-14);%,'Refine',-1);%odeset('RelTol',1e-3);%,'AbsTol',1e-7);%
                    solver_options=run_obj.set_solver_options(odeset_struct);
                    
% % %                     % run the simulation starting from half empty hillslope
% % %                     percentage_loaded=0.5;    
% % %                     % run transient simulation 
% % %                     run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
                    % run the simulation starting from the steady state condition
                    percentage_loaded=0;
                    recharge_averaged=1e3*24*3600*source_terms.recharge_mean; % recharge averaged in mm/d
                    state_values_initial=obj.prerun_steady_state(hs1D,recharge_averaged,ratio_P_R,'empty');
                    presteadystate_percentage_loaded=-2; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
                    % run transient simulation 
                    run_obj=run_obj.run_simulation(hs1D,source_terms,presteadystate_percentage_loaded,solver_options,ratio_P_R,state_values_initial,'empty');
                    
                    [x_S1,w_1,d1_2,angle1,x_Q1,f1,k1_2]=get_resampled_variables(run_obj.boussinesq_simulation.discretization);
                    slope_angle_top=interpn(x,slope_angle,x_Q1);
                    toc
                    
                    Q_temp=run_obj.simulation_results.compute_river_flow;%compute_seepage_total;
                    Q_out=Q_temp(range_);
                    
                    Q_temp2=run_obj.simulation_results.compute_river_flow_with_rooting(k1_2,slope_angle_top,soil_coef);
                    Q_out2=Q_temp2(range_);

                    
                    Q_out=Q_out+run_obj.boussinesq_simulation.source_terms.recharge_chronicle(range_)*x_S1(1)*w_1(1);
                    Q_out2=Q_out2+run_obj.boussinesq_simulation.source_terms.recharge_chronicle(range_)*x_S1(1)*w_1(1);
%                     toc
                    
                    load(file_path);
                    
                    if(length(Q_real)==length(Q_out))
                        residual2=Q_out2-Q_real;
                        residual2=nansum(residual2.^2)/nansum((Q_real-nanmean(Q_real)).^2);
                        
                        residual=Q_out-Q_real;
                        residual=nansum(residual.^2)/nansum((Q_real-nanmean(Q_real)).^2);
                    else
                        residual2=nan;
                        residual=nan;
                    end
                else
                    Q_out2=nan;
                    
                end
%             DSi_out=sum(DSi,1);      
        end
        
        
        function run_simulation_temp(f1,k1,d1,file_path)
            obj=simulation_set(file_path);
            obj=obj.instantiate_all_inputs_directory;
            
            
            hydro_loc=obj.hydrologic_inputs_directory{1};
            morpho_loc=obj.morphologic_inputs_directory{1};
            
            % 2/ read input files
            M=obj.read_input_file(morpho_loc);
            %#JM change after test
            x=M(:,1); w=M(:,2); slope_angle=M(:,3);
            z_top=cumtrapz(x,slope_angle);
            z_bottom=z_top-d1;
            k=k1*ones(size(x));
            f=f1*ones(size(x));
            d=z_top-z_bottom;
            
            hs1D=hillslope1D;
            hs1D=hs1D.set_properties(1,f,k);
            
            hs1D=hs1D.set_spatial_parameters(x,w,slope_angle,d);
            
            [M,input_type]=obj.read_input_file(hydro_loc);
            t=M(:,1)*3600*24;
            recharge_chronicle=(M(:,2))';
            
            ratio_P_R=1;
            source_terms=source('data_based');
            [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
            ratio_P_R=1;%0.38;
            
            % 3/ create a runs object and run the simulation
            run_obj=runs;
            % set the solver options default or assigned in parameters via an odeset structure
            % specify Refine options for real infiltrations chronicle because for accuracy you need
            % to force matlab ode15s to compute where you know sthg is happening
            odeset_struct=odeset('RelTol',1e-10);%2.5e-14);%,'Refine',-1);
            solver_options=run_obj.set_solver_options(odeset_struct);
            
            % run the simulation starting from half empty hillslope
            tic;
            percentage_loaded=0;
            recharge_averaged=1e3*24*3600*source_terms.recharge_mean; % recharge averaged in mm/d
            state_values_initial=obj.prerun_steady_state(hs1D,recharge_averaged,ratio_P_R);
            presteadystate_percentage_loaded=-2; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
            % run the simulation starting from the steady state condition
            toc; 
            tic;
            run_obj=run_obj.run_simulation(hs1D,source_terms,presteadystate_percentage_loaded,solver_options,ratio_P_R,state_values_initial);
            toc;
%             state_values_initial=prerun_steady_state(obj,hs1D,recharge_averaged,ratio_P_R)
%             run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
%             Q_temp=run_obj.simulation_results.compute_river_flow;%compute_seepage_total;

        end
        
        function [run_bedrock,run_regolith,run_soil]=simulation_3compartments(f_bed,k_bed,f_reg,k_reg,alpha_slope_reg,f_soil,k_soil,d_soil_init,file_path)
                occurence_slash=strfind(file_path,'\');
                folder_root=strcat(file_path(1:occurence_slash(end-1)),file_path(occurence_slash(end)+1:end-4));
                
                
                obj=simulation_set(folder_root);
                obj=obj.instantiate_all_inputs_directory;
                
                
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{1};
                
                % 1/ set structure file
                M=obj.read_input_file(morpho_loc);
                x=M(:,1); width=M(:,2); slope_angle=M(:,3); z=M(:,4);
                %#JM change 1
                z_land_surface=z(1)+cumtrapz(x,slope_angle);
                slope_bottom_soil=[0.2*slope_angle(1:5);slope_angle(6:end)];%slope_angle;%(linspace(0.95,1,length(x)))'.*slope_angle;%
                z_bottom_soil=z(1)+cumtrapz(x,slope_bottom_soil)-d_soil_init;
                
                slope_bottom_regolith=[0.16*slope_angle(1:11);slope_angle(12:end)];%((linspace(0.1^3,1^3,length(x))).^0.333)'.*slope_bottom_soil;%(mean(slope_bottom_soil)-0.04); % instead of 0.04 -> to get the 8m regolith depth on average on the hillslope
                z_bottom_regolith=z_land_surface(1)+cumtrapz(x,slope_bottom_regolith.*ones(size(x)))-alpha_slope_reg; % chose the slope such as in average z_land_surface-z_top=8m
                
                slope_bottom_bedrock=0*slope_angle;%(mean(slope_bottom_regolith)-0.1)*ones(size(x));%
                z_bottom_bedrock=z_bottom_regolith(1)+cumtrapz(x,slope_bottom_bedrock);
                d_bed=z_bottom_regolith-z_bottom_bedrock;
                d_reg=z_bottom_soil-z_bottom_regolith;
                d_soil=z_land_surface-z_bottom_soil;
                
                slope_angle=[0.4*slope_angle(1:5);slope_angle(6:end)];
                z_land_surface=z(1)+cumtrapz(x,slope_angle);
% % %                 width=trapz(M(:,1),M(:,2))/M(end,1)*ones(size(M(:,1)));
% %                 
% % %                 z_land_surface=z(1)+cumtrapz(x,slope_angle);
% % %                 slope_bottom_soil=slope_angle;
% % %                 z_bottom_soil=z_land_surface-d_soil_init;
% % %                 d_soil=z_land_surface-z_bottom_soil;
% % %                 slope_bottom_regolith=slope_angle; 
% % %                 z_bottom_regolith=z_bottom_soil-alpha_slope_reg;
% % %                 d_reg=z_bottom_soil-z_bottom_regolith;
% % %                 slope_bottom_bedrock=(0)*ones(size(x));
% % %                 z_bottom_bedrock=z_bottom_regolith(1)+cumtrapz(x,slope_bottom_bedrock)-2.7707;
% % %                 d_bed=z_bottom_regolith-z_bottom_bedrock;
% %                 
% %                 z_land_surface=z(1)+cumtrapz(x,slope_angle);
% %                 slope_bottom_soil=(linspace(0.9,1,length(x)))'.*slope_angle;
% %                 z_bottom_soil=z(1)+cumtrapz(x,slope_bottom_soil);
% %                 d_soil=z_land_surface-z_bottom_soil;
% % % %                 d_soil=d_soil(end:-1:1);
% % % %                 z_bottom_soil=z_land_surface-d_soil;
% % % %                 slope_bottom_soil=atan((z_bottom_soil(2:end)-z_bottom_soil(1:end-1))./(x(2:end)-x(1:end-1)));
% % % %                 slope_bottom_soil=[slope_bottom_soil;slope_bottom_soil(end)];
% %                 
% %                 slope_bottom_regolith=(linspace(0.85,1,length(x)))'.*slope_angle;%slope_bottom_soil; 
% %                 z_bottom_regolith=z_bottom_soil(1)+cumtrapz(x,slope_bottom_regolith);
% %                 d_reg=z_bottom_soil-z_bottom_regolith;
% % % %                 d_reg=d_reg(end:-1:1);
% % % %                 z_bottom_regolith=z_bottom_soil-d_reg;
% % % %                 slope_bottom_regolith=atan((z_bottom_regolith(2:end)-z_bottom_regolith(1:end-1))./(x(2:end)-x(1:end-1)));
% % % %                 slope_bottom_regolith=[slope_bottom_regolith;slope_bottom_regolith(end)];
% %                 
% %                 slope_bottom_bedrock=(linspace(0,1,length(x)))'.*slope_angle; %slope_bottom_bedrock=(linspace(0.4,1,length(x)))'.*slope_bottom_regolith; 
% %                 z_bottom_bedrock=z_bottom_regolith(1)+cumtrapz(x,slope_bottom_bedrock)-1;%-5.271;%%2.7707;
% %                 d_bed=z_bottom_regolith-z_bottom_bedrock;
% %                 
% %                 %                 slope_bottom_regolith=(trapz(x,width.*slope_angle)/trapz(x,width)/alpha_slope_reg)*ones(size(x));
% %                 %                 z_bottom_regolith=z_bottom_soil(1)+cumtrapz(x,slope_bottom_regolith);
% %                 %                 d_reg=z_bottom_soil-z_bottom_regolith;
                
                figure; hold on
                plot(x,z_bottom_bedrock)
                plot(x,z_bottom_regolith)
                plot(x,z_bottom_soil)
                plot(x,z_land_surface)
                
                % 2/ run the boussinesq simulations
                tic
                run_bedrock=obj.run_simulation_from_struct(x,f_bed,k_bed,width,slope_bottom_bedrock,d_bed);
                toc
                run_regolith=obj.run_simulation_from_struct(x,f_reg,k_reg,width,slope_bottom_regolith,d_reg,run_bedrock);
                toc
                run_soil=obj.run_simulation_from_struct(x,f_soil,k_soil,width,slope_bottom_soil,d_soil,run_regolith);
                toc
                
                simulation_set.plot_storage_structure(run_bedrock,run_regolith,run_soil,x,z_land_surface);
                
                [DPSA_soil,RF_soil]=compute_DPSA_RF(run_soil.simulation_results,run_soil.boussinesq_simulation);
                [~,RF_reg]=compute_DPSA_RF(run_regolith.simulation_results,run_regolith.boussinesq_simulation);
                [~,RF_bed]=compute_DPSA_RF(run_bedrock.simulation_results,run_bedrock.boussinesq_simulation);
                RF_bed(RF_bed<0)=0;
                t=datetime(datestr(run_bedrock.simulation_results.t/(24*3600)));
                figure; hold on
                plot(t,RF_bed)
                plot(t,RF_bed+RF_reg)
                plot(t,RF_bed+RF_reg+RF_soil)
                plot(t,RF_bed+RF_reg+RF_soil+DPSA_soil)
        end
        
        function plot_storage_structure(run_bedrock,run_regolith,run_soil,x,z_land_surface)
            [x_S_bed,w_bed,d_bed,angle_bed,x_Q_bed,f_bed,k_bed]=get_resampled_variables(run_bedrock.boussinesq_simulation.discretization);
            [x_S_reg,w_reg,d_reg,angle_reg,x_Q_reg,f_reg,k_reg]=get_resampled_variables(run_regolith.boussinesq_simulation.discretization);
            [x_S_soil,w_soil,d_soil,angle_soil,x_Q_soil,f_soil,k_soil]=get_resampled_variables(run_soil.boussinesq_simulation.discretization);
            
            h_bed=bsxfun(@rdivide,run_bedrock.simulation_results.S,f_bed.*w_bed);
            h_reg=bsxfun(@rdivide,run_regolith.simulation_results.S,f_reg.*w_reg);
            h_soil=bsxfun(@rdivide,run_soil.simulation_results.S,f_soil.*w_soil);
            
            Volume_soil=trapz(x_S_soil,run_soil.simulation_results.S);
            time_step_wettest=find(Volume_soil==max(Volume_soil));
            Volume_bed=trapz(x_S_bed,run_bedrock.simulation_results.S);
            time_step_driest=find(Volume_bed==min(Volume_bed));
            Volume_reg=trapz(x_S_reg,run_regolith.simulation_results.S);
            if(length(time_step_wettest)>1)
                time_step_wettest=find(Volume_reg==max(Volume_reg));
            end
            if(length(time_step_driest)>1)
                time_step_driest=find(Volume_reg==min(Volume_reg));
            end
            
            angle_soil=interpn(x_Q_soil,angle_soil,x_S_soil);
            angle_reg=interpn(x_Q_reg,angle_reg,x_S_reg);
            angle_bed=interpn(x_Q_bed,angle_bed,x_S_bed);
            
            z_land_surface=interpn(x,z_land_surface,x_S_soil);
            z_bottom_soil=z_land_surface(1)-d_soil(1)+cumtrapz(x_S_soil,angle_soil);
            z_bottom_reg=z_bottom_soil(1)-d_reg(1)+cumtrapz(x_S_reg,angle_reg);
            z_bottom_bed=z_bottom_reg(1)-d_bed(1)+cumtrapz(x_S_bed,angle_bed);
            
            
            figure; hold on
            
            plot(x_S_bed,z_bottom_bed+h_bed,'Color',[0, 0.75, 0.75])
            plot(x_S_reg,z_bottom_reg+h_reg,'Color',[0, 0.75, 0.75])%,'Color',[ 0    0.4470    0.7410])
            plot(x_S_soil,z_bottom_soil+h_soil,'Color',[0, 0.75, 0.75])
            
%             plot(x_S_bed,z_bottom_bed+h_bed(:,time_step_driest),'Color',[1, 0, 0],'LineWidth',2)
%             plot(x_S_reg,z_bottom_reg+h_reg(:,time_step_driest),'Color',[1, 0, 0],'LineWidth',2)%,'Color',[ 0    0.4470    0.7410])
%             plot(x_S_soil,z_bottom_soil+h_soil(:,time_step_driest),'Color',[1, 0, 0],'LineWidth',2)%,'Color',[0.3010    0.7450    0.9330])
%             
%             plot(x_S_bed,z_bottom_bed+h_bed(:,time_step_wettest),'Color',[ 0         0    1.0000],'LineWidth',2)
%             plot(x_S_reg,z_bottom_reg+h_reg(:,time_step_wettest),'Color',[ 0         0    1.0000],'LineWidth',2)%,'Color',[ 0    0.4470    0.7410])
%             plot(x_S_soil,z_bottom_soil+h_soil(:,time_step_wettest),'Color',[ 0         0    1.0000],'LineWidth',2)%,'Color',[0.3010    0.7450    0.9330])
            
            plot(x_S_soil,z_land_surface,'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880])%,'Color',[0.4660    0.6740    0.1880])
            plot(x_S_soil,z_bottom_soil,'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250])%,'Color',[0.2500    0.2500    0.2500])
            plot(x_S_reg,z_bottom_reg,'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980])%,'Color',[0.2500    0.2500    0.2500])
            plot(x_S_bed,z_bottom_bed,'LineWidth',2,'Color',[0.25, 0.25, 0.25])%,'Color',[0.2500    0.2500    0.2500])
            
%             plot(x_S_bed,z_bottom_bed+h_bed(:,(h_bed(8,:)==max(h_bed(8,:)))),'Color',[ 0         0    1.0000])
%             plot(x_S_reg,z_bottom_reg+h_reg(:,(h_reg(8,:)==max(h_reg(8,:)))),'Color',[ 0    0.4470    0.7410])
%             plot(x_S_soil,z_bottom_soil+h_soil(:,(h_soil(8,:)==max(h_soil(8,:)))),'Color',[0.3010    0.7450    0.9330])
        end
        
        function [run_bedrock,run_soil]=simulation_2compartments_bis(f_bed,k_bed,f_soil,k_soil,d_soil_init,file_path)
                occurence_slash=strfind(file_path,'\');
                folder_root=strcat(file_path(1:occurence_slash(end-1)),file_path(occurence_slash(end)+1:end-4));
                
                
                obj=simulation_set(folder_root);
                obj=obj.instantiate_all_inputs_directory;
                
                
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{1};
                
                % 1/ set structure file
                M=obj.read_input_file(morpho_loc);
                x=M(:,1); width=M(:,2); slope_angle=M(:,3); z=M(:,4);
                
                z_land_surface=z(1)+cumtrapz(x,slope_angle);
                slope_bottom_soil=slope_angle;%([(linspace(0,0.1,6)),(linspace(0.1,1,length(x)-6)).^0.2])'.*slope_angle;%(linspace(0.7,1,length(x)))'.*slope_angle;%
                z_bottom_soil=z_land_surface(1)+cumtrapz(x,slope_bottom_soil)-d_soil_init+2;
%                 z_bottom_soil=(z_land_surface).^(0.9);
%                 slope_bottom_soil=(z_bottom_soil(2:end)-z_bottom_soil(1:end-1))./(x(2:end)-x(1:end-1));
%                 slope_bottom_soil=[slope_bottom_soil;slope_bottom_soil(end)];
                d_soil=z_land_surface-z_bottom_soil;
               
                slope_bottom_bedrock=slope_angle;%zeros(size(slope_angle));%(linspace(0.7,1,length(x)))'.*slope_angle;%(trapz(x,width.*slope_angle)/trapz(x,width)/1.3)*ones(size(x));%(0)*ones(size(x));
                z_bottom_bedrock=z_bottom_soil(1)+cumtrapz(x,slope_bottom_bedrock)-d_soil_init+2;%-8+d_soil_init;
                d_bed=z_bottom_soil-z_bottom_bedrock;
                
                figure; hold on
                plot(x,z_bottom_bedrock)
                plot(x,z_bottom_soil)
                plot(x,z_land_surface)
                
                % 2/ run the boussinesq simulations
                tic
                run_bed1=obj.run_simulation_from_struct(x,f_soil,k_soil,width,slope_bottom_soil,d_soil);
                toc
                run_bedrock=obj.run_simulation_from_struct(x,f_bed,k_bed,width,slope_bottom_bedrock,d_bed);
                toc
                run_soil=obj.run_simulation_from_struct(x,f_soil,k_soil,width,slope_bottom_soil,d_soil,run_bedrock);
                toc
                
                [DPSA_soil,RF_soil]=compute_DPSA_RF(run_soil.simulation_results,run_soil.boussinesq_simulation);
                [~,RF_bed]=compute_DPSA_RF(run_bedrock.simulation_results,run_bedrock.boussinesq_simulation);
%                 t=datetime(datestr(run_bedrock.simulation_results.t/(24*3600)));
%                 figure; hold on
%                 plot(t,RF_bed)
%                 plot(t,RF_bed+RF_soil)
%                 plot(t,RF_bed+RF_soil+DPSA_soil)
        end
        
        function simulation_2compartments(k_bed,f_bed,k_soil,f_soil,file_path)
            soil_coef=1;
            bound_river_soil='empty';
            [~,~,Q_month,~,run_obj]=simulation_set.run_simulation_rooting(k_bed,soil_coef,file_path,f_bed);
            [x_S1,w_1,d1_2,angle1,x_Q1,f1,k1_2]=get_resampled_variables(run_obj.boussinesq_simulation.discretization);
           
            % 2nd option
            [DPSA_bed,RF_bed,DPSA_spat_bed,RF_spat_bed]=compute_DPSA_RF(run_obj.simulation_results,run_obj.boussinesq_simulation);

            
            qs1=DPSA_spat_bed;
            dx_Q=run_obj.simulation_results.x_Q(2:end)-run_obj.simulation_results.x_Q(1:end-1);
%             qs1(1,:)=run_obj.boussinesq_simulation.source_terms.recharge_chronicle*dx_Q(1)*w_1(1);
%             if(sum(qs1(:)<0)>0)
%                 Min_value_replaced=min(qs1(qs1<0));
%                 Number_values_replaced=sum(qs1(:)<0);
%                 fprintf(strcat('Warning: \n ',num2str(Number_values_replaced),' values have been replaced that were negative \n Values are between 0 and ',num2str(Min_value_replaced),'\n'));
%                 qs1(qs1<0)=0;
%             end
            qs1=bsxfun(@rdivide,qs1,w_1.*dx_Q);
            time_2006_2=time_properties(run_obj.simulation_results.t(1),run_obj.simulation_results.t(end),length(run_obj.simulation_results.t),'sec');
            seep=source('data_based');
            seep.time=time_2006_2;
            seep.recharge_chronicle=qs1;
            seep.recharge_mean=mean(seep.recharge_chronicle,2);
            tic
            occurence_slash=strfind(file_path,'\');
            folder_root=strcat(file_path(1:occurence_slash(end-1)),file_path(occurence_slash(end)+1:end-4));
            obj2=simulation_set(folder_root);
            
            x=x_Q1;
            z_bottom=cumtrapz(x,angle1);
            z_bottom_soil=z_bottom+[d1_2;2*d1_2(end)-d1_2(end-1)];
            z_landsurface=z_bottom_soil+1.5;
            d_soil=z_landsurface-z_bottom_soil;
            slope_angle2=(z_bottom_soil(2:end)-z_bottom_soil(1:end-1))./(x(2:end)-x(1:end-1));
            slope_angle2=[slope_angle2;slope_angle2(end)];
            w=[w_1;w_1(end)];

            k=k_soil*ones(size(x));
            f=f_soil*ones(size(x));
            hs1D=hillslope1D;
            hs1D=hs1D.set_properties(1,f,k);
            hs1D=hs1D.set_spatial_parameters(x,w,slope_angle2,d_soil);
            
            % run steady state simulation
            percentage_loaded=0;
            recharge_averaged=1e3*24*3600*seep.recharge_mean; % recharge averaged in mm/d
            ratio_P_R=1;
            state_values_initial=obj2.prerun_steady_state(hs1D,recharge_averaged,ratio_P_R,bound_river_soil);
            presteadystate_percentage_loaded=-2; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
            
            % run transient simulation
            run_obj2=runs;
            odeset_struct=odeset('RelTol',1e-5,'AbsTol',1e-6,'MaxStep',3600*24);%odeset('RelTol',2.5e-14,'AbsTol',1e-12);%1e-10);%,'Refine',-1);
            solver_options=run_obj2.set_solver_options(odeset_struct);
            run_obj2=run_obj2.run_simulation(hs1D,seep,presteadystate_percentage_loaded,solver_options,ratio_P_R,state_values_initial,bound_river_soil);
            toc
            
            t=datetime(datestr(run_obj.simulation_results.t/(24*3600)));
            load('H:\Data_Guillec_Daily\input_output_Guillec.mat');
            [DPSA_soil,RF_soil,DPSA_spat_soil,RF_spat_soil]=compute_DPSA_RF(run_obj2.simulation_results,run_obj2.boussinesq_simulation);
            
            
%             GW_bed=-run_obj.simulation_results.Q(2,:);
%             Flow_bed=GW_bed+RF_bed;
            Flow_bed=RF_bed;
%             GW_reg=-run_obj2.simulation_results.Q(2,:);
%             Flow_reg=GW_reg+RF_soil;%-RF_bed;
            Flow_reg=RF_soil;%-RF_bed;
            Flow_soil=DPSA_soil;%+run_obj2.boussinesq_simulation.source_terms.recharge_chronicle(1,:)*(x_Q1(2)-x_Q1(1))*w_1(1);
            
            figure; hold on
            plot(time_input(1:1500),Q_output(1:1500))
            plot(t,Q_month)
            plot(t,Flow_bed+Flow_reg+Flow_soil)
        end
        
        function [Q_out,DSi_out]=run_sim_transport(k1,d1,d2)
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\PF';
%             folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Synthetic';
%              folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Douffine2';
             folder_root='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\Guillec2';
            obj=simulation_set(folder_root);
            obj=obj.instantiate_all_inputs_directory;
            
            for i=1:length(obj.morphologic_inputs_directory)               
                % locations of different inputs file
                hydro_loc=obj.hydrologic_inputs_directory{1};
                morpho_loc=obj.morphologic_inputs_directory{i};
                
                % 2/ read input files
                M=obj.read_input_file(morpho_loc);
                %#JM change after test
                x=M(:,1); w=M(:,2); slope_angle=M(:,3);
                z_top=cumtrapz(x,slope_angle);
                z_bottom1=z_top(1)-d1; z_bottom2=z_top(end)-d2;
                
                %z_bottom=(z_bottom2-z_bottom1)/(z_top(end)-z_top(1))*(z_top-z_top(1))+z_bottom1;
                N=length(x);
                z_bottom=(linspace(z_bottom1,z_bottom2,N))';
                slope_angle_bottom=(z_bottom(2:end)-z_bottom(1:end-1))./(x(2:end)-x(1:end-1));
                slope_angle_bottom=[slope_angle_bottom;slope_angle(end)*slope_angle_bottom(end)/slope_angle(end-1)];
                d=z_top-z_bottom;

%                 slope_angle_bottom=zeros(size(x));
%                 z_bottom=-20*ones(size(x));
%                 d=z_top-z_bottom;
                
                f=0.2*ones(size(x));
                k=k1*ones(size(x));
                hs1D=hillslope1D;
                hs1D=hs1D.set_properties(i,f,k);
                hs1D=hs1D.set_spatial_parameters(x,w,slope_angle_bottom,d);
                
                [M,input_type]=obj.read_input_file(hydro_loc);
                t=M(:,1);
                recharge_chronicle=(M(:,2))';
                ratio_P_R=1.0417;%1.2432;%/0.38;
                source_terms=source('data_based');
                [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                ratio_P_R=1;%0.38;
                
                % 3/ create a runs object and run the simulation
                run_obj=runs;
                % set the solver options default or assigned in parameters via an odeset structure
                % specify Refine options for real infiltrations chronicle because for accuracy you need 
                % to force matlab ode15s to compute where you know sthg is happening
                odeset_struct=odeset('RelTol',1e-10);%2.5e-14);%,'Refine',-1);
                solver_options=run_obj.set_solver_options(odeset_struct);

                % run the simulation starting from half empty hillslope
                percentage_loaded=0;
                run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options,ratio_P_R);
                Q_temp=run_obj.simulation_results.compute_seepage_total;
                Q(i,:)=Q_temp(1531:1704);
                run_obj.simulation_results.save_velocity_field(fullfile(folder_root));
                [~,DSi(i,:)]=transport.test(run_obj,folder_root);
            end
%             Q_out=sum(Q,1);
%             DSi_out=sum(DSi,1); 
        end
        
        function [residual_tot,k,d,time_aqui,DPSA_prop]=cost_function(x,file_path)
            if(nargin<2)
                file_path='C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Douffine2.mat';
            end
            load(file_path);
%             f=x(1); k=x(2);
%             Q_out=simulation_set.run_simulation(f,k);
            %f=x(1); k=x(2); d1=x(3); d2=x(4);
%             k=x(1); d1=x(2); d2=x(3);
%             [Q_out,DSi_out]=simulation_set.run_sim_transport(k,d1,d2);
            if(length(x)==2)
                k=x(1); d=x(2);
            else
                k=x(1); d=0;
            end
            [Q_out,time_aqui,DPSA_prop]=simulation_set.run_simulation3(k,d,file_path);%,nanmean(Q_real));
            
            residual=Q_out-Q_real;
%             residual2=DSi_out-DSi_real;
            residual=nansum(residual.^2)/nansum((Q_real-nanmean(Q_real)).^2);
%             residual2=nansum(residual2.^2)/nansum((DSi_real-nanmean(DSi_real)).^2);
            residual_tot=residual;%residual2;%residual+residual2;
%             plot(t_real,Q_out);
            fprintf(strcat('residual:',num2str(residual_tot),'\t k:',num2str(k),'\n'));
        end
        
        function [residual_tot,k,d]=cost_function2(x,file_path)
            if(nargin<2)
                load('C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Ris.mat');
            else
                load(file_path);
            end
            k=x(1); d=x(2); f=x(3);
            Q_out=simulation_set.run_simulation3(k,d,f);
            
            
            residual=Q_out-Q_real;
%             residual2=DSi_out-DSi_real;
            residual=nansum(residual.^2)/nansum((Q_real-nanmean(Q_real)).^2);
%             residual2=nansum(residual2.^2)/nansum((DSi_real-nanmean(DSi_real)).^2);
            residual_tot=residual;%residual2;%residual+residual2;
            plot(t_real,Q_out);
        end
        
        function [res,k,d]=systematic_test
%             file_path={'C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Douffine2.mat','C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Guillec2.mat',...
%                 'C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Dourduff.mat','C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Penze.mat',...
%                 'C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Dossen.mat','C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Ris.mat'};

            file_path={'C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Ris.mat','C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Douffine2.mat'};
            
%             f_init=[0.1,0.2,0.3,0.4];
            k_init=[0.0036,0.0072,0.01,0.036,0.072,0.1,0.36,0.72,1,5,10];
            d_init=[0,0.05,0.1,0.2,0.3,0.4,0.6,0.8];%10,20,30,40,60,80,100];
            x=allcomb(k_init,d_init);
            size_x=size(x);
            res=nan(size_x(1),length(file_path));
            k=nan(size_x(1),length(file_path));
            d=nan(size_x(1),length(file_path));
%             f=nan(size_x(1),length(file_path));
            time_aqui=nan(size_x(1),length(file_path));
            DPSA_prop=nan(size_x(1),length(file_path));
            for j=1:2 %1:length(file_path)
                load(file_path{j});
                figure; hold on
                for i=1:size_x(1)
                    [res(i,j),k(i,j),d(i,j),time_aqui(i,j),DPSA_prop(i,j)]=simulation_set.cost_function(x(i,:),file_path{j});
                    %                 [res(i),k(i),d(i)]=simulation_set.cost_function2(x(i,:));
                    k_str=k(i,j); d_str=d(i,j); time_aqui_str=time_aqui(i,j); DPSA_prop_str=DPSA_prop(i,j); res_str=res(i,j);
                    fprintf(strcat(num2str(j),'\t',num2str(i),'\t k:',num2str(k_str),'\t alpha:',num2str(d_str),'\t residual:',num2str(res_str),'\t time_aqui:',num2str(time_aqui_str),'\t DPSA_prop:',num2str(DPSA_prop_str),'\t','\n'));
                end
                if(j==1) 
                    j2=6; 
                elseif(j==2) 
                    j2=1;
                end
                save(strcat('C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\temp2_',num2str(j2),'.mat'),'-v7.3');
                plot(t_real,Q_real,'kx-');
            end
        end
        
        function [k_opt,d_opt,residual]=optimize_simu
%             cost_func=@(x) simulation_set.cost_function(x);
%             f_init=0.3;
%             k_init=0.07;
%             x_init=[f_init,k_init];
%             x_bounds_min=[0.01,0.01];
%             x_bounds_max=[0.4,3];
            cost_func=@(x) simulation_set.cost_function(x);
            f_init=0.2;
            d1_init=40;
            d2_init=50;
            k_init=0.1;
            x_init=[k_init,d1_init];%,d2_init];
            x_bounds_min=[0.0036,0.5];%,0.5];
            x_bounds_max=[3.6,50];%,200];
%             x_init=[k_init,d1_init,d2_init];
%             x_bounds_min=[0.01,0.5,0.5];
%             x_bounds_max=[3,25,70];
            
%             
%             options=optimoptions(@lsqnonlin,'Display','iter');
%             [x_opt,residual]=lsqnonlin(cost_func,x_init,x_bounds_min,x_bounds_max,options);
%             options=saoptimset('Display','iter');
%             [x_opt,residual]=simulannealbnd(cost_func,x_init,x_bounds_min,x_bounds_max,options);
            options=optimoptions(@fmincon,'Display','iter','PlotFcn', @optimplotfval,'TolX',1e-3,'TolFun',1e-2);
            problem = createOptimProblem('fmincon','objective',cost_func,'x0',x_init,'lb',x_bounds_min,'ub',x_bounds_max,'options',options);
%             ms = MultiStart;
%             ms.UseParallel=true;
            gs = GlobalSearch;
            gs = GlobalSearch(gs,'TolX',1e-3,'TolFun',1e-2,'MaxTime',32400,'Display','iter');
            gs.NumTrialPoints=400;
            gs.NumStageOnePoints=200;
            tic
%             [x_opt,residual] = run(ms,problem,10);
            [x_opt,residual] = run(gs,problem);
%             [x_opt,residual]=fmincon(cost_func,x_init,[],[],[],[],x_bounds_min,x_bounds_max,[],options);
            toc
%             options=psoptimset('Display','iter','Cache','on','TimeLimit',3600,'PlotFcn', {@psplotbestf,@psplotfuncount,@psplotbestx});
%             [x_opt,residual,exitflag,output_inf]=patternsearch(cost_func,x_init,[],[],[],[],x_bounds_min,x_bounds_max,[],options);
%             options=optimoptions('particleswarm','Display','iter','MaxTime',3600,'PlotFcn', @psplotbestf);
%             [x_opt,residual,exitflag,output_inf]=particleswarm(cost_func,2,x_bounds_min,x_bounds_max,options);
%             f_opt=x_opt(1);
            k_opt=x_opt(1);
            d_opt=x_opt(2);
            
        end
        
        function [k_opt,residual]=optimize_simu2(x_init,file_path)
            cost_func=@(x) simulation_set.cost_function(x,file_path);
%             x_init=0.01;
            x_bounds_min=0.0036;%,0.5];
            x_bounds_max=3.6;%,200];

            options=optimoptions(@lsqnonlin,'Display','iter');
            [x_opt,residual]=lsqnonlin(cost_func,x_init,x_bounds_min,x_bounds_max,options);
%             tic
%             [x_opt,residual]=fmincon(cost_func,x_init,[],[],[],[],x_bounds_min,x_bounds_max,[],options);
%             toc
            k_opt=x_opt(1);
            
        end
        
        function [res,k,d]=systematic_optimization
            file_path={'C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Dourduff.mat','C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Penze.mat',...
                'C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Dossen.mat','C:\Users\Jean\Documents\ProjectDSi\BV_ecoflux\RealData\Ris.mat'};

            
%             f_init=[0.1,0.2,0.3,0.4];
            k_init=[0.1,0.036,0.036,0.036];
            
            k_opt=nan(1,4);
            residual=nan(1,4);
            
            for i=1:length(k_init)
                [k_opt(i),residual(i)]=simulation_set.optimize_simu2(k_init(i),file_path{i});
            end
        end
    end
end