
% add hs1D files to the directory
addpath(genpath('MyDirectory/hs1D_beta'));


%% Topographic analysis
fprintf('------------ Start the Topographic Analysis ------------');
DEM_filepath='MNT_PF_5m.tif'; % DEM file path to be analyzed
outlet_coordinates=[]; % coordinates of the outlet from where the cathcment will be extracted
critical_drainage_area=4000; % minimum area to initiate a stream [in number of pixels] 

% get the watershed object extracted from topographic analysis and the
% equivalent hillslope (1D summary of the catchment with two functions: the
% width and the aquifer depth function). See Troch et al. 2003, Bogaart et al. 2006
[watershed_topo,hillslope_equiv]=watershed.test(DEM_filepath,outlet_coordinates,critical_drainage_area);

% Some figures
% Catchment map
figure; hold on
imageschs(watershed_topo.DEM);
colorbar;
plot(watershed_topo.S.x,watershed_topo.S.y,'w.')
plot(watershed_topo.outlet(1),watershed_topo.outlet(2),'ro','Marker')
% equivalent hillslope plot
figure; hold on
subplot(1,2,1)
plot(hillslope_equiv.hsB.x,hillslope_equiv.hsB.w) 
subplot(1,2,2)
plot(hillslope_equiv.hsB.x,hillslope_equiv.hsB.z) 


%% Groundwater flow simulation
fprintf('------------ Start the Groundwater Flow Simulation ------------');
k=1e-5; % hydraulic conductivity [m/s]
phi_d=0.05;% drainable porosity [-]
Nx=100; % number of discretized elements
river_depth=2; % thickness of the aquifer at the river
aquifer_type='flat'; % two options: flat aquifer or sloped aquifer
% run the model for a given set of value. S is the groundwater storage
% [m^2], Q is the Boussinesq groundwater flow [m^3/s] and t is the time
% time series of length Nt [s]. S and Q are summarized as matrixes of size Nx*Nt.
[run_obj,t,S,Q]=simulation_set.run_simulation(phi_d,k*3600,river_depth,aquifer_type,Nx);
[DPSA_,RF_]=compute_DPSA_RF(run_obj.simulation_results,run_obj.boussinesq_simulation);

% Some figures


%% Groundwater transport simulation
phi_tot=0.2; % total porosity
[transp,~,t_out,transit_times,travel_distances,weights_,t_inj,~,DPSA_part,GW_part]=transport_2D_par_temp.transport_main(run_obj,phi_tot);
t_min=0.1; % in yrs
t_max=12; % in yrs
time_support=[0,logspace(log10(t_min/10),log10(t_max*10),1000)]*24*365*3600;%
ttds_=ttds.retrieve_ttds(t_out,transit_times,weights_,travel_distances,x_fin,time_support);

% Some figures
