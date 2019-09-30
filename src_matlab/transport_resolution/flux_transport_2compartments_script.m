%% Flux
f_deep=0.05;
k_deep=0.1;
f_soil=0.2;
k_soil=3;
% % run the hs1D flux simu
run_deep=simulation_set.run_simulation_rooting(k_deep,30/3600,'H:\Users\marcais\ProjectDSi\RealData\Guillec2.mat',f_deep);
run_deep_daily=simulation_set.run_simulation_rooting(k_deep,30/3600,'H:\Users\marcais\ProjectDSi\RealData\Guillec2.mat',f_deep);
run_soil=run_deep_daily.run_simulation_from_struct(x,f_soil,k_soil,w,slope_angle,2*ones(size(x)),run_deep_daily);

%% Transport
% GW transport
time_support=[linspace(0.25,95,95*4+1)*3600,linspace(4,364,361)*24*3600,linspace(1,100,100)*24*3600*365];
% run the particle tracking scheme transport for the deep part
[obj_deep,t_out_groundwater,transit_times_groundwater,distance_groundwater,weights_groundwater]=transport_2D_par.transport_main(run_deep);
% delete particles that directly fall on saturated areas
[t_out_groundwater1,transit_times_groundwater1,distance_groundwater1,weights_groundwater1]=delete_DPSA_particles(obj_deep,t_out_groundwater,transit_times_groundwater,distance_groundwater,weights_groundwater);
% compute the ttds for the deep part
ttds_deep_monthly=ttds.retrieve_ttds(t_out_groundwater1,transit_times_groundwater1,weights_groundwater1,distance_groundwater1,time_support);

% Soil transport
% run the particle tracking scheme transport for the shallow part
[obj_soil,t_out_soil,transit_times_soil,distance_soil,weights_soil]=transport_2D_par.transport_main(run_soil);
% compute the ttds for the shallow part
ttds_soil=ttds.retrieve_ttds(t_out_soil,transit_times_soil,weights_soil,distance_soil,time_support);

% merge different compartments
[~,~,~,RF_spat_deep]=compute_DPSA_RF(run_deep_daily.simulation_results,run_deep_daily.boussinesq_simulation);
Q_out_deep=sum(RF_spat_deep);
Q_out_deep(Q_out_deep<0)=0;

[DPSA_shallow,RF_shallow]=compute_DPSA_RF(run_soil.simulation_results,run_soil.boussinesq_simulation);
Q_out_shallow=DPSA_shallow+RF_shallow;
Q_out_shallow(Q_out_shallow<0)=0;

if(isequal(ttds_deep_monthly.sampling_time,ttds_soil.sampling_time))
    ttds_total=ttds;
    ttds_total=ttds_total.instantiate_ttds(time_support);
    ttds_total=ttds_total.merge_two_ttds_simulations(ttds_deep_monthly,ttds_soil,Q_out_deep(1:end-1)',Q_out_shallow(1:end-1)');
else
    % reinterpolation first
    ttds_deep_daily=ttds;
    ttds_deep_daily=ttds_deep_daily.instantiate_ttds(time_support);
    ttds_deep_daily=ttds_deep_daily.reinterpolate(ttds_deep_monthly,ttds_soil.sampling_time);
    ttds_deep_daily=ttds_deep_daily.compute_ttds_moments;
    ttds_deep_daily=ttds_deep_daily.compute_youngwaterfraction;
    % computation of the total ttds
    ttds_total=ttds;
    ttds_total=ttds_total.instantiate_ttds(time_support);
    ttds_total=ttds_total.merge_two_ttds_simulations(ttds_deep_daily,ttds_soil,Q_out_deep(1:end-1)',Q_out_shallow(1:end-1)');
end



function [t_out_groundwater1,transit_times_groundwater1,distance_groundwater1,weights_groundwater1]=delete_DPSA_particles...
                                        (obj_deep,t_out_groundwater,transit_times_groundwater,distance_groundwater,weights_groundwater)
    transit_times_groundwater1=transit_times_groundwater(~isnan(t_out_groundwater));
    weights_groundwater1=weights_groundwater(~isnan(t_out_groundwater));
    distance_groundwater1=distance_groundwater(~isnan(t_out_groundwater));
    DPSA_deep=[obj_deep.DPSA;obj_deep.DPSA(obj_deep.DPSA>0 & obj_deep.DPSA<1)];
    DPSA_deep=DPSA_deep(~isnan(t_out_groundwater));
    t_out_groundwater1=t_out_groundwater(~isnan(t_out_groundwater));

    transit_times_groundwater1=transit_times_groundwater1(DPSA_deep==0);
    weights_groundwater1=weights_groundwater1(DPSA_deep==0);
    distance_groundwater1=distance_groundwater1(DPSA_deep==0);
    t_out_groundwater1=t_out_groundwater1(DPSA_deep==0);
    DPSA_deep=DPSA_deep(DPSA_deep==0);
end


