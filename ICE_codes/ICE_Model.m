clear variables
% High contrast colors
% First triad
teal   = '#45DEC7';
pink   = '#DE459E';
yellow = '#DED645';
% Second triad
orange = '#DE7045';
green  = '#7FDE45';
purple = '#8545DE';
figpath = "/home/cgarcia/Documents/ICE_Modelling/ICE_doc/figures/";

% -------------------------------------------------------------------------
% Modeling the internal combustion engine from dynamo tests
% -------------------------------------------------------------------------

% 1. import the complete dataset for the first figure
ice_data_table = ICE_import( ...
	"ICE_data_2013_11_18@15_22_02.csv", ...
	[2 44808]);%3402 % 20s 4002
% 1.1 change the torque to Nm
% ice_data_table.Torque = ice_data_table.Torque./10;

% 2. plot the data, u: absorber RPM. y: engine RPM & Torque
data_fig = figure(1);
clf
colororder({teal, orange})
data_tile = tiledlayout("vertical","TileSpacing","tight","Padding","tight");
nexttile
hold on
yyaxis left
plot(ice_data_table.BoardTime,ice_data_table.EngineRPM_A,"Color",teal)
ylabel("Engine [RPM]")
yyaxis right
plot(ice_data_table.BoardTime,ice_data_table.Torque,"Color",orange)
ylabel("Torque [Nm]","Interpreter","latex")
xlabel("time [s]", "Interpreter","latex")
nexttile
plot(ice_data_table.BoardTime,ice_data_table.AbsorberRPM_C,"k")
xlabel("time [s]", "Interpreter","latex")
ylabel("Absorber [RPM]")
set(gcf,"PaperPosition",[0,0,10,13])
saveas(data_fig,strcat(figpath, "data_fig.fig"))
saveas(data_fig,strcat(figpath, "data_fig.eps"), "epsc")


% 3. Create the structure array of normalized data for the training of 
% the algorithm resampled.s
%%
% 3.1 Create the structure for the algorithm. Resampled
resample_period = 10;
resample_index = 1:resample_period:size(ice_data_table.BoardTime,1);%
edmd_data = struct(...
	'y', [ice_data_table.EngineRPM_A(resample_index), ice_data_table.Torque(resample_index)], ...
	'u', ice_data_table.AbsorberRPM_C(resample_index), ...
	't', ice_data_table.BoardTime(resample_index));
% 3.1 Create the pqEDMD wrapper object 
ice_pqEDMD = pqEDMDm( ...
	p = [1 2 3], ...
	q = [0.5 0.7 1 1.5 2], ...
	observable = @legendreObservable, ...
	dyn_dcp = @(dec_obs, dec_data)sidDecomposition(2, 3, dec_obs, dec_data));

% % 3.2 Fit the data. It returns an array of decompostions
dcps = ice_pqEDMD.fit(edmd_data);
% 
% %%
% % 3.3 Calculate the error of all the decompositions
err = arrayfun(@(dcpi)dcpi.abs_error(edmd_data), dcps);
% %%
% % 3.4 Extract the best decompostion
[er_best, in_best] = min(err);
dcp = dcps(in_best);

%%
% 3.5 Calculate the approximation of the dynamics
appx = dcp.pred_from_test(edmd_data);

% 3.6 Plot the approximation against the data
appx_fig = figure(2);
clf
colororder({teal, orange})
appx_tile = tiledlayout("vertical","TileSpacing","tight","Padding","tight");
nexttile
hold on 
yyaxis left
plot(edmd_data.t, edmd_data.y(:,1),"Color",teal,LineWidth=1.5)
plot(edmd_data.t, appx.y(:,1),'-.',"Color",pink,LineWidth=1.5)
ylabel("Engine [RPM]","Interpreter","latex")
yyaxis right
plot(edmd_data.t, edmd_data.y(:,2),'Color',orange,LineWidth=1.5)
plot(edmd_data.t, appx.y(:,2),'-.', "Color",green,LineWidth=1.5)
ylabel("Torque [Nm]","Interpreter","latex")
legend('RPM', 'RPM approximation','Torque','Torque approximation',"Location",'northoutside',"Interpreter","latex")
nexttile
plot(edmd_data.t, edmd_data.u,"k",LineWidth=1.5)
ylabel("Absorber [RPM]", "Interpreter","latex")
xlabel("time [s]")
set(gcf,"PaperPosition",[0,0,10,13])
saveas(appx_fig,strcat(figpath, "appx_fig.fig"))
saveas(appx_fig,strcat(figpath, "appx_fig.eps"), "epsc")
%%
ratio_fig = figure(3);
clf
tiledlayout("vertical","TileSpacing","tight","Padding","tight")
nexttile
plot(edmd_data.y(:,1),edmd_data.u,"Color",teal,LineWidth=1.5)
xlabel("Engine RPM", "Interpreter","latex")
ylabel("Absorber RPM", "Interpreter","latex")
nexttile
% 3.1 Calculate the gearbox ratio
ratio = edmd_data.y(:,1)./edmd_data.u;
plot(edmd_data.t,ratio,"Color",pink,LineWidth=1.5)
ylabel("Engine RPM/Absorber RPM","Interpreter","latex")
xlabel("time [s]","Interpreter","latex")
set(gcf,"PaperPosition",[0,0,10,13])
saveas(ratio_fig,strcat(figpath, "ratio_fig.fig"))
saveas(ratio_fig,strcat(figpath, "ratio_fig.eps"), "epsc")