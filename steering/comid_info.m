function info = comid_info(comid)

    % area of each subarea % DK: We don't use this
    subarea = {comid}';
    %area_km2 = [3469 1291 3570 2423 3424 3240 26099]';

    % weighting to use when spatially averaging the low frequency component of precipitation
    % (ie weight by area but ignore region G - see paper, section 3.5)
    %weighting_lowfreq = [3469 1291 3570 2423 3424 3240 0]';
	weighting_lowfreq = [1]';
    %SubAreaDetails = table(subarea, area_km2, weighting_lowfreq);
	SubAreaDetails = table(subarea, weighting_lowfreq);
	% DK: single weighting as each comid is a single subarea by itself.  


	% DK: Rainfall-runoff model related data are wiped out.


    % other parameters
    %pars.NumSubareas = size(RepCatchDetails, 1);
    pars.NumSubareas = 1; % DK: 1 comid = 1 subarea, always.
    %pars.WaterYearStart_clim = 'January'; % basis of water years for climate ('January' means calendar years are used)
	pars.WaterYearStart_clim = 'October'; % DK: California...
    %pars.WaterYearStart_flow = 'March';   % basis of water years for streamflow
	pars.WaterYearStart_flow = 'October'; % DK: California...
    %pars.StochRepLen_yrs  =  3000;        % length of stochastic replicates, in years
    pars.StochRepLen_yrs  =  1000;        % DK: According to Fowler et al (2022), it is ideal to have at least 1000 years? I don't remember that well.
    pars.Streamflow_BoxCox_Lambda = 0.79; % used during perturbation of the rainfall-runoff relationship (Section 2.4.5 / Table 1 / Section 3.7 / Supp Mat Section S7)

    % store all the above in structure 'info'
    %info = struct('WapabaParSets', WapabaParSets, 'SubAreaDetails', SubAreaDetails, 'RepCatchDetails', RepCatchDetails, 'FlowConversionFactors', FlowConversionFactors, 'pars', pars);
	info = struct('SubAreaDetails', SubAreaDetails, 'pars', pars);
    info.SubareaList = subarea;

end
