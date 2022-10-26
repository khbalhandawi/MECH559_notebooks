function Param = bb_Param
    
    % parameters = [8,-4,-3,-3,-3,3,2,0.8]; % <---- with these parameters, the constraints are inactive
    % parameters = [4,-4,-3,-3,-3,3,2,0.8]; % <---- with these parameters, the constraints are active
    % load surrogate_model.mat dmodel
    % 
    % save_data = true;
    % use_surrogate = false;
    % use_gradients = false;
    % Param = {save_data,use_surrogate,use_gradients,parameters,dmodel}; 

    global P
    Param = P;

return