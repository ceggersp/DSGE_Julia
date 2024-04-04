function IRF_compute(transition_matrix, policy_matrix, static_matrix, n_var, T, size_of_the_shock = 1);

    n_states = length(transition_matrix[:,1]);
    n_control = length(policy_matrix[:,1]);
    n_static = length(static_matrix[:,1]);
    impulse = zeros(n_states, 1);
    impulse[n_var,1] = 1
    response_matrix = zeros(1, n_control + n_states + n_static);
    for t = 0:T
        local transition_resp = (transition_matrix^t)*impulse;
        local policy_resp = policy_matrix*transition_resp;
        local static_resp = static_matrix*vcat(transition_resp, policy_resp);
        local response = hcat(transition_resp', policy_resp');
        local response = hcat(response, static_resp');
        response_matrix = vcat(response_matrix, response);
    end
    
    return size_of_the_shock*response_matrix[2:end,:]

end