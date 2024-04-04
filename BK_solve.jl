using LinearAlgebra
function BK_solve(A, B, en_index, ex_index, c_index, s_index)

    dyn_index = vcat(en_index, ex_index, c_index);
    n_c = length(c_index);
    n_en = length(en_index);
    n_ex = length(ex_index);
    n_s = length(s_index);
    n_dyn = length(dyn_index);
    n_total = n_dyn + n_s;
    println("Linear difference-equations system with:")
    println(string(n_c)*" forward-looking variable(s).")
    println(string(n_en)*" endogenous backward-looking variable(s).")
    println(string(n_ex)*" exogenous backward-looking variable(s).")
    println(string(n_s)*" endogenous static variable(s).")

    A_DD = A[dyn_index, dyn_index];
    A_DS = A[dyn_index, s_index];
    A_SD = A[s_index, dyn_index];
    A_SS = A[s_index, s_index];
    B_DD = B[dyn_index, dyn_index];
    B_DS = B[dyn_index, s_index];
    B_SD = B[s_index, dyn_index];
    B_SS = B[s_index, s_index];

    static_matrix = inv(I(n_s) - A_SS)*A_SD;
    A_mod = A_DD + A_DS*static_matrix;
    B_mod = B_DD + B_DS*static_matrix;
    W_mod = inv(I(n_dyn) - A_mod)*B_mod;

    Λ = eigvals(W_mod);
    abs_Λ = abs.(Λ);
    unstable_Λ = abs_Λ[abs_Λ .> 1];
    stable_Λ = abs_Λ[abs_Λ .<= 1];
    n_unstable_eigs = length(unstable_Λ);
    n_stable_eigs = length(stable_Λ);
    n_all_eigs = length(Λ);
    if n_unstable_eigs > n_c
        println("Blanchard-Kahn conditions not satisfied: there are MORE ("*string(n_unstable_eigs)*") eigen-values outside the unit-circle than forward looking variables.")
        println("The system has no stable solution.")
    elseif n_unstable_eigs < n_c
        println("Blanchard-Kahn conditions not satisfied: there are LESS ("*string(n_unstable_eigs)*") eigen-values outside the unit-circle than forward looking variables.")
        println("Infinitely many solutions exist.")
    else
        println("Blanchard-Kahn conditions satisfied: there are as many ("*string(n_unstable_eigs)*") eigen-values outside the unit-circle as forward looking variables.")
        println("A unique and stable solution exists.")
    end

    Λ_u = Diagonal(unstable_Λ);
    Λ_s = Diagonal(stable_Λ);
    Λ_matrix = Diagonal(Λ);

    J = eigvecs(W_mod);
    Jinv = inv(J);
    Jinv_ss = Jinv[collect(1:n_stable_eigs), collect(1:n_stable_eigs)];
    Jinv_us = Jinv[collect(n_stable_eigs+1:n_all_eigs), collect(1:n_stable_eigs)];
    Jinv_su = Jinv[collect(1:n_stable_eigs), collect(n_stable_eigs+1:n_all_eigs)];
    Jinv_uu = Jinv[collect(n_stable_eigs+1:n_all_eigs), collect(n_stable_eigs+1:n_all_eigs)];

    A_ss = W_mod[collect(1:n_stable_eigs), collect(1:n_stable_eigs)];
    A_us = W_mod[collect(n_stable_eigs+1:n_all_eigs), collect(1:n_stable_eigs)];
    A_su = W_mod[collect(1:n_stable_eigs), collect(n_stable_eigs+1:n_all_eigs)];
    A_uu = W_mod[collect(n_stable_eigs+1:n_all_eigs), collect(n_stable_eigs+1:n_all_eigs)];

    policy_matrix = -(Jinv_uu^(-1))*Jinv_us;
    transition_matrix = A_ss + A_su*policy_matrix;
    static_matrix = inv(I(n_s)-A[s_index, s_index])*A[s_index, dyn_index];
    return policy_matrix, transition_matrix, static_matrix

end