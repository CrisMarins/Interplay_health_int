using StatsBase, LinearAlgebra, CSV, DataFrames
using BenchmarkTools, Dates, Random
#using PyPlot

module various_functions
using StatsBase, Random
export activity_sampling, initial_states, rnd_activity_reduction, rnd_mask_reduction, overlap_states
export thresh, quello_principale, sample_neighbours, act_thresh, ab_mask_thresh, ab_mask_reductioin,ab_activity_reduction

function activity_sampling(eps, N_samp, α)
    u_vec = rand(N_samp)
    α1 = 1 - α
    a = (((1-eps^α1)*u_vec)  .+  eps^α1).^(1/α1)
    return a
end

function initial_states(N_samp, i0, f, v, vax_vec)
    S, I, R, R_v = N_samp, 0, 0, 0
    susceptible = ones(N_samp)
    recovered = zeros(N_samp)
    rec_vaccine = zeros(N_samp)
    infected = zeros(N_samp)
    idx = zeros(N_samp)#[i for i in 1:N_samp]
    for i in 1:N_samp
        idx[i] = copy(i)
        if vax_vec[i]==1 && rand()<= v
            rec_vaccine[i] = 1
            susceptible[i] = 0
            idx[i] = 0
        end
    end
    idx_S = filter!(x->x ≠ 0, idx )
    #println("index S")
    #println(idx_S)
    if i0<1
        i0 = Int(round(i0*N_samp))
    end
    idx_I = sample(idx_S, i0, replace = false)
    #println("index I")
    #println(idx_I)
    idx_I =copy(Int.(idx_I))
    infected[idx_I] .= 1
    susceptible[idx_I] .= 0
    S = sum(susceptible)
    I = sum(infected)
    R = sum(recovered)
    R_v = sum(rec_vaccine)
    return S, I, R, R_v, susceptible, infected, recovered, rec_vaccine
end

function overlap_states(N_samp, ψ, γ,α_m, adoption, overlap)
    S_red = ones(N_samp)
    I_red = ones(N_samp)
    vax_vec = zeros(N_samp)
    mask_red_s = ones(N_samp)
    mask_red_i = ones(N_samp)
    total_overlap = Int(floor(adoption*overlap*N_samp))
    total_excedence = Int(floor(adoption*N_samp)) - total_overlap
    if total_overlap + 2*total_excedence <= N_samp
        indices = shuffle(1:N_samp)
        overlap_vec = indices[1:total_overlap]
        exc1_vec = indices[(total_overlap+1 ): (total_overlap+total_excedence)]
        exc2_vec = indices[(total_excedence+total_overlap+1):(total_excedence*2+total_overlap)]
        behav1 = copy(overlap_vec)
        behav2 = copy(overlap_vec)
        append!(behav1, exc1_vec)
        append!(behav2, exc2_vec)
        mask_red_i[behav2] .= α_m
        mask_red_s = copy(mask_red_i)
        vax_vec[behav1] .= 1
        #S_red[behav1] .= ψ
        #I_red[behav1] .= γ
    end
    return mask_red_s, mask_red_i, vax_vec
end

#Ricorda \psi è la riduzione dei S, w frazione S
#Ricorda \gamma riduzione I, p frazione I
function rnd_activity_reduction(ψ, w, γ, p, N_samp)
    S_red = ones(N_samp)
    I_red = ones(N_samp)
    for i in 1:N_samp
        if rand()<w
            S_red[i] = copy(ψ)
        end
        if rand()<p
            I_red[i] = copy(γ)
        end
    end
    return S_red, I_red
end

function rnd_mask_reduction(N_samp, r_s, r_i, α_m)
    mask_red_s = ones(N_samp) 
    mask_red_i = ones(N_samp) 
    for i in 1:N_samp
        if rand()<r_s
            mask_red_s[i] =copy(α_m)
        end
        if rand()<r_i
            mask_red_i[i] =copy(α_m)
        end
    end
    return mask_red_s, mask_red_i
end

function act_thresh(eps, α, u) #u is the perecentage 
    if u == 0
        return eps
    elseif u == 1
        return 1
    end

    α1 = 1 - α
    x = (u *(1-eps^α1) + eps^α1).^(1/α1)
    return x
end

function ab_activity_reduction(ψ, w, γ, p, eps, N_samp, α, activity_sam)
    S_red = ones(N_samp)
    I_red = ones(N_samp)
    t_w_u = act_thresh(eps, α, w)
    t_p_u = act_thresh(eps, α, p)
    for i in 1:N_samp
        
        if activity_sam[i] <= t_w_u
            S_red[i] = copy(ψ)
        end


        if activity_sam[i] <= t_p_u
            I_red[i] = copy(γ)
        end
    end
    return S_red, I_red

end

function ab_mask_reductioin(N_samp, r_s, r_i, α_m, eps, α, pure_act, w)
    mask_red_s = ones(N_samp) 
    mask_red_i = ones(N_samp)
    t_r_s = act_thresh(eps, α, w)
    t_r_i = act_thresh(eps, α, w)
    for i in 1:N_samp
        
        if pure_act[i] >= t_r_s
            if rand() <= r_s
                mask_red_s[i]  = copy(α_m)
            end
        end

        if pure_act[i] >= t_r_i
            if rand() <= r_i 
                mask_red_i[i]  = copy(α_m)
            end
        end
    end
    return mask_red_s, mask_red_i
end

function ab_mask_thresh(mu, eps, α, m, r_s, r_i, α_m)
    t_r_s = act_thresh(eps, α, 1 - r_s)
    t_r_i = act_thresh(eps, α, 1 - r_i)
    n_rev_r_s = (1-t_r_s^(1-α))/(1 - eps^(1-α))
    n_rev_r_i = (1-t_r_i^(1-α))/(1 - eps^(1-α))
    mean_a = (1 - α)*(1- eps^(2-α))/((2-α)* (1-eps^(1-α)))
    mean_a2 = (1 - α)*(1- eps^(3-α))/((3-α)* (1-eps^(1-α)))
    mean_a_s_rev = (1 - α)*(1 - t_r_s^(2-α))/((2-α)* (1-eps^(1-α)))
    mean_a_i_rev = (1 - α)*(1 - t_r_i^(2-α))/((2-α)* (1-eps^(1-α)))
    mean2_a_s_rev = (1 - α)*(1 - t_r_s^(3-α))/((3-α)* (1-eps^(1-α)))
    mean2_a_i_rev = (1 - α)*(1 - t_r_i^(3-α))/((3-α)* (1-eps^(1-α)))
    
    n_rev = sort([n_rev_r_i, n_rev_r_s])
    m_a = sort([mean_a_i_rev, mean_a_s_rev])
    m2_a = sort([mean2_a_i_rev, mean2_a_s_rev])

    f_rev_n = 1 + α_m * n_rev[2] + α_m *(α_m -1)* n_rev[1]
    f_rev_a = mean_a + α_m * m_a[1] + α_m *(α_m -1)* m_a[2]
    f_rev_a2 = mean_a2 + α_m * m2_a[1] + α_m *(α_m -1)* m2_a[2]

    λ_t = mu* (m*(f_rev_a + sqrt(f_rev_n * f_rev_a2)))^(-1)
    return λ_t, mean_a
end

function sample_neighbours(index, m, N)
    vec =[1:index-1; index+1:N]
    return sample(vec, m, replace= false)
end

function thresh(mu, eps, α, m, v, f, α_m, r_s, r_i, ψ, p, γ, w)
    α_s = 1- r_s*(1- α_m)
    α_i = 1- r_i*(1- α_m)
    ψ_w = 1 - w*(1- ψ)
    γ_p = 1- p*(1- γ)
    mean_a = (1 - α)*(1- eps^(2-α))/((2-α)* (1-eps^(1-α)))
    mean_a2 = (1 - α)*(1- eps^(3-α))/((3-α)* (1-eps^(1-α)))
    λ_t = 2*mu*((m*(1-v*f)*α_s*α_i)*(mean_a*(ψ_w + γ_p)+ sqrt(mean_a^2 * (ψ_w - γ_p)^2 + 4* ψ_w * γ_p* mean_a2)))^(-1)
    return λ_t, mean_a
end


function simulate(i0, N, mu, m,activities_s , activities_i,mask_s, mask_i, lambda, f ,v, S, I, R, R_v, susceptible, infected, recovered, rec_vaccine)
    #lambda = LinRange(0.0, 1.0, points) lo voglio fare fuori
    
    R_inf = [i0*N]
    R_vec = []
    I_vec = []
    
    for k in lambda[1:end]

        #S, I, R, R_v, susceptible, infected, recovered, rec_vaccine = initial_states(i0,N, f, v)
        compartments = Dict( "S"=> [copy(S)], "I"=>[copy(I)], "R"=>[copy(R)], "R_v"=>[copy(R_v)])
        infected_sim = copy(infected)
        susceptible_sim = copy(susceptible)
        recovered_sim = copy(recovered)
        S_sim = copy(S)
        I_sim = copy(I)
        R_sim = copy(R)
        incidence = []
        spotty_vec =  []
        println("$k")
        pro = 0.5
        #while compartments["I"][end] > 0
        for ti in 0:1000
            new_I, new_R = [], []
            spotty_vec = rand(N)
            for i in 1:N
                if susceptible_sim[i] == 1
                    if rand()<activities_s[i]
                        neighbours = sample_neighbours(i, m, N)
                        for neighbour in neighbours
                            if infected_sim[neighbour]==1
                                if rand()<k*mask_s[i]*mask_i[neighbour]
                                    append!(new_I, copy(i))
                                end
                            end
                        end
                    end

                elseif infected_sim[i] == 1
                    if rand()<activities_i[i]
                        neighbours=sample_neighbours(i,m,N)
                        for neighbour in neighbours
                            if susceptible_sim[neighbour] == 1
                                if rand()<k*mask_s[neighbour]*mask_i[i]
                                    append!(new_I, copy(neighbour))
                                end
                            end
                        end
                    end
                    if rand()< mu
                        append!(new_R, copy(i))
                    end

                end

            end
            for j in new_I
                if infected_sim[j] == 0
                    infected_sim[j]= 1
                    susceptible_sim[j] = 0
                    S_sim = copy(S_sim)-1
                    I_sim = copy(I_sim) +1
                end
            end
            for j in new_R
                if infected_sim[j]==1
                    infected_sim[j] = 0
                    I_sim = copy(I_sim)-1
                    R_sim = copy(R_sim)+1
                end
            end
            append!(incidence, length(new_I))
            append!(compartments["S"], copy(S_sim))
            append!(compartments["I"], copy(I_sim))
            append!(compartments["R"], copy(R_sim))
        end
        append!(R_inf, copy(compartments["R"][end]))
        #return incidence, compartments
        #append!(R_vecs, compartments["R"])
        R_vec = copy(compartments["R"])
        I_vec = copy(compartments["I"])
        #println(typeof(R_vec[1]))
        #println(length(R_vec[1]))
        #println(R_vec)
        return R_vec, I_vec
    end
    
    #println("ueue")
end

function quello_principale(sim, eps, N_samp, α, i0, f, v, ψ, w, γ, p, r_s, r_i, α_m, mu, m, lambda, adoption,overlap )
    RR = zeros(sim, length(lambda))
    RR[1,:] = copy(lambda)
    R_vecs = Vector{Vector{Float64}}()
    I_vecs = Vector{Vector{Float64}}()
    mask_s = ones(N_samp)
    mask_i= ones(N_samp)
    S_red = ones(N_samp)
    I_red = ones(N_samp)
    #overlap = 0
    #adoption = 0.2
    for n in 1:sim
        println("overlap = $overlap , adoption = $adoption , mask")
        println("sim numero n $n")
        pure_act = activity_sampling(eps, N_samp, α)
        #S_red, I_red = ab_activity_reduction(ψ, w, γ, p, eps, N_samp, α, pure_act)
        #S_red, I_red = rnd_activity_reduction(ψ, w, γ, p, N_samp)
        mask_s, mask_i, vax_vec = overlap_states(N_samp, ψ, γ,α_m, adoption, overlap)
        #S_red, I_red , vax_vec = overlap_states(N_samp, ψ, γ,α_m, adoption, overlap)
        S, I, R, R_v, susceptible, infected, recovered, rec_vaccine = initial_states(N_samp, i0, f, v, vax_vec)
        activities_s = pure_act .* S_red 
        activities_i =  pure_act .* I_red
        #mask_s, mask_i = ab_mask_reductioin(N_samp, r_s, r_i, α_m, eps, α, pure_act, w)
        R_vec, I_vec = simulate(i0, N_samp, mu, m, activities_s , activities_i ,mask_s, mask_i, lambda, f, v, S, I, R, R_v, susceptible, infected, recovered, rec_vaccine)
        #println(typeof(R_vec))
        push!(R_vecs, R_vec)
        #println(typeof(R_vecs))
        push!(I_vecs, I_vec)
    end
    return R_vecs, I_vecs
end
end

using .various_functions

#df = CSV.read("/home/cmarinelli/Simulations/parameters.csv", DataFrame, decimal = ',')hai fatto finpo a 6
#hai fatto finpo a 6
df = CSV.read("/home/cmarinelli/Simulations/Parameters/parameters_casestudy_mask1.csv", DataFrame, decimal = ',')
println(df)
#tempo_inizio = 0
param = Dict(pairs(eachcol(df)))
nite = String[]
for key in keys(param)@eval  const 
    $key = param[$(QuoteNode(key))][1]
    val = param[key][1]
    str_key = "$(key) ="
    str_val = "$(val) "
    stro = str_key * str_val
    push!(nite, stro)
end
nute = join(nite)
a = typeof(thresh)
println(a)
#λ_t, media_a = ab_mask_thresh(mu, eps, α, m, r_s, r_i, α_m)
#lambda = round.([LinRange(0, 0.7*λ_t, 4)[1:end-1]; LinRange(0.7* λ_t, 1.3*λ_t,7); LinRange(1.3*λ_t,1,4)[2:end]], digits=3)
lambda = [8]/10
#cucu = fill("#", length(lambda))
#cucu[1] = nute
vi = 0.7
adoption = 0
#overlaps = [ 0.5, 0.75, 1]
#overlaps = [0, 0.25, 0.5, 0.75, 1]
#overlaps = [ 0.75, 1]
overlaps = 1
#R0 = copy(lambda) ./ λ_t
println(lambda)
#println(R0)
b = ARGS[1]
for overlap in overlaps
    tempo_inizio = now()
    R_vecs, I_vecs = quello_principale(sim, eps, N_samp, α, i0, f, vi, ψ, w, γ, p, r_s, r_i, α_m, mu, m, lambda, adoption, overlap)
    #println(R_vecs)
    #println(typeof(R_vecs))
    #RRdf = DataFrame(RR,Symbol.(R0))
    az = "_cs_$(b)_overlap_$(Int(floor(adoption*100)))_$(Int(floor(overlap*100)))"
    filename_r = "datafile_r$az"
    filename_i = "datafile_i$az"

    #println(filename)
    tempo_fine = now()
    tempo = tempo_fine - tempo_inizio
    df[!, "tempo"] .= tempo 


    dfr = DataFrame(R_vecs, :auto)
    dfi = DataFrame(I_vecs, :auto)

    CSV.write("/home/cmarinelli/Simulations/data_simulation/overlap_data/$filename_i.csv", dfi, decimal = '.')
    CSV.write("/home/cmarinelli/Simulations/data_simulation/overlap_data/$filename_r.csv", dfr, decimal = '.')

    CSV.write("/home/cmarinelli/Simulations/data_simulation/overlap_data/param_$az.csv", df, decimal = '.')
end