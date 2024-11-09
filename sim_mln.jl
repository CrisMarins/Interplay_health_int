using StatsBase, LinearAlgebra, CSV, DataFrames
using BenchmarkTools, Dates
#using PyPlot

module various_functions
using StatsBase
export activity_sampling, initial_states, rnd_activity_reduction, rnd_mask_reduction
export thresh, quello_principale, sample_neighbours

function activity_sampling(eps, N_samp, α)
    u_vec = rand(N_samp)
    α1 = 1 - α
    a = (((1-eps^α1)*u_vec)  .+  eps^α1).^(1/α1)
    return a
end

function initial_states(N_samp, i0, f, v)
    S, I, R, R_v = N_samp, 0, 0, 0
    susceptible = ones(N_samp)
    recovered = zeros(N_samp)
    rec_vaccine = zeros(N_samp)
    infected = zeros(N_samp)
    idx = zeros(N_samp)#[i for i in 1:N_samp]
    for i in 1:N_samp
        idx[i] = copy(i)
        if rand()<=f && rand()<= v
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


function simulate(i0, N, mu, m, activities_s, activities_i, mask_s, mask_i, lambda, f ,v, S, I, R, R_v, susceptible, infected, recovered, rec_vaccine)
    #lambda = LinRange(0.0, 1.0, points) lo voglio fare fuori
    
    R_inf = [i0*N]
    
    for k in lambda[2:end]

        #S, I, R, R_v, susceptible, infected, recovered, rec_vaccine = initial_states(i0,N, f, v)
        compartments = Dict( "S"=> [copy(S)], "I"=>[copy(I)], "R"=>[copy(R)], "R_v"=>[copy(R_v)], "N" =>[N] )
        infected_sim = copy(infected)
        susceptible_sim = copy(susceptible)
        recovered_sim = copy(recovered)
        S_sim = copy(S)
        I_sim = copy(I)
        R_sim = copy(R)
        incidence = []
        N_sim = copy(S) + copy(I) + copy(R) + copy(R_v)
        println("$k")
        while compartments["I"][end] > 0
            new_I, new_R = [], []

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
            N_sim = copy(S_sim) + copy(I_sim) + copy(R_sim) + copy(R_v)
            if N_sim != N
                println("aiuto errore")
            end
            append!(incidence, length(new_I))
            append!(compartments["S"], copy(S_sim))
            append!(compartments["I"], copy(I_sim))
            append!(compartments["R"], copy(R_sim))
            append!(compartments["N"], copy(N_sim))
        end
        append!(R_inf, copy(compartments["R"][end]))
        #return incidence, compartments
    end
    return R_inf
    #println("ueue")
end

function quello_principale(sim, eps, N_samp, α, i0, f, v, ψ, w, γ, p, r_s, r_i, α_m, mu, m, lambda)
    RR = zeros(sim+1, length(lambda))
    RR[1,:] = copy(lambda)
    for n in 1:sim
        println("sim numero n $n")
        pure_act = activity_sampling(eps, N_samp, α)
        S_red, I_red =rnd_activity_reduction(ψ, w, γ, p, N_samp)
        S, I, R, R_v, susceptible, infected, recovered, rec_vaccine = initial_states(N_samp, i0, f, v)
        activities_s = pure_act .* S_red 
        activities_i =  pure_act .* I_red
        mask_s, mask_i = rnd_mask_reduction(N_samp, r_s, r_i, α_m)
        RR[n+1, :]=simulate(i0, N_samp, mu, m, activities_s, activities_i ,mask_s, mask_i, lambda, f, v, S, I, R, R_v, susceptible, infected, recovered, rec_vaccine)
    end
    return RR
end
end

using .various_functions

df = CSV.read("/home/cmarinelli/Simulations/parameters.csv", DataFrame, decimal = ',')
println(df)
tempo_inizio = now()
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
λ_t, media_a = thresh(mu, eps, α, m, v, f, α_m, r_s, r_i, ψ, p, γ, w)
#lambda = round.([LinRange(0, 0.7*λ_t, 4)[1:end-1]; LinRange(0.7* λ_t, 1.3*λ_t,7); LinRange(1.3*λ_t,1,4)[2:end]], digits=3)
lambda = vec([0 1])
cucu = fill("#", length(lambda))
cucu[1] = nute
R0 = copy(lambda) ./ λ_t
println(R0)
RR = quello_principale(sim, eps, N_samp, α, i0, f, v, ψ, w, γ, p, r_s, r_i, α_m, mu, m, lambda)

RRdf = DataFrame(RR,Symbol.(R0))
a = ARGS[1]
filename = "datafile$a"
#println(filename)
tempo_fine = now()
tempo = tempo_fine - tempo_inizio
df[!, "tempo"] .= tempo 
CSV.write("/home/cmarinelli/Simulations/data_simulation/$filename.csv", RRdf, decimal = '.')
CSV.write("/home/cmarinelli/Simulations/data_simulation/param_$a.csv", df, decimal = '.')
