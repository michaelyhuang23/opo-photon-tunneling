using DifferentialEquations
using Plots
using PolynomialRoots
using JSON
using LsqFit


function calc_time(b_val, lambda_val, g_val)
    function parametrize(alpha, beta)
        u = asin(g_val * alpha / sqrt(lambda_val)) + asin(g_val * beta / sqrt(lambda_val))
        v = asin(g_val * alpha / sqrt(lambda_val)) - asin(g_val * beta / sqrt(lambda_val))
        return [u, v]
    end

    function sol1_time()
        function V(u, v)
            term1 = lambda_val * (cos(u) - cos(v)) - 2 * log(cos(u) + cos(v))
            term2 = 2 * g_val * b_val / sqrt(lambda_val) * log(cos(u) + cos(v))
            term3 = -2 * g_val * b_val / sqrt(lambda_val) * log(2 + 2 * sin((u + v) / 2) + 2 * sin((u - v) / 2) + cos(v) - cos(u))
            return term1 + term2 + term3
        end
        function Vuu(alpha)
            return (lambda_val - g_val^2 * b_val * alpha) / (lambda_val - g_val^2 * alpha^2) - lambda_val + 2 * g_val^2 * alpha^2
        end
        function Vvv(alpha)
            return (lambda_val - g_val^2 * b_val * alpha) / (lambda_val - g_val^2 * alpha^2) + lambda_val
        end
        alpha_roots = roots([-b_val, -(lambda_val - 1), 0, g_val^2])
        alpha_roots = sort(map(x -> real(x), alpha_roots))
        u_roots = map(x -> parametrize(x, x)[1], alpha_roots)
        Vs = V(u_roots[2], 0)
        Va = V(u_roots[1], 0)
        Vb = V(u_roots[3], 0)

        term1 = 2 * pi * exp((Vs - Va) / g_val^2) / sqrt(Vuu(alpha_roots[1]) * Vvv(alpha_roots[1]))
        term2 = 2 * pi * exp((Vs - Vb) / g_val^2) / sqrt(Vuu(alpha_roots[3]) * Vvv(alpha_roots[3]))
        term3 = sqrt(Vvv(alpha_roots[2]) / (-Vuu(alpha_roots[2])))
        term4 = exp((Va - Vb) / g_val^2)
        term5 = sqrt(Vuu(alpha_roots[1]) * Vvv(alpha_roots[1]) / Vuu(alpha_roots[3]) / Vvv(alpha_roots[3]))

        return (term1 + term2) * term3 / (term4 * term5 + 1)
    end
    return sol1_time()
end

function find_b_vals(lambda_val, g_val)
    left = BigFloat(0.0)
    right = BigFloat(4.0 * 10^10)
    possible_b_vals = []
    while right - left > BigFloat(1e-9)
        mid = (left + right) / 2
        try
            time = calc_time(Float64(mid), lambda_val, g_val)
            if time < 400
                right = mid
                push!(possible_b_vals, (Float64(mid), time))
            else
                left = mid
            end
        catch e
            right = mid
        end
        if length(possible_b_vals) > 2
            break
        end
    end
    return possible_b_vals
end

cases = []
lambda_val = 2.0
b_val = 0.0
for g_val in range(0.4, 1.0, 20)
    try
        time = calc_time(b_val, lambda_val, g_val)
        if time < 400
            push!(cases, Dict("g" => g_val, "lambda" => lambda_val, "b" => b_val, "analytic_time" => time))
            println("possible params: g: $g_val, lambda: $lambda_val, b: $b_val, time: $time")
        end
    catch e
        right = mid
    end
end


function sim_time(b_val, lambda_val, g_val, t_sim_end)
    function drift!(dalpha, alpha, p, t)
        dalpha[1] = b_val - alpha[1] + alpha[2] * (lambda_val - g_val^2 * alpha[1]^2)
        dalpha[2] = b_val - alpha[2] + alpha[1] * (lambda_val - g_val^2 * alpha[2]^2)
    end

    function noise!(dalpha, alpha, p, t)
        dalpha[1, 1] = sqrt(lambda_val - g_val^2 * alpha[1]^2)
        dalpha[2, 2] = sqrt(lambda_val - g_val^2 * alpha[2]^2)
        dalpha[1, 2] = 0.0
        dalpha[2, 1] = 0.0
    end

    t_span = (0.0, t_sim_end)
    alpha_0 = [Complex(-sqrt(lambda_val - 1) / g_val, 0.0), Complex(-sqrt(lambda_val - 1) / g_val, 0.0)]
    alpha_f = [Complex(sqrt(lambda_val - 1) / g_val, 0.0), Complex(sqrt(lambda_val - 1) / g_val, 0.0)]
    W = RealWienerProcess!(0.0, zeros(2), zeros(2))
    problem = SDEProblem(drift!, noise!, alpha_0, t_span, noise_rate_prototype=complex(zeros(2, 2)), noise=W)
    sol = solve(problem)
    ensamble_prob = EnsembleProblem(problem)

    # function condition(alpha, t, integrator)
    #     return real(alpha[1] + conj(alpha[2]))/2 > 0.0
    # end

    # function affect!(integrator)
    #     terminate!(integrator)
    # end

    # cb = ContinuousCallback(condition, affect!)
    sol = solve(ensamble_prob, EM(), dt=0.003, trajectories=1000, EnsembleThreads(), adaptive=false)
    successful_trajs = [traj for traj in sol if Symbol(traj.retcode) == :Success || Symbol(traj.retcode) == :Terminated]
    mean_amps = zeros(length(successful_trajs[1].t))
    mean_ts = zeros(length(successful_trajs[1].t))
    for traj in successful_trajs
        amps = map(x -> -real(x[1] + conj(x[2])), traj.u)
        mean_amps .+= amps
        mean_ts .+= traj.t
    end
    mean_amps ./= length(successful_trajs)
    mean_ts ./= length(successful_trajs)
    decay_model(x, p) = p[1] .* exp.(-p[2] .* x)
    fit = curve_fit(decay_model, mean_ts, mean_amps, [1.0, 1.0])
    fit_time = 1 / coef(fit)[2]
    end_time = 1.0#sum(end_times) / length(end_times)
    println("g: $g_val, lambda: $lambda_val, b: $b_val")
    println("success count: $(length(successful_trajs)), end times: $end_time, fit time: $fit_time")
    return end_time, fit_time
end


success_cases = []
for case in cases
    try
        end_time, fit_time = sim_time(case["b"], case["lambda"], case["g"], case["analytic_time"] * 4)
        case["sim_end_time"] = end_time
        case["sim_fit_time"] = fit_time
        push!(success_cases, case)
    catch e
        println(e)
    end
end

println("number of successes: ", size(success_cases))
open("output.json", "w") do file
    JSON.print(file, success_cases)
end
