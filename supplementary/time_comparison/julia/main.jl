using HTTP: request

"https://gist.githubusercontent.com/vankesteren/96207abcd16ecd01a2491bcbec12c73f/raw/1b59af6962a1107db5873eba59054acc3f9a8aac/Adamopt.jl" |>
  url -> request("GET", url) |>
  res -> String(res.body) |>
  str -> include_string(Main, str)

using .Adamopt

using SpecialFunctions
using Statistics
using CSV
using DataFrames
using Random
using Profile

function pressure_pre_calc(g_x, g_y, g_z, p_0, f0, c0, lambda, k, a, array_x, array_y, array_z)
    transducer_cont = ComplexF64[]
    for tr = 1:length(array_x)
        direction = 1
        if array_z[tr] > 0
            direction = -1
        end
        r_px = g_x - array_x[tr];
        r_py = g_y - array_y[tr];
        r_pz = g_z - array_z[tr];

        R = sqrt(r_px^2 + r_py^2 + r_pz^2)

        dotproduct = r_pz*direction;
        theta = acos(dotproduct/(R*sqrt(direction^2)));

        if theta == 0
            theta += eps(0.0);
        end

        D = 2 * besselj(1, k*a*sin(theta))/(k*a*sin(theta));
        push!(transducer_cont, (p_0/R) * D * exp(1im*(k*R)))
    end
    return transducer_cont
end

function pressure_calc(transducer_cont, phases)
    P = exp.(1im.*(phases)).*transducer_cont
    return sum(P)
end

function loss_fn(phases; transducer_matrix = transducer_matrix, tar_amp = tar_amp)
    loss_store = 0
    for N = 1:length(tar_amp)
        err_l = abs(pressure_calc(transducer_matrix[:,N], phases)) - tar_amp[N]
        loss_store += (err_l^2)
    end
    return loss_store
end

function finite_difference(phases; transducer_matrix = transducer_matrix, delta = delta, tar_amp = tar_amp)
    original_value = loss_fn(phases; transducer_matrix, tar_amp)
    finite_results = Float64[]
    for ii = 1:length(phases)
        try_phase = copy(phases)
        try_phase[ii] = try_phase[ii] + delta
        new_value = loss_fn(try_phase; transducer_matrix, tar_amp)
        push!(finite_results,(new_value - original_value) / delta)
    end
    return finite_results
end

function optimization_part(nn, geometry_array, amplitudes_array, length_trans, number_target, p_0, f0, c0, lambda, k, a, array_x, array_y, array_z, delta)
    tar_x = convert(Array, geometry_array[nn, 2:3:end])
    tar_y = convert(Array, geometry_array[nn, 3:3:end])
    tar_z = convert(Array, geometry_array[nn, 4:3:end])
    tar_amp = convert(Array, amplitudes_array[nn, 2:end])

    phases = rand(length_trans, 1);
    dopt   = Adam(phases, loss_fn, finite_difference)
    dopt.a = 0.1
    dopt.b1 = 0.9
    dopt.b2 = 0.999
    dopt.eps = 1e-8
    transducer_matrix = complex(zeros(length_trans, number_target))

    for N = 1:length(tar_x)
        transducer_matrix[:,N] = pressure_pre_calc(tar_x[N], tar_y[N], tar_z[N], p_0, f0, c0, lambda, k, a, array_x, array_y, array_z)
    end

    for i = 1:150
        step!(dopt; transducer_matrix = transducer_matrix, delta = delta, tar_amp = tar_amp)
    end
    return dopt, transducer_matrix, tar_amp
end

function main()
    sample_number = 1000
    array_store = CSV.File("tpos_14x14.csv", header=false, delim=',') |> DataFrame

    for cn = 1:5
        number_target = 2^cn

        amplitudes_array = CSV.File(string("amplitudes_", lpad(number_target,3,"0"), ".csv"), header=false, delim=',') |> DataFrame
        geometry_array = CSV.File(string("geometries_", lpad(number_target,3,"0"), ".csv"), header=false, delim=',') |> DataFrame
        array_x = array_store[!, 1];
        array_y = array_store[!, 2];
        array_z = array_store[!, 3];
        length_trans = length(array_x);

        phase_store = zeros(sample_number, length_trans+1)
        time_store = zeros(sample_number, length_trans+1)

        p_0 = 1.98
        f0 = 40e03
        c0 = 346.0
        lambda = c0/f0
        k = 2*pi*f0/c0
        a = 5e-03
        delta = 0.01*(pi/180)

        for nn = 1:sample_number 
            measured_time = @elapsed begin
                dopt, transducer_matrix, tar_amp = optimization_part(nn, geometry_array, amplitudes_array, length_trans, number_target, p_0, f0, c0, lambda, k, a, array_x, array_y, array_z, delta)
            end
            print(string("Step: ", dopt.t, " | Loss: ", dopt.loss(dopt.theta; transducer_matrix, tar_amp), "\n"))
            phase_store[nn, 1] = nn
            phase_store[nn, 2:end] = dopt.theta
            print(string("Ellapsed time: ", measured_time, "\n"))
            time_store[nn] = measured_time
        end
        CSV.write(string("phases_", lpad(number_target,3,"0"), ".csv"),  DataFrame(phase_store), writeheader=false)
        CSV.write(string("timed_", lpad(number_target,3,"0"), ".csv"),  DataFrame(time_store), writeheader=false)
    end
end

main()
