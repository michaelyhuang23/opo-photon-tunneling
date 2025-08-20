using Ferrite, SparseArrays, LinearAlgebra, PolynomialRoots, Arpack, Plots, StaticArrays, Statistics


# constants
b_val = 0.5
lambda_val = 1.8
g_val = 0.1
div_val = 1

function parametrize(alpha, beta)
    u = asin(g_val * alpha / sqrt(lambda_val)) + asin(g_val * beta / sqrt(lambda_val))
    v = asin(g_val * alpha / sqrt(lambda_val)) - asin(g_val * beta / sqrt(lambda_val))
    return [u, v]
end

alpha_roots = roots([-b_val, -(lambda_val - 1), 0, g_val^2])
alpha_roots = sort(map(x -> real(x), alpha_roots))
u_roots = map(x -> parametrize(x, x)[1], alpha_roots)

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

function V_grad(u, v)
    du = tan((u + v) / 2) + tan((u - v) / 2) - g_val * b_val / sqrt(lambda_val) * (sec((u + v) / 2) + sec((u - v) / 2)) - lambda_val * sin(u)
    dv = tan((u + v) / 2) - tan((u - v) / 2) - g_val * b_val / sqrt(lambda_val) * (sec((u + v) / 2) - sec((u - v) / 2)) + lambda_val * sin(v)
    return [du, dv]
end

function beta_A_calc()
    Va = V(u_roots[1], 0)
    Vb = V(u_roots[3], 0)
    term1 = exp((Va - Vb) / g_val^2)
    term2 = sqrt(Vuu(alpha_roots[1]) * Vvv(alpha_roots[1]) / Vuu(alpha_roots[3]) / Vvv(alpha_roots[3]))
    return term1 * term2
end


beta_B = -1.0
beta_A = beta_A_calc()

function div_J(u, v)
    radius = 0.005
    if (u - u_roots[1])^2 + v^2 < radius^2
        return div_val
    elseif (u - u_roots[3])^2 + v^2 < radius^2
        return -div_val
    else
        return 0
    end
end

println("beta_A: ", beta_A)
println("beta_B: ", beta_B)

# create grid and basis functions
grid = generate_grid(Quadrilateral, (400, 400), Vec(-2.0, -1.0), Vec(2.0, -1.0), Vec(2.0, 1.0), Vec(-2.0, 1.0))
dim = 2
polynomials = Lagrange{dim,RefCube,1}()
quadrature_rules = QuadratureRule{dim,RefCube}(2)
quadrature_rules_face = QuadratureRule{dim - 1,RefCube}(2)
cell_values = CellScalarValues(quadrature_rules, polynomials)
facevalues = FaceScalarValues(quadrature_rules_face, polynomials);

# add beta field
dof_handler = DofHandler(grid)
add!(dof_handler, :beta, 1)
close!(dof_handler)

# add dirichlet boundary condition
constraint_handler = ConstraintHandler(dof_handler)
# addvertexset!(grid, "LL", x -> x[1] ≈ -2.0 && x[2] ≈ -1.0)
#l_boundary_condition = Dirichlet(:beta, getfaceset(grid, "left"), (x, t) -> beta_A)
#r_boundary_condition = Dirichlet(:beta, getfaceset(grid, "right"), (x, t) -> beta_B)
#add!(constraint_handler, l_boundary_condition)
#add!(constraint_handler, r_boundary_condition)
close!(constraint_handler)


function assemble_element!(Ke::Matrix, Fe::Vector, cellvalues::CellScalarValues, cell::CellCache, cellcount::Int, facevalues::FaceScalarValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    # Reset to 0
    fill!(Ke, 0)
    fill!(Fe, 0)
    cell_coords = getcoordinates(cell)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dOmega = getdetJdV(cellvalues, q_point)
        point_coord = spatial_coordinate(cellvalues, q_point, cell_coords)
        # Loop over test basis
        for i in 1:n_basefuncs
            g_phi = shape_value(cellvalues, q_point, i)
            g_phi_grad = shape_gradient(cellvalues, q_point, i)
            # Add contribution to Fe
            Fe[i] += g_phi * div_J(point_coord[1], point_coord[2]) * exp(V(point_coord[1], point_coord[2]) / g_val^2) * dOmega
            # Loop over trial basis
            for j in 1:n_basefuncs
                beta_phi_grad = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += g_phi * dot(beta_phi_grad, V_grad(point_coord[1], point_coord[2])) * dOmega
                Ke[i, j] += g_val^2 * dot(g_phi_grad, beta_phi_grad) * dOmega
            end
        end
    end
    # add von neumann condition
    for face in 1:nfaces(cell)
        if onboundary(cell, face) #&&
            # ((cellcount, face) in getfaceset(grid, "bottom") ||
            #  (cellcount, face) in getfaceset(grid, "top"))
            reinit!(facevalues, cell, face)
            for q_point in 1:getnquadpoints(facevalues)
                dGamma = getdetJdV(facevalues, q_point)
                point_coord = spatial_coordinate(facevalues, q_point, cell_coords)
                n = getnormal(facevalues, q_point)
                # Loop over test basis
                for i in 1:n_basefuncs
                    g_phi = shape_value(facevalues, q_point, i)
                    # Loop over trial basis
                    for j in 1:n_basefuncs
                        beta_phi_grad = shape_gradient(facevalues, q_point, j)
                        # Ke[i, j] -= g_phi * dot(beta_phi_grad, n) * dGamma
                        # comment out to impose von neummann condition
                    end
                end
            end
        end
    end
    return Ke
end

function assemble_global(cellvalues::CellScalarValues, facevalues::FaceScalarValues, K::SparseMatrixCSC, dh::DofHandler)
    # Allocate the element stiffness matrix and element force vector
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    Fe = zeros(n_basefuncs)

    f = zeros(ndofs(dof_handler))
    # Create an assembler
    assembler = start_assemble(K, f)
    # Loop over all cels
    for (cellcount, cell) in enumerate(CellIterator(dh))
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        # Compute element contribution
        assemble_element!(Ke, Fe, cellvalues, cell, cellcount, facevalues)
        # Assemble Ke and fe into K and f
        assemble!(assembler, celldofs(cell), Ke, Fe)
    end
    return K, f
end

K = create_sparsity_pattern(dof_handler)  # create structure
K, f = assemble_global(cell_values, facevalues, K, dof_handler)

apply!(K, f, constraint_handler)

beta_sol = K \ f

# beta_sol = (beta_sol .- mean(beta_sol)) ./ std(beta_sol)

#beta_sol = clamp.(beta_sol, 0.8, 1)

vtk_grid("fem_sim", dof_handler) do vtk
    vtk_point_data(vtk, dof_handler, beta_sol)
end



