norm(x::Float64) = x
dot(x::Float64, y::Float64) = x*y

function conjugate_gradients(f::Function,
                             g::Function,
                             x0::Any,
                             finite_difference::Float64,
                             tolerance::Float64)

    x_old = x0
    x_new = x0
    y_old = Inf
    y_new = f(x_new)
    g_old = 0.0*x0
    g_new = g(x_new)
    direction_old = 0.0*x0
    direction_new = 0.0*x0

    println("f: $y_new")

    max_iterations = 10
    i = 0

    max_step = 0.2

    while abs(y_new - y_old) > tolerance && i <= max_iterations
        x_old = x_new

        a = abs(dot(g_new, g_old))
        b = sum(g_old.^2)
        if a < 0.5 *b
            gamma = dot(-g_new, g_old - g_new) / b
        else
            gamma = 0.0
        end
        
        direction_new = -g_new + gamma * direction_old
        direction_old = direction_new
        direction_norm = direction_new / norm(direction_new)

        x_step = x_old + direction_norm * finite_difference
        g_step = g(x_step)

        df = dot(-g_new, direction_norm) - dot(-g_step, direction_norm)
        curvature = df / finite_difference

        if curvature > 0.0
            step_size = -g_new / curvature
        else
            step_size = max_step
        end

        x_new += step_size .* direction_norm

        g_old = g_new
        g_new = g(x_new)

        y_old = y_new
        y_new = f(x_new)

        i += 1

        println("f_new: $y_new x_new: $x_new")
    end

    return x_new
end
