function gauss_rule(::Quad4, order::Integer = 2)
    if order != 2
        throw(ArgumentError("Only 2x2 Gauss quadrature is available for Quad4."))
    end
    a = inv(sqrt(3))
    return [(-a, -a), (a, -a), (a, a), (-a, a)], ones(4)
end

function gauss_rule(::Tri3, order::Integer = 1)
    if order != 1
        throw(ArgumentError("Only one-point quadrature is available for Tri3."))
    end
    return [(1 / 3, 1 / 3)], [0.5]
end
