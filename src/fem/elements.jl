abstract type AbstractElement end

struct Quad4 <: AbstractElement end
struct Tri3 <: AbstractElement end

function shape_functions(::Quad4, xi, eta)
    return 0.25 .* [
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ]
end

function shape_gradients(::Quad4, xi, eta)
    return 0.25 .* [
        -(1 - eta)  -(1 - xi)
         (1 - eta)  -(1 + xi)
         (1 + eta)   (1 + xi)
        -(1 + eta)   (1 - xi)
    ]
end

shape_functions(::Tri3, xi, eta) = [1 - xi - eta, xi, eta]

shape_gradients(::Tri3, xi, eta) = [
    -1 -1
     1  0
     0  1
]
