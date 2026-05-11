function write_results(path, state; step = 0)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# PhaseFieldFracture placeholder output")
        println(io, "step = ", step)
        println(io, "state = ", state)
    end
    return path
end
