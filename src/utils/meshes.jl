function make_l_shape_mesh(; h = 0.05)
    return (; kind = :l_shape, h)
end

function make_ct_specimen_mesh(; h = 0.05)
    return (; kind = :ct_specimen, h)
end
