function solve_block_system(Kuu, Kud, Kdu, Kdd, ru, rd)
    K = [Kuu Kud; Kdu Kdd]
    r = [ru; rd]
    return K \ r
end
