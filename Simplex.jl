
using JuMP
using Cbc
using JuMP
using PyPlot
using Optim
using GLPKMathProgInterface

debug = true

function printdbe(str::String)
    if (debug == true)
        println(str)
    end
end

function def_prob1()
    c = [4; 3]
    A = [2 1; 1 2]
    b = [4; 4]
    return c, A, b
end

function def_prob2()
    c = [ 1; 1]

    A = [0.5 -1; -4 1]
    b = [0.5; 1.0]
    return c, A, b
end

c, A, b = def_prob1()

num_folgas = size(c)[1]
c = [c;zeros(num_folgas,1)]
A = hcat(A, eye(num_folgas))

x = [0 ; 0]
n = size(A)[:2]
m = length(b)
#n_m = (length(c)-length(b))
n_m = n-m

#z = [0]
#z = [z; c[1:n_m]]

ind_base = [(m+1):n...]
printdbe("ind_base: $(ind_base)")
ind_nobase = [1:m...]
printdbe("ind_nobase: $(ind_nobase)")
#println("z: ", z)

println("\n\n\n")
custo_reduzido_nao_otimo = false
while(!custo_reduzido_nao_otimo)
    B = A[:, ind_base]
    printdbe("B: $(B)")
    N = A[:, ind_nobase]
    printdbe("N: $(N)")
    cB = c[ind_base]
    printdbe("cB: $(cB)")
    cN = c[ind_nobase]
    printdbe("cN: $(cN)")

    xB = B\b
    println("xB: ", xB)
    dB = -B\N
    println("dB: ", dB)

    z = cB'*xB
    println("z: ", z)
    y = B'\cB
    println("y: ", y)
    cr = (cN' - y'*N)'
    println("custo reduzido (cr): ", cr)
    #se não há cr(i) não negativo
    if (length(find(cred.>=0)) == 0)
        custo_reduzido_nao_otimo = true
    else
        ind_maior_crj = indmax(cr)
        printdb("ind_maior_crj: ", ind_maior_crj)
        r = xB./dB[:,j]
        println("r: ", r)
        ind = find(r.>=0)
        r[ind] = NaN
        println("r: ", r)
        i = indmin(abs.(r))

        aux = ind_base[i]
        println("aux: ", aux)
        ind_nobase[j]
        ind_base[i] =ind_nobase[j]
        println("ind_base[j]: ", ind_base[j])
        ind_nobase[j] = aux
        println("ind_nobase[j]: ", ind_nobase[j])
    end
    println("\n");
end
