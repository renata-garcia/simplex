using JuMP
using PyPlot
using Optim
using GLPKMathProgInterface

debug = true
epsilon = 1e-9
max_iter = 1000
fname = "simplex.log"

function printdbe(str::String)
    if (debug == true)
        println(str)
    end
end

function printlog(filestream, str::String)
    println(str)
    logWriteFile(filestream, (string(str, '\n')))
end

function logCreateFile()
    filestream = open(fname, "w")
    return filestream
end

function logWriteFile(filestream, str)
    write(filestream, str)
end

function logCloseFile(filestream)
    close(filestream)
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

function def_prob3()
    A = [3 3; 1 2; -1 -1]
    b = [4; 4; -1]
    c = [4;3]
    return c, A, b
end

function simplexFaseI(c, A, b)

    f =  logCreateFile()
    num_folgas = size(b)[1]
    println("num_folgas: $(num_folgas)")
    c = [c;zeros(num_folgas,1)]
    A = hcat(A, eye(num_folgas))

    x = [0 ; 0]
    n = size(A)[:2]
    m = length(b)
    n_m = n-m

    ind_base = [(m+1):n...]
    ind_nobase = [1:m...]

    if any(b .> 0)
        printlog(f,"\nSimplex Fase II")
        printdbe("c: $(c)")
        printdbe("A: $(A)")
        printdbe("b: $(b)")
    else
        i = indmin(b)
        ind_base[i] = n + m + 1
        ind_nobase = vcat(nobase, n + i)

        c = [c;-1]
        A = hcat(A, -ones(m))

        r1,r2,r3, ind_base_2,ind_nobase_2 = SimplexFaseII(f, c, A, b, ind_base, ind_nobase)

        printlog(f,"\nSimplex Fase I")

        i = indmax(ind_nobase_2)
        ind_nobase_2 = deleteat!(ind_nobase_2,i)

        x,z,status = SimplexFaseII(f, c, A, b,ind_base_2,ind_nobase_2)
        printlog(f,"\nSimplex Fase II")
    end

    a = simplexFaseII(f, c, A, b, ind_base, ind_nobase)
end

function simplexFaseII(f, c, A, b, ind_base, ind_nobase)

    for i=1:max_iter
        printlog(f,"\niteração: $(i)")
        printlog(f,"base: $(ind_base)")
        printlog(f,"no base: $(ind_nobase)")

        B = A[:, ind_base]
        printdbe("B: $(B)")
        N = A[:, ind_nobase]
        printdbe("N: $(N)")
        cB = c[ind_base]
        printdbe("cB: $(cB)")
        cN = c[ind_nobase]
        printdbe("cN: $(cN)")

        xB = B\b
        printdbe("xB: $(xB)")
        xw = zeros(size(A)[:2])
        xw[ind_base] = xB
        printlog(f, "xw: $(xw)")

        z = cB'*xB
        printlog(f,"z: $(z)")
        y = B'\cB
        printdbe("y: $(y)")
        cr = (cN' - y'*N)'
        printdbe("custo reduzido (cr): $(cr)")

        if all(cr .<= epsilon)
            z = cB' * xB
            println("Custo Total = ", z)
            println("Variáveis Básicas / Valores = ", ind_base, "/", xB)
            status=1
            return xB,z,status,ind_base,ind_nobase
        end

        cmax, j = findmax(cr)
        printdbe("cmax: $(cmax)")
        dB = -B\N[:,j]
        printdbe("dB: $(dB)")

        dall = zeros(length(c))
        dall[ind_nobase[j]] = 1
        dall[ind_base] = dB
        printdbe("dall: $(dall)")

        if all(dall .>= 0)
            printlog(f, "Problema Irrestrito")
            printlog(f, "Custo total infinito")
            printlog(f, "Direção extrema $(dall)")
            return Inf, dall, -1
        end
        #se não há cr(i) não negativo
        if (length(find(cr.>= epsilon)) == 0)
            return z, xw, 1
        else
            ind_maior_crj = indmax(cr)
            printdbe("ind_maior_crj: $(ind_maior_crj)")
            dB = -B\N
            r = xB./dB[:,ind_maior_crj]
            printdbe("r: $(r)")
            ind = find(r.>=0)
            r[ind] = NaN
            printdbe("r: $(r)")
            i = indmin(abs.(r))

            aux = ind_base[i]
            printdbe("aux: $(aux)")
            ind_base[i] = ind_nobase[ind_maior_crj]
            printdbe("ind_base[ind_maior_crj]: $(ind_base[ind_maior_crj])")
            ind_nobase[ind_maior_crj] = aux
            printdbe("ind_nobase[ind_maior_crj]: $(ind_nobase[ind_maior_crj])")
        end
    end
    logCloseFile()
    return z, x, 1
end

c, A, b = def_prob3()
ct, d, status = simplexFaseI(c, A, b)
printdbe("ct: $(ct)")
printdbe("d: $(d)")
printdbe("status: $(status)")
