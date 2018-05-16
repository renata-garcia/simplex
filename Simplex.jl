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
    logWriteFile(filestream, str)
end

function logCreateFile()
    filestream = open(fname, "w")
    return filestream
end

function logWriteFile(filestream, str)
    pwrite(filestream, str)
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

function simplexFaseI(c,A,b)
    f = logCreateFile()
    num_folgas = size(c)[1]
    c = [c;zeros(num_folgas,1)]
    A = hcat(A, eye(num_folgas))

    x = [0 ; 0]
    n = size(A)[:2]
    m = length(b)
    n_m = n-m

    ind_base = [(m+1):n...]
    ind_nobase = [1:m...]

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
        dB = -B\N
        printdbe("dB: $(dB)")

        z = cB'*xB
        printlog(f,"z: $(z)")
        y = B'\cB
        printdbe("y: $(y)")
        cr = (cN' - y'*N)'
        printdbe("custo reduzido (cr): $(cr)")
        #se não há cr(i) não negativo
        if (length(find(cr.>= epsilon)) == 0)
            return true
        else
            ind_maior_crj = indmax(cr)
            printdbe("ind_maior_crj: $(ind_maior_crj)")
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

    return false
end

c, A, b = def_prob1()
ret = simplexFaseI(c, A, b)
