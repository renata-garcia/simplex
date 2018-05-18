using JuMP
using PyPlot
using Optim

####################################################
###Implementação do SIMPLEX - 18/05/2018############
###Renata Garcia Oliveira - Pós - 1712535###########
###renata.garcia.eng@gmail.com######################

debug = true
epsilon = 1e-9
max_iter = 1000
filestream = open("simplex.log", "w")

function printdbe(str::String)
    if (debug == true)
        println(str)
    end
end

function printlog(str::String)
    println(str)
    logWriteFile(filestream, (string(str, '\n')))
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
    A = [2 2; 1 2; -1 -1]
    b = [4; 4; -1]
    c = [4;3]
    return c, A, b
end

function simplex(c, A, b)

    num_folgas = size(b)[1]
    c = [c;zeros(num_folgas,1)]
    A = hcat(A, eye(num_folgas))

    x = [0 ; 0]
    n = size(A)[:2]
    m = length(b)
    n_m = n-m

    ind_base = [(m+1):n...]
    ind_nobase = [1:m...]

    if all(b .> 0)
        printlog("\nSimplex Fase II")
        printdbe("c: $(c)")
        printdbe("A: $(A)")
        printdbe("b: $(b)")
	x,z,status = simplexFaseII(c, A, b, ind_base, ind_nobase)
    else
    	printlog("\nSimplex Fase I")
        printdbe("A: $(A)")
        printdbe("c: $(c)")
        printdbe("b: $(b)")
        ind_base,ind_nobase,status = simplexFaseI(c, A, b, ind_base, ind_nobase)
	if status == -2
            printdbe("Problema Inviável")
            return -2
        else
            printlog("\nSimplex Fase II")
            A = A[:,vcat(ind_base,ind_nobase)]
            c = c[vcat(ind_base,ind_nobase)]
	    printdbe("base: $(ind_base)")
	    printdbe("nobase: $(ind_nobase)")
	    x,z,status = simplexFaseII(c,A,b,ind_base,ind_nobase)
        end
    end    
end

function simplexFaseI(c, A, b, ind_base, ind_nobase)

    A = hcat(A,-1*ones(1,size(A)[1])')   
    c = vcat(zeros(length(c)),-1)

    m = size(A)[1] + 1
    n = size(A)[2] - m

    ind_nobase = vcat([i for i in 1:n],size(A)[2])
    ind_base  = [i for i in (n+1):(n+m-1)]

    B = A[:,ind_base]   
    N = A[:,ind_nobase] 
    xB = B \ b     

    j          = length(ind_nobase) 
    k          = indmin(xB)    
    r_ind_base  = ind_base[k]        
    r_ind_nobase = ind_nobase[j]      
    ind_base[k]      = r_ind_nobase   
    ind_nobase[j]      = r_ind_base  

    printlog("B: $(B)")
    printlog("xB: $(xB)")

    for i = 1:max_iter

        printdbe("Iteracao #$(i)")

        B = A[:,ind_base]
        N = A[:,ind_nobase]
        xB = B \ b       
        dB = -B \ N
        cB = c[ind_base] 
        cN = c[ind_nobase]
        cR = (cN' + cB'*dB)'
        j  = indmax(cR)     

        x      = zeros(n+m)
        x[ind_base]  = xB

        printdbe("x     = $(x)")
        printdbe("ind_base  = $(ind_base)")
        printdbe("ind_nobase = $(ind_nobase)")
        printdbe("")

        if all(cR .<= epsilon)
            ind_base
            ind_nobase = ind_nobase[ind_nobase .!= size(A)[2]]
            status = 1
            return (ind_base, ind_nobase, status)
        end

        if all(dB[:,j] .>= epsilon)
            status = -2
            return (status);
        end

        r  = xB./(dB[:,j])
        k  = indmax(r[r .< 0])

        auxB  = ind_base[k]
        ind_base[k] = ind_nobase[j]
        ind_nobase[j] = auxB
    end	
end

function simplexFaseII(c, A, b, ind_base, ind_nobase)

    for i=1:max_iter
        printlog("\niteração: $(i)")
        printlog("base: $(ind_base)")
        printlog("no base: $(ind_nobase)")

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
        printlog("xw: $(xw)")

        z = cB'*xB
        printlog("z: $(z)")
        y = B'\cB
        printdbe("y: $(y)")
        cr = (cN' - y'*N)'
        printdbe("custo reduzido (cr): $(cr)")

        if all(cr .<= epsilon)
            z = cB' * xB
            printdbe("Custo Total = $(z)")
            printdbe("Variáveis Básicas / Valores = $(ind_base) / $(xB)")
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
            printlog("Problema Irrestrito")
            printlog("Custo total infinito")
            printlog("Direção extrema $(dall)")
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

printlog("\n\n######SIMPLEX PROB 1########")
c, A, b = def_prob1()
ct, d, status = simplex(c, A, b)
printdbe("ct: $(ct)")
printdbe("d: $(d)")
printdbe("status: $(status)")


printlog("\n\n######SIMPLEX PROB 2########")
c, A, b = def_prob2()
ct, d, status = simplex(c, A, b)
printdbe("ct: $(ct)")
printdbe("d: $(d)")
printdbe("status: $(status)")


printlog("\n\n######SIMPLEX PROB 3########")
c, A, b = def_prob3()
ct, d, status = simplex(c, A, b)
printdbe("ct: $(ct)")
printdbe("d: $(d)")
printdbe("status: $(status)")


