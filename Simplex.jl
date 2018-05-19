using JuMP
using PyPlot
using Optim

####################################################
###Implementação do SIMPLEX - 18/05/2018############
###Renata Garcia Oliveira - Pós - 1712535###########
###renata.garcia.eng@gmail.com######################

debug = false
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
        A1, b, c1, status = simplexFaseI(c, A, b)
	if status == 1
            printlog("\nSimplex Fase II")

	    num_folgas = size(b)[1]
	    c1 = [c;zeros(num_folgas,1)]
	    A1 = hcat(A, eye(num_folgas))

	    x = [0 ; 0]
	    n = size(A)[:2]
	    m = length(b)
	    n_m = n-m

	    ind_base = [(m+1):(n+1)...]
	    ind_nobase = [1:m...]
	    printdbe("ind_base: $(ind_base)")
	    printdbe("ind_nobase: $(ind_nobase)")
            #A = A[:,vcat(ind_base,ind_nobase)]
            #c = c[vcat(ind_base,ind_nobase)]
	    #printdbe("base: $(ind_base)")
	    #printdbe("nobase: $(ind_nobase)")
	    x,z,status = simplexFaseII(c1,A1,b,ind_base,ind_nobase)
	else
            printdbe("Problema Inviável")
            return -2
        end
    end    
end

function simplexFaseI(c, A, b)
    m, n = size(A)
    nv = n - m
    
    Aw = -1 * ones(m, n + 1)
    Aw[:,1:n] = A
    bw = b
    cw = zeros(n+1)
    cw[n+1] = 1
    
    x = zeros(n+1)
    bidx = [i for i in (nv+1):(n)]
    nidx = [i for i in 1:(n+1) if !(i in bidx)]
    
    status = 3 
    it = 0
    while status > 1
        it += 1

        B = Aw[:, bidx]
        N = Aw[:, nidx]
        
        xn = zeros(length(nidx))
        d = B \ b
        db = B \ N

        xb = d - db * xn
        x[bidx] = xb
        x[nidx] = xn

        cn = cw[nidx]
        cb = cw[bidx]

        cr = cn' - cb' * db
        z = cb' * xb + (cr) * xn

        if it == 1
            nvbidx = length(nidx)
        else
            nvbidx = indmin(cr)
        end

        nvnidx = indmin(xb ./ abs.(db[:,nvbidx]))
        
        if minimum(cr) >= 0 && it > 1
            status = 1

            x = zeros(n+1)
            x[bidx] = xb

            printdbe("it: $(it)")
	    printdbe("x: $(x)")
	    printdbe("bidx: $(bidx)")
	    printdbe("nidx: $(nidx)")
	    printdbe("z: $(z)")
	    printdbe("status: $(status)")
            break
        end
       
        old_bidx = deepcopy(bidx)
        bidx[nvnidx] = nidx[nvbidx] 
        nidx[nvbidx] = old_bidx[nvnidx]
    end

    orig_bidx = [i for i in bidx if i!=(n+1)]
    orig_nidx = [i for i in nidx if i!=(n+1)]

    A1 = zeros(size(A))
    A1[:, 1:length(orig_nidx)] = A[:,orig_nidx] 
    A1[:, (length(orig_nidx)+1):n] = A[:, orig_bidx] 
    c1 = zeros(n)
    c1[1:length(orig_nidx)] = c[orig_nidx]
    c1[(length(orig_nidx)+1):end] = c[orig_bidx]

    printdbe("A1: $(A1)")
    printdbe("b: $(b)")
    printdbe("c1: $(c1)")

    return A1, b, c1, status
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


