#Autor: Paweł Okrutny 221478

include("blocksys.jl")

using blocksys


action = 7


#Wylicz wektor x do pliku
#Standardowa eliminacja Gaussa
if (action == 0)
	A = readA("A.txt")
	b = readB("b.txt")
	#printMatrix(A[1])
	x = GaussElimination(A[1],A[2], b)
	#printMatrix(A[1])
	arrayToFile("x.txt", x)
end

#Wylicz wektor x do pliku
#Eliminacja Gaussa bez wyboru elementu
if (action == 1)
	A = readA("A.txt")
	b = readB("b.txt")
	#printMatrix(A[1])
	x = SpecificGaussElimination(A[1],A[2],A[3], b)
	#printMatrix(A[1])
	arrayToFile("x.txt", x)
end

#Wylicz wektor x do pliku
#Eliminacja Gaussa z czesciowym wyborem elementu glownego
if (action == 2)
	A = readA("A.txt")
	b = readB("b.txt")
	#printMatrix(A[1])
	x = SpecificGaussEliminationPivot(A[1],A[2],A[3], b)
	#printMatrix(A[1])
	arrayToFile("x.txt", x)
end

#Wylicz blad wzgledny do pliku
#Eliminacja Gaussa bez wyboru elementu
if (action == 3)
	A = readA("A.txt")
	x = ones(A[2])
	b = bFromAx(A[1],A[2],A[3], x)
	x1 = SpecificGaussElimination(A[1],A[2],A[3], b)
	relativeError(x1, x, "error.txt")
end

#Wylicz blad wzgledny do pliku
#Eliminacja Gaussa z czesciowym wyborem elementu glownego
if (action == 4)
	A = readA("A.txt")
	x = ones(A[2])
	b = bFromAx(A[1],A[2],A[3], x)
	x1 = SpecificGaussEliminationPivot(A[1],A[2],A[3], b)
	relativeError(x1, x, "error.txt")
end

#Wylicz blad wzgledny do pliku
#Standardowa eliminacja Gaussa
if (action == 5)
	A = readA("A.txt")
	x = ones(A[2])
	b = bFromAx(A[1],A[2],A[3], x)
	x1 = GaussElimination(A[1],A[2], b)
	relativeError(x1, x, "error.txt")
end

#Wylicz blad wzgledny do pliku
#Eliminacja Gaussa wbudowana w jezyku Julia
if (action == 6)
	A = readA("A.txt")
	x = ones(A[2])
	b = bFromAx(A[1],A[2],A[3], x)
	x1 = A[1]\b
	relativeError(x1, x, "error.txt")
end

#Porownanie czasow dzialania
if (action == 7)
	A = readA("A.txt")
	b = readB("b.txt")
	#SpecificGaussElimination(A[1],A[2],A[3], b)
	A = readA("A.txt")
	b = readB("b.txt")
	#SpecificGaussEliminationPivot(A[1],A[2],A[3], b)
	A = readA("A.txt")
	b = readB("b.txt")
	#GaussElimination(A[1],A[2], b)
	A[1]\b


	A = readA("A.txt")
	b = readB("b.txt")
	println("Standardowa implementacja eliminacji Gaussa")
	#@time GaussElimination(A[1],A[2], b)
	println()

	A = readA("A.txt")
	b = readB("b.txt")
	println("Eliminacja Gaussa bez wyboru elementu glownego")
	@time SpecificGaussElimination(A[1],A[2],A[3], b)
	println()
	
	A = readA("A.txt")
	b = readB("b.txt")
	println("Eliminacja Gaussa z czesciowym wyborem elementu glownego")
	@time SpecificGaussEliminationPivot(A[1],A[2],A[3], b)
	println()
	
	A = readA("A.txt")
	b = readB("b.txt")
	println("Eliminacja Gaussa wbudowana w jezyku Julia")
	@time A[1]\b

end

