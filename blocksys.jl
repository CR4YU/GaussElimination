#Autor: Paweł Okrutny 221478

module blocksys

	export 
		readA, 
		readB,
		GaussElimination,
		SpecificGaussElimination,
		SpecificGaussEliminationPivot,
		arrayToFile,
		bFromAx,
		relativeError,
		printMatrix

	function SpecificGaussEliminationPivot(A::SparseMatrixCSC, n::Int64, l::Int64, b::Array)
	#Funkcja realizujaca algorytm eliminacji Gaussa z czesciowym wyborem elementu glownego
	#Algorytm uwzglednia specyficzna postac macierzy
	#Input:
	#     A - macierz
	#     n - rozmiar macierzy A
	#     l - rozmiar blokow w macierzy A
	#     b - wektor prawych stron
	#Return:
	#     wektor x jako rozwiazanie ukladu Ax = b

		for k=1 : n-1					#petla zewnetrzna, poruszanie sie po skosie

			if (div(k,l)*l + l < n)		#ograniczenie dla petli (unikanie zerowych elementow)
				z = div(k,l)*l + l
			else
				z = n
			end

			if(k + 2*l < n)				#ograniczenie dla petli (unianie zerowych elementow)
				maxrange = k + 2*l
			else
				maxrange = n
			end

			max = abs(A[k,k])
        	id = k
			for i=k+1 : z
          		if (abs(A[i,k]) > max)	#szukanie maxa w kolumnie
            		max = abs(A[i,k])
            		id = i
          		end
        	end
        	if id != k
        		b[id], b[k] = b[k], b[id]	#zamiana wartosci w wektorze b
        		for j=k : maxrange
        			A[k,j],A[id,j] = A[id,j],A[k,j]	#zamiana wierszy w macierzy
        		end
        	end

			for i=k+1 : z
				m = A[i,k]/A[k,k]
				for j=k : maxrange
					A[i,j] = A[i,j] - m*A[k,j]	#eliminacja
				end
				b[i] = b[i] - m*b[k]
			end

		end

		x = Array(Float64, n)
		for k=n : -1 : 1
			s = b[k]
			if (k>n-2*l)
				from = n
			else
				from = k + 2*l
			end
			for j=from : -1 : k+1
				s = s - A[k,j]*x[j]	#odczytywanie wartosci wektora x
			end
			x[k] = s/A[k,k]
		end
		return x
	end #of Gauss


	function SpecificGaussElimination(A::SparseMatrixCSC, n::Int64, l::Int64, b::Array)
	#Algorytm eliminacji Gaussa przystosowany do specyficznej postaci macierzy
	#Input:
	#     A - macierz
	#     n - rozmiar macierzy A
	#     l - rozmiar blokow w macierzy A
	#     b - wektor prawych stron
	#Return:
	#     wektor x jako rozwiazanie ukladu Ax = b
		for k=1 : n-1
			if (k+l < n)
				z = k+l
			else
				z = n
			end			
			if (div(k,l)*l + l < n)
					z2 = div(k,l)*l + l
				else
					z2 = n
				end
				if (A[k,k] == 0.0)
					println("Blad. Dzielenie przez zero")
					return [0]
				end
			for i=k+1 : z2
				m = A[i,k]/A[k,k]
				for j=k : z
					A[i,j] =A[i,j] - m*A[k,j]
				end
				b[i] = b[i] - m*b[k]
			end
		end
		x = Array(Float64, n)
		for k=n : -1 : 1
			s = b[k]
			if (k > n-l)
				point = n
			else
				point = k + l
			end
			for j=point : -1 : k+1
				s = s - A[k,j]*x[j]
			end
			x[k] = s/A[k,k]
		end
		return x
	end #of Gauss


	function GaussElimination(A::SparseMatrixCSC, n::Int64, b::Array)
	#Standardowa implementacja algorytmu eliminacji Gaussa
		for k=1 : n-1
			for i=k+1 : n
				m = A[i,k]/A[k,k]
				for j=k : n
					A[i,j] = A[i,j] - m*A[k,j]
				end
				b[i] = b[i] - m*b[k]
			end
		end
		x = Array(Float64, n)
		for k=n : -1 : 1
			s = b[k]
			for j=n : -1 : k+1
				s = s - A[k,j]*x[j]
			end
			x[k] = s/A[k,k]
		end
		return x
	end #of Gauss


	function readA(filename::String)
	#Czytanie macierzy A z zadanego pliku
		lines = readdlm(filename,' ')
		n = lines[1,1]
		l = lines[1,2]
		I = Int[]
		J = Int[]
		V = Float64[]
		for i=2 : (size(lines)[1])
			push!(I,lines[i,1])
			push!(J,lines[i,2])
			push!(V,lines[i,3])
		end
		A = sparse(I,J,V)
		return A, n, l
	end

	function readB(filename::String)
	#Czytanie wektora prawych stron b z pliku
		lines = readdlm(filename,' ')
		b = Float64[]
		for i=2 : (size(lines)[1])
			push!(b,lines[i,1])
		end
		return b
	end

	function arrayToFile(output::String, x::Array)
	#Zapisywanie wektora do zadanego pliku
		f = open(output, "w")
		for i=1 : length(x)
			write(f,"$(x[i])\n")
		end
	end

	function bFromAx(A::SparseMatrixCSC, n::Int64, l::Int64, x::Array)
	#Oblczanie Ax, zwraca wektor prawych stron
		b = Array(Float64, n)
		for i=1 : n
			s = 0.0
			if (i + l < n)
				to = i+l
			else
				to = n
			end
			if (i <= l)
				from = 1
			else
				from = div(i-1,l)*l
			end
			for k=from : to
				s += A[i,k] * x[k]
			end
			b[i] = s
		end
	return b
	end

	function relativeError(x1::Array, x2::Array, output::String)
	#Zapisuje do pliku blad wzgledny
		if(length(x1) != length(x2))
			return "error"
		end
		error = norm(abs(x1-x2))/norm(abs(x1))
		f = open(output, "w")
		write(f, "$error \n")
		for i=1 : length(x1)
			write(f,"$(x1[i])\n")
		end
	end

	function printMatrix(A::SparseMatrixCSC)
	#Wyswietla macierz A
		n = Int(sqrt(length(A)))
		for i=1 : n 
			print("    $i")
		end
		println("\n")
		for i=1 : n
			print("$i  ")
			for j=1 : n
				if i == j print("|") end
				@printf("%.1f  ", A[i,j])
			end
			println("\n")
		end
	end

end #of module